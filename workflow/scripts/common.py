"""Small shared helpers for the per-gene screen and the aggregators.

The code intentionally avoids heavy dependencies in the QC layer so that a
single gene can be inspected without a full Snakemake run.
"""
from __future__ import annotations

import gzip
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

STOP_CODONS = {"TAA", "TAG", "TGA"}
ALLOWED_NT = set("ACGTNRYKMSWBDHV-?.")


def open_text(path: str | os.PathLike):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_fasta(path: str | os.PathLike) -> Dict[str, str]:
    """Read a FASTA file into an ordered dict; headers are split on whitespace."""
    seqs: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                if not name:
                    raise ValueError(f"Empty FASTA header in {path}")
                if name in seqs:
                    raise ValueError(f"Duplicate FASTA id {name!r} in {path}")
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before first header in {path}")
                chunks.append(line)
    if name is not None:
        seqs[name] = "".join(chunks).upper()
    if not seqs:
        raise ValueError(f"No sequences in {path}")
    return seqs


def write_fasta(path: str | os.PathLike, seqs: Mapping[str, str], order: Optional[Sequence[str]] = None) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(order) if order is not None else list(seqs)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wt") as fh:
        for key in keys:
            fh.write(f">{key}\n")
            seq = seqs[key]
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")
    os.replace(tmp, path)


def codon_qc(seqs: Mapping[str, str]) -> Dict[str, float | int]:
    """Basic checks that must pass before HyPhy is allowed to see an alignment."""
    if len(seqs) < 2:
        raise ValueError("At least two sequences are required")
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError(f"Sequences have different lengths: {sorted(lengths)}")
    n_sites = lengths.pop()
    if n_sites == 0:
        raise ValueError("Empty alignment")
    if n_sites % 3 != 0:
        raise ValueError(f"Alignment length {n_sites} is not a multiple of 3")

    n_codons = n_sites // 3
    total_chars = n_sites * len(seqs)
    gap_chars = 0
    n_chars = 0
    invalid_chars = 0
    internal_stops = 0
    terminal_stops = 0
    for seq in seqs.values():
        gap_chars += seq.count("-") + seq.count(".")
        n_chars += seq.count("N") + seq.count("n")
        invalid_chars += sum(1 for ch in seq if ch not in ALLOWED_NT)
        for i in range(n_codons):
            codon = seq[3 * i:3 * i + 3]
            if codon in STOP_CODONS:
                if i == n_codons - 1:
                    terminal_stops += 1
                else:
                    internal_stops += 1

    variable_sites = 0
    variable_codons = 0
    for i in range(n_codons):
        column = [s[3 * i:3 * i + 3] for s in seqs.values()]
        if len({c for c in column if set(c) <= set("ACGT")}) > 1:
            variable_codons += 1
        for j in range(3):
            nts = {s[3 * i + j] for s in seqs.values() if s[3 * i + j] not in "-.N?"}
            if len(nts) > 1:
                variable_sites += 1

    return {
        "n_taxa": len(seqs),
        "n_codons": n_codons,
        "n_sites": n_sites,
        "gap_fraction": gap_chars / total_chars if total_chars else math.nan,
        "n_fraction": n_chars / total_chars if total_chars else math.nan,
        "invalid_fraction": invalid_chars / total_chars if total_chars else math.nan,
        "internal_stops": internal_stops,
        "terminal_stops": terminal_stops,
        "variable_sites": variable_sites,
        "variable_codons": variable_codons,
    }


def bh_fdr(pvalues: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg q-values; NaNs are returned unchanged."""
    p = np.asarray(list(pvalues), dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0.0, 1.0)
    q[ok] = out
    return q


def parse_rate_distributions(rates) -> Tuple[float, float, int]:
    """Return (max omega, proportion of sites with omega>1, number of classes).

    HyPhy writes a list of [omega, proportion] pairs for aBSREL.  Some versions
    wrap the list in a dict, so both shapes are accepted.
    """
    if rates is None:
        return math.nan, 0.0, 0
    pairs: List[Tuple[float, float]] = []
    if isinstance(rates, dict):
        if "omega" in rates and "proportion" in rates:
            pairs.append((float(rates["omega"]), float(rates["proportion"])))
        else:
            for value in rates.values():
                if isinstance(value, (list, tuple)) and len(value) >= 2:
                    pairs.append((float(value[0]), float(value[1])))
    elif isinstance(rates, (list, tuple)):
        for value in rates:
            if isinstance(value, (list, tuple)) and len(value) >= 2:
                pairs.append((float(value[0]), float(value[1])))
    pairs = [(w, p) for w, p in pairs if math.isfinite(w) and math.isfinite(p)]
    if not pairs:
        return math.nan, 0.0, 0
    omega_max = max(w for w, _ in pairs)
    prop_gt1 = sum(p for w, p in pairs if w > 1.0)
    return omega_max, prop_gt1, len(pairs)


def write_tsv(df: pd.DataFrame, path: str | os.PathLike) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, path)


def write_json(obj, path: str | os.PathLike) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "wt") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)
    os.replace(tmp, path)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def as_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y", "bat"}
