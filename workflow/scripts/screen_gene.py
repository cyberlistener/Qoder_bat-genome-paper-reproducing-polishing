#!/usr/bin/env python3
"""Run one gene through QC -> tree pruning -> aBSREL -> tidy branch table.

This is the per-gene unit of work in the Snakefile.  It is deliberately
idempotent and writes a status JSON even when a gene is skipped, so a run over
17k genes does not stop because one alignment is malformed.

Important: the species tree must have branch lengths in substitutions/site.
aBSREL is a branch-site model; coalescent units from raw ASTRAL output are not
the same quantity.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Tuple

import dendropy
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (  # noqa: E402
    codon_qc,
    now_utc,
    parse_rate_distributions,
    read_fasta,
    write_fasta,
    write_json,
    write_tsv,
)

BRANCH_COLUMNS = [
    "gene_id",
    "branch_id",
    "node_label",
    "tested",
    "corrected_p",
    "uncorrected_p",
    "lrt",
    "baseline_omega",
    "rate_classes",
    "omega_max",
    "prop_sites_omega_gt1",
    "selected",
    "mapping_ok",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gene", required=True)
    p.add_argument("--alignment", required=True)
    p.add_argument("--tree", required=True, help="Labeled species tree (substitutions/site)")
    p.add_argument("--species-metadata", required=True)
    p.add_argument("--out-branch", required=True)
    p.add_argument("--out-status", required=True)
    p.add_argument("--min-taxa-fraction", type=float, default=0.6)
    p.add_argument("--min-taxa-absrel", type=int, default=4)
    p.add_argument("--min-codons", type=int, default=0)
    p.add_argument("--allow-internal-stops", action="store_true")
    p.add_argument("--max-invalid-fraction", type=float, default=0.05)
    p.add_argument("--p-threshold", type=float, default=0.05)
    p.add_argument("--branches", default="All")
    p.add_argument("--genetic-code", default="Universal")
    p.add_argument("--hyphy", default="hyphy")
    p.add_argument("--timeout-sec", type=int, default=7200)
    p.add_argument("--keep-temp", action="store_true")
    return p.parse_args()


def load_tree(path: str) -> dendropy.Tree:
    return dendropy.Tree.get(
        path=path,
        schema="newick",
        preserve_underscores=True,
        rooting="default-rooted",
    )


def prune_tree(full_tree: dendropy.Tree, taxa: List[str]) -> dendropy.Tree:
    # Re-parsing the Newick is cheap at this scale and guarantees that the
    # pruned copy has its own taxon namespace.
    pruned = dendropy.Tree.get(
        data=full_tree.as_string(schema="newick"),
        schema="newick",
        preserve_underscores=True,
        rooting="default-rooted",
    )
    pruned.retain_taxa_with_labels(taxa)
    pruned.suppress_unifurcations()
    return pruned


def label_pruned_tree(pruned: dendropy.Tree, full_tree: dendropy.Tree) -> Dict[str, str]:
    """Give every pruned node a unique label and map it to a full-tree branch id.

    The map is by descendant leaf set -> MRCA in the full tree.  This is stable
    when a gene lacks some species: the remaining branch still points at the
    same full-tree node.
    """
    leaf_labels = {leaf.taxon.label for leaf in pruned.leaf_node_iter()}
    used = set(leaf_labels)
    for node in pruned.preorder_node_iter():
        if node.is_leaf():
            node.label = node.taxon.label
            continue
        idx = 0
        while True:
            idx += 1
            candidate = f"n{idx:05d}"
            if candidate not in used:
                break
        node.label = candidate
        used.add(candidate)

    label_to_full: Dict[str, str] = {}
    full_ids = {}
    for node in pruned.preorder_node_iter():
        if node.parent_node is None:
            continue
        desc = [leaf.taxon.label for leaf in node.leaf_iter()]
        mrca = full_tree.mrca(taxon_labels=desc)
        if mrca is None or mrca.label is None:
            raise ValueError(f"cannot map pruned branch {node.label!r} back to the full tree")
        full_id = str(mrca.label)
        if node.label in label_to_full and label_to_full[node.label] != full_id:
            raise ValueError(f"label collision for {node.label!r}")
        label_to_full[str(node.label)] = full_id
        full_ids.setdefault(full_id, []).append(str(node.label))

    duplicates = {k: v for k, v in full_ids.items() if len(v) > 1}
    if duplicates:
        raise ValueError(f"multiple pruned branches map to one full-tree branch: {duplicates}")
    return label_to_full


def run_absrel(args: argparse.Namespace, aln: Path, tree: Path, out_json: Path, log_path: Path) -> Tuple[bool, str]:
    cmd = [
        args.hyphy,
        "absrel",
        "--alignment", str(aln),
        "--tree", str(tree),
        "--branches", args.branches,
        "--output", str(out_json),
        "--code", args.genetic_code,
    ]
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    started = time.time()
    with open(log_path, "wt") as log:
        log.write("CMD: " + " ".join(cmd) + "\n")
        log.flush()
        try:
            proc = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=args.timeout_sec,
                env=env,
                check=False,
            )
        except subprocess.TimeoutExpired:
            return False, f"aBSREL timed out after {args.timeout_sec}s"
        except FileNotFoundError:
            return False, f"aBSREL executable not found: {args.hyphy}"
    if proc.returncode != 0:
        return False, f"aBSREL exited with status {proc.returncode}"
    if not out_json.exists() or out_json.stat().st_size == 0:
        return False, "aBSREL did not create a non-empty JSON output"
    return True, f"aBSREL finished in {time.time() - started:.1f}s"


def parse_absrel(json_path: Path, label_to_full: Dict[str, str], gene: str, p_threshold: float):
    with open(json_path, "rt") as fh:
        j = json.load(fh)

    ba = j.get("branch attributes", {})
    if isinstance(ba, dict) and "0" in ba:
        ba = ba["0"]
    tested_map = j.get("tested", {})
    if isinstance(tested_map, dict) and "0" in tested_map:
        tested_map = tested_map["0"]
    if not isinstance(ba, dict) or not ba:
        raise ValueError("aBSREL JSON has no 'branch attributes' object")

    rows = []
    mapping_failures = []
    for key, attrs in ba.items():
        if not isinstance(attrs, dict):
            continue
        label = str(attrs.get("original name") or key)
        full_id = label_to_full.get(label)
        mapping_ok = full_id is not None
        if full_id is None:
            full_id = label_to_full.get(str(key))
            mapping_ok = full_id is not None
        if full_id is None:
            full_id = f"UNMAPPED::{label}"
            mapping_failures.append(label)

        tested_value = tested_map.get(key, tested_map.get(label, "background")) if isinstance(tested_map, dict) else None
        cp = attrs.get("Corrected P-value")
        up = attrs.get("Uncorrected P-value")
        cp = float(cp) if cp not in (None, "", "NA") else math.nan
        up = float(up) if up not in (None, "", "NA") else math.nan
        tested = (str(tested_value).lower() == "test") if tested_value is not None else np.isfinite(cp)
        omega_max, prop_gt1, n_classes = parse_rate_distributions(attrs.get("Rate Distributions"))
        selected = bool(
            tested
            and np.isfinite(cp)
            and cp < p_threshold
            and np.isfinite(omega_max)
            and omega_max > 1.0
        )
        rows.append(
            {
                "gene_id": gene,
                "branch_id": full_id,
                "node_label": label,
                "tested": tested,
                "corrected_p": cp,
                "uncorrected_p": up,
                "lrt": attrs.get("LRT"),
                "baseline_omega": attrs.get("Baseline MG94xREV omega ratio"),
                "rate_classes": attrs.get("Rate classes", n_classes),
                "omega_max": omega_max,
                "prop_sites_omega_gt1": prop_gt1,
                "selected": selected,
                "mapping_ok": mapping_ok,
            }
        )
    df = pd.DataFrame(rows, columns=BRANCH_COLUMNS)
    version = str(j.get("analysis", {}).get("version", ""))
    return df, version, mapping_failures


def empty_branch_table() -> pd.DataFrame:
    return pd.DataFrame(columns=BRANCH_COLUMNS)


def main() -> int:
    args = parse_args()
    started = time.time()
    gene = args.gene
    status: Dict[str, object] = {
        "gene_id": gene,
        "status": "error",
        "reason": "",
        "qc": {},
        "n_branches_total": 0,
        "n_branches_tested": 0,
        "n_branches_selected": 0,
        "mapping_failures": [],
        "hyphy_version": "",
        "started_utc": now_utc(),
    }

    def finish_status() -> None:
        status["runtime_sec"] = round(time.time() - started, 2)
        status["finished_utc"] = now_utc()
        write_json(status, args.out_status)

    try:
        seqs_all = read_fasta(args.alignment)
        full_tree = load_tree(args.tree)
        tree_taxa = {leaf.taxon.label for leaf in full_tree.leaf_node_iter()}
        meta = pd.read_csv(args.species_metadata, sep="\t", dtype=str)
        if "taxon" not in meta.columns:
            raise ValueError("species metadata has no 'taxon' column")
        metadata_taxa = set(meta["taxon"].astype(str))
        if not tree_taxa <= metadata_taxa:
            missing = sorted(tree_taxa - metadata_taxa)
            raise ValueError(f"tree taxa missing from species metadata: {missing[:10]}")

        common = [t for t in seqs_all if t in tree_taxa]
        dropped_aln = sorted(set(seqs_all) - set(common))
        missing_tree = sorted(tree_taxa - set(common))
        status["qc"] = {
            "alignment_taxa": len(seqs_all),
            "tree_taxa": len(tree_taxa),
            "retained_taxa": len(common),
            "dropped_alignment_taxa": dropped_aln,
            "tree_taxa_absent_from_alignment": missing_tree,
        }

        min_taxa = max(args.min_taxa_absrel, int(math.ceil(args.min_taxa_fraction * len(tree_taxa))))
        if len(common) < min_taxa:
            status.update(
                status="skipped",
                reason=f"only {len(common)} taxa after intersecting alignment and tree; need >= {min_taxa}",
            )
            write_tsv(empty_branch_table(), args.out_branch)
            finish_status()
            return 0

        seqs = {t: seqs_all[t] for t in common}
        qc = codon_qc(seqs)
        status["qc"].update(qc)
        if qc["n_codons"] < args.min_codons:
            status.update(status="skipped", reason=f"alignment has {qc['n_codons']} codons; need >= {args.min_codons}")
            write_tsv(empty_branch_table(), args.out_branch)
            finish_status()
            return 0
        if qc["internal_stops"] and not args.allow_internal_stops:
            status.update(
                status="skipped",
                reason=f"{qc['internal_stops']} internal stop codons; mask them or pass --allow-internal-stops",
            )
            write_tsv(empty_branch_table(), args.out_branch)
            finish_status()
            return 0
        if qc["invalid_fraction"] > args.max_invalid_fraction:
            status.update(
                status="skipped",
                reason=f"invalid characters fraction {qc['invalid_fraction']:.3f} > {args.max_invalid_fraction}",
            )
            write_tsv(empty_branch_table(), args.out_branch)
            finish_status()
            return 0

        tempdir = Path(tempfile.mkdtemp(prefix=f"absrel_{gene}_"))
        try:
            pruned = prune_tree(full_tree, common)
            if len([1 for _ in pruned.leaf_node_iter()]) != len(common):
                raise ValueError("pruning did not retain the expected taxa")
            label_to_full = label_pruned_tree(pruned, full_tree)
            aln_path = tempdir / f"{gene}.fa"
            tree_path = tempdir / f"{gene}.nwk"
            json_path = tempdir / f"{gene}.ABSREL.json"
            hyphy_log = tempdir / f"{gene}.hyphy.log"
            write_fasta(aln_path, seqs, order=common)
            pruned.write(path=tree_path, schema="newick", suppress_rooting=True)

            ok, message = run_absrel(args, aln_path, tree_path, json_path, hyphy_log)
            if not ok:
                tail = ""
                if hyphy_log.exists():
                    tail = "\n".join(hyphy_log.read_text(errors="replace").splitlines()[-20:])
                status.update(status="error", reason=f"{message}\n{tail}")
                write_tsv(empty_branch_table(), args.out_branch)
                finish_status()
                return 0

            df, version, failures = parse_absrel(json_path, label_to_full, gene, args.p_threshold)
            write_tsv(df, args.out_branch)
            status.update(
                status="ok",
                reason=message,
                n_branches_total=int(len(df)),
                n_branches_tested=int(df["tested"].sum()),
                n_branches_selected=int(df["selected"].sum()),
                mapping_failures=failures,
                hyphy_version=version,
            )
            if args.keep_temp:
                print(f"[info] kept temp dir: {tempdir}", file=sys.stderr)
        finally:
            if not args.keep_temp:
                shutil.rmtree(tempdir, ignore_errors=True)

    except Exception as exc:  # keep the run alive and make the failure inspectable
        status.update(status="error", reason=f"{type(exc).__name__}: {exc}")
        write_tsv(empty_branch_table(), args.out_branch)
        finish_status()
        return 0

    finish_status()
    print(
        f"[{status['status']}] {gene}: {status.get('n_branches_tested', 0)} tested, "
        f"{status.get('n_branches_selected', 0)} selected",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
