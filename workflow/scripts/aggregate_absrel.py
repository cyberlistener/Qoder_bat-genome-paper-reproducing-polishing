#!/usr/bin/env python3
"""Aggregate per-gene aBSREL tables into branch-, gene- and run-level tables.

The paper corrected for multiple branches *within* a gene (aBSREL's BH over
branches).  This script keeps that primary definition and adds an OPTIONAL
gene-level BH column (`q_gene`, `selected_any_fdr`) as a sensitivity analysis.
Do not silently treat the optional column as the paper's primary result.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import bh_fdr, now_utc, write_json, write_tsv  # noqa: E402

META_COLS = [
    "branch_id",
    "original_label",
    "branch_type",
    "branch_length",
    "n_descendant_taxa",
    "order",
    "taxon_group",
    "is_order_stem",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--branch-dir", required=True)
    p.add_argument("--status-dir", required=True)
    p.add_argument("--branch-metadata", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--p-threshold", type=float, default=0.05)
    p.add_argument("--fdr", type=float, default=0.05)
    return p.parse_args()


def read_statuses(status_dir: str) -> pd.DataFrame:
    rows = []
    for path in sorted(glob.glob(os.path.join(status_dir, "*.status.json"))):
        try:
            with open(path, "rt") as fh:
                obj = json.load(fh)
            if not isinstance(obj, dict):
                raise ValueError("status JSON is not an object")
        except Exception as exc:
            obj = {
                "gene_id": Path(path).name.replace(".status.json", ""),
                "status": "error",
                "reason": f"unreadable status: {exc}",
            }
        obj.setdefault("gene_id", Path(path).name.replace(".status.json", ""))
        obj.setdefault("status", "error")
        obj["status_file"] = path
        rows.append(obj)
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["gene_id", "status", "reason", "runtime_sec", "hyphy_version"])
    out["status"] = out["status"].fillna("error").astype(str)
    return out


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.branch_dir, "*.branch.tsv")))
    if not files:
        raise SystemExit(f"no *.branch.tsv files in {args.branch_dir}")

    frames: List[pd.DataFrame] = []
    for path in files:
        try:
            df = pd.read_csv(
                path,
                sep="\t",
                dtype={"gene_id": "string", "branch_id": "string", "node_label": "string"},
            )
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        frames.append(df)
    if not frames:
        raise SystemExit("all per-gene tables are empty")
    all_df = pd.concat(frames, ignore_index=True)
    del frames

    for col in ["gene_id", "branch_id", "node_label"]:
        all_df[col] = all_df[col].astype("string")
    for col in ["tested", "selected", "mapping_ok"]:
        all_df[col] = all_df[col].astype(str).str.lower().isin(["true", "t", "1", "yes"])
    all_df["corrected_p"] = pd.to_numeric(all_df["corrected_p"], errors="coerce")
    all_df["omega_max"] = pd.to_numeric(all_df["omega_max"], errors="coerce")

    meta = pd.read_csv(args.branch_metadata, sep="\t", dtype={"branch_id": "string"})
    missing_meta = [c for c in META_COLS if c not in meta.columns]
    if missing_meta:
        raise SystemExit(f"branch metadata is missing columns: {missing_meta}")
    meta_small = meta[META_COLS].copy()
    all_df = all_df.merge(meta_small, on="branch_id", how="left", validate="many_to_one")
    all_df["mapping_ok"] = all_df["mapping_ok"].fillna(False).astype(bool)
    all_df["order"] = all_df["order"].fillna("UNMAPPED")
    all_df["taxon_group"] = all_df["taxon_group"].fillna("UNMAPPED")

    tested = all_df[all_df["tested"]].copy()
    tested["neg_log10_p"] = tested["corrected_p"].map(
        lambda p: -math.log10(p) if pd.notna(p) and p > 0 else (300.0 if pd.notna(p) and p == 0 else np.nan)
    )

    # Gene-level table -----------------------------------------------------
    gene_rows = []
    for gene, g in all_df.groupby("gene_id", sort=False):
        sel = g["selected"]
        selected_orders = sorted(
            {str(o) for o in g.loc[sel, "order"].dropna().unique() if str(o) not in {"multiple", "mixed", "UNMAPPED"}}
        )
        gene_rows.append(
            {
                "gene_id": gene,
                "n_branches_total": int(len(g)),
                "n_branches_tested": int(g["tested"].sum()),
                "n_branches_selected": int(sel.sum()),
                "min_corrected_p": g.loc[g["tested"], "corrected_p"].min(),
                "max_omega": g["omega_max"].max(),
                "selected_any": bool(sel.any()),
                "selected_orders": ";".join(selected_orders),
                "selected_bat": bool((sel & (g["taxon_group"] == "bat")).any()),
                "selected_non_bat": bool((sel & (g["taxon_group"] == "non-bat")).any()),
                "selected_mixed": bool((sel & (g["taxon_group"] == "mixed")).any()),
                "mapping_ok": bool(g["mapping_ok"].all()),
            }
        )
    gene = pd.DataFrame(gene_rows)
    gene["q_gene"] = bh_fdr(gene["min_corrected_p"])
    gene["selected_any_fdr"] = gene["q_gene"] < args.fdr
    gene = gene.sort_values(["selected_any", "min_corrected_p", "gene_id"], ascending=[False, True, True], kind="stable")

    # Branch-level table ---------------------------------------------------
    branch_rows = []
    for branch, g in all_df.groupby("branch_id", sort=False):
        n_tested = int(g["tested"].sum())
        n_selected = int(g["selected"].sum())
        branch_rows.append(
            {
                "branch_id": branch,
                "n_genes_total": int(len(g)),
                "n_genes_tested": n_tested,
                "n_genes_selected": n_selected,
                "selected_fraction": (n_selected / n_tested) if n_tested else np.nan,
                "min_corrected_p": g.loc[g["tested"], "corrected_p"].min(),
                "max_omega": g["omega_max"].max(),
            }
        )
    branch = pd.DataFrame(branch_rows).merge(meta_small, on="branch_id", how="left", validate="one_to_one")
    branch = branch.sort_values(["order", "branch_type", "branch_id"], kind="stable")

    # Write outputs --------------------------------------------------------
    branch_tests_path = outdir / "branch_tests.tsv.gz"
    tmp = branch_tests_path.with_name(branch_tests_path.name + ".tmp")
    tested.to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, branch_tests_path)

    gene_path = outdir / "gene_summary.tsv.gz"
    tmp = gene_path.with_name(gene_path.name + ".tmp")
    gene.to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, gene_path)

    write_tsv(branch, outdir / "branch_summary.tsv")

    statuses = read_statuses(args.status_dir)
    errors = statuses[statuses["status"].isin(["error", "skipped"])].copy()
    if not errors.empty:
        keep = [c for c in ["gene_id", "status", "reason", "runtime_sec", "hyphy_version"] if c in errors.columns]
        write_tsv(errors[keep].sort_values(["status", "gene_id"]), outdir / "errors.tsv")

    unmapped = all_df[~all_df["mapping_ok"]]
    if not unmapped.empty:
        write_tsv(
            unmapped[["gene_id", "branch_id", "node_label", "tested", "corrected_p", "selected"]],
            outdir / "unmapped_branches.tsv",
        )

    status_counts: Dict[str, int] = statuses["status"].value_counts().to_dict() if not statuses.empty else {}
    versions = sorted({str(v) for v in statuses.get("hyphy_version", pd.Series(dtype=str)).dropna() if str(v)})
    run_summary = {
        "generated_utc": now_utc(),
        "p_threshold": args.p_threshold,
        "gene_fdr": args.fdr,
        "n_gene_tables": len(files),
        "n_genes_in_status": int(len(statuses)),
        "status_counts": status_counts,
        "hyphy_versions": versions,
        "n_branch_rows": int(len(all_df)),
        "n_branch_tests": int(len(tested)),
        "n_selected_tests": int(tested["selected"].sum()),
        "n_genes_selected_any": int(gene["selected_any"].sum()),
        "n_genes_selected_any_fdr": int(gene["selected_any_fdr"].sum()),
        "n_mapping_failures": int((~all_df["mapping_ok"]).sum()),
        "genes": int(len(gene)),
        "branches": int(len(branch)),
    }
    write_json(run_summary, outdir / "run_summary.json")
    print(json.dumps(run_summary, indent=2), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
