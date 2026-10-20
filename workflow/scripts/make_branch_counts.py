#!/usr/bin/env python3
"""Count positively selected genes in one GO term per branch.

Output is the input table for the negative-binomial regression.  In addition
to the paper's branch length, it exposes the number of GO genes actually tested
on each branch (`n_genes_tested`).  That column can be used as a log offset to
separate "more opportunity because more genes were testable" from "more
selection".
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--branch-tests", required=True)
    p.add_argument("--branch-summary", default=None)
    p.add_argument("--branch-metadata", required=True)
    p.add_argument("--gene2go", required=True)
    p.add_argument("--go-id", required=True)
    p.add_argument("--out", required=True)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    go_id = args.go_id.replace("_", ":", 1) if args.go_id.startswith("GO_") else args.go_id

    tests = pd.read_csv(args.branch_tests, sep="\t", dtype={"gene_id": "string", "branch_id": "string"})
    meta = pd.read_csv(args.branch_metadata, sep="\t", dtype={"branch_id": "string"})
    g2g = pd.read_csv(args.gene2go, sep="\t", dtype=str)

    for col in ["gene_id", "go_id"]:
        if col not in g2g.columns:
            raise SystemExit(f"gene2go is missing column {col!r}")
    g2g["go_id"] = g2g["go_id"].str.replace("_", ":", n=1)
    go_genes = set(g2g.loc[g2g["go_id"] == go_id, "gene_id"].dropna().astype(str))
    if not go_genes:
        print(f"[warn] no genes found for {go_id} in {args.gene2go}", file=sys.stderr)

    tests["tested"] = tests["tested"].astype(str).str.lower().isin(["true", "t", "1", "yes"])
    tests["selected"] = tests["selected"].astype(str).str.lower().isin(["true", "t", "1", "yes"])
    sub = tests[tests["gene_id"].isin(go_genes) & tests["tested"]].copy()

    n_tested = sub.groupby("branch_id")["gene_id"].nunique()
    n_selected = sub[sub["selected"]].groupby("branch_id")["gene_id"].nunique()

    out = meta.copy()
    out["n_genes_tested"] = out["branch_id"].map(n_tested).fillna(0).astype(int)
    out["count"] = out["branch_id"].map(n_selected).fillna(0).astype(int)
    out["fraction"] = np.where(out["n_genes_tested"] > 0, out["count"] / out["n_genes_tested"], np.nan)

    if args.branch_summary:
        bs = pd.read_csv(args.branch_summary, sep="\t", dtype={"branch_id": "string"})
        if {"branch_id", "n_genes_tested"}.issubset(bs.columns):
            out = out.merge(
                bs[["branch_id", "n_genes_tested"]].rename(columns={"n_genes_tested": "n_genes_tested_all"}),
                on="branch_id",
                how="left",
                validate="one_to_one",
            )
    if "n_genes_tested_all" not in out.columns:
        all_tested = tests[tests["tested"]].groupby("branch_id")["gene_id"].nunique()
        out["n_genes_tested_all"] = out["branch_id"].map(all_tested).fillna(0).astype(int)

    out["label"] = out["original_label"].fillna("").astype(str)
    out.loc[out["label"].eq(""), "label"] = out.loc[out["label"].eq(""), "branch_id"]
    out["taxon"] = out["taxon_group"]
    keep = [
        "branch_id", "label", "order", "taxon", "taxon_group", "branch_type",
        "branch_length", "n_descendant_taxa", "is_order_stem",
        "n_genes_tested", "n_genes_tested_all", "count", "fraction",
    ]
    keep = [c for c in keep if c in out.columns]
    out = out[keep].sort_values(["order", "branch_type", "branch_id"], kind="stable")
    write_tsv(out, args.out)
    print(
        f"[ok] {go_id}: {len(go_genes)} GO genes, {int(out['count'].sum())} selected branch-gene pairs, "
        f"{len(out)} branches -> {args.out}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
