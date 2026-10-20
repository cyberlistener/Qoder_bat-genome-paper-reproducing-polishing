#!/usr/bin/env python3
"""Label every internal node of the species tree and emit branch metadata.

Why this exists
---------------
aBSREL returns results keyed by tree node labels.  If internal nodes are
unlabelled (or if two nodes share a label), per-gene results cannot be joined
back to branches in a stable way.  This script assigns deterministic, unique
labels, keeps the original labels for display, and annotates each branch with
the set of descendant orders / bat status.

The input tree must carry branch lengths in substitutions per site.  A tree
whose branch lengths are coalescent units (e.g. raw ASTRAL output) or millions
of years must not be passed to aBSREL without conversion.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Set

import dendropy
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import as_bool, write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree", required=True, help="Input Newick tree (branch lengths in substitutions/site)")
    p.add_argument("--species-metadata", required=True, help="TSV with taxon, order, is_bat")
    p.add_argument("--out-tree", required=True, help="Output Newick with unique node labels")
    p.add_argument("--out-metadata", required=True, help="Output branch metadata TSV")
    p.add_argument("--allow-zero-lengths", action="store_true", help="Do not fail on zero-length branches")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    meta = pd.read_csv(args.species_metadata, sep="\t", dtype=str)
    required = {"taxon", "order", "is_bat"}
    missing = required - set(meta.columns)
    if missing:
        raise SystemExit(f"species metadata is missing columns: {sorted(missing)}")
    meta["taxon"] = meta["taxon"].astype(str)
    meta["order"] = meta["order"].astype(str)
    meta["is_bat"] = meta["is_bat"].map(as_bool)
    if meta["taxon"].duplicated().any():
        dup = meta.loc[meta["taxon"].duplicated(), "taxon"].tolist()
        raise SystemExit(f"duplicate taxa in species metadata: {dup[:10]}")
    meta_by_taxon: Dict[str, dict] = meta.set_index("taxon").to_dict("index")

    tree = dendropy.Tree.get(
        path=args.tree,
        schema="newick",
        preserve_underscores=True,
        rooting="default-rooted",
    )
    for node in tree.preorder_node_iter():
        if node.is_leaf() and (node.taxon is None or not node.taxon.label):
            raise SystemExit("tree contains an unnamed leaf")

    leaf_labels = [leaf.taxon.label for leaf in tree.leaf_node_iter()]
    missing_meta = sorted(set(leaf_labels) - set(meta_by_taxon))
    if missing_meta:
        raise SystemExit(
            "tree leaves missing from species metadata: " + ", ".join(missing_meta[:20])
            + (" ..." if len(missing_meta) > 20 else "")
        )
    extra_meta = sorted(set(meta_by_taxon) - set(leaf_labels))
    if extra_meta:
        print(f"[warn] metadata taxa absent from tree: {extra_meta[:20]}", file=sys.stderr)

    no_length = [n for n in tree.preorder_node_iter() if n.parent_node is not None and n.edge_length is None]
    if no_length:
        raise SystemExit(
            "tree has branches without lengths; aBSREL needs substitutions/site branch lengths. "
            "Use an IQ-TREE/RAxML tree, or convert the units explicitly."
        )
    zero_length = [n for n in tree.preorder_node_iter() if n.parent_node is not None and float(n.edge_length) <= 0]
    if zero_length and not args.allow_zero_lengths:
        print(f"[warn] {len(zero_length)} zero-length branches; they will be flagged but not removed", file=sys.stderr)

    used = set()
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            node.label = node.taxon.label
        else:
            node.label = None
        if node.label:
            used.add(node.label)

    idx = 0
    for node in tree.preorder_node_iter():
        if node.is_leaf():
            continue
        while True:
            idx += 1
            candidate = f"n{idx:05d}"
            if candidate not in used:
                break
        node.label = candidate
        used.add(candidate)

    rows = []
    for node in tree.preorder_node_iter():
        parent = node.parent_node
        if parent is None:
            continue
        leaves = [leaf.taxon.label for leaf in node.leaf_iter()]
        orders: Set[str] = {str(meta_by_taxon[t]["order"]) for t in leaves}
        bats: Set[bool] = {bool(meta_by_taxon[t]["is_bat"]) for t in leaves}
        order = next(iter(orders)) if len(orders) == 1 else "multiple"
        if bats == {True}:
            group = "bat"
        elif bats == {False}:
            group = "non-bat"
        else:
            group = "mixed"

        is_stem = False
        if len(orders) == 1 and parent is not None:
            parent_orders = {str(meta_by_taxon[leaf.taxon.label]["order"]) for leaf in parent.leaf_iter()}
            is_stem = len(parent_orders) > 1

        rows.append(
            {
                "branch_id": str(node.label),
                "parent_id": str(parent.label) if parent.label else "",
                "branch_type": "terminal" if node.is_leaf() else "internal",
                "branch_length": float(node.edge_length),
                "n_descendant_taxa": len(leaves),
                "order": order,
                "taxon_group": group,
                "is_order_stem": bool(is_stem),
                "descendant_taxa": ",".join(sorted(leaves)),
            }
        )

    # Recover original internal labels by descendant leaf set; this is stable
    # even though the tree has just been relabelled.
    original_tree = dendropy.Tree.get(
        path=args.tree, schema="newick", preserve_underscores=True, rooting="default-rooted"
    )
    # Match by descendant leaf set, which is stable under relabelling.
    original_by_taxa = {}
    for node in original_tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        key = tuple(sorted(leaf.taxon.label for leaf in node.leaf_iter()))
        original_by_taxa[key] = node.label
    for row in rows:
        key = tuple(sorted(row["descendant_taxa"].split(",")))
        row["original_label"] = original_by_taxa.get(key) or ""

    out = pd.DataFrame(rows)
    out = out.sort_values(["branch_type", "order", "branch_id"], kind="stable")
    write_tsv(out, args.out_metadata)
    Path(args.out_tree).parent.mkdir(parents=True, exist_ok=True)
    tree.write(path=args.out_tree, schema="newick", suppress_rooting=True)
    print(f"[ok] {len(out)} branches; wrote {args.out_tree} and {args.out_metadata}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
