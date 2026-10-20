# Toy fixtures (synthetic)

All files here are **synthetic**. They exist to test file formats, tree pruning, branch mapping, JSON parsing and the R scripts. They are not sampled from the paper and must not be used for biological inference.

## Files

```
examples/toy/
├── make_toy_data.py                 # deterministic generator
├── alignments/
│   ├── toy_gene1.fa                 # 6 taxa, 30 codons
│   ├── toy_gene2.fa                 # 5 taxa (art1 absent), 20 codons
│   └── toy_gene3.fa                 # 6 taxa, 20 codons
└── resources/
    ├── species_tree.nwk             # 6 taxa, substitutions/site
    ├── species_metadata.tsv
    ├── gene_map.tsv                 # toy gene -> real human Ensembl ID
    ├── gene2go.tsv
    ├── branch_counts_toy.tsv        # hand-made regression input (9 branches)
    └── zover_records.tsv            # 24 records, 8 studies
```

The tree:

```
(((bat1:0.05,bat2:0.08):0.12,(pri1:0.06,(rod1:0.03,rod2:0.04):0.09):0.15):0.10,art1:0.20);
```

The root is a trifurcation, so the labeled tree has 9 branches. `branch_counts_toy.tsv` uses the labels that `make_branch_metadata.py` assigns in preorder (`n00001`, `n00002`, `n00003`); if you change the tree, regenerate/inspect the labels before using this table.

## Run the pipeline on toy data

```bash
# optional: regenerate
python examples/toy/make_toy_data.py

snakemake -s workflow/Snakefile \
  --configfile workflow/config/config.toy.yaml \
  --cores 2 --use-conda all

snakemake -s workflow/Snakefile \
  --configfile workflow/config/config.toy.yaml \
  --cores 2 --use-conda downstream
```

The toy config uses short chains and `min_branches: 8`; the regression will run on 9 branches, which is enough for code-path testing only. `min_foreground: 1` lets the enrichment script run with two GO:0002376 genes; g:Profiler may still return no terms.

## Run the R scripts directly

```bash
Rscript R/count2branches_optimized.R \
  --input examples/toy/resources/branch_counts_toy.tsv \
  --outdir results/toy_regression \
  --quick --min-branches 8

Rscript R/fisher_zover_optimized.R \
  --input examples/toy/resources/zover_records.tsv \
  --outdir results/toy_zover \
  --n-boot 200
```

The branch counts table is hand-made and does **not** come from the toy aBSREL run; it only exercises the regression code. The ZOVER table is also hand-made; its Fisher/GEE numbers are not meaningful.

## Why `toy_gene2` has only 5 taxa

It deliberately omits `art1` so the pipeline must:

1. intersect alignment taxa with tree taxa;
2. prune `art1` from the tree;
3. suppress unifurcations;
4. map the remaining pruned branches back to the full tree by descendant leaf sets.

If branch mapping is wrong, `toy_gene2` is the fixture that should expose it.
