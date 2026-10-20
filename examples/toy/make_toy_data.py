#!/usr/bin/env python3
"""Generate the tiny synthetic fixtures under examples/toy/.

These files exist only to exercise file formats and plumbing.  They are not
sampled from the paper and must never be used for biological inference.
Run from the repository root:

    python examples/toy/make_toy_data.py
"""
from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
ALN = HERE / "alignments"
RES = HERE / "resources"
ALN.mkdir(parents=True, exist_ok=True)
RES.mkdir(parents=True, exist_ok=True)

BASE = [
    "ATG", "GCT", "GAT", "AAA", "TTC", "CAA", "GAA", "CTG", "ATT", "GTT",
    "TCT", "ACT", "CCT", "CGT", "GGT", "CAT", "TAT", "TGT", "TGG", "GCA",
    "GCC", "GCG", "AAG", "GAC", "TTT", "CAG", "GAG", "CTC", "ATC", "GTC",
]
STOPS = {"TAA", "TAG", "TGA"}


def seq(n_codons: int, changes: dict[int, str]) -> str:
    codons = BASE[:n_codons].copy()
    for pos, codon in changes.items():
        codons[pos - 1] = codon
    out = "".join(codons)
    assert len(out) == n_codons * 3
    assert not any(out[i:i + 3] in STOPS for i in range(0, len(out), 3))
    return out


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("wt") as fh:
        for name, s in records:
            fh.write(f">{name}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")


# Gene 1: six taxa, 30 codons, mostly synonymous variation.
write_fasta(ALN / "toy_gene1.fa", [
    ("bat1", seq(30, {})),
    ("bat2", seq(30, {10: "GTA", 14: "CGC", 20: "GCC", 25: "TTC"})),
    ("pri1", seq(30, {3: "GAC", 8: "CTA", 12: "ACC", 18: "TGC"})),
    ("rod1", seq(30, {4: "AAG", 9: "ATC", 15: "GGC", 21: "GCA"})),
    ("rod2", seq(30, {5: "TTT", 11: "TCC", 16: "CAC", 22: "GCA"})),
    ("art1", seq(30, {6: "CAG", 13: "CCC", 17: "TAC", 23: "AAA"})),
])

# Gene 2: art1 is absent; the pipeline must prune it from the tree.
write_fasta(ALN / "toy_gene2.fa", [
    ("bat1", seq(20, {})),
    ("bat2", seq(20, {7: "GAG", 14: "CGC"})),
    ("pri1", seq(20, {2: "GCC", 9: "ATC", 15: "GGC"})),
    ("rod1", seq(20, {4: "AAG", 11: "TCC", 18: "TGC"})),
    ("rod2", seq(20, {5: "TTT", 12: "ACC", 20: "GCG"})),
])

# Gene 3: six taxa, 20 codons.
write_fasta(ALN / "toy_gene3.fa", [
    ("bat1", seq(20, {})),
    ("bat2", seq(20, {3: "GAC", 10: "GTC"})),
    ("pri1", seq(20, {6: "CAG", 13: "CCA"})),
    ("rod1", seq(20, {8: "CTC", 16: "CAC"})),
    ("rod2", seq(20, {9: "ATA", 17: "TAC"})),
    ("art1", seq(20, {1: "ATA", 14: "CGA"})),
])

# Branch lengths are substitutions/site.  The root is a trifurcation.
(RES / "species_tree.nwk").write_text(
    "(((bat1:0.05,bat2:0.08):0.12,(pri1:0.06,(rod1:0.03,rod2:0.04):0.09):0.15):0.10,"
    "art1:0.20);\n"
)

(RES / "species_metadata.tsv").write_text(
    "taxon\torder\tis_bat\n"
    "bat1\tChiroptera\tTRUE\n"
    "bat2\tChiroptera\tTRUE\n"
    "pri1\tPrimates\tFALSE\n"
    "rod1\tRodentia\tFALSE\n"
    "rod2\tRodentia\tFALSE\n"
    "art1\tArtiodactyla\tFALSE\n"
)

# Real human Ensembl IDs are used so that g:Profiler can resolve them.
(RES / "gene_map.tsv").write_text(
    "gene_id\tensembl_gene_id\n"
    "toy_gene1\tENSG00000141510\n"
    "toy_gene2\tENSG00000012048\n"
    "toy_gene3\tENSG00000157764\n"
)

(RES / "gene2go.tsv").write_text(
    "gene_id\tgo_id\n"
    "toy_gene1\tGO:0002376\n"
    "toy_gene2\tGO:0002376\n"
    "toy_gene3\tGO:0006955\n"
)

# A hand-made regression input matching the labeled toy tree.  The first row
# is the chiropteran stem; the third internal row is the rodent stem.
(RES / "branch_counts_toy.tsv").write_text(
    "branch_id\tlabel\torder\ttaxon\ttaxon_group\tbranch_type\tbranch_length\t"
    "n_descendant_taxa\tis_order_stem\tn_genes_tested\tn_genes_tested_all\tcount\tfraction\n"
    "n00001\tChiroptera_stem\tChiroptera\tbat\tbat\tinternal\t0.12\t2\tTRUE\t2\t3\t2\t1.0\n"
    "bat1\tbat1\tChiroptera\tbat\tbat\tterminal\t0.05\t1\tFALSE\t2\t3\t1\t0.5\n"
    "bat2\tbat2\tChiroptera\tbat\tbat\tterminal\t0.08\t1\tFALSE\t2\t3\t1\t0.5\n"
    "n00002\tEuarchontoglires_stem\tmultiple\tmixed\tmixed\tinternal\t0.15\t3\tFALSE\t2\t3\t0\t0.0\n"
    "pri1\tpri1\tPrimates\tnon-bat\tnon-bat\tterminal\t0.06\t1\tFALSE\t2\t3\t0\t0.0\n"
    "n00003\tRodentia_stem\tRodentia\tnon-bat\tnon-bat\tinternal\t0.09\t2\tTRUE\t2\t3\t1\t0.5\n"
    "rod1\trod1\tRodentia\tnon-bat\tnon-bat\tterminal\t0.03\t1\tFALSE\t2\t3\t0\t0.0\n"
    "rod2\trod2\tRodentia\tnon-bat\tnon-bat\tterminal\t0.04\t1\tFALSE\t2\t3\t1\t0.5\n"
    "art1\tart1\tArtiodactyla\tnon-bat\tnon-bat\tterminal\t0.20\t1\tFALSE\t2\t3\t0\t0.0\n"
)

# Record-level toy ZOVER table.  Study ids are the clusters.
zover_rows = [
    ("R01", "S1", "bat", "Rhinolophidae", 1),
    ("R02", "S1", "bat", "Rhinolophidae", 1),
    ("R03", "S1", "bat", "Rhinolophidae", 0),
    ("R04", "S2", "bat", "Rhinolophidae", 1),
    ("R05", "S2", "bat", "Rhinolophidae", 0),
    ("R06", "S2", "bat", "Rhinolophidae", 0),
    ("R07", "S3", "bat", "Hipposideridae", 1),
    ("R08", "S3", "bat", "Hipposideridae", 0),
    ("R09", "S3", "bat", "Hipposideridae", 0),
    ("R10", "S4", "bat", "Pteropodidae", 0),
    ("R11", "S4", "bat", "Pteropodidae", 0),
    ("R12", "S4", "bat", "Pteropodidae", 0),
    ("R13", "S5", "rodent", "Muridae", 0),
    ("R14", "S5", "rodent", "Muridae", 0),
    ("R15", "S5", "rodent", "Muridae", 0),
    ("R16", "S6", "rodent", "Muridae", 0),
    ("R17", "S6", "rodent", "Muridae", 0),
    ("R18", "S6", "rodent", "Muridae", 0),
    ("R19", "S7", "rodent", "Cricetidae", 1),
    ("R20", "S7", "rodent", "Cricetidae", 0),
    ("R21", "S7", "rodent", "Cricetidae", 0),
    ("R22", "S8", "rodent", "Cricetidae", 0),
    ("R23", "S8", "rodent", "Cricetidae", 0),
    ("R24", "S8", "rodent", "Cricetidae", 0),
]
with (RES / "zover_records.tsv").open("wt") as fh:
    fh.write("record_id\tstudy_id\thost_group\thost_family\tis_coronavirus\n")
    for row in zover_rows:
        fh.write("\t".join(map(str, row)) + "\n")

print("wrote toy fixtures under", HERE)
