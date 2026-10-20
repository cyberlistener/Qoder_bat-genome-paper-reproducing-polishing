#!/usr/bin/env Rscript
# Background-aware functional enrichment for per-order selected genes.
#
# Why this replaces the original one-liner:
#   The paper's snippet used gost(..., significant = TRUE) without a custom
#   background.  g:Profiler then treats all annotated human genes as the
#   statistical domain, even though only the 17,130 screened genes could ever
#   be selected.  That inflates enrichment for terms that are simply common in
#   the screened set and makes P values from different orders hard to compare.
#
# This script runs each order twice:
#   1. custom background = genes with a tested branch in that order;
#   2. default background = g:Profiler's annotated human domain (closest to
#      the published analysis, kept for comparison).
#
# Usage:
#   Rscript R/enrichment_gprofiler_optimized.R \
#     --branch-tests results/aggregate/branch_tests.tsv.gz \
#     --gene-summary results/aggregate/gene_summary.tsv.gz \
#     --branch-metadata results/branch_metadata.tsv \
#     --gene-map resources/gene_map.tsv \
#     --outdir results/enrichment

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(gprofiler2)
})

option_list <- list(
  make_option("--branch-tests", type = "character", help = "aggregated branch_tests.tsv.gz"),
  make_option("--gene-summary", type = "character", default = NULL,
              help = "gene_summary.tsv.gz (used only for a sanity count)"),
  make_option("--branch-metadata", type = "character", default = NULL,
              help = "branch metadata (order/taxon_group are already in branch_tests)"),
  make_option("--gene-map", type = "character", help = "TSV: gene_id, ensembl_gene_id"),
  make_option("--outdir", type = "character", default = "results/enrichment"),
  make_option("--organism", type = "character", default = "hsapiens"),
  make_option("--sources", type = "character", default = "GO:BP,GO:CC,GO:MF,KEGG,REAC,WP"),
  make_option("--fdr", type = "double", default = 0.05),
  make_option("--min-foreground", type = "integer", default = 10)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$`branch-tests`)) stop("--branch-tests is required", call. = FALSE)
if (is.null(opt$`gene-map`)) stop("--gene-map is required", call. = FALSE)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

tests <- read_tsv(opt$`branch-tests`, show_col_types = FALSE)
needed <- c("gene_id", "branch_id", "selected", "tested", "order", "taxon_group")
miss <- setdiff(needed, names(tests))
if (length(miss)) stop("branch tests missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
tests <- tests %>%
  mutate(
    selected = tolower(as.character(selected)) %in% c("true", "t", "1", "yes"),
    tested = tolower(as.character(tested)) %in% c("true", "t", "1", "yes"),
    order = as.character(order),
    taxon_group = as.character(taxon_group)
  ) %>%
  filter(tested, !is.na(gene_id))

gene_map <- read_tsv(opt$`gene-map`, show_col_types = FALSE)
if (!all(c("gene_id", "ensembl_gene_id") %in% names(gene_map))) {
  stop("--gene-map must have columns gene_id and ensembl_gene_id", call. = FALSE)
}
gene_map <- gene_map %>%
  mutate(
    gene_id = as.character(gene_id),
    ensembl_gene_id = sub("\\..*$", "", as.character(ensembl_gene_id))
  ) %>%
  distinct(gene_id, .keep_all = TRUE)

tests <- tests %>% left_join(gene_map, by = "gene_id")
if (all(is.na(tests$ensembl_gene_id))) {
  warning("No gene IDs mapped to Ensembl; falling back to gene_id. g:Profiler may not recognise transcript IDs.")
  tests$ensembl_gene_id <- tests$gene_id
}
tests <- tests %>% filter(!is.na(ensembl_gene_id), ensembl_gene_id != "")

orders <- sort(setdiff(unique(tests$order), c("multiple", "mixed", "UNMAPPED", "", NA)))
if (!length(orders)) stop("no order labels found in branch_tests", call. = FALSE)
sources <- trimws(strsplit(opt$sources, ",", fixed = TRUE)[[1]])

flatten_list_cols <- function(df) {
  if (is.null(df) || !nrow(df)) return(df)
  for (nm in names(df)) {
    if (is.list(df[[nm]])) {
      df[[nm]] <- vapply(df[[nm]], function(x) paste(x, collapse = ";"), character(1))
    }
  }
  df
}

run_gost <- function(query, background, use_custom) {
  args <- list(
    query = query,
    organism = opt$organism,
    sources = sources,
    correction_method = "fdr",
    significant = FALSE,
    user_threshold = 1,
    evcodes = TRUE
  )
  if (use_custom) {
    args$domain_scope <- "custom"
    args$custom_bg <- background
  }
  do.call(gost, args)
}

annotate_result <- function(res, order_label, use_custom) {
  if (is.null(res) || is.null(res$result) || !nrow(res$result)) return(NULL)
  df <- flatten_list_cols(as.data.frame(res$result, stringsAsFactors = FALSE))
  df$order <- order_label
  df$background <- if (use_custom) "custom_tested_in_order" else "gprofiler_annotated"
  df$p_adjust_method <- "gprofiler_fdr"
  df$significant <- !is.na(df$p_value) & df$p_value < opt$fdr
  df$neg_log10_p <- ifelse(!is.na(df$p_value) & df$p_value > 0, -log10(df$p_value), NA_real_)
  if (all(c("intersection_size", "query_size", "term_size", "effective_domain_size") %in% names(df))) {
    df$fold_enrichment <- (df$intersection_size / df$query_size) /
      (df$term_size / df$effective_domain_size)
  }
  df
}

custom_rows <- list()
default_rows <- list()
size_rows <- list()

for (o in orders) {
  d <- tests %>% filter(order == o)
  fg <- unique(d$ensembl_gene_id[d$selected])
  bg <- unique(d$ensembl_gene_id)
  status <- "ok"
  if (length(fg) < opt$min_foreground) {
    status <- "skipped_too_few_foreground"
    size_rows[[o]] <- data.frame(
      order = o, n_foreground = length(fg), n_background = length(bg),
      n_terms_custom = 0L, n_terms_default = 0L, status = status,
      stringsAsFactors = FALSE
    )
    next
  }

  res_custom <- tryCatch(run_gost(fg, bg, TRUE), error = function(e) {
    warning("custom gost failed for ", o, ": ", conditionMessage(e))
    NULL
  })
  res_default <- tryCatch(run_gost(fg, bg, FALSE), error = function(e) {
    warning("default gost failed for ", o, ": ", conditionMessage(e))
    NULL
  })

  a <- annotate_result(res_custom, o, TRUE)
  b <- annotate_result(res_default, o, FALSE)
  if (!is.null(a)) custom_rows[[o]] <- a else status <- "custom_no_terms"
  if (!is.null(b)) default_rows[[o]] <- b
  size_rows[[o]] <- data.frame(
    order = o, n_foreground = length(fg), n_background = length(bg),
    n_terms_custom = if (is.null(a)) 0L else nrow(a),
    n_terms_default = if (is.null(b)) 0L else nrow(b),
    status = status,
    stringsAsFactors = FALSE
  )
}

bind_or_empty <- function(rows, schema) {
  if (!length(rows)) {
    out <- as.data.frame(setNames(replicate(length(schema), character(0), simplify = FALSE), schema),
                         stringsAsFactors = FALSE)
    return(out)
  }
  bind_rows(rows)
}

enrichment_schema <- c(
  "query", "significant", "p_value", "term_size", "query_size",
  "intersection_size", "effective_domain_size", "source", "term_id",
  "term_name", "order", "background", "p_adjust_method",
  "neg_log10_p", "fold_enrichment"
)

custom <- bind_or_empty(custom_rows, enrichment_schema)
default <- bind_or_empty(default_rows, enrichment_schema)
sizes <- bind_or_empty(size_rows, c("order", "n_foreground", "n_background",
                                    "n_terms_custom", "n_terms_default", "status"))

write_tsv(custom, file.path(opt$outdir, "gprofiler_custom.tsv"))
write_tsv(default, file.path(opt$outdir, "gprofiler_default.tsv"))
write_tsv(sizes, file.path(opt$outdir, "foreground_sizes.tsv"))
capture.output(sessionInfo(), file = file.path(opt$outdir, "sessionInfo.txt"))

message("custom terms: ", nrow(custom), "; default terms: ", nrow(default),
        "; orders: ", length(orders))
