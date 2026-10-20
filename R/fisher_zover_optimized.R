#!/usr/bin/env Rscript
# ZOVER-style coronavirus detection analysis, with a clustering sensitivity.
#
# The published analysis uses a two-tailed Fisher exact test on a 2 x 2 table
# of viral records.  That treats every record as independent.  Records from the
# same metagenomic study (and often the same host species/sample) are not
# independent, so the Fisher P value is best read as a descriptive association.
# This script therefore reports:
#   1. the paper-style Fisher test;
#   2. a study-clustered GEE logistic model (geepack);
#   3. a study-level cluster bootstrap CI for the odds ratio.
#
# Input TSV (one row per viral record):
#   record_id       optional unique record id
#   study_id        study / publication / dataset id (cluster)
#   host_group      bat or rodent
#   host_family     bat family (Rhinolophidae, Hipposideridae, ...); optional
#   is_coronavirus  1/0 (or TRUE/FALSE)
#
# Usage:
#   Rscript R/fisher_zover_optimized.R \
#     --input examples/toy/resources/zover_records.tsv \
#     --outdir results/zover

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(geepack)
})

option_list <- list(
  make_option("--input", type = "character", help = "record-level TSV"),
  make_option("--outdir", type = "character", default = "results/zover"),
  make_option("--n-boot", type = "integer", default = 1000),
  make_option("--seed", type = "integer", default = 20250101)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("--input is required", call. = FALSE)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
set.seed(opt$seed)

as_flag <- function(x) {
  if (is.logical(x)) return(x)
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

dat0 <- read_tsv(opt$input, show_col_types = FALSE)
required <- c("study_id", "host_group", "is_coronavirus")
missing <- setdiff(required, names(dat0))
if (length(missing)) {
  stop("input is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}
if (!"record_id" %in% names(dat0)) dat0$record_id <- as.character(seq_len(nrow(dat0)))
if (!"host_family" %in% names(dat0)) dat0$host_family <- NA_character_

dat <- dat0 %>%
  transmute(
    record_id = as.character(record_id),
    study_id = as.character(study_id),
    host_group = tolower(trimws(as.character(host_group))),
    host_family = trimws(as.character(host_family)),
    coronavirus = as_flag(is_coronavirus)
  ) %>%
  mutate(
    host_group = case_when(
      host_group %in% c("chiroptera", "chiropteran", "bat") ~ "bat",
      host_group %in% c("rodentia", "rodent") ~ "rodent",
      TRUE ~ host_group
    )
  ) %>%
  filter(!is.na(study_id), study_id != "", !is.na(host_group), !is.na(coronavirus))

if (!nrow(dat)) stop("no usable records", call. = FALSE)

# Normalise common spelling variants.
dat <- dat %>%
  mutate(
    host_family = case_when(
      grepl("^rhinoloph", host_family, ignore.case = TRUE) ~ "Rhinolophidae",
      grepl("^hipposider", host_family, ignore.case = TRUE) ~ "Hipposideridae",
      host_family %in% c("", "NA") ~ NA_character_,
      TRUE ~ host_family
    )
  )

or_haldane <- function(tab) {
  # tab: rows = unexposed/exposed, columns = non-case/case
  a <- tab[2, 2]; b <- tab[2, 1]; c <- tab[1, 2]; d <- tab[1, 1]
  ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
}

make_tab <- function(d, exposure) {
  tab <- table(factor(exposure, levels = c(FALSE, TRUE)),
               factor(d$coronavirus, levels = c(FALSE, TRUE)))
  tab
}

fisher_row <- function(d, exposure, comparison) {
  tab <- make_tab(d, exposure)
  ft <- tryCatch(fisher.test(tab, alternative = "two.sided"), error = function(e) NULL)
  data.frame(
    comparison = comparison,
    n_records = nrow(d),
    n_studies = n_distinct(d$study_id),
    n_cases = sum(d$coronavirus),
    n_exposed = sum(exposure),
    odds_ratio = if (is.null(ft)) NA_real_ else unname(ft$estimate),
    ci_low = if (is.null(ft)) NA_real_ else unname(ft$conf.int[1]),
    ci_high = if (is.null(ft)) NA_real_ else unname(ft$conf.int[2]),
    p_value = if (is.null(ft)) NA_real_ else ft$p.value,
    method = "two-sided Fisher exact; records treated as independent",
    stringsAsFactors = FALSE
  )
}

gee_row <- function(d, exposure, comparison) {
  d <- d %>% mutate(.exposed = as.integer(exposure), .study = factor(study_id))
  if (nlevels(d$.study) < 2 || length(unique(d$.exposed)) < 2) {
    return(data.frame(
      comparison = comparison, n_records = nrow(d), n_studies = n_distinct(d$study_id),
      beta = NA_real_, robust_se = NA_real_, odds_ratio = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
      method = "GEE not estimable (one study or no exposure variation)",
      stringsAsFactors = FALSE
    ))
  }
  fit <- tryCatch(
    geeglm(coronavirus ~ .exposed, id = .study, data = d,
           family = binomial(link = "logit"), corstr = "exchangeable"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(
      comparison = comparison, n_records = nrow(d), n_studies = n_distinct(d$study_id),
      beta = NA_real_, robust_se = NA_real_, odds_ratio = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_,
      method = "GEE failed to converge",
      stringsAsFactors = FALSE
    ))
  }
  co <- coef(summary(fit))
  beta <- co[".exposed", "Estimate"]
  se <- co[".exposed", "Std.err"]
  p <- co[".exposed", "Pr(>|W|)"]
  data.frame(
    comparison = comparison,
    n_records = nrow(d),
    n_studies = n_distinct(d$study_id),
    beta = beta,
    robust_se = se,
    odds_ratio = exp(beta),
    ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se),
    p_value = p,
    method = "GEE logistic, exchangeable working correlation, cluster = study_id",
    stringsAsFactors = FALSE
  )
}

boot_row <- function(d, exposure, comparison, n_boot) {
  studies <- unique(d$study_id)
  if (length(studies) < 3) {
    return(data.frame(comparison = comparison, n_boot = 0L,
                      boot_or_low = NA_real_, boot_or_high = NA_real_,
                      method = "cluster bootstrap skipped (<3 studies)",
                      stringsAsFactors = FALSE))
  }
  vals <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    pick <- sample(studies, length(studies), replace = TRUE)
    idx <- unlist(lapply(pick, function(s) which(d$study_id == s)), use.names = FALSE)
    b <- d[idx, , drop = FALSE]
    e <- exposure[idx]
    vals[i] <- or_haldane(make_tab(b, e))
  }
  vals <- vals[is.finite(vals)]
  data.frame(
    comparison = comparison,
    n_boot = length(vals),
    boot_or_low = if (length(vals)) unname(quantile(vals, 0.025)) else NA_real_,
    boot_or_high = if (length(vals)) unname(quantile(vals, 0.975)) else NA_real_,
    method = "study-level cluster bootstrap percentile CI (Haldane-corrected OR)",
    stringsAsFactors = FALSE
  )
}

comparisons <- list(
  bats_vs_rodents = list(
    data = dat %>% filter(host_group %in% c("bat", "rodent")),
    exposure = function(d) d$host_group == "bat"
  ),
  rhinolophidae_vs_other_bats = list(
    data = dat %>% filter(host_group == "bat", !is.na(host_family)),
    exposure = function(d) d$host_family == "Rhinolophidae"
  ),
  hipposideridae_vs_other_bats = list(
    data = dat %>% filter(host_group == "bat", !is.na(host_family)),
    exposure = function(d) d$host_family == "Hipposideridae"
  )
)

fisher_rows <- list()
gee_rows <- list()
boot_rows <- list()
study_rows <- list()

for (nm in names(comparisons)) {
  d <- comparisons[[nm]]$data
  if (!nrow(d)) next
  exposure <- comparisons[[nm]]$exposure(d)
  fisher_rows[[nm]] <- fisher_row(d, exposure, nm)
  gee_rows[[nm]] <- gee_row(d, exposure, nm)
  boot_rows[[nm]] <- boot_row(d, exposure, nm, opt$n_boot)
  study_rows[[nm]] <- d %>%
    group_by(study_id) %>%
    summarise(
      n_records = n(),
      n_cases = sum(coronavirus),
      n_bat = sum(host_group == "bat"),
      n_rhinolophidae = sum(host_family == "Rhinolophidae", na.rm = TRUE),
      n_hipposideridae = sum(host_family == "Hipposideridae", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(comparison = nm, .before = 1)
}

empty_row <- function(cols) as.data.frame(setNames(replicate(length(cols), logical(0), simplify = FALSE), cols))
fisher <- if (length(fisher_rows)) bind_rows(fisher_rows) else empty_row(c(
  "comparison", "n_records", "n_studies", "n_cases", "n_exposed", "odds_ratio",
  "ci_low", "ci_high", "p_value", "method"
))
if (nrow(fisher)) {
  fisher$p_bh <- p.adjust(fisher$p_value, method = "BH")
} else {
  fisher$p_bh <- numeric(0)
}

empty_df <- function(cols) {
  as.data.frame(setNames(replicate(length(cols), character(0), simplify = FALSE), cols),
                stringsAsFactors = FALSE)
}
gee <- if (length(gee_rows)) bind_rows(gee_rows) else empty_df(c(
  "comparison", "n_records", "n_studies", "beta", "robust_se", "odds_ratio",
  "ci_low", "ci_high", "p_value", "method"
))
boot <- if (length(boot_rows)) bind_rows(boot_rows) else empty_df(c(
  "comparison", "n_boot", "boot_or_low", "boot_or_high", "method"
))
study <- if (length(study_rows)) bind_rows(study_rows) else empty_df(c(
  "comparison", "study_id", "n_records", "n_cases", "n_bat",
  "n_rhinolophidae", "n_hipposideridae"
))

write_tsv(fisher, file.path(opt$outdir, "zover_fisher.tsv"))
write_tsv(gee, file.path(opt$outdir, "zover_gee.tsv"))
write_tsv(boot, file.path(opt$outdir, "zover_bootstrap.tsv"))
write_tsv(study, file.path(opt$outdir, "zover_study_summary.tsv"))
capture.output(sessionInfo(), file = file.path(opt$outdir, "sessionInfo.txt"))

message("wrote ZOVER results to ", opt$outdir)
message("Fisher P values are descriptive; read the GEE/bootstrap columns for study-clustered uncertainty.")
