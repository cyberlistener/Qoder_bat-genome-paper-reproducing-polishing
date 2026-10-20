#!/usr/bin/env Rscript
# Bayesian negative-binomial regression of per-branch selected-gene counts.
#
# This is a cleaned-up, corrected version of the authors' count2branches_v2.R.
# The original script labelled `count ~ x:taxon` as "two intercepts, two
# slopes", but that formula has ONE intercept and TWO slopes.  The full
# two-intercept/two-slope model is `count ~ x * taxon`.  Both are fitted here
# and labelled correctly.
#
# It also adds: fixed seed, explicit priors, Rhat/ESS/divergence checks,
# WAIC + LOO, posterior predictive check, an optional log(opportunity) offset,
# and sensitivity models that separate terminal from ancestral order-stem
# branches.
#
# Usage:
#   Rscript R/count2branches_optimized.R \
#     --input results/aggregate/branch_counts_GO_0002376.tsv \
#     --outdir results/regression

suppressPackageStartupMessages({
  library(brms)
  library(loo)
  library(optparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(posterior)
})

option_list <- list(
  make_option("--input", type = "character", help = "branch counts TSV"),
  make_option("--outdir", type = "character", default = "results/regression"),
  make_option("--branch-length-col", type = "character", default = "branch_length"),
  make_option("--seed", type = "integer", default = 20250101),
  make_option("--chains", type = "integer", default = 4),
  make_option("--iter", type = "integer", default = 4000),
  make_option("--warmup", type = "integer", default = 2000),
  make_option("--min-branches", type = "integer", default = 20,
              help = "refuse to fit with fewer usable branches (toy data may lower this)"),
  make_option("--quick", action = "store_true", default = FALSE,
              help = "short chains for plumbing tests only; do not use for inference")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("--input is required", call. = FALSE)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (opt$quick) {
  opt$chains <- 2L
  opt$iter <- 1000L
  opt$warmup <- 500L
}
if (opt$warmup >= opt$iter) stop("--warmup must be smaller than --iter", call. = FALSE)
options(mc.cores = opt$chains)
set.seed(opt$seed)

as_flag <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

dat0 <- read_tsv(opt$input, show_col_types = FALSE)
required <- c("branch_id", "count", "taxon", "branch_type", "is_order_stem",
              "n_genes_tested", opt$branch_length_col)
missing <- setdiff(required, names(dat0))
if (length(missing)) {
  stop("input is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

dat <- dat0 %>%
  mutate(
    bl = as.numeric(.data[[opt$branch_length_col]]),
    count = as.integer(count),
    n_genes_tested = as.integer(n_genes_tested),
    taxon_raw = as.character(taxon),
    taxon_model = factor(if_else(taxon_raw == "bat", "bat", "non-bat"),
                         levels = c("non-bat", "bat")),
    branch_type = factor(as.character(branch_type), levels = c("terminal", "internal")),
    is_order_stem = as_flag(is_order_stem),
    taxon_group = if ("taxon_group" %in% names(.)) as.character(taxon_group) else NA_character_,
    label = if ("label" %in% names(.)) as.character(label) else as.character(branch_id),
    order = if ("order" %in% names(.)) as.character(order) else NA_character_
  ) %>%
  filter(!is.na(bl), bl > 0, !is.na(count), count >= 0, !is.na(taxon_model)) %>%
  droplevels()

if (nrow(dat) < opt$`min-branches`) {
  stop("fewer than ", opt$`min-branches`, " usable branches; refusing to fit", call. = FALSE)
}
if (nlevels(droplevels(dat$taxon_model)) < 2) {
  stop("input has only one taxon level; the bat-vs-non-bat contrast is not estimable", call. = FALSE)
}
write_tsv(dat, file.path(opt$outdir, "data_used.tsv"))

priors <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(gamma(2, 0.1), class = "shape")
)

model_specs <- list(
  m0 = list(formula = bf(count ~ sqrt(bl)),
            text = "count ~ sqrt(bl)  [one intercept, one slope]"),
  m1 = list(formula = bf(count ~ sqrt(bl) + taxon_model),
            text = "count ~ sqrt(bl) + taxon  [two intercepts, one slope]"),
  m2 = list(formula = bf(count ~ sqrt(bl) * taxon_model),
            text = "count ~ sqrt(bl) * taxon  [two intercepts, two slopes]"),
  m3 = list(formula = bf(count ~ sqrt(bl):taxon_model),
            text = "count ~ sqrt(bl):taxon  [one intercept, two slopes]")
)

fit_one <- function(spec, data, name) {
  message("fitting ", name, ": ", spec$text)
  fit <- brm(
    formula = spec$formula,
    data = data,
    family = negbinomial(link = "log"),
    prior = priors,
    chains = opt$chains,
    cores = opt$chains,
    iter = opt$iter,
    warmup = opt$warmup,
    seed = opt$seed,
    control = list(adapt_delta = 0.95, max_treedepth = 12),
    save_pars = save_pars(all = TRUE),
    refresh = 0,
    silent = 2
  )
  fit <- add_criterion(fit, c("waic", "loo"), force_save = TRUE)
  saveRDS(fit, file.path(opt$outdir, paste0("fit_", name, ".rds")))
  fit
}

fits_main <- lapply(names(model_specs), function(nm) fit_one(model_specs[[nm]], dat, nm))
names(fits_main) <- names(model_specs)

# Sensitivity models ---------------------------------------------------------
dat_offset <- dat %>% filter(n_genes_tested > 0)
dat_no_stem <- dat %>% filter(!is_order_stem)
dat_terminal <- dat %>% filter(branch_type == "terminal") %>% droplevels()
dat_excl_mixed <- if (any(dat$taxon_group == "mixed", na.rm = TRUE)) {
  dat %>% filter(is.na(taxon_group) | taxon_group != "mixed") %>% droplevels()
} else {
  dat[0, , drop = FALSE]
}

valid_sens <- function(d, needs_branch_type = FALSE) {
  if (nrow(d) < opt$`min-branches`) return(FALSE)
  if (nlevels(droplevels(d$taxon_model)) < 2) return(FALSE)
  if (needs_branch_type && nlevels(droplevels(d$branch_type)) < 2) return(FALSE)
  TRUE
}

sens_specs <- list()
if (valid_sens(dat_offset)) {
  sens_specs$m1_offset <- list(
    formula = bf(count ~ sqrt(bl) + taxon_model + offset(log(n_genes_tested))),
    text = "count ~ sqrt(bl) + taxon + offset(log(n_genes_tested))",
    data = dat_offset
  )
}
if (valid_sens(dat, needs_branch_type = TRUE)) {
  sens_specs$m1_branch_type <- list(
    formula = bf(count ~ sqrt(bl) + taxon_model + branch_type),
    text = "count ~ sqrt(bl) + taxon + branch_type",
    data = dat
  )
}
if (valid_sens(dat_no_stem)) {
  sens_specs$m1_no_order_stem <- list(
    formula = bf(count ~ sqrt(bl) + taxon_model),
    text = "count ~ sqrt(bl) + taxon; order-stem branches removed",
    data = dat_no_stem
  )
}
if (valid_sens(dat_terminal)) {
  sens_specs$m1_terminal_only <- list(
    formula = bf(count ~ sqrt(bl) + taxon_model),
    text = "count ~ sqrt(bl) + taxon; terminal branches only",
    data = dat_terminal
  )
}
if (valid_sens(dat_excl_mixed)) {
  sens_specs$m1_exclude_mixed <- list(
    formula = bf(count ~ sqrt(bl) + taxon_model),
    text = "count ~ sqrt(bl) + taxon; mixed-descendant branches removed",
    data = dat_excl_mixed
  )
}
sens_specs$m1_log_bl <- list(
  formula = bf(count ~ log(bl) + taxon_model),
  text = "count ~ log(bl) + taxon; Box-Cox-style alternative to sqrt",
  data = dat
)

fit_sens_one <- function(spec, name) {
  tryCatch(
    fit_one(spec, spec$data, name),
    error = function(e) {
      message("skipping sensitivity model ", name, ": ", conditionMessage(e))
      NULL
    }
  )
}

fits_sens <- lapply(names(sens_specs), function(nm) fit_sens_one(sens_specs[[nm]], nm))
names(fits_sens) <- names(sens_specs)
fits_sens <- fits_sens[!vapply(fits_sens, is.null, logical(1))]


# Diagnostics and model comparison ------------------------------------------
diag_one <- function(fit) {
  draws <- as_draws_df(fit)
  rhat_vals <- rhat(draws)
  ess_vals <- ess_bulk(draws)
  np <- nuts_params(fit)
  div <- if (nrow(np)) sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE) else NA_real_
  data.frame(
    rhat_max = max(rhat_vals, na.rm = TRUE),
    ess_bulk_min = min(ess_vals, na.rm = TRUE),
    divergences = div
  )
}

criterion_row <- function(fit, name, text) {
  d <- diag_one(fit)
  waic_est <- fit$criteria$waic$estimates
  loo_est <- fit$criteria$loo$estimates
  r2 <- tryCatch(median(bayes_R2(fit)[, "Estimate"], na.rm = TRUE), error = function(e) NA_real_)
  data.frame(
    model = name,
    formula = text,
    n_obs = nrow(fit$data),
    waic = if (is.null(waic_est)) NA_real_ else unname(waic_est["waic", "Estimate"]),
    waic_se = if (is.null(waic_est)) NA_real_ else unname(waic_est["waic", "SE"]),
    looic = if (is.null(loo_est)) NA_real_ else unname(loo_est["looic", "Estimate"]),
    looic_se = if (is.null(loo_est)) NA_real_ else unname(loo_est["looic", "SE"]),
    elpd_loo = if (is.null(loo_est)) NA_real_ else unname(loo_est["elpd_loo", "Estimate"]),
    p_loo = if (is.null(loo_est)) NA_real_ else unname(loo_est["p_loo", "Estimate"]),
    bayes_R2_median = r2,
    rhat_max = d$rhat_max,
    ess_bulk_min = d$ess_bulk_min,
    divergences = d$divergences,
    stringsAsFactors = FALSE
  )
}

main_table <- bind_rows(lapply(names(fits_main), function(nm) {
  criterion_row(fits_main[[nm]], nm, model_specs[[nm]]$text)
}))

loos <- setNames(lapply(fits_main, function(f) f$criteria$loo), names(fits_main))
comp_df <- tryCatch({
  comp <- loo_compare(loos)
  out <- as.data.frame(comp)
  out$model <- rownames(out)
  out[, c("model", "elpd_diff", "se_diff")]
}, error = function(e) {
  warning("LOO comparison failed: ", conditionMessage(e), call. = FALSE)
  data.frame(model = names(fits_main), elpd_diff = NA_real_, se_diff = NA_real_)
})
main_table <- left_join(main_table, comp_df, by = "model")
main_table$rank <- rank(-main_table$elpd_loo, ties.method = "min")
write_tsv(main_table, file.path(opt$outdir, "model_comparison_main.tsv"))

if (length(fits_sens)) {
  sens_table <- bind_rows(lapply(names(fits_sens), function(nm) {
    criterion_row(fits_sens[[nm]], nm, sens_specs[[nm]]$text)
  }))
} else {
  sens_table <- data.frame(
    model = character(), formula = character(), n_obs = integer(),
    waic = numeric(), waic_se = numeric(), looic = numeric(), looic_se = numeric(),
    elpd_loo = numeric(), p_loo = numeric(), bayes_R2_median = numeric(),
    rhat_max = numeric(), ess_bulk_min = numeric(), divergences = numeric(),
    stringsAsFactors = FALSE
  )
}
sens_table$note <- "Models use different rows (offset/terminal/order-stem/mixed filters); LOO/WAIC are descriptive, not a head-to-head test."
write_tsv(sens_table, file.path(opt$outdir, "model_comparison_sensitivity.tsv"))

# The paper's key directional claim: bat intercept > non-bat intercept.
hyp <- hypothesis(fits_main$m1, "taxon_modelbat > 0")
hyp_df <- as.data.frame(hyp$hypothesis)
write_tsv(hyp_df, file.path(opt$outdir, "hypothesis_taxon.tsv"))

# Observed vs expected and fitted curve --------------------------------------
pred <- fitted(fits_main$m1, newdata = dat, summary = TRUE)
dat_out <- dat %>%
  mutate(
    expected = pred[, "Estimate"],
    expected_lower = pred[, "Q2.5"],
    expected_upper = pred[, "Q97.5"],
    residual = count - expected
  ) %>%
  select(branch_id, label, order, taxon_raw, branch_type, is_order_stem,
         bl, n_genes_tested, count, expected, expected_lower, expected_upper, residual)
write_tsv(dat_out, file.path(opt$outdir, "branch_observed_expected.tsv"))

grid <- expand.grid(
  bl = seq(min(dat$bl), max(dat$bl), length.out = 200),
  taxon_model = factor(c("non-bat", "bat"), levels = c("non-bat", "bat"))
)
grid_pred <- fitted(fits_main$m1, newdata = grid, summary = TRUE)
grid$expected <- grid_pred[, "Estimate"]
grid$expected_lower <- grid_pred[, "Q2.5"]
grid$expected_upper <- grid_pred[, "Q97.5"]
write_tsv(grid, file.path(opt$outdir, "fitted_curve.tsv"))

p <- ggplot() +
  geom_ribbon(data = grid, aes(x = bl, ymin = expected_lower, ymax = expected_upper,
                               fill = taxon_model), alpha = 0.2) +
  geom_line(data = grid, aes(x = bl, y = expected, colour = taxon_model), linewidth = 0.8) +
  geom_point(data = dat, aes(x = bl, y = count, colour = taxon_model,
                             shape = branch_type), size = 2, alpha = 0.75) +
  labs(x = "Branch length (substitutions/site)", y = "Selected genes in GO term",
       colour = "Taxon", fill = "Taxon", shape = "Branch type") +
  theme_minimal(base_size = 12)
ggsave(file.path(opt$outdir, "regression_curve.png"), p, width = 7, height = 5, dpi = 150)

png(file.path(opt$outdir, "pp_check_m1.png"), width = 1200, height = 800, res = 150)
print(pp_check(fits_main$m1, ndraws = 100) + theme_minimal(base_size = 12))
dev.off()

# Diagnostics text -----------------------------------------------------------
diag_lines <- c(
  paste("R version:", R.version.string),
  paste("brms:", as.character(packageVersion("brms"))),
  paste("seed:", opt$seed),
  paste("chains/iter/warmup:", opt$chains, opt$iter, opt$warmup),
  paste("usable branches:", nrow(dat)),
  paste("main model max Rhat:", signif(max(main_table$rhat_max, na.rm = TRUE), 4)),
  paste("main model min bulk ESS:", signif(min(main_table$ess_bulk_min, na.rm = TRUE), 4)),
  paste("main model divergences:", paste(main_table$model, main_table$divergences, collapse = "; ")),
  "",
  "Interpretation note:",
  "  m3 was called 'two intercepts, two slopes' in the original code, but",
  "  count ~ sqrt(bl):taxon has one intercept and two slopes.  m2 is the",
  "  actual two-intercept/two-slope model.  See docs/03_findings_zh.md.",
  "  The optional offset model uses n_genes_tested (GO genes with a tested",
  "  branch), which the original analysis did not include."
)
writeLines(diag_lines, file.path(opt$outdir, "diagnostics.txt"))
capture.output(sessionInfo(), file = file.path(opt$outdir, "sessionInfo.txt"))

message("done: ", opt$outdir)
