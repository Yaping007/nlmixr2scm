# ============================================================================
# plot_combined_figures.R  --  PsN/NONMEM SCM: advan4 + ODE COMBINED figures
# ----------------------------------------------------------------------------
# Draws the four headline figures with BOTH PsN structural platforms side by
# side -- NONMEM analytic (ADVAN4) and NONMEM general-ODE (ADVAN13) -- exactly
# the way the nlmixr2 linCmt|ode (and focei|vae) comparisons were rendered:
# the two platforms are separate `structure` values faceted within ONE figure.
#
# The two PsN sweeps live in SEPARATE aggregated dirs:
#   advan4 : output/psn_scm_full0727_aggregated1  (structure = "advan4")
#   ode    : output/psn_scm_ode0729_aggregated    (structure = "ode")
# Their CSV schemas are identical, so we ROW-BIND each scm_*.csv into one
# combined dir; the `structure` column then carries the platform, and the
# existing facet-by-structure logic in every fig_* function does the rest.
#
# Run from the repo root:
#   Rscript script/psn_scm/plot_combined_figures.R
# ============================================================================

suppressPackageStartupMessages({library(readr); library(dplyr)})

source("script/viz/fig_diag_rates.R")
source("script/viz/fig_power.R")
source("script/viz/fig_covsel_heatmap.R")
source("script/viz/fig_error_metrics.R")

adv <- "output/psn_scm_full0727_aggregated1"  # advan4 (analytic; fuller parse)
ode <- "output/psn_scm_ode0729_aggregated"    # ODE (ADVAN13)
cmb <- "output/psn_scm_combined_aggregated"   # merged advan4 + ode
out <- "output/figures/psn_scm_combined"

# ---- merge: row-bind every scm_*.csv present in BOTH dirs ------------------
# Read all columns as character so per-file type guesses (e.g. an all-NA column
# guessed <chr> in one dir but <dbl> in the other) can't block the row-bind; the
# downstream fig_* re-read re-guesses types normally from the written text.
dir.create(cmb, showWarnings = FALSE, recursive = TRUE)
csvs <- intersect(list.files(adv, "\\.csv$"), list.files(ode, "\\.csv$"))
for (f in csvs) {
  a <- readr::read_csv(file.path(adv, f), col_types = readr::cols(.default = "c"))
  o <- readr::read_csv(file.path(ode, f), col_types = readr::cols(.default = "c"))
  readr::write_csv(dplyr::bind_rows(a, o), file.path(cmb, f), na = "")
}
message("merged ", length(csvs), " CSVs -> ", cmb,
        "  (structures: ",
        paste(unique(readr::read_csv(file.path(cmb, "scm_power.csv"),
                                     show_col_types = FALSE)$structure),
              collapse = ", "), ")")

# ---- 1. Convergence / run-quality rates (advan4 | ODE facets) --------------
# Keep only the two headline gates: Converged (%) and Cond# < 1000 (%).
p_conv <- fig_diag_rates(
  cmb, csv_name = "scm_diag_rates.csv",
  metric = c("Converged", "CNBelowCutoff"),
  structure = NULL,                    # NULL -> both platforms as facet columns
  estimator = "nonmem_scm", outer_opt = "focei",
  title_prefix = "PsN SCM", save = TRUE, out_dir = out)

# ---- 2. Power curves (advan4 | ODE facets) --------------------------------
p_pow <- fig_power(
  cmb, csv_name = "scm_power.csv",
  metric = c("Power", "PowerCN", "PowerMinSuc"),
  structure = NULL,
  estimator = "nonmem_scm", outer_opt = "focei",
  title_prefix = "PsN SCM", save = TRUE, out_dir = out)

# ---- 3. Covariate-selection heatmap (advan4 + ODE side by side) -----------
p_sel <- fig_covsel_heatmap_scm(
  agg_dir = cmb, csv_name = "scm_covsel_by_covar.csv",
  structure = "both",                  # "both" -> N(rows) x structure(cols) grid
  estimator = "nonmem_scm", outer_opt = "focei",
  save = TRUE, out_dir = out)

# ---- 4. Estimation accuracy (RMRSE / MARE), scn16, both platforms ---------
p_est <- fig_error_metrics(
  cmb, scenario = 16, structure = NULL,
  estimator = "nonmem_scm", outer_opt = "focei",
  metrics = c("RMRSE", "MARE"),
  csv_all = "scm_estim_all.csv", csv_cond = "scm_estim_cond.csv",
  save = TRUE, out_dir = out)

message("\nAll four COMBINED (advan4 + ODE) figures written to ", out)
