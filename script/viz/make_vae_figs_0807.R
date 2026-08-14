# ============================================================================
# make_vae_figs_0807.R
# Regenerate the four nlmixr2-VAE covariate-selection scorecard-companion
# figures for the est-7.0.3 run (aggregated 2026-08-07), mirroring
# make_scm_figs_0805.R but pointed at the VAE aggregated folder.
#
#   1. convergence / cov-step run-quality curves   (fig_diag_rates)
#   2. power curves (Power / PowerCN / PowerMinSuc) (fig_power)
#   3. covariate-selection error-pattern heatmaps   (fig_covsel_heatmap)
#   4. estimation accuracy (RMRSE / MARE)           (fig_error_metrics)
#
# The VAE aggregated CSVs carry NO estimator / outer_opt columns (unlike SCM),
# so the estimator= / outer_opt= knobs are left at their NULL defaults.
#
# Usage:  source("script/viz/make_vae_figs_0807.R")
# ============================================================================

AGG_DIR <- "output/vae_covsel_full0807_est703_aggregated"
FIG_DIR <- "output/figures/nlmixr2_vae_0807"
PREFIX  <- "nlmixr2 VAE (est 7.0.3)"

source("script/viz/fig_diag_rates.R")
source("script/viz/fig_power.R")
source("script/viz/fig_covsel_heatmap.R")
source("script/viz/fig_error_metrics.R")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. convergence + cov-step run quality --------------------------------
fig_diag_rates(
  agg_dir      = AGG_DIR,
  metric       = c("Converged", "CNBelowCutoff"),
  save         = TRUE, out_dir = FIG_DIR,
  csv_name     = "vae_diag_rates.csv",
  title_prefix = PREFIX
)

# ---- 2. power curves ------------------------------------------------------
fig_power(
  agg_dir      = AGG_DIR,
  metric       = c("Power", "PowerCN", "PowerMinSuc"),
  save         = TRUE, out_dir = FIG_DIR,
  csv_name     = "vae_power.csv",
  title_prefix = PREFIX
)

# ---- 3. selection-pattern heatmaps ----------------------------------------
# (a) all N, both structures side-by-side
fig_covsel_heatmap(
  agg_dir   = AGG_DIR, sample_N = NULL, structure = "both",
  save      = TRUE, out_dir = FIG_DIR,
  csv_name  = "vae_covsel_by_covar.csv",
  title_prefix = PREFIX
)
# (b) N80 linCmt zoom (single structure close-up)
fig_covsel_heatmap(
  agg_dir   = AGG_DIR, sample_N = 80, structure = "linCmt", layout = "slide",
  save      = TRUE, out_dir = FIG_DIR,
  csv_name  = "vae_covsel_by_covar.csv",
  title_prefix = PREFIX
)

# ---- 4. estimation accuracy (RMRSE / MARE) --------------------------------
fig_error_metrics(
  agg_dir   = AGG_DIR, scenario = 16,
  metrics   = c("RMRSE", "MARE"),
  save      = TRUE, out_dir = FIG_DIR,
  csv_all   = "vae_estim_all.csv",
  csv_cond  = "vae_estim_cond.csv"
)

message("\nAll VAE figures written to ", FIG_DIR)
print(list.files(FIG_DIR, pattern = "\\.png$"))
