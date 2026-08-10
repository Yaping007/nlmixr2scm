# ============================================================================
# make_scm_figs_0805.R
# Regenerate the four nlmixr2-SCM scorecard-companion figures for the
# latest-CRAN-package run, mirroring output/figures/nlmixr2_scm_0730 but
# pointed at the new aggregated folder.
#
#   1. convergence / cov-step run-quality curves   (fig_diag_rates)
#   2. power curves (Power / PowerCN / PowerMinSuc) (fig_power)
#   3. covariate-selection error-pattern heatmaps   (fig_covsel_heatmap_scm)
#   4. estimation accuracy (RMRSE / MARE)           (fig_error_metrics)
#
# Usage:  source("script/viz/make_scm_figs_0805.R")
# ============================================================================

AGG_DIR <- "output/scm_focei_bobyqa_703est_0805_aggregated"
FIG_DIR <- "output/figures/nlmixr2_scm_0805"
EST     <- "focei"
OUT     <- "bobyqa"
PREFIX  <- "nlmixr2 SCM (CRAN)"

source("script/viz/fig_diag_rates.R")
source("script/viz/fig_power.R")
source("script/viz/fig_covsel_heatmap.R")
source("script/viz/fig_error_metrics.R")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. convergence + cov-step run quality --------------------------------
fig_diag_rates(
  agg_dir      = AGG_DIR,
  metric       = c("Converged", "CNBelowCutoff"),
  estimator    = EST, outer_opt = OUT,
  save         = TRUE, out_dir = FIG_DIR,
  csv_name     = "scm_diag_rates.csv",
  title_prefix = PREFIX
)

# ---- 2. power curves ------------------------------------------------------
fig_power(
  agg_dir      = AGG_DIR,
  metric       = c("Power", "PowerCN", "PowerMinSuc"),
  estimator    = EST, outer_opt = OUT,
  save         = TRUE, out_dir = FIG_DIR,
  csv_name     = "scm_power.csv",
  title_prefix = PREFIX
)

# ---- 3. selection-pattern heatmaps ----------------------------------------
# (a) all N, both structures side-by-side
fig_covsel_heatmap_scm(
  agg_dir   = AGG_DIR, sample_N = NULL, structure = "both",
  estimator = EST, outer_opt = OUT,
  save      = TRUE, out_dir = FIG_DIR,
  title_prefix = PREFIX
)
# (b) N80 linCmt zoom (single structure close-up)
fig_covsel_heatmap_scm(
  agg_dir   = AGG_DIR, sample_N = 80, structure = "linCmt",
  estimator = EST, outer_opt = OUT, layout = "slide",
  save      = TRUE, out_dir = FIG_DIR,
  title_prefix = PREFIX
)

# ---- 4. estimation accuracy (RMRSE / MARE) --------------------------------
fig_error_metrics(
  agg_dir   = AGG_DIR, scenario = 16,
  estimator = EST, outer_opt = OUT,
  metrics   = c("RMRSE", "MARE"),
  save      = TRUE, out_dir = FIG_DIR,
  csv_all   = "scm_estim_all.csv",
  csv_cond  = "scm_estim_cond.csv"
)

message("\nAll SCM figures written to ", FIG_DIR)
print(list.files(FIG_DIR, pattern = "\\.png$"))
