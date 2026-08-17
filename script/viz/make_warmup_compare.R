# ============================================================================
# make_warmup_compare.R
# Contrast the SCM warm-start (Mechanism B, profileInitOnStall) ON vs OFF.
#
#   warmup-on  : output/scm_focei_bobyqa_703est_0805_aggregated  (package
#                default -> profileInitOnStall = TRUE)
#   warmup-off : output/scm_profile_off_aggregated               (profileInit
#                warm-start disabled)
#
# Both runs are focei / bobyqa.  The profile-off run only carries the
# N = 80 / linCmt cell, so ALL comparisons are hard-restricted to N = 80,
# structure = linCmt for an apples-to-apples contrast.
#
# Figures (output/figures/warmup_compare/):
#   paired/   on_ / off_ convergence + power (single-source functions)
#   covsel/   signed-difference heatmap (on - off) + per-scenario agreement
#   error/    three-way dots (on vs off) + signed difference, MARE + RMRSE
#
# Usage:  source("script/viz/make_warmup_compare.R")
# ============================================================================

source("script/viz/fig_diag_rates.R")
source("script/viz/fig_power.R")
source("script/viz/fig_covsel_compare.R")
source("script/viz/fig_error_compare.R")

ON_DIR  <- "output/scm_focei_bobyqa_703est_0805_aggregated"
OFF_DIR <- "output/scm_profile_off_aggregated"
FIG     <- "output/figures/warmup_compare"

EST <- "focei"; OUT <- "bobyqa"
SN  <- 80;      STR <- "linCmt"

## register palette colours for the two warm-up conditions so the reused
## compare functions (keyed on method name) render with real colours.
.PAL_METHOD_E[["warmup-on"]]  <- "#0F8B8D"   # teal
.PAL_METHOD_E[["warmup-off"]] <- "#C1666B"   # muted red
.PAL_METHOD_C[["warmup-on"]]  <- "#0F8B8D"
.PAL_METHOD_C[["warmup-off"]] <- "#C1666B"

## ---- 1. paired convergence + power (single-source; on then off) ------------
paired_dir <- file.path(FIG, "paired")

for (cond in c("on", "off")) {
  dir_c  <- if (cond == "on") ON_DIR else OFF_DIR
  prefix <- if (cond == "on") "warmup ON" else "warmup OFF"

  fig_diag_rates(
    agg_dir = dir_c, metric = c("Converged", "CNBelowCutoff"),
    sample_N = SN, structure = STR, estimator = EST, outer_opt = OUT,
    save = TRUE, out_dir = file.path(paired_dir, cond),
    csv_name = "scm_diag_rates.csv", title_prefix = prefix)

  fig_power(
    agg_dir = dir_c, metric = c("Power", "PowerCN", "PowerMinSuc"),
    sample_N = SN, structure = STR, estimator = EST, outer_opt = OUT,
    save = TRUE, out_dir = file.path(paired_dir, cond),
    csv_name = "scm_power.csv", title_prefix = prefix)
}

## ---- 2. covariate-selection difference (on - off) --------------------------
cov_src <- list(
  "warmup-on"  = list(path = file.path(ON_DIR,  "scm_covsel_by_covar.csv"),
                      estimator = EST, outer_opt = OUT),
  "warmup-off" = list(path = file.path(OFF_DIR, "scm_covsel_by_covar.csv"),
                      estimator = EST, outer_opt = OUT)
)

fig_covsel_diff(sources = cov_src, reference = "warmup-off",
                methods = "warmup-on", structure = STR, sample_N = SN,
                save = TRUE, out_dir = file.path(FIG, "covsel"))
fig_covsel_agreement(sources = cov_src, reference = "warmup-off",
                     methods = "warmup-on", structure = STR, sample_N = SN,
                     save = TRUE, out_dir = file.path(FIG, "covsel"))

## ---- 3. estimation-accuracy comparison (on vs off) -------------------------
err_src <- list(
  "warmup-on"  = list(dir = ON_DIR,  csv_all = "scm_estim_all.csv",
                      csv_cond = "scm_estim_cond.csv",
                      estimator = EST, outer_opt = OUT),
  "warmup-off" = list(dir = OFF_DIR, csv_all = "scm_estim_all.csv",
                      csv_cond = "scm_estim_cond.csv",
                      estimator = EST, outer_opt = OUT)
)

for (met in c("MARE", "RMRSE")) {
  fig_error_threeway(
    sources = err_src, methods = c("warmup-on", "warmup-off"),
    scenario = 16, structure = STR, sample_N = SN, metric = met, x_cap = 70,
    save = TRUE, out_dir = file.path(FIG, "error"))
}
fig_error_threeway_diff(
  sources = err_src, reference = "warmup-off", methods = "warmup-on",
  scenario = 16, structure = STR, sample_N = SN, metric = "MARE",
  save = TRUE, out_dir = file.path(FIG, "error"))

message("\nAll warm-up on/off comparison figures written under ", FIG)
