# ============================================================================
# make_composite_figs_est703.R
# Three-method comparison composite (nlmixr2-SCM / PsN-SCM / nlmixr2-VAE) for
# the est-7.0.3 rerun, mirroring output/figures/Composite_graph_Est702 but
# pointed at the new aggregated folders.
#
#   nlmixr2-SCM : output/scm_focei_bobyqa_703est_0805_aggregated  (focei/bobyqa)
#   PsN-SCM     : output/psn_scm_combined_aggregated              (unchanged;
#                 NONMEM does not depend on the nlmixr2est version)
#   nlmixr2-VAE : output/vae_covsel_full0807_est703_aggregated    (est 7.0.3)
#
# The three-way figure functions deliberately fold linCmt/advan4 -> "Analytic"
# and ode -> "ODE" so the mixed nlmixr2 + PsN panels share one structure axis;
# that folding is KEPT here (per request).
#
# Reproduces the four Est702 sub-folders:
#   threeway_N80/    power + diagnostic three-way curves
#   covsel_compare/  covariate-selection difference + agreement
#   error_compare/   estimation-accuracy three-way + difference
#   runtime_compare/ wall-time by-method + three-way + grid
#
# Usage:  source("script/viz/make_composite_figs_est703.R")
# ============================================================================

source("script/viz/fig_power_threeway.R")
source("script/viz/fig_diag_threeway.R")
source("script/viz/fig_covsel_compare.R")
source("script/viz/fig_covsel_heatmap.R")
source("script/viz/fig_error_compare.R")
source("script/viz/fig_runtime_compare.R")

## ---- est-7.0.3 folders ------------------------------------------------------
SCM_DIR <- "output/scm_focei_bobyqa_703est_0805_aggregated"
PSN_DIR <- "output/psn_scm_combined_aggregated"
VAE_DIR <- "output/vae_covsel_full0807_est703_aggregated"

COMP <- "output/figures/Composite_garphEst703"   # existing target folder

## helper: build a per-method source list with swapped CSV paths --------------
.mk_src <- function(scm_path, psn_path, vae_path, extra_scm = list(),
                    extra_psn = list(), extra_vae = list()) {
  list(
    "nlmixr2-SCM" = c(list(path = scm_path, estimator = "focei",
                           outer_opt = "bobyqa"), extra_scm),
    "PsN-SCM"     = c(list(path = psn_path, estimator = "nonmem_scm",
                           outer_opt = "focei"), extra_psn),
    "nlmixr2-VAE" = c(list(path = vae_path, estimator = NULL,
                           outer_opt = NULL), extra_vae)
  )
}

## ---- 1. power three-way -----------------------------------------------------
pow_src <- .mk_src(
  file.path(SCM_DIR, "scm_power.csv"),
  file.path(PSN_DIR, "scm_power.csv"),
  file.path(VAE_DIR, "vae_power.csv"))

fig_power_threeway_N80(sources = pow_src, sample_N = 80, metric = "Power",
                       save = TRUE, out_dir = file.path(COMP, "threeway_N80"))
fig_power_threeway_grid(sources = pow_src, save = TRUE,
                        out_dir = file.path(COMP, "threeway_N80"))

## ---- 2. diagnostic three-way ------------------------------------------------
diag_src <- .mk_src(
  file.path(SCM_DIR, "scm_diag_rates.csv"),
  file.path(PSN_DIR, "scm_diag_rates.csv"),
  file.path(VAE_DIR, "vae_diag_rates.csv"))

fig_diag_threeway_grid(sources = diag_src, save = TRUE,
                       out_dir = file.path(COMP, "threeway_N80"))

## ---- 3. covariate-selection comparison --------------------------------------
cov_src <- .mk_src(
  file.path(SCM_DIR, "scm_covsel_by_covar.csv"),
  file.path(PSN_DIR, "scm_covsel_by_covar.csv"),
  file.path(VAE_DIR, "vae_covsel_by_covar.csv"))

fig_covsel_diff(sources = cov_src, structure = "linCmt", sample_N = 80,
                save = TRUE, out_dir = file.path(COMP, "covsel_compare"))
fig_covsel_agreement(sources = cov_src, structure = "linCmt", sample_N = 80,
                     save = TRUE, out_dir = file.path(COMP, "covsel_compare"))

## ---- 3b. per-method covariate-selection PATTERN heatmaps (CL only) ---------
## CL-focused, Analytic|ODE side-by-side, N = 40 vs 300 to show the FP/FN
## sample-size contrast.  These use the per-method heatmap (NOT the diff), so
## they carry the enlarged in-tile numbers / axis + strip fonts / panel gap.
fig_covsel_heatmap_scm(
  agg_dir   = SCM_DIR, structure = "both", params = "CL", sample_N = c(40, 300),
  estimator = "focei", outer_opt = "bobyqa",
  save = TRUE, out_dir = file.path(COMP, "covsel_pattern"))
fig_covsel_heatmap(
  agg_dir   = VAE_DIR, structure = "both", params = "CL", sample_N = c(40, 300),
  save = TRUE, out_dir = file.path(COMP, "covsel_pattern"))

## ---- 4. estimation-accuracy comparison --------------------------------------
err_src <- list(
  "nlmixr2-SCM" = list(dir = SCM_DIR, csv_all = "scm_estim_all.csv",
                       csv_cond = "scm_estim_cond.csv",
                       estimator = "focei", outer_opt = "bobyqa"),
  "PsN-SCM"     = list(dir = PSN_DIR, csv_all = "scm_estim_all.csv",
                       csv_cond = "scm_estim_cond.csv",
                       estimator = "nonmem_scm", outer_opt = "focei"),
  "nlmixr2-VAE" = list(dir = VAE_DIR, csv_all = "vae_estim_all.csv",
                       csv_cond = "vae_estim_cond.csv",
                       estimator = NULL, outer_opt = NULL)
)

fig_error_threeway(sources = err_src, scenario = 16, structure = "linCmt", x_cap = 70,
                   sample_N = 80, metric = "MARE", save = TRUE,
                   out_dir = file.path(COMP, "error_compare"))
fig_error_threeway(sources = err_src, scenario = 16, structure = "linCmt", x_cap = 70,
                   sample_N = 80, metric = "RMRSE", save = TRUE,
                   out_dir = file.path(COMP, "error_compare"))

fig_error_threeway_diff(sources = err_src, scenario = 16, structure = "linCmt",
                        sample_N = 80, metric = "MARE", save = TRUE,
                        out_dir = file.path(COMP, "error_compare"))

## ---- 5. runtime comparison --------------------------------------------------
rt_src <- .mk_src(
  file.path(SCM_DIR, "scm_diag_long.csv"),
  file.path(PSN_DIR, "scm_diag_long.csv"),
  file.path(VAE_DIR, "vae_diag_long.csv"),
  extra_scm = list(wall = "wall_total_sec"),
  extra_psn = list(wall = "wall_total_sec"),
  extra_vae = list(wall = "fit_runtime_sec"))

for (m in names(rt_src))
  fig_runtime_bymethod(method = m, sources = rt_src, save = TRUE,
                       out_dir = file.path(COMP, "runtime_compare"))
fig_runtime_threeway(sources = rt_src, sample_N = 80, save = TRUE,
                     out_dir = file.path(COMP, "runtime_compare"))
fig_runtime_grid(sources = rt_src, log_y = TRUE, save = TRUE,
                 out_dir = file.path(COMP, "runtime_compare"))
fig_runtime_grid(sources = rt_src, log_y = FALSE, save = TRUE,
                 out_dir = file.path(COMP, "runtime_compare"))

message("\nAll composite (est 7.0.3) figures written under ", COMP)
