# ============================================================================
# make_poweronly_N80_threeway.R
# ----------------------------------------------------------------------------
# Power-SHAPE-ONLY covariate space, N = 80, three methods overlaid:
#   nlmixr2-SCM (focei + bobyqa) | PsN-SCM (nonmem_scm) | nlmixr2-VAE
#
# Renders the four reference comparison figures (same layout as the
# competing-shape threeway figures) for the power-only runs:
#   (a) convergence + Cond# < 1000   -> fig_diag_threeway_grid
#   (b) selection power              -> fig_power_threeway_grid
#   (c) runtime (median, IQR band)   -> fig_runtime_threeway
#   (d) estimation accuracy (scn16)  -> fig_error_threeway
#
# INPUT TREES (power-only aggregated dirs)
#   nlmixr2-SCM : output/scm_poweronly_est703_08132026_aggregated  (focei/bobyqa)
#   PsN-SCM     : output/psn_scm_poweronly_advan4_aggregated  (Analytic)  +
#                 output/psn_scm_poweronly_ode_aggregated     (ODE)
#                 -> row-bound into ONE combined dir (PsN ships split by struct)
#   nlmixr2-VAE : output/vae_covsel_poweronly_N80_08162026_aggregated
#
# OUTPUT: output/figures/poweronly_N80/
#
# Usage:  source("script/viz/make_poweronly_N80_threeway.R")
# ============================================================================

suppressPackageStartupMessages({ library(readr); library(dplyr) })

## ---- 0. source the four figure families (defines loaders + fns) ------------
invisible(lapply(
  c("fig_diag_threeway.R", "fig_power_threeway.R",
    "fig_runtime_compare.R", "fig_error_compare.R"),
  function(f) source(file.path("script/viz", f))))

OUT <- "output/figures/poweronly_N80"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 1. merge the split PsN dirs (advan4 + ode) into ONE combined dir ------
PSN_A4  <- "output/psn_scm_poweronly_advan4_aggregated"
PSN_ODE <- "output/psn_scm_poweronly_ode_aggregated"
PSN_CMB <- "output/psn_scm_poweronly_combined_aggregated"
dir.create(PSN_CMB, showWarnings = FALSE, recursive = TRUE)

for (f in c("scm_power.csv", "scm_diag_rates.csv", "scm_diag_long.csv",
            "scm_estim_all.csv", "scm_estim_cond.csv")) {
  parts <- Filter(Negate(is.null), lapply(
    c(file.path(PSN_A4, f), file.path(PSN_ODE, f)),
    function(p) if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL))
  if (length(parts))
    write_csv(bind_rows(parts), file.path(PSN_CMB, f))
}

## ---- 2. power-only sources for each figure family --------------------------
SCM <- "output/scm_poweronly_est703_08132026_aggregated"
VAE <- "output/vae_covsel_poweronly_N80_08162026_aggregated"

.TW_SOURCES <- list(
  "nlmixr2-SCM" = list(path = file.path(SCM, "scm_power.csv"),
                       estimator = "focei", outer_opt = "bobyqa"),
  "PsN-SCM"     = list(path = file.path(PSN_CMB, "scm_power.csv"),
                       estimator = "nonmem_scm", outer_opt = "focei"),
  "nlmixr2-VAE" = list(path = file.path(VAE, "vae_power.csv"),
                       estimator = NULL, outer_opt = NULL))

.TWD_SOURCES <- list(
  "nlmixr2-SCM" = list(path = file.path(SCM, "scm_diag_rates.csv"),
                       estimator = "focei", outer_opt = "bobyqa"),
  "PsN-SCM"     = list(path = file.path(PSN_CMB, "scm_diag_rates.csv"),
                       estimator = "nonmem_scm", outer_opt = "focei"),
  "nlmixr2-VAE" = list(path = file.path(VAE, "vae_diag_rates.csv"),
                       estimator = NULL, outer_opt = NULL))

.RT_SOURCES <- list(
  "nlmixr2-SCM" = list(path = file.path(SCM, "scm_diag_long.csv"),
                       wall = "wall_total_sec", estimator = "focei", outer_opt = "bobyqa"),
  "PsN-SCM"     = list(path = file.path(PSN_CMB, "scm_diag_long.csv"),
                       wall = "wall_total_sec", estimator = "nonmem_scm", outer_opt = "focei"),
  "nlmixr2-VAE" = list(path = file.path(VAE, "vae_diag_long.csv"),
                       wall = "fit_runtime_sec", estimator = NULL, outer_opt = NULL))

.ERR_SOURCES <- list(
  "nlmixr2-SCM" = list(dir = SCM, csv_all = "scm_estim_all.csv",
                       csv_cond = "scm_estim_cond.csv",
                       estimator = "focei", outer_opt = "bobyqa"),
  "PsN-SCM"     = list(dir = PSN_CMB, csv_all = "scm_estim_all.csv",
                       csv_cond = "scm_estim_cond.csv",
                       estimator = "nonmem_scm", outer_opt = "focei"),
  "nlmixr2-VAE" = list(dir = VAE, csv_all = "vae_estim_all.csv",
                       csv_cond = "vae_estim_cond.csv",
                       estimator = NULL, outer_opt = NULL))

## ---- 3. render the four figures (single N=80 -> grid auto-restricts) --------
# (a) convergence + Cond# < 1000  (method cols x structure rows)
p_diag <- fig_diag_threeway_grid(.TWD_SOURCES, save = TRUE, out_dir = OUT)

# (b) power (Power / PowerCN / PowerMinSuc linetypes)
p_pow  <- fig_power_threeway_grid(.TW_SOURCES, save = TRUE, out_dir = OUT)

# (c) runtime (structure facets, methods coloured, IQR band)
p_rt   <- fig_runtime_threeway(.RT_SOURCES, sample_N = 80, save = TRUE, out_dir = OUT)

# (d) estimation accuracy (scenario 16, both structures)
p_err_a <- fig_error_threeway(.ERR_SOURCES, structure = "linCmt",
                              sample_N = 80, save = TRUE, out_dir = OUT)
p_err_o <- fig_error_threeway(.ERR_SOURCES, structure = "ode",
                              sample_N = 80, save = TRUE, out_dir = OUT)

message("Power-only N=80 threeway figures written to ", OUT)
