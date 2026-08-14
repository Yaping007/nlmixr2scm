# ============================================================================
# plot_ode_figures.R  --  PsN/NONMEM SCM (ODE) figure driver
# ----------------------------------------------------------------------------
# Renders the four headline figures for the NONMEM-ODE SCM sweep, using the
# SAME viz functions + calling convention as the advan4 run (README §7):
#   structure = "ode", estimator = "nonmem_scm", outer_opt = "focei".
#
# Reads:  output/psn_scm_ode0729_aggregated/scm_*.csv
# Writes: output/figures/psn_scm_ode/fig_*.{png,pdf}
#
# Run from the repo root:
#   Rscript script/psn_scm/plot_ode_figures.R
# ============================================================================

source("script/viz/fig_diag_rates.R")
source("script/viz/fig_power.R")
source("script/viz/fig_covsel_heatmap.R")
source("script/viz/fig_error_metrics.R")

agg <- "output/psn_scm_ode0729_aggregated"
out <- "output/figures/psn_scm_ode"

# ---- 1. Convergence / run-quality rates -----------------------------------
p_conv <- fig_diag_rates(
  agg, csv_name = "scm_diag_rates.csv",
  metric = c("Converged", "CNBelowCutoff", "ConvergedStrict"),
  structure = "ode", estimator = "nonmem_scm", outer_opt = "focei",
  title_prefix = "PsN SCM (ODE)", save = TRUE, out_dir = out)

# ---- 2. Power curves ------------------------------------------------------
p_pow <- fig_power(
  agg, csv_name = "scm_power.csv",
  metric = c("Power", "PowerCN", "PowerMinSuc"),
  structure = "ode", estimator = "nonmem_scm", outer_opt = "focei",
  title_prefix = "PsN SCM (ODE)", save = TRUE, out_dir = out)

# ---- 3. Covariate-selection pattern heatmap -------------------------------
p_sel <- fig_covsel_heatmap_scm(
  agg_dir = agg, csv_name = "scm_covsel_by_covar.csv",
  structure = "ode", estimator = "nonmem_scm", outer_opt = "focei",
  save = TRUE, out_dir = out)

# ---- 4. Estimation accuracy (RMRSE / MARE), scenario 16 -------------------
p_est <- fig_error_metrics(
  agg, scenario = 16, structure = "ode",
  estimator = "nonmem_scm", outer_opt = "focei",
  metrics = c("RMRSE", "MARE"),
  csv_all = "scm_estim_all.csv", csv_cond = "scm_estim_cond.csv",
  save = TRUE, out_dir = out)

message("\nAll four ODE figures written to ", out)
