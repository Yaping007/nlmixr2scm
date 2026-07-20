# ---- probe: optExpression TRUE vs FALSE on the ifoceif ODE build ------------
# Finding (nlmixr2est 6.2.0): the FOCEi analytic outer-gradient CSE pass in
# .foceiAnalyticAugModelDirs is ALREADY serial + chunked (vapply over 40-line
# chunks, no daemons). The old "7 chunks, 4 daemons" parallel deadlock was
# fixed upstream by that chunked rewrite. So rxThreads does NOT gate this pass;
# the only real question is the build-vs-execute tradeoff:
#
#   optExpression=FALSE : cheap build, redundant expressions re-evaluated every
#                         gradient/ODE step (our current setting, 307 s).
#   optExpression=TRUE  : longer build (find+substitute duplicate expressions
#                         across the huge augmented gradient model), but leaner
#                         per-iteration code.
#
# This probe times ifoceif with optExpression=TRUE and compares to the saved
# optExpression=FALSE fit (307 s). setTimeLimit caps a runaway build so a hang
# aborts instead of stalling the console.
suppressPackageStartupMessages({ library(nlmixr2); library(rxode2) })

sd <- "script"
source(file.path(sd, "refit_helpers.R"),      chdir = FALSE)
source(file.path(sd, "true_model_factory.R"), chdir = FALSE)
source(file.path(sd, "estimator_factory.R"),  chdir = FALSE)

N <- 80L; scn <- 16L; ds <- 1L; bnd <- "narrow"; opt <- "lbfgsb3c"

master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scn))
sim_slice <- load_scenario_dataset(master_rds, scn, ds)
ds_nm <- to_nm_dataset(sim_slice) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))
true_mod <- make_true_model(scn, boundary = bnd, structure = "ode")

# Start from the exact irlsfoceif control, then flip optExpression back ON.
cfg <- make_est_control("irlsfoceif", outer_opt = opt, tier = "final")
ctrl_optT <- cfg$ctrl
ctrl_optT$optExpression <- TRUE   # the one knob under test

cat("rxThreads (unchanged):", rxode2::getRxThreads(), "\n")
cat(">>> ifoceif with optExpression=TRUE ...\n")

BUILD_CAP_SEC <- 900   # abort if the build/fit runs away (>15 min)
t0 <- Sys.time()
fit_optT <- tryCatch({
  setTimeLimit(elapsed = BUILD_CAP_SEC, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  nlmixr2(true_mod, ds_nm, est = "ifoceif", control = ctrl_optT)
}, error = function(e) { message("optExpression=TRUE aborted: ",
                                 conditionMessage(e)); NULL })
t_optT <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

cat("\n================ optExpression: TRUE vs FALSE ================\n")
if (!is.null(fit_optT)) {
  cat(sprintf("%-22s %10s %12s %14s\n", "config", "time_s", "objf", "grad(header)"))
  cat(sprintf("%-22s %10.1f %12.4f %14s\n", "optExpression=TRUE",
              t_optT, fit_optT$objf,
              if (grepl("analytic", paste(utils::capture.output(print(fit_optT)),
                                          collapse=" "))) "analytic" else "FD"))
  cat(sprintf("%-22s %10s %12.4f %14s\n", "optExpression=FALSE (saved)",
              "307.0", -1591.4887, "analytic"))
  cat(sprintf("\nspeedup (FALSE->TRUE): %+.1f s (%.0f%%)\n",
              307.0 - t_optT, 100 * (307.0 - t_optT) / 307.0))
  cat("objf match:", abs(fit_optT$objf - (-1591.4887)) < 1e-2, "\n")
  cat("$time setup:", round(fit_optT$time$setup, 1),
      " optimize:", round(fit_optT$time$optimize, 1),
      " covariance:", round(fit_optT$time$covariance, 1), "\n")
} else {
  cat("optExpression=TRUE did NOT complete within", BUILD_CAP_SEC, "s.\n")
  cat("=> keep optExpression=FALSE (current setting).\n")
}
