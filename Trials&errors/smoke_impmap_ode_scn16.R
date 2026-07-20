# ==============================================================================
# smoke_impmap_ode_scn16.R
# ------------------------------------------------------------------------------
# Phase 0.2 gate: does importance-sampling EM (impmap / imp) actually FIT the
# bounded-theta ODE model and return a GENUINE objective function value?
#
# This is the test Gaussian quadrature failed for SAEM.  We compare, on
# scn16 / ds01 / N80, the ODE-form true model fit with:
#     est = "saem"    (reference: no native OFV, borrows focei fallback)
#     est = "imp"     (importance sampling EM, no MAP proposal)
#     est = "impmap"  (importance sampling EM, MAP proposal)  <- primary
#
# It also resolves the NAMING question: we print each fit's $est string and
# whether it carries a real $objf, so we know whether impmap is a standalone
# estimator (-> cell name impmap_NA) or a saem+IS pipeline (-> saem_impmap).
#
# Runs against the DEV scratch library (~/R-nlmixr2-dev, nlmixr2est 6.2.0),
# leaving the stable main library untouched.
#
# Usage:
#   Rscript script/smoke_impmap_ode_scn16.R
# ==============================================================================

## ---- Use the dev scratch library ------------------------------------------
scratch <- path.expand("~/R-nlmixr2-dev")
if (dir.exists(scratch)) .libPaths(c(scratch, .libPaths()))
cat("nlmixr2est lib:", find.package("nlmixr2est"), "\n")

suppressPackageStartupMessages({
  library(nlmixr2)
  library(nlmixr2est)
  library(dplyr)
})
cat("nlmixr2est version:", as.character(packageVersion("nlmixr2est")), "\n")
stopifnot("impmapControl" %in% getNamespaceExports("nlmixr2est"))

## ---- Sources: helpers + ODE factory ---------------------------------------
.script_dir <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)))),
  error = function(e) file.path(getwd(), "script")
)
source(file.path(.script_dir, "refit_helpers.R"),      chdir = FALSE)  # PsN_scenarios, load_scenario_dataset, to_nm_dataset, diagnose_fit
source(file.path(.script_dir, "true_model_factory.R"), chdir = FALSE)  # make_true_model_ode
source(file.path(.script_dir, "scm_bench_helpers.R"),  chdir = FALSE)  # base_2cmt_oral_ode (for reference)

## ---- Data: scn16 / ds01 / N80, NM-shaped ODE-ready ------------------------
N <- 80L; scenario_id <- 16L; dataset_id <- 1L
master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scenario_id))
stopifnot(file.exists(master_rds))

sim_slice <- load_scenario_dataset(master_rds, scenario_id, dataset_id)
ds_nm <- to_nm_dataset(sim_slice) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))

cat(sprintf("data: %d rows, %d subjects\n", nrow(ds_nm), dplyr::n_distinct(ds_nm$ID)))

## ---- ODE true model (bounded-theta, narrow) -------------------------------
mod_ode <- make_true_model_ode(scenario_id, "narrow")
cat("model structure attr:", attr(mod_ode, "structure"), "\n")

## ---- Fit helper: time + capture, never abort the whole run ----------------
out_dir <- file.path("output", "smoke", "impmap_ode_scn16_ds01_N80")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fit_one <- function(tag, est, control) {
  cat(sprintf("\n===== %-8s (est = \"%s\") =====\n", tag, est))
  t0 <- Sys.time()
  fit <- tryCatch(
    nlmixr2(mod_ode, ds_nm, est = est, control = control),
    error = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(fit)) {
    return(data.frame(tag = tag, est = est, fit_est = NA, runtime_s = dt,
                      objf = NA_real_, cond_num = NA_real_, status = "FAILED"))
  }
  saveRDS(fit, file.path(out_dir, sprintf("fit_%s.rds", tag)))
  print(fit)                                       # raw nlmixr2 fit to console
  objf <- tryCatch(as.numeric(fit$objf),               error = function(e) NA_real_)
  cn   <- tryCatch(as.numeric(fit$conditionNumberCor), error = function(e) NA_real_)
  data.frame(tag = tag, est = est,
             fit_est   = tryCatch(as.character(fit$est), error = function(e) NA),
             runtime_s = round(dt, 1),
             objf      = objf, cond_num = cn, status = "OK")
}

## ---- Run: saem reference, then imp / impmap -------------------------------
res <- dplyr::bind_rows(
  fit_one("saem",   "saem",   saemControl(print = 0)),
  fit_one("imp",    "imp",    impControl(print = 0)),
  fit_one("impmap", "impmap", impmapControl(print = 0))
)

cat("\n\n================ SUMMARY ================\n")
print(res, row.names = FALSE)
saveRDS(res, file.path(out_dir, "summary.rds"))
cat("\nsaved to:", out_dir, "\n")
