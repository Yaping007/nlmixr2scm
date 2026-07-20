# ============================================================================
# smoke_irlsfoceif_vae_scn16.R
# ----------------------------------------------------------------------------
# One-shot SMOKE TEST for the two estimators currently EXCLUDED from the
# benchmark grid: `irlsfoceif` and `vae`.  Runs BOTH:
#     (a) true-model REFIT   (no SCM search)
#     (b) SCM search         (runSCM)
# on a single cell:  N = 80, scenario = 16, dataset = 1, boundary = "narrow".
#
# Why this script exists
# ----------------------
# `is_valid_combo()` in estimator_factory.R deliberately blocks both
# estimators, so run_bench_cell() / bench_refit_one() refuse to run them.
# This script calls nlmixr2()/runSCM() directly, BYPASSING that gate, and
# captures whatever breaks (error message + traceback + fit object) so the
# problematic output can be forwarded to the manager for advice.
#
# Known, expected pain points (documented in estimator_factory.R header):
#   * irlsfoceif : fits, but its nlmixr2est class does NOT populate
#                  fit$conditionNumberCor; the cov2cor+kappa fallback may
#                  also fail -> never satisfies the PMx-strict convergence
#                  flag.  Captured, not "fixed".
#   * vae        : make_est_control() has no `vae` branch (would stop()), so
#                  we build vaeControl() by hand.  runSCM() is expected to
#                  ERROR because the vae-returned ui renames cl/vc ->
#                  lTVCL/lTVVc and the search cannot locate the parameters.
#                  Wrapped in tryCatch; the error is the deliverable.
#
# ALL estimator settings below are taken VERBATIM from the production
# pipeline (estimator_factory.R + PerformanceEvaluation_scm_bench.R +
# bench_refit_estimators.R) so the smoke output reflects the real thing.
#
# Output:
#   output/smoke/scn16_ds01_N80_irlsfoceif_vae.rds   (full capture bundle)
#   console echo of a compact status table
#
# Run:
#   Rscript script/smoke_irlsfoceif_vae_scn16.R
#   # or source() interactively
# ============================================================================

suppressPackageStartupMessages({
  library(nlmixr2)
  library(nlmixr2est)
  library(rxode2)
  library(dplyr)
  library(tibble)
  library(nlmixr2scm)
})

# ---- Locate helper scripts (repo-root relative) ----------------------------
.script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(fa)) return(normalizePath(dirname(fa[1]), mustWork = FALSE))
  if (dir.exists("script")) return("script")
  getwd()
})()

source(file.path(.script_dir, "refit_helpers.R"))       # to_nm_dataset,
                                                         # load_scenario_dataset,
                                                         # diagnose_fit_table3,
                                                         # PsN_scenarios, true_params
source(file.path(.script_dir, "true_model_factory.R"))  # make_true_model
source(file.path(.script_dir, "scm_bench_helpers.R"))   # base_2cmt_oral_linCmt,
                                                         # runSCM_traced,
                                                         # package_scm_result,
                                                         # diagnose_fit,
                                                         # extract_params_long,
                                                         # rel_err_one
source(file.path(.script_dir, "estimator_factory.R"))   # make_est_control,
                                                         # seed_for_dataset, seed_all

# ---- Module constants needed by the driver code (normally set in the -------
# ---- PerformanceEvaluation_scm_bench.R driver, not in the helpers) ---------
`%||%`           <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x
DOSE_MG          <- 100
scm_bench_vars   <- c("cl", "vc")
scm_bench_covars <- c("BW", "CrCL", "BMI")
scm_bench_cats   <- c("SEX", "RACE")
scm_bench_shapes <- c("power", "lin")

# ---- Cell coordinates ------------------------------------------------------
N        <- 80L
scen     <- 16L
ds       <- 1L
boundary <- "narrow"

sim_path         <- sprintf("Inputdataset/sim_obs_N%d/sim_obs_scenario_%02d.rds", N, scen)
true_params_path <- "Inputdataset/true_params_long.rds"

stopifnot(file.exists(sim_path), file.exists(true_params_path))

# ---- Capture helper: run an expression, keep result OR full error info -----
.capture <- function(expr) {
  tb <- NULL
  withCallingHandlers(
    tryCatch(
      list(ok = TRUE, value = force(expr), error = NULL, traceback = NULL),
      error = function(e) list(ok = FALSE, value = NULL,
                               error = conditionMessage(e),
                               traceback = tb)
    ),
    error = function(e) {
      tb <<- paste(utils::limitedLabels(sys.calls()), collapse = "\n")
    }
  )
}

# ---- Control builder: bypasses is_valid_combo() ----------------------------
# vae is NOT a make_est_control() branch (it would stop("Unknown estimator")),
# so we hand-build vaeControl() with the same tier semantics the driver uses:
#   screen -> covMethod = "" , calcTables = FALSE
#   final  -> covMethod = "r,s", calcTables = TRUE
make_ctrl <- function(est, opt, tier) {
  if (est == "vae") {
    cov <- if (tier == "final") "r,s" else ""
    return(nlmixr2est::vaeControl(covMethod = cov,
                                  calcTables = (tier == "final"),
                                  print = 0))
  }
  # irlsfoceif routes through the real factory (returns a foceiControl)
  make_est_control(est, outer_opt = opt, tier = tier)$ctrl
}

# ---- Data: two shapes -------------------------------------------------------
# SCM path keeps the pre-dose obs row (matches PerformanceEvaluation_scm_bench);
# Refit path drops EVID==0 & TIME==0 (matches bench_refit_estimators).
sim_all     <- readRDS(sim_path)
true_params <- readRDS(true_params_path)

ds_scm <- to_nm_dataset(sim_all) |>
  dplyr::filter(DATASET == ds) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))

ds_refit <- ds_scm |>
  dplyr::filter(!(EVID == 0L & TIME == 0))

true_mod    <- make_true_model(scen, boundary = boundary)
bounds_spec <- attr(true_mod, "bounds_spec")

# ---- (a) REFIT the true model ----------------------------------------------
run_refit <- function(est, opt) {
  seed <- seed_for_dataset(ds); seed_all(seed, est)     # 1001
  ctrl <- make_ctrl(est, opt, "final")

  cap <- .capture({
    t0  <- Sys.time()
    fit <- nlmixr2(true_mod, ds_refit, est = est, control = ctrl)
    list(
      fit     = fit,
      runtime = as.numeric(difftime(Sys.time(), t0, units = "secs")),
      diag_t3 = tryCatch(diagnose_fit_table3(fit, bounds_spec = bounds_spec),
                         error = function(e) conditionMessage(e)),
      diag    = tryCatch(diagnose_fit(fit),          error = function(e) conditionMessage(e)),
      est_long= tryCatch(extract_params_long(fit),   error = function(e) conditionMessage(e)),
      rel_err = tryCatch(rel_err_one(extract_params_long(fit), true_params, scen),
                         error = function(e) conditionMessage(e))
    )
  })
  cap$estimator <- est; cap$outer_opt <- opt; cap$mode <- "refit"
  cap
}

# ---- (b) SCM search --------------------------------------------------------
run_scm <- function(est, opt) {
  seed <- seed_for_dataset(ds); seed_all(seed, est)     # 1001
  screen <- make_ctrl(est, opt, "screen")
  final  <- make_ctrl(est, opt, "final")

  cap <- .capture({
    # Base fit (screening control), exactly as the driver does it
    t0   <- Sys.time()
    base <- nlmixr2(base_2cmt_oral_linCmt, ds_scm, est = est, control = screen)
    t_base <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # SCM search (workers = 1 for a quick local smoke; driver uses 3)
    scm <- runSCM_traced(
      label      = sprintf("SMOKE_N%d_scn%02d_%s_%s_ds%02d",
                           N, scen, est, opt %||% "NA", ds),
      data       = ds_scm,
      fit        = base,
      varsVec    = scm_bench_vars,
      covarsVec  = scm_bench_covars,
      catvarsVec = scm_bench_cats,
      shapes     = scm_bench_shapes,
      searchType = "scm",
      control    = screen,
      saveModels = FALSE,
      workers    = 1L,
      print      = 0,
      maxRetries = 0L,
      confirm    = FALSE
    )
    t_scm <- as.numeric(attr(scm, "elapsed_s"))

    # Post-SCM covariance refit strategy (vae -> in-fit cov, no refit)
    final_ctrl <- if (est == "vae") NULL else final
    packaged <- package_scm_result(
      label       = sprintf("SMOKE_%s_%s", est, opt %||% "NA"),
      scm_res     = scm,
      runtime_sec = t_scm,
      true_long   = true_params,
      scenario_id = scen,
      final_ctrl  = final_ctrl,
      final_est   = est
    )

    list(base_runtime = t_base, scm_runtime = t_scm,
         summaryTable = scm$summaryTable, packaged = packaged)
  })
  cap$estimator <- est; cap$outer_opt <- opt; cap$mode <- "scm"
  cap
}

# ---- Drive both estimators, both modes -------------------------------------
specs <- list(
  list(est = "irlsfoceif", opt = "lbfgsb3c"),   # practical optimizer for IRLS
  list(est = "vae",        opt = NA_character_)  # vae has no outer_opt
)

results <- list()
for (s in specs) {
  tag <- sprintf("%s_%s", s$est, s$opt %||% "NA")
  message(sprintf("\n================ SMOKE: %s ================", tag))

  message(sprintf("  [refit] %s ...", tag))
  results[[paste0(tag, "__refit")]] <- run_refit(s$est, s$opt)

  message(sprintf("  [scm]   %s ...", tag))
  results[[paste0(tag, "__scm")]]   <- run_scm(s$est, s$opt)
}

# ---- Compact status table for the console ----------------------------------
status_tbl <- dplyr::bind_rows(lapply(names(results), function(k) {
  r <- results[[k]]
  tibble::tibble(
    cell   = k,
    est    = r$estimator,
    opt    = r$outer_opt %||% "NA",
    mode   = r$mode,
    ok     = r$ok,
    error  = r$error %||% NA_character_
  )
}))

# ---- Save full bundle for the manager --------------------------------------
out_dir <- "output/smoke"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_rds <- file.path(out_dir, "scn16_ds01_N80_irlsfoceif_vae.rds")

bundle <- list(
  meta = list(
    N = N, scenario = scen, dataset = ds, boundary = boundary,
    specs = specs, timestamp = Sys.time(),
    sim_path = sim_path, true_params_path = true_params_path
  ),
  results     = results,
  status      = status_tbl,
  sessionInfo = utils::sessionInfo()
)
saveRDS(bundle, out_rds)

message("\n================ SMOKE SUMMARY ================")
print(status_tbl)
message(sprintf("\nFull capture bundle saved to: %s", out_rds))
