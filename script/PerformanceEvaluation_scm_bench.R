# ============================================================================
# PerformanceEvaluation_scm_bench.R
# ----------------------------------------------------------------------------
# Per-dataset SCM driver for the estimator x optimizer benchmark.
# Mirrors PerformanceEvaluation_sc_n_ds_HPCE.R but is parameterised over
#   --N, --scenario, --dataset, --estimator, --outer_opt
# and writes to a unified output tree:
#   output/scm_bench/N{NN}/scn{SS}/{est}_{opt}/res_ds{DDD}.rds
#
# Estimator handling:
#   * focei/foceif/irlsfoceif : screening + tight-tol refit both use the
#                                selected (est, outer_opt) from the factory.
#   * saem                    : screening uses saemControl; the tight-tol
#                                covariance refit uses foceif+lbfgsb3c
#                                (refit_estimator="foceif" logged).
#   * vae                     : screening uses vaeControl; final refit is
#                                skipped (VAE covariance is computed in-fit
#                                via vaeControl(covMethod="r,s")).  If the
#                                r,s attempt errors, retries once with "".
#
# CLI:
#   Rscript script/PerformanceEvaluation_scm_bench.R \
#       --N 80 --scenario 16 --dataset 1 \
#       --estimator foceif --outer_opt lbfgsb3c \
#       [--force_rerun] [--force_repackage]
#
# Input data layout (unchanged from PerformanceEvaluation_sc_n_ds_HPCE.R):
#   Inputdataset/sim_obs_N{NN}/sim_obs_scenario_{SS}.rds
#   Inputdataset/true_params_long.rds
# ============================================================================

suppressPackageStartupMessages({
  library(nlmixr2scm)
  library(nlmixr2)
  library(nlmixr2est)
  library(rxode2)
  library(dplyr)
  library(tibble)
})

# ---- Locate helper scripts (repo-root relative) ----------------------------
.script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg)) return(normalizePath(dirname(file_arg[1]), mustWork = FALSE))
  if (!is.null(sys.frames()) && length(sys.frames())) {
    fr <- sys.frame(1L)
    if (!is.null(fr$ofile)) return(normalizePath(dirname(fr$ofile), mustWork = FALSE))
  }
  "script"
})()

source(file.path(.script_dir, "scm_bench_helpers.R"))
source(file.path(.script_dir, "estimator_factory.R"))

# ---- CLI parser ------------------------------------------------------------
parse_args <- function(argv) {
  opts <- list(N = 80L, scenario = 16L, dataset = 1L,
               estimator = "focei", outer_opt = NA_character_,
               force_rerun = FALSE, force_repackage = FALSE,
               out_root = "output/scm_bench",
               input_root = "Inputdataset",
               true_params_path = "Inputdataset/true_params_long.rds")
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    val <- function() { i <<- i + 1L; argv[i] }
    switch(a,
      "--N"                = { opts$N        <- as.integer(val()) },
      "--scenario"         = { opts$scenario <- as.integer(val()) },
      "--dataset"          = { opts$dataset  <- as.integer(val()) },
      "--estimator"        = { opts$estimator<- val() },
      "--outer_opt"        = { v <- val(); opts$outer_opt <- if (v %in% c("NA","na","")) NA_character_ else v },
      "--force_rerun"      = { opts$force_rerun     <- TRUE },
      "--force_repackage"  = { opts$force_repackage <- TRUE },
      "--out_root"         = { opts$out_root <- val() },
      "--input_root"       = { opts$input_root <- val() },
      "--true_params"      = { opts$true_params_path <- val() },
      stop(sprintf("Unknown arg: %s", a))
    )
    i <- i + 1L
  }
  opts
}

# ---- Cell path helper ------------------------------------------------------
cell_dir <- function(out_root, N, scenario, estimator, outer_opt) {
  opt_tag <- if (is.na(outer_opt) || nchar(outer_opt) == 0L) "NA" else outer_opt
  file.path(out_root, sprintf("N%d", N),
            sprintf("scn%02d", scenario),
            paste(estimator, opt_tag, sep = "_"))
}

# ---- SCM covariate axes (scenario-16 superset used for all scenarios) -----
DOSE_MG          <- 100
scm_bench_vars   <- c("cl", "vc")
scm_bench_covars <- c("BW", "CrCL", "BMI")
scm_bench_cats   <- c("SEX", "RACE")
scm_bench_shapes <- c("power", "lin")

# ---- VAE screening wrapper (with graceful fallback) -----------------------
# Returns list(fit=..., covMethod_used=..., status=...)
fit_vae_screen <- function(ui, data, screen_ctrl) {
  status <- "vae_ok"
  fit <- tryCatch(
    nlmixr2(ui, data, est = "vae", control = screen_ctrl),
    error = function(e) {
      msg <- conditionMessage(e)
      warning("VAE fit with covMethod='r,s' failed: ", msg, call. = FALSE)
      NULL
    }
  )
  if (is.null(fit)) {
    # Retry with covMethod=""
    fallback_ctrl <- screen_ctrl
    fallback_ctrl$covMethod <- ""
    status <- "vae_fallback_no_cov"
    fit <- tryCatch(
      nlmixr2(ui, data, est = "vae", control = fallback_ctrl),
      error = function(e) {
        warning("VAE fallback (covMethod='') also failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(fit)) status <- "vae_failed"
  }
  list(fit = fit, status = status)
}

# ---- Per-cell driver -------------------------------------------------------
run_bench_cell <- function(opts) {
  # Validate combo
  if (!is_valid_combo(opts$estimator, opts$outer_opt)) {
    stop(sprintf("Invalid (est, outer_opt) combo: (%s, %s). See valid_combos().",
                 opts$estimator, opts$outer_opt))
  }

  # Paths
  save_dir  <- cell_dir(opts$out_root, opts$N, opts$scenario,
                        opts$estimator, opts$outer_opt)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  ds_tag    <- sprintf("ds%02d", opts$dataset)
  res_path  <- file.path(save_dir, sprintf("res_%s.rds", ds_tag))
  scm_path  <- file.path(save_dir, sprintf("scm_%s.rds", ds_tag))
  err_path  <- file.path(save_dir, sprintf("%s_ERROR.txt", ds_tag))

  # Tier 1: cached packaged result
  if (!opts$force_rerun && !opts$force_repackage && file.exists(res_path)) {
    message(sprintf(">>> [%s|%s_%s|N%d|scn%02d] tier-1 cached, skipping",
                    ds_tag, opts$estimator,
                    opts$outer_opt %||% "NA", opts$N, opts$scenario))
    return(invisible(readRDS(res_path)))
  }

  # Input data
  sim_path <- file.path(opts$input_root, sprintf("sim_obs_N%d", opts$N),
                       sprintf("sim_obs_scenario_%02d.rds", opts$scenario))
  if (!file.exists(sim_path)) stop("Sim RDS not found: ", sim_path)
  if (!file.exists(opts$true_params_path))
    stop("true_params RDS not found: ", opts$true_params_path)

  sim_all     <- readRDS(sim_path)
  true_params <- readRDS(opts$true_params_path)

  ds_i <- to_nm_dataset(sim_all) %>%
    dplyr::filter(DATASET == opts$dataset) %>%
    dplyr::select(-SCENARIO, -DATASET) %>%
    dplyr::mutate(ID = as.integer(ID),
                  SEX = as.integer(SEX),
                  RACE = as.integer(RACE))

  # Seed for reproducibility
  seed <- seed_for_dataset(opts$dataset)
  seed_all(seed, opts$estimator)

  # Build controls: screen (no cov, no tables) for base fit + every SCM
  # candidate LRT; final (r,s cov + tables) for the single post-SCM refit
  # (see package_scm_result()). Both tiers share the tuned sigdig/derivEps/
  # ODE tols per (est, outer_opt).
  screen_bundle <- make_est_control(opts$estimator, opts$outer_opt, "screen")
  final_bundle  <- make_est_control(opts$estimator, opts$outer_opt, "final")

  # Tier 2: cached SCM result -> repackage only
  if (!opts$force_rerun && file.exists(scm_path)) {
    message(sprintf(">>> [%s] tier-2 repackaging from %s", ds_tag, basename(scm_path)))
    return(invisible(tryCatch({
      scm_res     <- readRDS(scm_path)
      t_scm_sec   <- suppressWarnings(as.numeric(attr(scm_res, "elapsed_s")))
      final_ctrl_for_refit <- if (opts$estimator == "saem") {
        make_saem_refit_control()
      } else if (opts$estimator == "vae") {
        NULL  # VAE cov done in-fit; no separate refit
      } else {
        final_bundle$ctrl
      }
      refit_est <- if (opts$estimator == "saem") "foceif" else opts$estimator
      test_i <- package_scm_result(
        label       = sprintf("N%d_scn%02d_%s_%s_%s",
                              opts$N, opts$scenario, opts$estimator,
                              opts$outer_opt %||% "NA", ds_tag),
        scm_res     = scm_res,
        runtime_sec = t_scm_sec,
        true_long   = true_params,
        scenario_id = opts$scenario,
        final_ctrl  = final_ctrl_for_refit,
        final_est   = refit_est
      )
      out <- assemble_output(opts, test_i, seed, t_base_sec = NA_real_,
                             t_scm_sec = t_scm_sec, resumed = TRUE)
      saveRDS(out, res_path)
      if (file.exists(err_path)) unlink(err_path)
      out
    }, error = function(e) {
      msg <- conditionMessage(e)
      writeLines(c(format(Sys.time()), "tier-2 repackage failed:", msg), err_path)
      warning(sprintf("[%s] tier-2 FAILED: %s", ds_tag, msg))
      NULL
    })))
  }

  # Tier 3: full pipeline
  message(sprintf("\n>>> [%s|%s_%s|N%d|scn%02d] tier-3 START %s",
                  ds_tag, opts$estimator, opts$outer_opt %||% "NA",
                  opts$N, opts$scenario, format(Sys.time(), "%H:%M:%S")))

  tryCatch({
    # Base fit (uses screening control for the selected estimator)
    if (opts$estimator == "vae") {
      t_base <- system.time(vae_res <- fit_vae_screen(base_2cmt_oral_linCmt,
                                                      ds_i, screen_bundle$ctrl))
      fit_base_i <- vae_res$fit
      vae_status <- vae_res$status
      if (is.null(fit_base_i)) stop("VAE base fit failed after fallback: ", vae_status)
    } else {
      t_base <- system.time(
        fit_base_i <- nlmixr2(base_2cmt_oral_linCmt, ds_i,
                              est = opts$estimator, control = screen_bundle$ctrl)
      )
      vae_status <- NA_character_
    }
    t_base_sec <- as.numeric(t_base["elapsed"])

    # SCM screening
    scm_res <- runSCM_traced(
      label      = sprintf("N%d_scn%02d_%s_%s_%s",
                           opts$N, opts$scenario, opts$estimator,
                           opts$outer_opt %||% "NA", ds_tag),
      data       = ds_i,
      fit        = fit_base_i,
      varsVec    = scm_bench_vars,
      covarsVec  = scm_bench_covars,
      catvarsVec = scm_bench_cats,
      shapes     = scm_bench_shapes,
      # 2026-07-13: non-zero SCM candidate inits (A) --------------------
      # runSCM defaults every candidate slope to init=0. At theta=0 the
      # log-additive covariate contribution is identically zero, so the
      # finite-difference gradient of the outer OFV w.r.t. that theta is
      # ~0 modulo FOCEi inner-ETA noise. Gradient optimizers (nlminb,
      # lbfgsb3c) then quit with "false convergence (8)" before the LRT
      # can see any real improvement -> real covariates get rejected ->
      # Power collapses. Bobyqa is derivative-free and unaffected. Match
      # PsN convention + our bench_refit setting: 0.5 for continuous,
      # log(1.5) ~ 0.405 for categorical.
      inits      = list(
        power = list(est = 0.5,       lower = -5, upper = 5),
        lin   = list(est = 0.5,       lower = -5, upper = 5),
        cat   = list(est = log(1.5),  lower = -5, upper = 5)
      ),
      searchType = "scm",
      control    = screen_bundle$ctrl,
      saveModels = FALSE,
      workers    = 3L,
      print      = 0,
      maxRetries = 0L,
      confirm    = FALSE
    )
    t_scm_sec <- as.numeric(attr(scm_res, "elapsed_s"))
    saveRDS(scm_res, scm_path)

    # Tight-tol covariance refit strategy
    final_ctrl_for_refit <- if (opts$estimator == "saem") {
      make_saem_refit_control()
    } else if (opts$estimator == "vae") {
      NULL
    } else {
      final_bundle$ctrl
    }
    refit_est <- if (opts$estimator == "saem") "foceif" else opts$estimator

    test_i <- package_scm_result(
      label       = sprintf("N%d_scn%02d_%s_%s_%s",
                            opts$N, opts$scenario, opts$estimator,
                            opts$outer_opt %||% "NA", ds_tag),
      scm_res     = scm_res,
      runtime_sec = t_scm_sec,
      true_long   = true_params,
      scenario_id = opts$scenario,
      final_ctrl  = final_ctrl_for_refit,
      final_est   = refit_est
    )
    test_i$vae_status <- vae_status

    out <- assemble_output(opts, test_i, seed, t_base_sec, t_scm_sec,
                           resumed = FALSE)
    saveRDS(out, res_path)
    if (file.exists(err_path)) unlink(err_path)

    message(sprintf("<<< [%s] DONE base=%.1fs scm=%.1fs refit=%.1fs cov=%s",
                    ds_tag, t_base_sec, t_scm_sec,
                    test_i$refit_runtime_sec %||% NA,
                    if (isTRUE(test_i$cov_done)) "ok" else "FAIL"))
    invisible(out)
  }, error = function(e) {
    msg <- conditionMessage(e)
    writeLines(c(format(Sys.time()), msg), err_path)
    warning(sprintf("[%s] FAILED: %s (see %s)", ds_tag, msg, err_path),
            call. = FALSE, immediate. = TRUE)
    invisible(NULL)
  })
}

# ---- Output assembly (uniform 5-key schema) --------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

assemble_output <- function(opts, test_i, seed, t_base_sec, t_scm_sec, resumed) {
  list(
    # 5 coordinate keys
    sample_N    = opts$N,
    scenario_id = opts$scenario,
    dataset_id  = opts$dataset,
    estimator   = opts$estimator,
    outer_opt   = opts$outer_opt,
    # Provenance
    seed        = seed,
    resumed     = resumed,
    timestamp   = Sys.time(),
    # Timings
    t_base_sec        = t_base_sec,
    t_scm_sec         = t_scm_sec,
    t_refit_sec       = test_i$refit_runtime_sec,
    t_total_sec       = sum(c(t_base_sec, t_scm_sec, test_i$refit_runtime_sec),
                            na.rm = TRUE),
    # Packaged SCM result
    test        = test_i
  )
}

# ---- CLI entry -------------------------------------------------------------
if (!interactive()) {
  argv <- commandArgs(trailingOnly = TRUE)
  opts <- parse_args(argv)
  message(sprintf("[scm_bench] N=%d scn=%d ds=%d est=%s opt=%s",
                  opts$N, opts$scenario, opts$dataset,
                  opts$estimator, opts$outer_opt %||% "NA"))
  t_run <- system.time(res <- run_bench_cell(opts))
  message(sprintf("[scm_bench] finished in %.1f min",
                  t_run["elapsed"] / 60))
}
