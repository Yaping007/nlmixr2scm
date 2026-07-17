# ============================================================================
# PerformanceEvaluation_scm_bench.R
# ----------------------------------------------------------------------------
# Per-dataset SCM driver for the estimator x optimizer x structure benchmark.
# Parameterised over
#   --N, --scenario, --dataset, --estimator, --outer_opt, --structure
# and writes SCHEMA-2.1 records (via package_scm_schema21) to:
#   output/scm_bench/N{NN}/scn{SS}_{structure}/{est}_{opt}/res_ds{DDD}.rds
#   (+ .fit.rds refit winner, + .meta.json manifest)
#
# Model structure (--structure):
#   * linCmt : analytic 2-cmt oral (base_2cmt_oral_linCmt). foceif/irlsfoceif
#              degrade to focei here (no ODE to interaction-linearise) but are
#              still run to demonstrate 'prefer linCmt when possible'.
#   * ode    : explicit 2-cmt oral ODE (base_2cmt_oral_ode); unlocks analytic
#              gradients for foceif/irlsfoceif.
#
# Estimator handling (runSCM-based two-tier: screen search + tight-tol refit):
#   * focei/foceif/irlsfoceif : screening + tight-tol covariance refit both use
#                                the selected (est, outer_opt) from the factory.
#   * vae                     : screening uses vaeControl; final refit is
#                                skipped (VAE covariance is computed in-fit
#                                via vaeControl(covMethod="r,s")).  If the
#                                r,s attempt errors, retries once with "".
#
# CLI:
#   Rscript script/PerformanceEvaluation_scm_bench.R \
#       --N 80 --scenario 16 --dataset 1 \
#       --estimator foceif --outer_opt lbfgsb3c --structure ode \
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
source(file.path(.script_dir, "refit_helpers.R"))    # diagnose_fit_table3, %||%
source(file.path(.script_dir, "output_schema.R"))    # assemble_common, write_fit_sidecar, SCHEMA_VERSION

# ---- CLI parser ------------------------------------------------------------
parse_args <- function(argv) {
  opts <- list(N = 80L, scenario = 16L, dataset = 1L,
               estimator = "focei", outer_opt = NA_character_,
               structure = "linCmt",
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
      "--structure"        = { opts$structure <- val() },
      "--force_rerun"      = { opts$force_rerun     <- TRUE },
      "--force_repackage"  = { opts$force_repackage <- TRUE },
      "--out_root"         = { opts$out_root <- val() },
      "--input_root"       = { opts$input_root <- val() },
      "--true_params"      = { opts$true_params_path <- val() },
      stop(sprintf("Unknown arg: %s", a))
    )
    i <- i + 1L
  }
  if (!opts$structure %in% c("linCmt", "ode"))
    stop(sprintf("--structure must be 'linCmt' or 'ode', got '%s'", opts$structure))
  opts
}

# ---- Cell path helper ------------------------------------------------------
cell_dir <- function(out_root, N, scenario, structure, estimator, outer_opt) {
  opt_tag <- if (is.na(outer_opt) || nchar(outer_opt) == 0L) "NA" else outer_opt
  file.path(out_root, sprintf("N%d", N),
            sprintf("scn%02d_%s", scenario, structure),
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

  # Parallelism policy (fair + reproducible):
  #   runSCM forks `workers` child processes for SCM candidate LRTs. rxode2's
  #   ODE solver uses OpenMP. fork() over a live OpenMP thread pool is UNDEFINED
  #   behaviour -- it silently corrupts the ODE solve in the children, producing
  #   non-deterministic, wrong fits (frozen-at-init OFVs, junk LRT deltas).
  #   Therefore we PIN rxThreads=1: OpenMP off, so the SCM fork is clean and
  #   every cell is bitwise reproducible. Parallelism comes from the fork
  #   (workers) and, across datasets, the LSF array width.
  scm_workers <- 3L
  rx_threads  <- 1L
  rxode2::setRxThreads(rx_threads)
  Sys.setenv(OMP_NUM_THREADS = "1")   # also quiet BLAS/OpenMP in child procs

  # Paths
  save_dir  <- cell_dir(opts$out_root, opts$N, opts$scenario,
                        opts$structure, opts$estimator, opts$outer_opt)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  ds_tag    <- sprintf("ds%03d", opts$dataset)
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
  # candidate LRT; final (r,s cov + tables) for the single post-SCM tight-tol
  # covariance refit (two-tier runSCM approach). Both tiers share the tuned
  # sigdig/derivEps/ODE tols per (est, outer_opt).
  screen_bundle <- make_est_control(opts$estimator, opts$outer_opt, "screen")
  final_bundle  <- make_est_control(opts$estimator, opts$outer_opt, "final")

  # Base (covariate-free) model matching the requested structure; carries
  # attr(,"structure") so assemble_common()->fit_model_type() sets model_type.
  base_mod <- switch(opts$structure,
    "linCmt" = base_2cmt_oral_linCmt,
    "ode"    = base_2cmt_oral_ode,
    stop("Unknown structure: ", opts$structure))

  # Tier 2: cached SCM result -> repackage only
  if (!opts$force_rerun && file.exists(scm_path)) {
    message(sprintf(">>> [%s] tier-2 repackaging from %s", ds_tag, basename(scm_path)))
    return(invisible(tryCatch({
      scm_res     <- readRDS(scm_path)
      t_scm_sec   <- suppressWarnings(as.numeric(attr(scm_res, "elapsed_s")))
      # VAE covariance is done in-fit (no separate refit); classical estimators
      # get the tight-tol final refit inside package_scm_schema21.
      refit_ctrl_for_cell <- if (opts$estimator == "vae") NULL else final_bundle$ctrl
      package_scm_schema21(
        scm_res         = scm_res,
        true_mod        = base_mod,
        scenario_id     = opts$scenario,
        true_params     = true_params,
        estimator       = opts$estimator,
        outer_opt       = opts$outer_opt,
        refit_ctrl      = refit_ctrl_for_cell,
        refit_estimator = opts$estimator,
        runtimes        = list(base_sec = NA_real_, scm_sec = t_scm_sec),
        identity        = list(sample_N = opts$N, dataset_id = opts$dataset,
                               scm_workers = scm_workers, rx_threads = rx_threads),
        structure       = opts$structure,
        out_rds         = res_path
      )
      if (file.exists(err_path)) unlink(err_path)
      readRDS(res_path)
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
      t_base <- system.time(vae_res <- fit_vae_screen(base_mod,
                                                      ds_i, screen_bundle$ctrl))
      fit_base_i <- vae_res$fit
      vae_status <- vae_res$status
      if (is.null(fit_base_i)) stop("VAE base fit failed after fallback: ", vae_status)
    } else {
      t_base <- system.time(
        fit_base_i <- nlmixr2(base_mod, ds_i,
                              est = opts$estimator, control = screen_bundle$ctrl)
      )
      vae_status <- NA_character_
    }
    t_base_sec <- as.numeric(t_base["elapsed"])
    # CPU-seconds for the base fit (self + forked-child); cpu/wall = speedup.
    t_base_cpu <- unname(sum(t_base[c("user.self", "sys.self",
                                      "user.child", "sys.child")], na.rm = TRUE))

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
      # 2026-07-13: reverted non-zero SCM candidate inits (A). Non-zero
      # inits shift the LRT null: H0(theta=0) vs H1(theta=0.5) creates a
      # large structural OFV drop at the init itself, so LRT selects
      # every candidate without the optimizer moving. Keep runSCM's
      # default init=0 so LRT reflects genuine slope evidence.
      searchType = "scm",
      control    = screen_bundle$ctrl,
      saveModels = FALSE,
      workers    = 3L,
      print      = 0,
      maxRetries = 0L,
      confirm    = FALSE
    )
    t_scm_sec <- as.numeric(attr(scm_res, "elapsed_s"))
    t_scm_cpu <- suppressWarnings(as.numeric(attr(scm_res, "cpu_s")))
    saveRDS(scm_res, scm_path)

    # Tight-tol covariance refit strategy: VAE cov is in-fit (skip refit);
    # classical estimators refit with the cell's own final settings inside
    # package_scm_schema21.
    refit_ctrl_for_cell <- if (opts$estimator == "vae") NULL else final_bundle$ctrl

    rec <- package_scm_schema21(
      scm_res         = scm_res,
      true_mod        = base_mod,
      scenario_id     = opts$scenario,
      true_params     = true_params,
      estimator       = opts$estimator,
      outer_opt       = opts$outer_opt,
      refit_ctrl      = refit_ctrl_for_cell,
      refit_estimator = opts$estimator,
      runtimes        = list(base_sec = t_base_sec, scm_sec = t_scm_sec,
                             base_cpu = t_base_cpu, scm_cpu = t_scm_cpu),
      identity        = list(sample_N = opts$N, dataset_id = opts$dataset,
                             scm_workers = scm_workers, rx_threads = rx_threads),
      structure       = opts$structure,
      out_rds         = res_path
    )
    if (file.exists(err_path)) unlink(err_path)

    message(sprintf("<<< [%s] DONE base=%.1fs scm=%.1fs refit=%.1fs cov=%s | wall=%.1fs cpu=%.1fs hog=%.2f",
                    ds_tag, t_base_sec, t_scm_sec,
                    rec$runtime$refit_sec %||% NA,
                    if (isTRUE(rec$scm$cov_done)) "ok" else "FAIL",
                    rec$runtime$total_sec %||% NA,
                    rec$cpu$total_sec %||% NA,
                    rec$cpu$hog_factor %||% NA))
    invisible(rec)
  }, error = function(e) {
    msg <- conditionMessage(e)
    writeLines(c(format(Sys.time()), msg), err_path)
    warning(sprintf("[%s] FAILED: %s (see %s)", ds_tag, msg, err_path),
            call. = FALSE, immediate. = TRUE)
    invisible(NULL)
  })
}

# ---- fallback operator (also defined in refit_helpers.R; kept for safety) --
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

# ---- CLI entry -------------------------------------------------------------
if (!interactive()) {
  argv <- commandArgs(trailingOnly = TRUE)
  opts <- parse_args(argv)
  message(sprintf("[scm_bench] N=%d scn=%d ds=%d est=%s opt=%s struct=%s",
                  opts$N, opts$scenario, opts$dataset,
                  opts$estimator, opts$outer_opt %||% "NA", opts$structure))
  t_run <- system.time(res <- run_bench_cell(opts))
  message(sprintf("[scm_bench] finished in %.1f min",
                  t_run["elapsed"] / 60))
}
