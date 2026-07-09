# ==============================================================================
# refit_one_dataset.R
# ------------------------------------------------------------------------------
# CLI-invokable driver that fits the "true" model to ONE (cohort, scenario,
# boundary, dataset) cell and saves a slim result RDS.  Designed to be called
# once per task in a HPCE LSF job array (10 pilot / 250 full tasks per array).
#
# Usage (CLI):
#   Rscript scripts/refit_one_dataset.R \
#     --cohort N300 --scenario 16 --boundary narrow --dataset 1 \
#     [--master_rds <override>] [--out_dir <override>] [--overwrite TRUE]
#
# Usage (interactive / testing):
#   source("scripts/refit_one_dataset.R")
#   refit_one_dataset(scenario_id = 16, boundary = "narrow",
#                     dataset_id = 1, cohort = "N300")
#
# Cohort -> defaults:
#   N40  -> simulated_virtual_dataset_eta_filtered_N40/sim_obs_all_scenarios.rds
#            outputs/refit_true_N40/scn<NN>_<bnd>/res_ds<DDD>.rds
#   N80  -> simulated_virtual_dataset_eta_filtered_N80/sim_obs_all_scenarios.rds
#            outputs/refit_true_N80/scn<NN>_<bnd>/res_ds<DDD>.rds
#   N300 -> simulated_virtual_dataset_eta_filtered/sim_obs_all_scenarios.rds
#            outputs/refit_true_N300/scn<NN>_<bnd>/res_ds<DDD>.rds
#
# Output:
#   <out_dir>/res_ds<DDD>.rds  -- slim named list, see refit_one_dataset() docs.
#   <out_dir>/res_ds<DDD>_ERROR.txt  -- written on tryCatch failure only.
#
# Cache behaviour: skips silently if res_ds<DDD>.rds already exists.
# ==============================================================================

## ---- Repo-root anchor + helpers ------------------------------------------
.this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  }
  ## Interactive fallback: prefer ./scripts/ if it holds the helpers,
  ## else current working dir.
  cand <- file.path(getwd(), "scripts", "refit_helpers.R")
  if (file.exists(cand)) return(dirname(cand))
  getwd()
}
.script_dir <- .this_script_dir()
source(file.path(.script_dir, "refit_helpers.R"), chdir = FALSE)
source(file.path(.script_dir, "true_model_factory.R"), chdir = FALSE)

suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
})

## ---- FOCEi control (tight tolerance, cov step enabled) --------------------
##   Matches scripts/PerformanceEvaluation04062026.R:1696 (`scm_focei_final`).
##   sigdig = 4 + covMethod = "r,s" -> tight-tol stationary point + a well-
##   conditioned R-matrix so parFixedDf and cond_num_cor populate.
##   NOTE: atol/rtol are ODE-solver tolerances, passed via rxControl.
refit_focei_control <- function() {
  nlmixr2est::foceiControl(
    sigdig             = 4,
    outerOpt           = "bobyqa",
    print              = 0,
    calcTables         = FALSE,      # no IPRED/CWRES tables needed here
    covMethod          = "r,s",
    stickyRecalcN      = 20,
    maxOuterIterations = 2000,
    maxInnerIterations = 2000,
    rxControl          = rxode2::rxControl(atol = 1e-8, rtol = 1e-6)
  )
}

## ---- Cohort -> (master RDS, out_dir) defaults ----------------------------
##   The three sim cohorts share a canonical directory naming convention.
##   Note: N=300 lives in the un-suffixed *_eta_filtered/ dir (historical).
.cohort_defaults <- function(cohort) {
  stopifnot(cohort %in% c("N40", "N80", "N300"))
  master_dir <- switch(cohort,
    N40  = "simulated_virtual_dataset_eta_filtered_N40",
    N80  = "simulated_virtual_dataset_eta_filtered_N80",
    N300 = "simulated_virtual_dataset_eta_filtered"
  )
  list(
    master_rds = file.path(master_dir, "sim_obs_all_scenarios.rds"),
    out_root   = file.path("outputs", paste0("refit_true_", cohort))
  )
}

## ---- Main driver ----------------------------------------------------------
refit_one_dataset <- function(scenario_id,
                              boundary,
                              dataset_id,
                              cohort     = "N300",
                              master_rds = NULL,
                              out_dir    = NULL,
                              overwrite  = FALSE) {

  stopifnot(is.numeric(scenario_id), length(scenario_id) == 1L,
            is.character(boundary),  length(boundary) == 1L,
            boundary %in% c("none", "wide", "narrow", "tight"),
            is.numeric(dataset_id),  length(dataset_id) == 1L,
            is.character(cohort),    length(cohort) == 1L,
            cohort %in% c("N40", "N80", "N300"))

  defaults <- .cohort_defaults(cohort)
  if (is.null(master_rds)) master_rds <- defaults$master_rds
  if (is.null(out_dir)) {
    out_dir <- file.path(defaults$out_root,
                         sprintf("scn%02d_%s", scenario_id, boundary))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ds_tag   <- sprintf("ds%03d", dataset_id)
  out_rds  <- file.path(out_dir, sprintf("res_%s.rds", ds_tag))
  out_err  <- file.path(out_dir, sprintf("res_%s_ERROR.txt", ds_tag))

  if (file.exists(out_rds) && !overwrite) {
    message(sprintf("[skip] cached %s", out_rds))
    return(invisible(readRDS(out_rds)))
  }

  cat(sprintf(
    "[refit] cohort=%-4s scenario=%02d boundary=%-6s dataset=%03d  ->  %s\n",
    cohort, scenario_id, boundary, dataset_id, out_rds
  ))

  ## ---- Wrap the whole fit in tryCatch so a single-cell failure doesn't
  ##      take down the job array.
  result <- tryCatch({
    ## 1. Load + slice the master RDS to one (scenario, dataset)
    sim_slice <- load_scenario_dataset(master_rds, scenario_id, dataset_id)
    if (nrow(sim_slice) == 0L) {
      stop("Empty slice for scenario=", scenario_id, ", dataset=", dataset_id)
    }
    ## Mirror the NM-dataset post-processing used in the original smoke tests
    ## (scripts/PerformanceEvaluation04062026.R:981):
    ##   - drop SCENARIO / DATASET (nlmixr2 chokes on unknown columns for cov)
    ##   - coerce ID, SEX, RACE to integer (needed for categorical covariates)
    ##   - drop the pre-dose observation at TIME == 0 (EVID = 0, cp_obs = 0)
    ##     which breaks the proportional-only error likelihood
    ##     (SD = prop.err * |cp| = 0 when cp = 0 causes log(0) at that row).
    ##     The EVID = 1 dose row at TIME = 0 is kept.
    ds_nm <- to_nm_dataset(sim_slice) |>
      dplyr::filter(!(EVID == 0L & TIME == 0)) |>
      dplyr::select(-SCENARIO, -DATASET) |>
      dplyr::mutate(
        ID   = as.integer(ID),
        SEX  = as.integer(SEX),
        RACE = as.integer(RACE)
      )

    ## 2. Build the scenario-and-boundary-aware true model
    true_mod <- make_true_model(scenario_id, boundary = boundary)

    ## 3. Fit
    t0  <- Sys.time()
    fit <- nlmixr2(true_mod, ds_nm,
                   est     = "focei",
                   control = refit_focei_control())
    runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    ## 4. Extract Table-3 diagnostics + estimates + relative errors
    bounds_spec <- attr(true_mod, "bounds_spec")
    diag_t3     <- diagnose_fit_table3(fit, bounds_spec = bounds_spec)
    final_est   <- extract_params_long(fit)
    rel_err     <- rel_err_one(final_est, true_params, scenario_id)

    ## parFixedDf (only populated when covMethod completed).
    parFixed <- tryCatch(as.data.frame(fit$parFixedDf),
                         error = function(e) NULL)

    list(
      cohort       = cohort,
      scenario_id  = as.integer(scenario_id),
      dataset_id   = as.integer(dataset_id),
      boundary     = boundary,
      final_est    = final_est,
      rel_err      = rel_err,
      diag         = diag_t3,
      parFixed     = parFixed,
      bounds_spec  = bounds_spec,
      cont_params  = attr(true_mod, "estimated_cont_params"),
      cat_params   = attr(true_mod, "categorical_params"),
      fn_text      = attr(true_mod, "fn_text"),
      runtime_sec  = runtime_sec,
      timestamp    = Sys.time(),
      status       = "ok"
    )
  }, error = function(e) {
    msg <- sprintf(
      "[FAIL] cohort=%s scenario=%02d boundary=%s dataset=%03d\nerror: %s\ncalltrace:\n%s\n",
      cohort, scenario_id, boundary, dataset_id, conditionMessage(e),
      paste(deparse(sys.calls()), collapse = "\n")
    )
    cat(msg, file = out_err)
    list(
      cohort      = cohort,
      scenario_id = as.integer(scenario_id),
      dataset_id  = as.integer(dataset_id),
      boundary    = boundary,
      status      = "error",
      error_msg   = conditionMessage(e),
      timestamp   = Sys.time()
    )
  })

  saveRDS(result, out_rds)
  invisible(result)
}

## ---- CLI dispatcher -------------------------------------------------------
##   Runs only when the script is invoked non-interactively AND
##   `commandArgs(trailingOnly = TRUE)` is non-empty (i.e. an actual Rscript
##   call with flags).  Sourcing from an interactive session or another
##   script leaves the dispatcher inert so `refit_one_dataset()` can be
##   called programmatically.
.parse_cli <- function(argv) {
  # Very small --key value parser (no third-party deps).
  out <- list()
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      # Support both "--k v" and "--k=v" forms
      if (grepl("=", key, fixed = TRUE)) {
        parts <- strsplit(key, "=", fixed = TRUE)[[1L]]
        out[[parts[1L]]] <- parts[2L]
        i <- i + 1L
      } else {
        if (i + 1L > length(argv) || startsWith(argv[i + 1L], "--")) {
          out[[key]] <- TRUE
          i <- i + 1L
        } else {
          out[[key]] <- argv[i + 1L]
          i <- i + 2L
        }
      }
    } else {
      i <- i + 1L
    }
  }
  out
}

.run_cli <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  if (length(argv) == 0L) return(invisible(NULL))
  args <- .parse_cli(argv)

  required <- c("scenario", "boundary", "dataset")
  missing_flags <- setdiff(required, names(args))
  if (length(missing_flags) > 0L) {
    stop("Missing required flag(s): ", paste(missing_flags, collapse = ", "),
         "\nUsage: Rscript refit_one_dataset.R --scenario <int> --boundary <none|wide|narrow|tight> --dataset <int> [--cohort <N40|N80|N300>] [--master_rds <path>] [--out_dir <path>]")
  }

  refit_one_dataset(
    scenario_id = as.integer(args$scenario),
    boundary    = as.character(args$boundary),
    dataset_id  = as.integer(args$dataset),
    cohort      = as.character(args$cohort %||% "N300"),
    master_rds  = args$master_rds,     # NULL -> derive from cohort
    out_dir     = args$out_dir,        # NULL -> derive from cohort
    overwrite   = isTRUE(as.logical(args$overwrite))
  )
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0L) {
  .run_cli()
}
