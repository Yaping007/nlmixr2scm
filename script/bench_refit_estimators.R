# ==============================================================================
# bench_refit_estimators.R
# ------------------------------------------------------------------------------
# CLI driver that fits the TRUE model (no SCM search) to ONE
# (N, scenario, dataset, estimator, outer_opt) cell and saves a slim RDS.
#
# Purpose: isolate estimator/optimizer behaviour from the SCM search itself.
# Used to tune SAEM (nBurn/nEm), IRLS-FOCEI (cov step, hess type), and
# lbfgsb3c-based fits, before pushing tuned controls back into the SCM bench.
#
# Layout mirrors refit_one_dataset.R but adds estimator/outer_opt keys.
#
# Usage (CLI):
#   Rscript script/bench_refit_estimators.R \
#     --N 80 --scenario 16 --dataset 1 \
#     --estimator foceif --outer_opt lbfgsb3c \
#     [--boundary narrow] [--out_dir <path>] [--overwrite TRUE]
#
# SAEM has no outer_opt; pass --outer_opt NA (or omit).  The driver honours
# needs_post_refit from make_est_control() and refits with foceif for cov.
#
# Output:
#   <out_dir>/res_ds<DDD>.rds
#   <out_dir>/res_ds<DDD>_ERROR.txt    (on tryCatch failure)
#   default out_dir: output/bench_refit_estimators/N<NN>/scn<SS>/<est>_<opt>/
#
# Cache: skip silently if res_ds<DDD>.rds exists (unless --overwrite TRUE).
# ==============================================================================

## ---- Repo-root anchor -----------------------------------------------------
.this_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa) == 1L) return(dirname(normalizePath(sub("^--file=", "", fa))))
  if (!is.null(sys.frame(1)$ofile))
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  cand <- file.path(getwd(), "script", "refit_helpers.R")
  if (file.exists(cand)) return(dirname(cand))
  getwd()
}
.script_dir <- .this_script_dir()
source(file.path(.script_dir, "refit_helpers.R"),      chdir = FALSE)
source(file.path(.script_dir, "true_model_factory.R"), chdir = FALSE)
source(file.path(.script_dir, "estimator_factory.R"),  chdir = FALSE)
source(file.path(.script_dir, "output_schema.R"),      chdir = FALSE)

suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

## ---- Resolve per-scenario simulated RDS -----------------------------------
## New layout (2026-07):
##   Inputdataset/sim_obs_N{NN}/sim_obs_scenario_{SS}.rds
## Each RDS holds all datasets for ONE scenario (SCENARIO/DATASET cols kept).
.resolve_sim_rds <- function(N, scenario_id,
                             input_root = "Inputdataset") {
  N <- as.integer(N)
  if (!N %in% c(40L, 80L, 300L)) stop("Unsupported N=", N)
  file.path(input_root,
            sprintf("sim_obs_N%d", N),
            sprintf("sim_obs_scenario_%02d.rds", as.integer(scenario_id)))
}

.cell_dir <- function(estimator, outer_opt) {
  opt <- if (is.na(outer_opt) || is.null(outer_opt)) "NA" else outer_opt
  paste(estimator, opt, sep = "_")
}

## ---- Main driver ----------------------------------------------------------
bench_refit_one <- function(N, scenario_id, dataset_id,
                            estimator, outer_opt = NA_character_,
                            boundary   = "narrow",
                            master_rds = NULL,
                            out_dir    = NULL,
                            overwrite  = FALSE) {

  stopifnot(is_valid_combo(estimator, outer_opt))
  stopifnot(boundary %in% c("none", "wide", "narrow", "tight"))

  N <- as.integer(N)
  if (is.null(master_rds)) master_rds <- .resolve_sim_rds(N, scenario_id)
  if (is.null(out_dir)) {
    out_dir <- file.path("output", "bench_refit_estimators",
                         sprintf("N%d",  N),
                         sprintf("scn%02d", scenario_id),
                         .cell_dir(estimator, outer_opt))
  }
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ds_tag  <- sprintf("ds%03d", dataset_id)
  out_rds <- file.path(out_dir, sprintf("res_%s.rds", ds_tag))
  out_err <- file.path(out_dir, sprintf("res_%s_ERROR.txt", ds_tag))

  if (file.exists(out_rds) && !overwrite) {
    message(sprintf("[skip] cached %s", out_rds))
    return(invisible(readRDS(out_rds)))
  }

  cat(sprintf(
    "[bench-refit] N=%d scn=%02d ds=%03d %s/%s boundary=%s -> %s\n",
    N, scenario_id, dataset_id, estimator,
    outer_opt %||% "NA", boundary, out_rds
  ))

  result <- tryCatch({
    ## 1. Load slice + NM-shape it
    sim_slice <- load_scenario_dataset(master_rds, scenario_id, dataset_id)
    if (nrow(sim_slice) == 0L)
      stop("Empty slice for scenario=", scenario_id, ", dataset=", dataset_id)

    ds_nm <- to_nm_dataset(sim_slice) |>
      dplyr::filter(!(EVID == 0L & TIME == 0)) |>
      dplyr::select(-SCENARIO, -DATASET) |>
      dplyr::mutate(
        ID   = as.integer(ID),
        SEX  = as.integer(SEX),
        RACE = as.integer(RACE)
      )

    ## 2. True model (ODE is canonical; see output_schema.R model_type tag)
    true_mod    <- make_true_model(scenario_id, boundary = boundary,
                                   structure = "ode")
    bounds_spec <- attr(true_mod, "bounds_spec")

    ## 3. Fit with the chosen (est, outer_opt) at final tier
    cfg <- make_est_control(estimator, outer_opt = outer_opt, tier = "final")

    ## Grid label maps to nlmixr2's est= dispatch string via nlmixr_est_name
    ## (foceif/irlsfoceif -> "focei"; the interaction/muModel behaviour is
    ## carried by the control from make_est_control, see estimator_factory.R).
    est_dispatch <- nlmixr_est_name(estimator)

    t0  <- Sys.time()
    fit <- nlmixr2(true_mod, ds_nm, est = est_dispatch, control = cfg$ctrl)
    t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    ## 4. Optional post-refit for SAEM (cov step delegated to focei+bobyqa)
    refit_estimator <- estimator
    t_refit <- NA_real_
    if (isTRUE(cfg$needs_post_refit)) {
      cat("[bench-refit]  ... post-SAEM focei cov refit (bobyqa)\n")
      t1 <- Sys.time()
      refit_ctrl <- make_saem_refit_control()
      fit_refit  <- tryCatch(
        ## est="focei" (FOCEi + interaction, bobyqa, FD) for the post-SAEM
        ## covariance refit. NOT foceif: under bobyqa the analytic gradient
        ## (foceif's only edge) downgrades to FD, so foceif==focei here.
        nlmixr2(true_mod, ds_nm, est = nlmixr_est_name("focei"),
                control = refit_ctrl),
        error = function(e) { message("post-refit failed: ", conditionMessage(e)); NULL }
      )
      t_refit <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
      if (!is.null(fit_refit)) {
        fit <- fit_refit
        refit_estimator <- "focei"
      }
    }

    ## 5. Diagnostics + estimates + rel err via the shared schema assembler.
    ##    assemble_common() stamps schema_version + model_type and builds the
    ##    estimator-agnostic core (objf, converged, cond_num_cor, diag,
    ##    diag_t3, rel_err, parFixed, cov, *_params, fn_text) and the canonical
    ##    fit_runtime_sec. Driver-specific identity keys and the refit runtime
    ##    are layered on top here.
    common <- assemble_common(
      fit         = fit,
      true_mod    = true_mod,
      scenario_id = scenario_id,
      true_params = true_params,
      runtime_sec = t_fit,
      status      = "ok",
      estimator   = refit_estimator
    )

    c(list(
      sample_N          = N,
      scenario_id       = as.integer(scenario_id),
      dataset_id        = as.integer(dataset_id),
      estimator         = estimator,
      outer_opt         = if (is.null(outer_opt)) NA_character_ else as.character(outer_opt),
      refit_estimator   = refit_estimator,
      boundary          = boundary,
      refit_runtime_sec = t_refit,
      timestamp         = Sys.time()
    ), common)
  }, error = function(e) {
    msg <- sprintf(
      "[FAIL] N=%d scn=%02d ds=%03d %s/%s boundary=%s\nerror: %s\ncalls:\n%s\n",
      N, scenario_id, dataset_id, estimator, outer_opt %||% "NA", boundary,
      conditionMessage(e),
      paste(deparse(sys.calls()), collapse = "\n"))
    cat(msg, file = out_err)
    list(
      schema_version = SCHEMA_VERSION,
      model_type  = "ode",
      sample_N    = N,
      scenario_id = as.integer(scenario_id),
      dataset_id  = as.integer(dataset_id),
      estimator   = estimator,
      outer_opt   = if (is.null(outer_opt)) NA_character_ else as.character(outer_opt),
      boundary    = boundary,
      status      = "error",
      error_msg   = conditionMessage(e),
      timestamp   = Sys.time()
    )
  })

  write_fit_sidecar(result, out_rds,
                    fit = if (identical(result$status, "ok") &&
                              exists("fit", inherits = FALSE)) fit else NULL)
  invisible(result)
}

## ---- CLI dispatcher -------------------------------------------------------
.parse_cli <- function(argv) {
  out <- list(); i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]
    if (startsWith(a, "--")) {
      key <- sub("^--", "", a)
      if (grepl("=", key, fixed = TRUE)) {
        parts <- strsplit(key, "=", fixed = TRUE)[[1L]]
        out[[parts[1L]]] <- parts[2L]; i <- i + 1L
      } else {
        if (i + 1L > length(argv) || startsWith(argv[i + 1L], "--")) {
          out[[key]] <- TRUE; i <- i + 1L
        } else {
          out[[key]] <- argv[i + 1L]; i <- i + 2L
        }
      }
    } else i <- i + 1L
  }
  out
}

.run_cli <- function() {
  argv <- commandArgs(trailingOnly = TRUE)
  if (length(argv) == 0L) return(invisible(NULL))
  a <- .parse_cli(argv)

  required <- c("N", "scenario", "dataset", "estimator")
  miss <- setdiff(required, names(a))
  if (length(miss))
    stop("Missing required flag(s): ", paste(miss, collapse = ", "),
         "\nUsage: Rscript script/bench_refit_estimators.R --N <int> --scenario <int> --dataset <int> --estimator <focei|foceif|irlsfoceif|saem> [--outer_opt <nlminb|lbfgsb3c|bobyqa|NA>] [--boundary <narrow|wide|none|tight>] [--out_dir <path>] [--overwrite TRUE]")

  outer_opt <- a$outer_opt
  if (is.null(outer_opt) || identical(outer_opt, "NA") || identical(outer_opt, TRUE))
    outer_opt <- NA_character_

  bench_refit_one(
    N           = as.integer(a$N),
    scenario_id = as.integer(a$scenario),
    dataset_id  = as.integer(a$dataset),
    estimator   = as.character(a$estimator),
    outer_opt   = outer_opt,
    boundary    = as.character(a$boundary %||% "narrow"),
    master_rds  = a$master_rds,
    out_dir     = a$out_dir,
    overwrite   = isTRUE(as.logical(a$overwrite))
  )
}

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0L) {
  .run_cli()
}
