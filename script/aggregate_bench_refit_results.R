# ==============================================================================
# aggregate_bench_refit_results.R
# ------------------------------------------------------------------------------
# Aggregate the true-model refit benchmark produced by bench_refit_estimators.R.
#
# Directory scanned:
#   output/bench_refit_estimators/N<NN>/scn<SS>/<est>_<opt>/res_ds<DDD>.rds
#
# Two output families:
#
#  A. Convergence / model stability
#     bench_refit_diag_long   (per fit; converged, cov_ok, cond#, timings, msg)
#     bench_refit_diag_rates  (per cell; %converged, %cov, med/max CN, med OFV)
#
#  B. Estimation precision
#     bench_refit_rse_long    (per (fit, parameter); rel_err vs truth)
#     bench_refit_rmse_cells  (per cell x parameter; N, MedRE, RMRSE, MedAbsRE)
#     bench_refit_rmse_cells_success  (same, restricted to converged fits)
#
# Usage:
#   Rscript script/aggregate_bench_refit_results.R [--root output] [--out_dir <path>]
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(readr)
  library(fs)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# ---- discovery ------------------------------------------------------------
discover_refit_files <- function(root = "output") {
  base <- file.path(root, "bench_refit_estimators")
  if (!dir.exists(base))
    stop("Directory not found: ", base)
  files <- list.files(base, pattern = "^res_ds\\d+\\.rds$",
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("No res_ds*.rds under ", base)

  # Walk 3 parents up from each file:
  # .../N<NN>/scn<SS>/<est>_<opt>/res_ds<DDD>.rds
  cell_dir <- basename(dirname(files))                    # <est>_<opt>
  scn_dir  <- basename(dirname(dirname(files)))           # scn<SS>
  n_dir    <- basename(dirname(dirname(dirname(files))))  # N<NN>

  tibble(
    file_path  = files,
    n_dir      = n_dir,
    scn_dir    = scn_dir,
    cell_dir   = cell_dir,
    dataset_id = as.integer(sub("^res_ds(\\d+)\\.rds$", "\\1", basename(files))),
    sample_N   = suppressWarnings(as.integer(sub("^N",   "", n_dir))),
    scenario   = suppressWarnings(as.integer(sub("^scn", "", scn_dir))),
    estimator  = sub("_[^_]+$", "", cell_dir),
    outer_opt  = sub("^[^_]+_",  "", cell_dir)
  ) |>
    select(sample_N, scenario, estimator, outer_opt, dataset_id, file_path)
}

# ---- unpack one RDS -------------------------------------------------------
.unpack_diag_one <- function(r, meta) {
  d <- r$diag %||% list()
  d3 <- r$diag_t3 %||% list()
  tibble(
    sample_N        = r$sample_N    %||% meta$sample_N,
    scenario        = r$scenario_id %||% meta$scenario,
    estimator       = r$estimator   %||% meta$estimator,
    outer_opt       = r$outer_opt   %||% meta$outer_opt,
    dataset_id      = r$dataset_id  %||% meta$dataset_id,
    boundary        = as.character(r$boundary %||% NA_character_),
    refit_estimator = as.character(r$refit_estimator %||% NA_character_),
    status          = as.character(r$status %||% NA_character_),
    # 2026-07-12: Recompute `converged` from numerical evidence rather than
    # trusting the saved flag. Older HPCE runs used a stale diagnose_fit()
    # that gated on fit$convergence==0, which lbfgsb3c fails ("code 8, false
    # convergence") even for perfectly good fits. Rule: finite OFV +
    # cov_ok + finite cond_num_cor. Falls back to saved flag if any input
    # is missing (e.g. status="fit_failed").
    converged       = {
      objf_ok <- is.finite(as.numeric(d$objf %||% NA_real_))
      cov_ok_ <- isTRUE(as.logical(d$cov_ok %||% NA))
      cn_ok   <- is.finite(as.numeric(d$cond_num_cor %||% NA_real_))
      recomputed <- objf_ok && cov_ok_ && cn_ok
      if (is.na(objf_ok) || is.na(cov_ok_)) as.logical(d$converged %||% NA)
      else recomputed
    },
    cn_below_cutoff = {
      # PMx-standard: cond_num_cor <= 1000 (NONMEM default warning threshold,
      # Khandelwal 2019 Fig 7). NA/non-finite -> FALSE (safe).
      cn <- suppressWarnings(as.numeric(d$cond_num_cor %||% NA_real_))
      isTRUE(is.finite(cn) && cn <= 1000)
    },
    cov_ok          = as.logical(d$cov_ok      %||% NA),
    objf            = as.numeric(d$objf        %||% NA_real_),
    cond_num_cor    = as.numeric(d$cond_num_cor %||% NA_real_),
    cond_num_cor_source = as.character(d$cond_num_cor_source %||% NA_character_),
    convergence_code    = suppressWarnings(as.integer(d$convergence_code %||% NA_integer_)),
    # Table-3 style extras (may be NULL for non-final tier or SAEM):
    min_suc         = as.logical(d3$min_suc  %||% NA),
    cov_step        = as.logical(d3$cov_step %||% NA),
    est_bnd         = as.logical(d3$est_bnd  %||% NA),
    phys_bnd        = as.logical(d3$phys_bnd %||% NA),
    # rnd_err / zero_grad are NA (not FALSE) for SAEM -- gradient-free, no
    # outer optimizer. est_family + cov_source let the report footnote that
    # SAEM cond numbers come from linFim, not the r/s sandwich.
    rnd_err         = as.logical(d3$rnd_err   %||% NA),
    zero_grad       = as.logical(d3$zero_grad %||% NA),
    est_family      = as.character(d3$est_family %||% NA_character_),
    cov_source      = as.character(d3$cov_source %||% NA_character_),
    diag_note       = as.character(d3$diag_note  %||% NA_character_),
    fit_runtime_sec   = as.numeric(r$fit_runtime_sec   %||% NA_real_),
    refit_runtime_sec = as.numeric(r$refit_runtime_sec %||% NA_real_),
    message         = as.character(d$message %||% NA_character_)
  )
}

.unpack_rse_one <- function(r, meta) {
  rel <- r$rel_err
  if (is.null(rel) || !nrow(rel)) return(tibble())
  tibble(
    sample_N   = r$sample_N    %||% meta$sample_N,
    scenario   = r$scenario_id %||% meta$scenario,
    estimator  = r$estimator   %||% meta$estimator,
    outer_opt  = r$outer_opt   %||% meta$outer_opt,
    dataset_id = r$dataset_id  %||% meta$dataset_id,
    boundary   = as.character(r$boundary %||% NA_character_),
    parameter  = rel$parameter,
    true_value = rel$true_value,
    estimate   = rel$estimate,
    rel_err    = rel$rel_err
  )
}

# ---- driver ---------------------------------------------------------------
load_all_refit_results <- function(root = "output", verbose = TRUE) {
  idx <- discover_refit_files(root)
  if (verbose) message(sprintf(
    "Loading %d RDS across %d cells",
    nrow(idx),
    dplyr::n_distinct(paste(idx$sample_N, idx$scenario, idx$estimator, idx$outer_opt))))

  diag_rows <- vector("list", nrow(idx))
  rse_rows  <- vector("list", nrow(idx))
  for (i in seq_len(nrow(idx))) {
    meta <- as.list(idx[i, ])
    r <- tryCatch(readRDS(meta$file_path), error = function(e) NULL)
    if (is.null(r)) {
      message("[skip] cannot read ", meta$file_path); next
    }
    diag_rows[[i]] <- .unpack_diag_one(r, meta)
    rse_rows[[i]]  <- .unpack_rse_one(r,  meta)
  }
  list(
    diag_long = dplyr::bind_rows(diag_rows),
    rse_long  = dplyr::bind_rows(rse_rows)
  )
}

# ---- convergence / stability ----------------------------------------------
compute_refit_diag_rates <- function(diag_long) {
  diag_long |>
    # Build PMx-standard strict convergence flag: numerical convergence AND
    # CN <= 1000 AND no bounded parameter at its active bound. NA-safe:
    # any missing input -> FALSE (so aggregation never breaks).
    dplyr::mutate(
      converged_strict = (converged %in% TRUE) &
                         (cn_below_cutoff %in% TRUE) &
                         !(est_bnd %in% TRUE)
    ) |>
    dplyr::group_by(sample_N, scenario, estimator, outer_opt, boundary) |>
    dplyr::summarise(
      n_total             = dplyr::n(),
      n_ok_status         = sum(status == "ok", na.rm = TRUE),
      Converged_pct       = 100 * mean(converged, na.rm = TRUE),
      ConvergedStrict_pct = 100 * mean(converged_strict),  # no NA in strict
      CovOk_pct           = 100 * mean(cov_ok,    na.rm = TRUE),
      CNBelowCutoff_pct   = 100 * mean(cn_below_cutoff),   # no NA
      MinSuc_pct          = 100 * mean(min_suc,   na.rm = TRUE),
      EstBnd_pct          = 100 * mean(est_bnd,   na.rm = TRUE),
      PhysBnd_pct         = 100 * mean(phys_bnd,  na.rm = TRUE),
      MedCN               = stats::median(cond_num_cor, na.rm = TRUE),
      MaxCN               = suppressWarnings(max(cond_num_cor, na.rm = TRUE)),
      P90CN               = as.numeric(stats::quantile(cond_num_cor, 0.9, na.rm = TRUE)),
      MedObjF             = stats::median(objf, na.rm = TRUE),
      MedFitSec           = stats::median(fit_runtime_sec,   na.rm = TRUE),
      MedRefitSec         = stats::median(refit_runtime_sec, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(dplyr::across(where(is.numeric), \(x) ifelse(is.finite(x), x, NA_real_)))
}

# ---- estimation precision -------------------------------------------------
compute_refit_rmse <- function(rse_long, restrict_success = NULL) {
  if (!is.null(restrict_success)) {
    key <- c("sample_N","scenario","estimator","outer_opt","boundary","dataset_id")
    # 2026-07-12: restrict to PMx-strict convergence (numerical + CN<=1000 +
    # no boundary hit), matching Khandelwal Fig 7 practice. NA-safe: any
    # missing input treated as FALSE so aggregation never breaks.
    ok <- restrict_success |>
      dplyr::mutate(
        converged_strict = (converged %in% TRUE) &
                           (cn_below_cutoff %in% TRUE) &
                           !(est_bnd %in% TRUE)
      ) |>
      dplyr::filter(converged_strict) |>
      dplyr::select(dplyr::all_of(key))
    rse_long <- rse_long |> dplyr::inner_join(ok, by = key)
  }
  rse_long |>
    dplyr::group_by(sample_N, scenario, estimator, outer_opt, boundary, parameter) |>
    dplyr::summarise(
      N          = dplyr::n(),
      MedRE_pct  = 100 * stats::median(rel_err, na.rm = TRUE),
      MeanRE_pct = 100 * mean(rel_err,          na.rm = TRUE),
      RMRSE_pct  = 100 * sqrt(mean(rel_err^2,   na.rm = TRUE)),
      MedAbsRE_pct = 100 * stats::median(abs(rel_err), na.rm = TRUE),
      P90AbsRE_pct = 100 * as.numeric(stats::quantile(abs(rel_err), 0.9, na.rm = TRUE)),
      MedEst     = stats::median(estimate, na.rm = TRUE),
      MedTrue    = stats::median(true_value, na.rm = TRUE),
      .groups = "drop"
    )
}

na_true <- function(x) !is.na(x) & isTRUE(x) | (is.logical(x) & !is.na(x) & x)

# ---- top-level driver -----------------------------------------------------
aggregate_bench_refit_run <- function(root    = "output",
                                      out_dir = file.path(root, "bench_refit_aggregated"),
                                      write_outputs = TRUE, verbose = TRUE) {
  loaded <- load_all_refit_results(root, verbose = verbose)

  diag_rates       <- compute_refit_diag_rates(loaded$diag_long)
  rmse_all         <- compute_refit_rmse(loaded$rse_long)
  rmse_success     <- compute_refit_rmse(loaded$rse_long,
                                         restrict_success = loaded$diag_long)

  if (write_outputs) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(loaded$diag_long, file.path(out_dir, "bench_refit_diag_long.csv"))
    readr::write_csv(loaded$rse_long,  file.path(out_dir, "bench_refit_rse_long.csv"))
    readr::write_csv(diag_rates,       file.path(out_dir, "bench_refit_diag_rates.csv"))
    readr::write_csv(rmse_all,         file.path(out_dir, "bench_refit_rmse_cells.csv"))
    readr::write_csv(rmse_success,     file.path(out_dir, "bench_refit_rmse_cells_success.csv"))
    saveRDS(list(diag_long = loaded$diag_long, rse_long = loaded$rse_long,
                 diag_rates = diag_rates,
                 rmse_all = rmse_all, rmse_success = rmse_success),
            file.path(out_dir, "bench_refit_aggregated.rds"))
    if (verbose) message("Wrote 5 CSVs + 1 RDS to ", out_dir)
  }

  invisible(list(
    diag_long = loaded$diag_long,
    rse_long  = loaded$rse_long,
    diag_rates = diag_rates,
    rmse_all = rmse_all,
    rmse_success = rmse_success
  ))
}

# ---- CLI ------------------------------------------------------------------
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

if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0L) {
  a <- .parse_cli(commandArgs(trailingOnly = TRUE))
  aggregate_bench_refit_run(
    root    = a$root    %||% "output",
    out_dir = a$out_dir %||% file.path(a$root %||% "output", "bench_refit_aggregated")
  )
}
