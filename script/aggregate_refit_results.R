# ==============================================================================
# aggregate_refit_results.R
# ------------------------------------------------------------------------------
# Load + aggregate the per-dataset RDS files produced by refit_one_dataset.R
# on HPCE for the Figure 3 / Table 3 refit robustness study.
#
# Produces FOUR long-format tibbles keyed on
# (cohort, scenario, boundary, [dataset_id, parameter]):
#
#   1. refit_diag_long   -- one row per (cohort, scenario, boundary, dataset)
#                            with 7 Table-3 diagnostic fields + PhysBnd +
#                            objf + cond_num_cor + runtime.  Answers Q3
#                            ("extract $diag items per dataset").
#
#   2. refit_rse_long    -- one row per (cohort, scenario, boundary, dataset,
#                            parameter): true_value, estimate, rel_err,
#                            sq_rel, abs_rel.  Answers Q1
#                            ("extract RSE per dataset per parameter").
#
#   3. refit_rmse_cells  -- one row per (cohort, scenario, boundary, parameter):
#                            RMRSE_pct, MARE_pct, MedRE_pct, MeanRE_pct, n_used.
#                            Answers Q2 ("aggregate RSE to RMSE per cell").
#                            Available in two flavours via status_filter:
#                              "all"          -- every fit contributes
#                              "success_only" -- min_suc & cov_step & !phys_bnd
#
#   4. refit_diag_rates  -- one row per (cohort, scenario, boundary) with
#                            % rates for each Table-3 diagnostic + median /
#                            max condition number + median runtime.  Ready
#                            to feed a Table-3 style wide pivot.
#
# Directory layout expected (matches scripts/hpce/refit_array.lsf output):
#
#   <root>/refit_true_<COHORT>/scn<NN>_<BND>/res_ds<DDD>.rds
#
# Usage (interactive):
#   source("scripts/aggregate_refit_results.R")
#   res <- aggregate_refit_run(root = "output")   # HPCE default
#   str(res, max.level = 1)
#
#   # inspect scenario 16 x N=300 diagnostic rates across boundaries
#   dplyr::filter(res$diag_rates, cohort == "N300", scenario == 16)
#
#   # inspect RMSE of the 4 covariate thetas across boundaries
#   dplyr::filter(res$rmse_cells, cohort == "N300", scenario == 16,
#                 parameter %in% c("CLBW", "CLcrCL", "VcBW", "VcSEX"))
#
# Usage (batch):
#   Rscript scripts/aggregate_refit_results.R [--root output]
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Helpers to unpack one RDS ---------------------------------------------
# refit_one_dataset() writes a named list with (relevant fields):
#   $cohort, $scenario_id, $dataset_id, $boundary, $status
#   $rel_err     -- 13-row tibble: parameter, true_value, estimate,
#                   abs_err, rel_err, rel_err_pct  (already computed)
#   $diag        -- 14-field named list (Table 3 mapping + PhysBnd + cov)
#   $runtime_sec, $timestamp
# For error results, only the meta keys + $status + $error_msg are present.

.unpack_diag <- function(r) {
  d <- r$diag %||% list()
  tibble::tibble(
    cohort       = r$cohort      %||% NA_character_,
    scenario     = r$scenario_id %||% NA_integer_,
    boundary     = r$boundary    %||% NA_character_,
    dataset_id   = r$dataset_id  %||% NA_integer_,
    status       = r$status      %||% NA_character_,
    converged    = as.logical(d$converged    %||% NA),
    min_suc      = as.logical(d$min_suc      %||% NA),
    cov_step     = as.logical(d$cov_step     %||% NA),
    rnd_err      = as.logical(d$rnd_err      %||% NA),
    zero_grad    = as.logical(d$zero_grad    %||% NA),
    est_bnd      = as.logical(d$est_bnd      %||% NA),
    phys_bnd     = as.logical(d$phys_bnd     %||% NA),
    objf         = as.numeric(d$objf         %||% NA_real_),
    cond_num_cor = as.numeric(d$cond_num_cor %||% NA_real_),
    bound_hits   = as.character(d$bound_hits %||% NA_character_),
    phys_hits    = as.character(d$phys_hits  %||% NA_character_),
    message      = as.character(d$message    %||% NA_character_),
    runtime_sec  = as.numeric(r$runtime_sec  %||% NA_real_),
    error_msg    = as.character(r$error_msg  %||% NA_character_)
  )
}

.unpack_rse <- function(r) {
  rel <- r$rel_err
  if (is.null(rel) || !nrow(rel)) return(tibble::tibble())
  tibble::tibble(
    cohort      = r$cohort      %||% NA_character_,
    scenario    = r$scenario_id %||% NA_integer_,
    boundary    = r$boundary    %||% NA_character_,
    dataset_id  = r$dataset_id  %||% NA_integer_,
    parameter   = rel$parameter,
    true_value  = rel$true_value,
    estimate    = rel$estimate,
    abs_err     = rel$abs_err,
    rel_err     = rel$rel_err,
    rel_err_pct = rel$rel_err_pct,
    sq_rel      = rel$rel_err^2,
    abs_rel     = abs(rel$rel_err)
  )
}

# Parse cohort / scenario / boundary from the file path so older RDS files
# missing those meta fields still get correctly grouped.
.meta_from_path <- function(path) {
  parts <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1]]
  cohort <- {
    m <- regmatches(parts, regexpr("^refit_true_(.*)$", parts))
    if (length(m)) sub("^refit_true_", "", m[1]) else NA_character_
  }
  leaf <- basename(dirname(path))                  # scn16_narrow
  sc   <- suppressWarnings(as.integer(sub("^scn(\\d+)_.*$", "\\1", leaf)))
  bnd  <- sub("^scn\\d+_", "", leaf)
  ds   <- suppressWarnings(as.integer(sub("^res_ds(\\d+)\\.rds$", "\\1",
                                          basename(path))))
  list(cohort = cohort, scenario = sc, boundary = bnd, dataset_id = ds)
}

load_one_refit_rds <- function(path) {
  r <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(r)) {
    message(sprintf("[skip] cannot read %s", path))
    return(list(diag = tibble::tibble(), rse = tibble::tibble()))
  }
  meta <- .meta_from_path(path)
  # Backfill any missing meta from the path
  r$cohort      <- r$cohort      %||% meta$cohort
  r$scenario_id <- r$scenario_id %||% meta$scenario
  r$boundary    <- r$boundary    %||% meta$boundary
  r$dataset_id  <- r$dataset_id  %||% meta$dataset_id
  list(diag = .unpack_diag(r), rse = .unpack_rse(r))
}

# ---- File discovery --------------------------------------------------------
discover_refit_files <- function(root       = "output",
                                 cohorts    = c("N40", "N80", "N300"),
                                 boundaries = c("none", "wide", "narrow"),
                                 scenarios  = 2:16) {
  purrr::map_dfr(cohorts, function(co) {
    co_dir <- file.path(root, paste0("refit_true_", co))
    if (!dir.exists(co_dir)) return(NULL)
    purrr::map_dfr(scenarios, function(sc) {
      purrr::map_dfr(boundaries, function(bnd) {
        sc_dir <- file.path(co_dir, sprintf("scn%02d_%s", sc, bnd))
        if (!dir.exists(sc_dir)) return(NULL)
        rds <- list.files(sc_dir, pattern = "^res_ds\\d+\\.rds$",
                          full.names = TRUE)
        if (!length(rds)) return(NULL)
        tibble::tibble(
          cohort     = co,
          scenario   = sc,
          boundary   = bnd,
          dataset_id = as.integer(sub("^res_ds(\\d+)\\.rds$", "\\1",
                                      basename(rds))),
          file_path  = rds,
          has_err    = file.exists(sub("\\.rds$", "_ERROR.txt", rds))
        )
      })
    })
  })
}

# ---- Bulk loader -----------------------------------------------------------
load_all_refit_results <- function(root = "output", verbose = TRUE, ...) {
  idx <- discover_refit_files(root = root, ...)
  if (!nrow(idx)) {
    warning("No refit RDS files found under ", root)
    return(list(index     = idx,
                diag_long = tibble::tibble(),
                rse_long  = tibble::tibble()))
  }
  if (verbose) {
    message(sprintf(
      "Loading %d RDS files across %d (cohort, scenario, boundary) cells",
      nrow(idx),
      dplyr::n_distinct(idx[, c("cohort", "scenario", "boundary")])
    ))
  }
  loaded <- purrr::map(idx$file_path, load_one_refit_rds)
  list(
    index     = idx,
    diag_long = purrr::map_dfr(loaded, "diag"),
    rse_long  = purrr::map_dfr(loaded, "rse")
  )
}

# ---- Aggregation: RSE -> RMSE per cell -------------------------------------
#   For each (cohort, scenario, boundary, parameter):
#     RMRSE_pct  = 100 * sqrt(mean(rel_err^2))         (Root Mean Rel Sq Err)
#     MARE_pct   = 100 * median(|rel_err|)             (Median Abs Rel Err)
#     MedRE_pct  = 100 * median(rel_err)               (bias, signed)
#     MeanRE_pct = 100 * mean(rel_err)                 (mean bias, signed)
#     n_used     = # of finite contributions
#
#   status_filter:
#     "all"          -> every fit's rel_err contributes (paper-style unfiltered)
#     "success_only" -> keep only fits with min_suc & cov_step & !phys_bnd
#                      (drops the degenerate MLE runs with TH_SEX_VC=371 etc)
#
#   Drops rows where true_value == 0 (parameter inactive in scenario) or
#   rel_err is NA / non-finite -- matches rel_err_one()'s NA convention.
compute_rmse_per_cell <- function(rse_long,
                                  diag_long     = NULL,
                                  status_filter = c("all", "success_only")) {
  status_filter <- match.arg(status_filter)

  df <- rse_long |>
    dplyr::filter(is.finite(true_value), true_value != 0,
                  is.finite(rel_err))

  if (status_filter == "success_only") {
    if (is.null(diag_long))
      stop("status_filter='success_only' requires diag_long argument")
    # Treat unknown flags (NA in older schema) as pass-through: min_suc=NA
    # -> keep, phys_bnd=NA -> keep. Only explicit FALSE / TRUE drop the fit.
    na_true  <- function(x) ifelse(is.na(x), TRUE,  x)
    na_false <- function(x) ifelse(is.na(x), FALSE, x)
    keep <- diag_long |>
      dplyr::filter(status == "ok",
                    na_true(min_suc),
                    na_true(cov_step),
                    !na_false(phys_bnd)) |>
      dplyr::select(cohort, scenario, boundary, dataset_id) |>
      dplyr::distinct()
    df <- dplyr::semi_join(df, keep,
                           by = c("cohort", "scenario", "boundary", "dataset_id"))
  }

  df |>
    dplyr::group_by(cohort, scenario, boundary, parameter) |>
    dplyr::summarise(
      true_value = dplyr::first(true_value),
      n_used     = dplyr::n(),
      MedRE_pct  = 100 * stats::median(rel_err),
      MeanRE_pct = 100 * mean(rel_err),
      MARE_pct   = 100 * stats::median(abs_rel),
      RMRSE_pct  = 100 * sqrt(mean(sq_rel)),
      .groups    = "drop"
    ) |>
    dplyr::arrange(cohort, scenario, boundary, parameter)
}

# ---- Aggregation: Table-3 diagnostic rates per cell ------------------------
compute_diag_rates_per_cell <- function(diag_long) {
  diag_long |>
    dplyr::group_by(cohort, scenario, boundary) |>
    dplyr::summarise(
      n_total        = dplyr::n(),
      n_status_ok    = sum(status == "ok", na.rm = TRUE),
      MinSuc_pct     = 100 * mean(min_suc,  na.rm = TRUE),
      CovStep_pct    = 100 * mean(cov_step, na.rm = TRUE),
      RndErr_pct     = 100 * mean(rnd_err,  na.rm = TRUE),
      PhysBnd_pct    = 100 * mean(phys_bnd, na.rm = TRUE),
      EstBnd_pct     = 100 * mean(est_bnd,  na.rm = TRUE),
      MedCN          = stats::median(cond_num_cor, na.rm = TRUE),
      StCN           = stats::sd(cond_num_cor,     na.rm = TRUE),
      MaxCN          = suppressWarnings(max(cond_num_cor, na.rm = TRUE)),
      MedObjF        = stats::median(objf, na.rm = TRUE),
      MedRuntime_min = stats::median(runtime_sec, na.rm = TRUE) / 60,
      .groups        = "drop"
    ) |>
    dplyr::mutate(MaxCN = ifelse(is.finite(MaxCN), MaxCN, NA_real_)) |>
    dplyr::arrange(cohort, scenario, boundary)
}

# ---- Wide-format renderers (Table 3 style pivot) --------------------------
# Table 3 style: rows = scenario, cols = boundary × metric
table3_diag_wide <- function(diag_rates, cohort_target = "N300",
                             metrics = c("MinSuc_pct", "CovStep_pct",
                                         "PhysBnd_pct", "MedCN")) {
  diag_rates |>
    dplyr::filter(cohort == cohort_target) |>
    dplyr::select(scenario, boundary, dplyr::all_of(metrics)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(metrics),
                        names_to = "metric", values_to = "value") |>
    tidyr::pivot_wider(names_from = c(boundary, metric),
                       values_from = value, names_sep = "_") |>
    dplyr::arrange(scenario)
}

# RMSE wide pivot for one parameter: rows = scenario, cols = boundary
table3_rmse_wide <- function(rmse_cells, cohort_target = "N300",
                             parameter_target = "CLBW") {
  rmse_cells |>
    dplyr::filter(cohort    == cohort_target,
                  parameter == parameter_target) |>
    dplyr::select(scenario, boundary, RMRSE_pct) |>
    tidyr::pivot_wider(names_from = boundary, values_from = RMRSE_pct) |>
    dplyr::arrange(scenario)
}

# ---- Top-level driver ------------------------------------------------------
aggregate_refit_run <- function(root          = "output",
                                out_dir       = file.path(root, "refit_aggregated"),
                                write_outputs = TRUE,
                                verbose       = TRUE,
                                ...) {
  loaded <- load_all_refit_results(root = root, verbose = verbose, ...)
  if (!nrow(loaded$diag_long)) return(loaded)

  diag_rates      <- compute_diag_rates_per_cell(loaded$diag_long)
  rmse_all        <- compute_rmse_per_cell(loaded$rse_long,
                                           status_filter = "all")
  rmse_success    <- compute_rmse_per_cell(loaded$rse_long,
                                           diag_long     = loaded$diag_long,
                                           status_filter = "success_only")

  if (write_outputs) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(loaded$index,   file.path(out_dir, "refit_file_index.csv"))
    readr::write_csv(loaded$diag_long, file.path(out_dir, "refit_diag_long.csv"))
    readr::write_csv(loaded$rse_long,  file.path(out_dir, "refit_rse_long.csv"))
    readr::write_csv(diag_rates,       file.path(out_dir, "refit_diag_rates.csv"))
    readr::write_csv(rmse_all,         file.path(out_dir, "refit_rmse_all.csv"))
    readr::write_csv(rmse_success,     file.path(out_dir, "refit_rmse_success.csv"))
    saveRDS(list(
      index         = loaded$index,
      diag_long     = loaded$diag_long,
      rse_long      = loaded$rse_long,
      diag_rates    = diag_rates,
      rmse_all      = rmse_all,
      rmse_success  = rmse_success,
      created_at    = Sys.time()
    ), file.path(out_dir, "refit_aggregated.rds"))
    if (verbose)
      message(sprintf("Wrote 6 CSVs + 1 RDS to %s/", out_dir))
  }

  invisible(list(
    index        = loaded$index,
    diag_long    = loaded$diag_long,
    rse_long     = loaded$rse_long,
    diag_rates   = diag_rates,
    rmse_all     = rmse_all,
    rmse_success = rmse_success
  ))
}

# ---- CLI dispatcher --------------------------------------------------------
# Run everything when invoked as: Rscript aggregate_refit_results.R [--root <dir>]
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0L) {
  argv <- commandArgs(trailingOnly = TRUE)
  root <- "output"
  i <- 1L
  while (i <= length(argv)) {
    if (argv[i] == "--root" && i + 1L <= length(argv)) {
      root <- argv[i + 1L]; i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  aggregate_refit_run(root = root)
}
