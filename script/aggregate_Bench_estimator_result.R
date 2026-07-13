# ==============================================================================
# aggregate_Bench_estimator_result.R
# ------------------------------------------------------------------------------
# Load + aggregate the per-dataset RDS files produced by
# PerformanceEvaluation_scm_bench.R on HPCE for the estimator x optimizer
# SCM benchmark, using the same convention as the pilot report in
# scenario16console_output_report.R (Power / PowerCN / PowerMinSuc / relpower).
#
# Directory layout expected:
#   <root>/scm_bench/N<NN>/scn<SS>/<est>_<opt>/res_ds<DD>.rds
# (opt = "NA" for saem)
#
# Per-dataset extraction (extract-per-ds, aggregate-per-cell):
#
#   -- diagnostics from r$test$diag:
#        converged, cov_ok, objf, cond_num_cor, message
#   -- timings from top-level: t_base_sec, t_scm_sec, t_refit_sec, t_total_sec
#   -- rel_err from r$test$rel_err: parameter/true_value/estimate/rel_err
#   -- SCM selection from r$test$selected: matched to truth (extract_true_relations
#      + match_selected_to_truth) to give n_true_hit, n_false_pos, exact_match
#
# Per-cell aggregation, keyed on (sample_N, scenario, estimator, outer_opt):
#
#   1. bench_diag_rates -- %Converged, %CovStep, Med/Max CN, MedObjF,
#                          Med runtime.
#
#   2. bench_rmse_all / bench_rmse_success -- per (cell, parameter):
#        RMRSE_pct = 100*sqrt(mean(rel_err^2))
#        MARE_pct  = 100*median(|rel_err|)
#        MedRE_pct, MeanRE_pct (bias), n_used
#      success_only := converged & cov_ok.
#
#   3. bench_power -- per cell: Power / PowerCN / PowerMinSuc
#        Power       = mean(exact_match)             over all datasets
#        PowerCN     = mean(exact_match | cond_num_cor < cn_cor_cut)
#                       -- denominator = # datasets with acceptable CN
#        PowerMinSuc = mean(exact_match | converged)
#                       -- denominator = # datasets with min success
#      "exact_match" = SCM selected exactly the true covariate-parameter
#      shape triples for that scenario (no misses, no false positives).
#
#   4. bench_relpower -- per (cell, k): fraction of datasets recovering
#      at least k of the n_true true covariates (k = 1..n_true).
#
#   5. bench_rmrse_cond -- per (cell, parameter) RMRSE conditional on
#      exact_match (SCM recovered the true structure).
#
# Usage:
#   source("script/aggregate_Bench_estimator_result.R")
#   tp  <- readRDS("Inputdataset/true_params_long.rds")
#   res <- aggregate_bench_run(root = "output", true_params = tp)
#
#   dplyr::filter(res$power,      sample_N == 80, scenario == 16)
#   dplyr::filter(res$relpower,   sample_N == 80, scenario == 16)
#   dplyr::filter(res$rmse_success, sample_N == 80, scenario == 16,
#                 parameter %in% c("TVCL","TVVc"))
#
# CLI:
#   Rscript script/aggregate_Bench_estimator_result.R \
#       --root output --true_params Inputdataset/true_params_long.rds
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Truth utilities (borrowed from scenario16console_output_report.R) -----
# All Khandelwal-2019 style benchmark covariates supported by the SCM grid.
# Maps a "parameter" name in true_params_long.rds to (var, covar, shape).
BENCH_SHAPE_MAP <- list(
  CLBW    = list(var = "cl", covar = "BW",   shape = "power"),
  CLcrCL  = list(var = "cl", covar = "CrCL", shape = "power"),
  CLCrCL  = list(var = "cl", covar = "CrCL", shape = "power"),
  CLBMI   = list(var = "cl", covar = "BMI",  shape = "power"),
  CLSEX   = list(var = "cl", covar = "SEX",  shape = "cat"),
  CLRACE  = list(var = "cl", covar = "RACE", shape = "cat"),
  VcBW    = list(var = "vc", covar = "BW",   shape = "power"),
  VccrCL  = list(var = "vc", covar = "CrCL", shape = "power"),
  VcCrCL  = list(var = "vc", covar = "CrCL", shape = "power"),
  VcBMI   = list(var = "vc", covar = "BMI",  shape = "power"),
  VcSEX   = list(var = "vc", covar = "SEX",  shape = "cat"),
  VcRACE  = list(var = "vc", covar = "RACE", shape = "cat")
)

.canon_shape_tbl <- function(rel) {
  if (is.null(rel) || !nrow(rel))
    return(tibble::tibble(var = character(), covar = character(),
                          shape = character()))
  rel |>
    dplyr::mutate(
      var   = as.character(var),
      covar = as.character(covar),
      shape = as.character(shape),
      shape = ifelse(grepl("^[0-9]+$", shape), "cat", shape)
    ) |>
    dplyr::distinct(var, covar, shape)
}

extract_true_relations <- function(true_long, scenario_id,
                                   shape_map = BENCH_SHAPE_MAP) {
  # Accept legacy column names
  tp <- true_long
  if (!"scenario" %in% names(tp) && "scenario_id" %in% names(tp))
    tp$scenario <- tp$scenario_id
  if (!"true_value" %in% names(tp) && "value" %in% names(tp))
    tp$true_value <- tp$value
  if (!"parameter" %in% names(tp) && "param" %in% names(tp))
    tp$parameter <- tp$param

  hits <- tp |>
    dplyr::filter(scenario == scenario_id,
                  parameter %in% names(shape_map),
                  !is.na(true_value), true_value != 0) |>
    dplyr::pull(parameter)
  if (!length(hits))
    return(tibble::tibble(var = character(), covar = character(),
                          shape = character()))
  purrr::map_dfr(hits, function(p) {
    m <- shape_map[[p]]
    tibble::tibble(var = m$var, covar = m$covar, shape = m$shape)
  })
}

match_selected_to_truth <- function(selected, true_rel) {
  if (is.null(selected) || !nrow(selected)) {
    return(list(n_true_hit = 0L, n_false_pos = 0L,
                exact_match = (nrow(true_rel) == 0L)))
  }
  sel <- .canon_shape_tbl(selected)
  tru <- .canon_shape_tbl(true_rel)
  true_hit  <- dplyr::inner_join(sel, tru, by = c("var","covar","shape"))
  false_pos <- dplyr::anti_join (sel, tru, by = c("var","covar","shape"))
  miss      <- dplyr::anti_join (tru, sel, by = c("var","covar","shape"))
  list(n_true_hit  = nrow(true_hit),
       n_false_pos = nrow(false_pos),
       exact_match = (nrow(false_pos) == 0L) && (nrow(miss) == 0L))
}

# ---- File discovery --------------------------------------------------------
discover_bench_files <- function(root = "output") {
  base <- file.path(root, "scm_bench")
  if (!dir.exists(base)) {
    warning("Directory not found: ", base); return(tibble::tibble())
  }
  rds <- list.files(base, pattern = "^res_ds\\d+\\.rds$",
                    recursive = TRUE, full.names = TRUE)
  if (!length(rds)) return(tibble::tibble())
  parts <- strsplit(gsub("\\\\", "/", rds), "/", fixed = TRUE)
  tibble::tibble(
    file_path  = rds,
    cell_dir   = vapply(parts, function(p) p[length(p) - 1L], character(1)),
    scn_dir    = vapply(parts, function(p) p[length(p) - 2L], character(1)),
    n_dir      = vapply(parts, function(p) p[length(p) - 3L], character(1)),
    dataset_id = as.integer(sub("^res_ds(\\d+)\\.rds$", "\\1", basename(rds)))
  ) |>
    dplyr::mutate(
      sample_N  = suppressWarnings(as.integer(sub("^N",   "", n_dir))),
      scenario  = suppressWarnings(as.integer(sub("^scn", "", scn_dir))),
      estimator = sub("_[^_]+$", "", cell_dir),
      outer_opt = sub("^[^_]+_",  "", cell_dir)
    ) |>
    dplyr::select(sample_N, scenario, estimator, outer_opt, dataset_id, file_path)
}

# ---- Unpack one RDS (diag + rse + match to truth) --------------------------
# Re-runs diagnose_fit() on the stored final_fit when the diag list looks
# "stale" (missing convergence_code, or converged=FALSE with a valid cov_done).
# This lets us rescue old RDS produced before the 2026-07 diagnose_fit rewrite
# without re-running any HPCE fits.
.rediagnose_if_stale <- function(t) {
  d <- t$diag %||% list()
  is_stale <- is.null(d$convergence_code) || is.null(d$cond_num_cor_source)
  fit <- t$final_fit
  if (is_stale && !is.null(fit) && exists("diagnose_fit", mode = "function")) {
    d <- tryCatch(diagnose_fit(fit), error = function(e) d)
  }
  d
}

.unpack_diag <- function(r, meta, hit) {
  t <- r$test %||% list()
  d <- .rediagnose_if_stale(t)
  cn_val <- suppressWarnings(as.numeric(d$cond_num_cor %||% NA_real_))
  tibble::tibble(
    sample_N        = r$sample_N    %||% meta$sample_N,
    scenario        = r$scenario_id %||% meta$scenario,
    estimator       = r$estimator   %||% meta$estimator,
    outer_opt       = r$outer_opt   %||% meta$outer_opt,
    dataset_id      = r$dataset_id  %||% meta$dataset_id,
    converged       = as.logical(d$converged %||% NA),
    # PMx-standard CN gate (Khandelwal 2019 Fig 7 / NONMEM default = 1000).
    # NA / non-finite -> FALSE so downstream aggregation cannot break.
    cn_below_cutoff = isTRUE(is.finite(cn_val) && cn_val <= 1000),
    cov_ok          = as.logical(t$cov_done  %||% d$cov_ok %||% NA),
    objf            = as.numeric(d$objf         %||% NA_real_),
    cond_num_cor    = cn_val,
    cond_num_cor_source = as.character(d$cond_num_cor_source %||% NA_character_),
    convergence_code    = suppressWarnings(as.integer(d$convergence_code %||% NA_integer_)),
    n_selected      = if (!is.null(t$selected)) nrow(t$selected) else NA_integer_,
    n_true_hit      = hit$n_true_hit,
    n_false_pos     = hit$n_false_pos,
    exact_match     = hit$exact_match,
    refit_estimator = as.character(t$refit_estimator %||% NA_character_),
    t_base_sec      = as.numeric(r$t_base_sec  %||% NA_real_),
    t_scm_sec       = as.numeric(r$t_scm_sec   %||% NA_real_),
    t_refit_sec     = as.numeric(r$t_refit_sec %||% NA_real_),
    t_total_sec     = as.numeric(r$t_total_sec %||% NA_real_),
    message         = as.character(d$message   %||% NA_character_)
  )
}

.unpack_rse <- function(r, meta) {
  rel <- r$test$rel_err
  if (is.null(rel) || !nrow(rel)) return(tibble::tibble())
  tibble::tibble(
    sample_N    = r$sample_N    %||% meta$sample_N,
    scenario    = r$scenario_id %||% meta$scenario,
    estimator   = r$estimator   %||% meta$estimator,
    outer_opt   = r$outer_opt   %||% meta$outer_opt,
    dataset_id  = r$dataset_id  %||% meta$dataset_id,
    parameter   = rel$parameter,
    true_value  = rel$true_value,
    estimate    = rel$estimate,
    rel_err     = rel$rel_err,
    rel_err_pct = rel$rel_err_pct,
    sq_rel      = rel$rel_err^2,
    abs_rel     = abs(rel$rel_err)
  )
}

load_one_bench_rds <- function(path, meta, true_rel_by_scn) {
  r <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(r)) {
    message(sprintf("[skip] cannot read %s", path))
    return(list(diag = tibble::tibble(), rse = tibble::tibble()))
  }
  scn_key <- as.character(r$scenario_id %||% meta$scenario)
  true_rel <- true_rel_by_scn[[scn_key]] %||%
              tibble::tibble(var = character(), covar = character(),
                             shape = character())
  hit <- match_selected_to_truth(r$test$selected, true_rel)
  list(diag = .unpack_diag(r, meta, hit),
       rse  = .unpack_rse(r, meta))
}

# ---- Bulk loader -----------------------------------------------------------
load_all_bench_results <- function(root = "output", true_params = NULL,
                                   verbose = TRUE) {
  idx <- discover_bench_files(root)
  if (!nrow(idx)) {
    warning("No bench RDS files found under ", root, "/scm_bench")
    return(list(index = idx, diag_long = tibble::tibble(),
                rse_long = tibble::tibble(), true_rel_by_scn = list()))
  }

  # Pre-compute truth per scenario (once) for fast per-file matching
  scns <- sort(unique(idx$scenario))
  true_rel_by_scn <- if (!is.null(true_params)) {
    setNames(lapply(scns, function(s)
      extract_true_relations(true_params, s)), as.character(scns))
  } else {
    setNames(lapply(scns, function(s)
      tibble::tibble(var = character(), covar = character(),
                     shape = character())), as.character(scns))
  }

  if (verbose) {
    message(sprintf(
      "Loading %d RDS across %d cells (N x scn x estimator x outer_opt)",
      nrow(idx),
      dplyr::n_distinct(idx[, c("sample_N","scenario","estimator","outer_opt")])
    ))
  }
  loaded <- purrr::map(seq_len(nrow(idx)), function(i) {
    load_one_bench_rds(idx$file_path[i], meta = as.list(idx[i, ]),
                       true_rel_by_scn = true_rel_by_scn)
  })
  list(
    index           = idx,
    diag_long       = purrr::map_dfr(loaded, "diag"),
    rse_long        = purrr::map_dfr(loaded, "rse"),
    true_rel_by_scn = true_rel_by_scn
  )
}

# ---- Per-dataset summary (one row per fit) --------------------------------
# Joins diag_long (already per-dataset) with a per-fit RMRSE roll-up that
# collapses ACROSS parameters within the fit. Each *_across_params_pct
# column is therefore a single scalar per fit, computed over the fit's
# identifiable parameter set (true_value != 0 & finite rel_err):
#
#   RMRSE_across_params_pct = 100 * sqrt(mean(rel_err^2))
#   MARE_across_params_pct  = 100 * median(|rel_err|)
#   MedRE_across_params_pct = 100 * median(rel_err)               (signed bias)
#   MaxAbsRE_across_params_pct = 100 * max(|rel_err|)             (worst param)
#
# Use as a coarse fit-quality score for ranking / outlier hunting. Do NOT
# compare directly across cells if n_params_used differs (e.g. SAEM drops
# unidentifiable params). For per-parameter comparisons use rse_long or
# the per-cell rmse_* tibbles.
compute_bench_per_dataset <- function(diag_long, rse_long) {
  rse_ds <- rse_long |>
    dplyr::filter(is.finite(true_value), true_value != 0,
                  is.finite(rel_err)) |>
    dplyr::group_by(sample_N, scenario, estimator, outer_opt, dataset_id) |>
    dplyr::summarise(
      n_params_used              = dplyr::n(),
      MedRE_across_params_pct    = 100 * stats::median(rel_err),
      MARE_across_params_pct     = 100 * stats::median(abs_rel),
      RMRSE_across_params_pct    = 100 * sqrt(mean(sq_rel)),
      MaxAbsRE_across_params_pct = 100 * max(abs_rel),
      .groups                    = "drop"
    )
  dplyr::left_join(diag_long, rse_ds,
    by = c("sample_N","scenario","estimator","outer_opt","dataset_id")) |>
    dplyr::arrange(sample_N, scenario, estimator, outer_opt, dataset_id)
}

# ---- Aggregation: convergence + cov step + timings per cell ---------------
compute_bench_diag_rates <- function(diag_long) {
  diag_long |>
    dplyr::mutate(
      converged_strict = (converged %in% TRUE) & (cn_below_cutoff %in% TRUE)
      # note: est_bnd not tracked in SCM diag_long (base model has no
      # bounded covariate thetas); strict flag omits that criterion here.
    ) |>
    dplyr::group_by(sample_N, scenario, estimator, outer_opt) |>
    dplyr::summarise(
      n_total             = dplyr::n(),
      Converged_pct       = 100 * mean(converged, na.rm = TRUE),
      ConvergedStrict_pct = 100 * mean(converged_strict),
      CovStep_pct         = 100 * mean(cov_ok,    na.rm = TRUE),
      CNBelowCutoff_pct   = 100 * mean(cn_below_cutoff),
      MedCN               = stats::median(cond_num_cor, na.rm = TRUE),
      MaxCN               = suppressWarnings(max(cond_num_cor, na.rm = TRUE)),
      MedObjF             = stats::median(objf, na.rm = TRUE),
      MedNSelected        = stats::median(n_selected, na.rm = TRUE),
      MedTotal_min        = stats::median(t_total_sec, na.rm = TRUE) / 60,
      MedBase_min         = stats::median(t_base_sec,  na.rm = TRUE) / 60,
      MedSCM_min          = stats::median(t_scm_sec,   na.rm = TRUE) / 60,
      .groups             = "drop"
    ) |>
    dplyr::mutate(MaxCN = ifelse(is.finite(MaxCN), MaxCN, NA_real_)) |>
    dplyr::arrange(sample_N, scenario, estimator, outer_opt)
}

# ---- Aggregation: RMRSE per (cell, parameter) -----------------------------
compute_bench_rmse <- function(rse_long, diag_long = NULL,
                               status_filter = c("all","success_only")) {
  status_filter <- match.arg(status_filter)
  df <- rse_long |>
    dplyr::filter(is.finite(true_value), true_value != 0,
                  is.finite(rel_err))
  if (status_filter == "success_only") {
    if (is.null(diag_long))
      stop("status_filter='success_only' requires diag_long")
    # PMx-strict: converged (numerical) AND CN <= 1000. NA-safe.
    keep <- diag_long |>
      dplyr::mutate(
        converged_strict = (converged %in% TRUE) & (cn_below_cutoff %in% TRUE)
      ) |>
      dplyr::filter(converged_strict) |>
      dplyr::select(sample_N, scenario, estimator, outer_opt, dataset_id) |>
      dplyr::distinct()
    df <- dplyr::semi_join(df, keep,
      by = c("sample_N","scenario","estimator","outer_opt","dataset_id"))
  }
  df |>
    dplyr::group_by(sample_N, scenario, estimator, outer_opt, parameter) |>
    dplyr::summarise(
      true_value = dplyr::first(true_value),
      n_used     = dplyr::n(),
      MedRE_pct  = 100 * stats::median(rel_err),
      MeanRE_pct = 100 * mean(rel_err),
      MARE_pct   = 100 * stats::median(abs_rel),
      RMRSE_pct  = 100 * sqrt(mean(sq_rel)),
      .groups    = "drop"
    ) |>
    dplyr::arrange(sample_N, scenario, estimator, outer_opt, parameter)
}

# ---- Aggregation: RMRSE conditional on exact_match ------------------------
compute_bench_rmse_cond <- function(rse_long, diag_long) {
  keep <- diag_long |>
    dplyr::filter(isTRUE(exact_match) | exact_match) |>
    dplyr::select(sample_N, scenario, estimator, outer_opt, dataset_id) |>
    dplyr::distinct()
  df <- rse_long |>
    dplyr::filter(is.finite(true_value), true_value != 0,
                  is.finite(rel_err)) |>
    dplyr::semi_join(keep,
      by = c("sample_N","scenario","estimator","outer_opt","dataset_id"))
  if (!nrow(df))
    return(tibble::tibble(sample_N=integer(), scenario=integer(),
                          estimator=character(), outer_opt=character(),
                          parameter=character(), n_used=integer(),
                          RMRSE_pct=numeric(), MARE_pct=numeric()))
  df |>
    dplyr::group_by(sample_N, scenario, estimator, outer_opt, parameter) |>
    dplyr::summarise(
      true_value = dplyr::first(true_value),
      n_used     = dplyr::n(),
      MedRE_pct  = 100 * stats::median(rel_err),
      MARE_pct   = 100 * stats::median(abs_rel),
      RMRSE_pct  = 100 * sqrt(mean(sq_rel)),
      .groups    = "drop"
    ) |>
    dplyr::arrange(sample_N, scenario, estimator, outer_opt, parameter)
}

# ---- Aggregation: Power / PowerCN / PowerMinSuc per cell ------------------
# Same convention as compute_power_block() in scenario16console_output_report.R:
#
#   Power       = # exact_match          / N
#   PowerCN     = # exact_match & ok_cn  / # ok_cn        (cond_num_cor < cut)
#   PowerMinSuc = # exact_match & ok_min / # ok_min       (converged)
#
compute_bench_power <- function(diag_long, cn_cor_cut = 1000) {
  diag_long |>
    dplyr::group_by(sample_N, scenario, estimator, outer_opt) |>
    dplyr::summarise(
      N              = dplyr::n(),
      n_exact        = sum(exact_match, na.rm = TRUE),
      n_ok_cn        = sum(cond_num_cor < cn_cor_cut & !is.na(cond_num_cor)),
      n_ok_min       = sum(converged, na.rm = TRUE),
      n_exact_cn     = sum(exact_match &
                             (cond_num_cor < cn_cor_cut &
                                !is.na(cond_num_cor)), na.rm = TRUE),
      n_exact_min    = sum(exact_match & converged, na.rm = TRUE),
      Power          = n_exact     / N,
      PowerCN        = ifelse(n_ok_cn  > 0, n_exact_cn  / n_ok_cn,  NA_real_),
      PowerMinSuc    = ifelse(n_ok_min > 0, n_exact_min / n_ok_min, NA_real_),
      .groups        = "drop"
    ) |>
    dplyr::arrange(sample_N, scenario, estimator, outer_opt)
}

# ---- Aggregation: relative power per (cell, k) ----------------------------
# For each cell: fraction of datasets recovering >= k of the true covariates,
# for k = 1..n_true(scenario).
compute_bench_relpower <- function(diag_long, true_rel_by_scn) {
  n_true_by_scn <- vapply(true_rel_by_scn, nrow, integer(1))
  diag_long |>
    dplyr::group_by(sample_N, scenario, estimator, outer_opt) |>
    dplyr::group_modify(function(df, key) {
      n_true <- as.integer(n_true_by_scn[as.character(key$scenario)] %||% 0L)
      if (n_true == 0L)
        return(tibble::tibble(k = integer(), n_at_least_k = integer(),
                              fraction_at_least_k = numeric(),
                              n_true = integer()))
      tibble::tibble(
        k                   = seq_len(n_true),
        n_at_least_k        = vapply(seq_len(n_true), function(k)
          sum(df$n_true_hit >= k, na.rm = TRUE), integer(1)),
        fraction_at_least_k = vapply(seq_len(n_true), function(k)
          mean(df$n_true_hit >= k, na.rm = TRUE), numeric(1)),
        n_true              = n_true
      )
    }) |>
    dplyr::ungroup() |>
    dplyr::arrange(sample_N, scenario, estimator, outer_opt, k)
}

# ---- Top-level driver ------------------------------------------------------
aggregate_bench_run <- function(root          = "output",
                                out_dir       = file.path(root, "scm_bench_aggregated"),
                                true_params   = NULL,
                                cn_cor_cut    = 1000,
                                write_outputs = TRUE,
                                verbose       = TRUE) {
  loaded <- load_all_bench_results(root = root, true_params = true_params,
                                   verbose = verbose)
  if (!nrow(loaded$diag_long)) return(loaded)

  diag_rates   <- compute_bench_diag_rates(loaded$diag_long)
  per_dataset  <- compute_bench_per_dataset(loaded$diag_long, loaded$rse_long)
  rmse_all     <- compute_bench_rmse(loaded$rse_long,  status_filter = "all")
  rmse_success <- compute_bench_rmse(loaded$rse_long, loaded$diag_long,
                                     status_filter = "success_only")
  rmse_cond    <- compute_bench_rmse_cond(loaded$rse_long, loaded$diag_long)
  power        <- compute_bench_power(loaded$diag_long, cn_cor_cut = cn_cor_cut)
  relpower     <- compute_bench_relpower(loaded$diag_long, loaded$true_rel_by_scn)

  if (write_outputs) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(loaded$index,     file.path(out_dir, "bench_file_index.csv"))
    readr::write_csv(loaded$diag_long, file.path(out_dir, "bench_diag_long.csv"))
    readr::write_csv(loaded$rse_long,  file.path(out_dir, "bench_rse_long.csv"))
    readr::write_csv(per_dataset,      file.path(out_dir, "bench_per_dataset.csv"))
    readr::write_csv(diag_rates,       file.path(out_dir, "bench_diag_rates.csv"))
    readr::write_csv(rmse_all,         file.path(out_dir, "bench_rmse_all.csv"))
    readr::write_csv(rmse_success,     file.path(out_dir, "bench_rmse_success.csv"))
    readr::write_csv(rmse_cond,        file.path(out_dir, "bench_rmse_cond.csv"))
    readr::write_csv(power,            file.path(out_dir, "bench_power.csv"))
    readr::write_csv(relpower,         file.path(out_dir, "bench_relpower.csv"))
    saveRDS(list(
      index         = loaded$index,
      diag_long     = loaded$diag_long,
      rse_long      = loaded$rse_long,
      per_dataset   = per_dataset,
      diag_rates    = diag_rates,
      rmse_all      = rmse_all,
      rmse_success  = rmse_success,
      rmse_cond     = rmse_cond,
      power         = power,
      relpower      = relpower,
      true_rel_by_scn = loaded$true_rel_by_scn,
      cn_cor_cut    = cn_cor_cut,
      created_at    = Sys.time()
    ), file.path(out_dir, "bench_aggregated.rds"))
    if (verbose)
      message(sprintf("Wrote 10 CSVs + 1 RDS to %s/", out_dir))
  }

  invisible(list(
    index        = loaded$index,
    diag_long    = loaded$diag_long,
    rse_long     = loaded$rse_long,
    per_dataset  = per_dataset,
    diag_rates   = diag_rates,
    rmse_all     = rmse_all,
    rmse_success = rmse_success,
    rmse_cond    = rmse_cond,
    power        = power,
    relpower     = relpower,
    true_rel_by_scn = loaded$true_rel_by_scn
  ))
}

# ---- CLI dispatcher --------------------------------------------------------
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0L) {
  argv <- commandArgs(trailingOnly = TRUE)
  root <- "output"; tp_path <- "Inputdataset/true_params_long.rds"
  i <- 1L
  while (i <= length(argv)) {
    if (argv[i] == "--root" && i + 1L <= length(argv)) {
      root <- argv[i + 1L]; i <- i + 2L
    } else if (argv[i] == "--true_params" && i + 1L <= length(argv)) {
      tp_path <- argv[i + 1L]; i <- i + 2L
    } else { i <- i + 1L }
  }
  tp <- if (file.exists(tp_path)) readRDS(tp_path) else NULL
  aggregate_bench_run(root = root, true_params = tp)
}
