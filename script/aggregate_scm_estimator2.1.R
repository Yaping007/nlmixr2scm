# ==============================================================================
# aggregate_scm_estimator2.1.R
# ------------------------------------------------------------------------------
# Aggregate the SCHEMA-2.1 SCM estimator x optimizer benchmark produced by
# PerformanceEvaluation_scm_bench.R (classical runSCM forward/backward LRT,
# two-tier control).  Sibling of aggregate_vae_covsel.R -- same operating-
# characteristics outputs, adapted to the SCM record shape.
#
# NOTE ON NAMING: the "2.1" suffix distinguishes this from the OLD reporting
# schema aggregators (e.g. aggregate_bench_refit_results.R).  It reads the flat
# schema-2.1 records (r$runtime / r$cpu / r$scm blocks) written by
# package_scm_schema21().
#
# Directory scanned:
#   output/scm_bench/N<N>/scn<SS>_<structure>/<est>_<opt>/res_ds<DDD>.rds
# (the extra <est>_<opt> level -- absent in the VAE tree -- is parsed into the
#  estimator/outer_opt grouping keys; the *.fit.rds sidecars are skipped.)
#
# Grouping cell = (sample_N, scenario, structure, estimator, outer_opt).
# The completed sweep is focei_bobyqa only, but the keys are kept general so
# additional estimator x optimizer cells aggregate without change.
#
# RECORD SHAPE (schema 2.1, from package_scm_schema21):
#   r$diag / r$diag_t3 : convergence + stability flags (same as VAE / refit)
#   r$rel_err          : tibble(parameter, true_value, estimate, rel_err, ...)
#   r$scm$selected     : tibble(var, covar, shape, theta_name, estimate)
#                        -- NOTE: SCM DOES search shape (unlike VAE); we still
#                           score selection on (var, covar) for comparability
#                           with the DGP truth and the VAE pilot.
#   r$scm$step_hist    : runSCM summaryTable (not aggregated here)
#   r$scm$cov_done     : logical, covariance of the refit winner succeeded
#   r$runtime          : list(base_sec, scm_sec, refit_sec, total_sec)  WALL
#   r$cpu              : list(base_sec, scm_sec, refit_sec, total_sec,
#                             hog_factor)                                CPU
#   r$scm_workers, r$rx_threads, r$sample_N, r$scenario_id, r$dataset_id, ...
#
# TRUE COVARIATE SET -- DERIVED, not stored.
#   Unlike the VAE driver, the SCM driver does NOT stamp r$covsel$true_set.
#   We reconstruct it per scenario from the PsN_scenarios indicators (embedded
#   below, verbatim from refit_helpers.R):
#     I_BW_CL   -> cl~BW (power)    I_CRCL_CL -> cl~CrCL (power)
#     I_BW_VC   -> vc~BW (power)    I_SEX_VC  -> vc~SEX  (cat)
#
# CONVERGENCE.
#   `converged` is RECOMPUTED from numerical evidence (finite OFV + cov_ok +
#   finite cond_num_cor), matching aggregate_bench_refit_results.R, so a stale
#   saved flag (e.g. lbfgsb3c "code 8 false convergence") cannot bias rates.
#
# TIMING.
#   $runtime (WALL) is the usable timing block; scm_sec is the fork phase.
#   $cpu / hog_factor from proc.time() do NOT capture the fork workers' CPU
#   (children are not reaped), so in-record hog_factor is an under-report
#   (~0.03-0.25) -- for the honest CPU / parallelism story use LSF CPU via
#   script/aggregate_lsf_cpu.R.  We surface the WALL phase medians here and
#   carry the (caveated) cpu/hog columns through untouched.
#
# REFERENCE CAVEAT (estimation metrics only, identical to the VAE aggregator):
#   covariate_beta rel_err is centring-invariant (clean); structural_intercept
#   (TVCL/TVVc) is reference-dependent (interpret with care).
#
# OUTPUTS (to output/scm_bench_aggregated/):
#   scm_file_index.csv      one row per discovered RDS
#   scm_diag_long.csv       per fit: convergence, stability, timing, selection
#   scm_rse_long.csv        per (fit, parameter): rel_err + param_class
#   scm_covsel_long.csv     per (fit, var, covar): in_true / in_scm / verdict
#   scm_diag_rates.csv      per cell: %converged, CN, WALL/CPU timing, hog
#   scm_estim_all.csv       per (cell, parameter): MedRE / MARE / RMRSE (all)
#   scm_estim_success.csv   same, strict-converged fits only
#   scm_estim_cond.csv      same, exact-match fits only
#   scm_power.csv           per cell: Power / PowerCN / PowerMinSuc
#   scm_relpower.csv        per (cell, k): fraction recovering >= k true covs
#   scm_covsel_by_covar.csv per (cell, var, covar): detection rate / TP-FP-FN
#   scm_bench_aggregated.rds  bundle of every tibble + params + created_at
#
# Usage:
#   Rscript "script/aggregate_scm_estimator2.1.R" [--root output] \
#           [--sub scm_bench] [--out_dir <path>]
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# ---- Per-scenario TRUE covariate set (embedded PsN_scenarios) --------------
# Verbatim indicator table from refit_helpers.R so the aggregator is
# self-contained (no heavy sourcing / side effects).
.PsN_scenarios <- data.frame(
  scenario  = 1:16,
  I_BW_CL   = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
  I_CRCL_CL = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
  I_BW_VC   = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
  I_SEX_VC  = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
)

.true_set_for <- function(scenario_id) {
  scn <- .PsN_scenarios[.PsN_scenarios$scenario == scenario_id, , drop = FALSE]
  if (nrow(scn) != 1L)
    return(tibble::tibble(var = character(), covar = character()))
  rows <- list()
  if (isTRUE(scn$I_BW_CL   == 1)) rows[[length(rows)+1L]] <- tibble::tibble(var="cl", covar="BW")
  if (isTRUE(scn$I_CRCL_CL == 1)) rows[[length(rows)+1L]] <- tibble::tibble(var="cl", covar="CrCL")
  if (isTRUE(scn$I_BW_VC   == 1)) rows[[length(rows)+1L]] <- tibble::tibble(var="vc", covar="BW")
  if (isTRUE(scn$I_SEX_VC  == 1)) rows[[length(rows)+1L]] <- tibble::tibble(var="vc", covar="SEX")
  if (length(rows) == 0L)
    return(tibble::tibble(var = character(), covar = character()))
  dplyr::bind_rows(rows)
}

# ---- Parameter classification (same rule as the VAE aggregator) -----------
.classify_param <- function(parameter) {
  dplyr::case_when(
    parameter %in% c("TVCL", "TVVc") ~ "structural_intercept",
    grepl("^TV", parameter) |
      parameter %in% c("var_CL", "var_Vc", "cov_VcCL", "ResErr") ~ "structural_other",
    TRUE ~ "covariate_beta"
  )
}

# ---- Selection scoring on (var, covar) ------------------------------------
.canon_pairs <- function(tbl) {
  if (is.null(tbl) || !nrow(tbl))
    return(tibble::tibble(var = character(), covar = character()))
  tbl |>
    dplyr::transmute(var = as.character(var), covar = as.character(covar)) |>
    dplyr::distinct()
}

match_selected_to_truth <- function(selected, true_set) {
  sel <- .canon_pairs(selected)
  tru <- .canon_pairs(true_set)
  if (!nrow(sel)) {
    return(list(n_true_hit = 0L, n_false_pos = 0L,
                exact_match = (nrow(tru) == 0L)))
  }
  true_hit  <- dplyr::inner_join(sel, tru, by = c("var", "covar"))
  false_pos <- dplyr::anti_join(sel, tru, by = c("var", "covar"))
  miss      <- dplyr::anti_join(tru, sel, by = c("var", "covar"))
  list(n_true_hit  = nrow(true_hit),
       n_false_pos = nrow(false_pos),
       exact_match = (nrow(false_pos) == 0L) && (nrow(miss) == 0L))
}

# ---- File discovery -------------------------------------------------------
discover_scm_files <- function(root = "output", sub = "scm_bench") {
  base <- file.path(root, sub)
  if (!dir.exists(base)) {
    warning("Directory not found: ", base); return(tibble::tibble())
  }
  rds <- list.files(base, pattern = "^res_ds\\d+\\.rds$",
                    recursive = TRUE, full.names = TRUE)
  rds <- rds[!grepl("\\.fit\\.rds$", rds)]     # drop the raw-fit sidecars
  if (!length(rds)) {
    warning("No res_ds*.rds under ", base); return(tibble::tibble())
  }
  parts    <- strsplit(gsub("\\\\", "/", rds), "/", fixed = TRUE)
  est_dir  <- vapply(parts, function(p) p[length(p) - 1L], character(1)) # <est>_<opt>
  cell_dir <- vapply(parts, function(p) p[length(p) - 2L], character(1)) # scn<SS>_<struct>
  n_dir    <- vapply(parts, function(p) p[length(p) - 3L], character(1)) # N<N>
  err_txt  <- sub("\\.rds$", "_ERROR.txt", rds)

  tibble::tibble(
    file_path  = rds,
    n_dir      = n_dir,
    cell_dir   = cell_dir,
    est_dir    = est_dir,
    dataset_id = as.integer(sub("^res_ds(\\d+)\\.rds$", "\\1", basename(rds))),
    has_error  = file.exists(err_txt)
  ) |>
    dplyr::mutate(
      sample_N  = suppressWarnings(as.integer(sub("^N", "", n_dir))),
      scenario  = suppressWarnings(as.integer(sub("^scn0*(\\d+)_.*$", "\\1", cell_dir))),
      structure = sub("^scn\\d+_", "", cell_dir),
      estimator = sub("_[^_]+$", "", est_dir),   # split on LAST underscore
      outer_opt = sub("^.*_",    "", est_dir)
    ) |>
    dplyr::select(sample_N, scenario, structure, estimator, outer_opt,
                  dataset_id, has_error, file_path) |>
    dplyr::arrange(sample_N, scenario, structure, estimator, outer_opt, dataset_id)
}

# ---- Unpack one RDS -------------------------------------------------------
.unpack_diag <- function(r, meta, hit) {
  d  <- r$diag    %||% list()
  d3 <- r$diag_t3 %||% list()
  pick <- function(k) d3[[k]] %||% d[[k]]
  cn_val <- suppressWarnings(as.numeric(pick("cond_num_cor") %||% NA_real_))
  rt  <- r$runtime %||% list()
  cpu <- r$cpu     %||% list()
  # Recompute `converged` from numerical evidence (finite OFV + cov_ok +
  # finite CN); fall back to the saved flag if inputs are missing.
  objf_val <- as.numeric(pick("objf") %||% NA_real_)
  cov_ok_v <- as.logical(pick("cov_ok") %||% NA)
  converged_recomp <- {
    objf_ok <- is.finite(objf_val)
    if (is.na(objf_ok) || is.na(cov_ok_v)) as.logical(pick("converged") %||% NA)
    else objf_ok && isTRUE(cov_ok_v) && is.finite(cn_val)
  }
  tibble::tibble(
    sample_N        = r$sample_N    %||% meta$sample_N,
    scenario        = r$scenario_id %||% meta$scenario,
    structure       = r$structure   %||% meta$structure,
    estimator       = r$estimator   %||% meta$estimator,
    outer_opt       = r$outer_opt   %||% meta$outer_opt,
    dataset_id      = r$dataset_id  %||% meta$dataset_id,
    status          = as.character(r$status %||% NA_character_),
    converged       = as.logical(converged_recomp),
    cn_below_cutoff = isTRUE(is.finite(cn_val) && cn_val <= 1000),
    cov_ok          = cov_ok_v,
    cov_done        = as.logical(r$scm$cov_done %||% NA),
    objf            = objf_val,
    cond_num_cor    = cn_val,
    cond_num_cor_source = as.character(pick("cond_num_cor_source") %||% NA_character_),
    convergence_code    = suppressWarnings(as.integer(pick("convergence_code") %||% NA_integer_)),
    min_suc         = as.logical(d3$min_suc  %||% NA),
    cov_step        = as.logical(d3$cov_step %||% NA),
    est_bnd         = as.logical(d3$est_bnd  %||% NA),
    phys_bnd        = as.logical(d3$phys_bnd %||% NA),
    est_family      = as.character(d3$est_family %||% NA_character_),
    cov_source      = as.character(d3$cov_source %||% NA_character_),
    n_selected      = if (!is.null(r$scm$selected)) nrow(r$scm$selected) else 0L,
    n_true          = nrow(.true_set_for(r$scenario_id %||% meta$scenario)),
    n_true_hit      = hit$n_true_hit,
    n_false_pos     = hit$n_false_pos,
    exact_match     = hit$exact_match,
    # WALL-clock phase timing (usable)
    wall_base_sec   = as.numeric(rt$base_sec  %||% NA_real_),
    wall_scm_sec    = as.numeric(rt$scm_sec   %||% NA_real_),
    wall_refit_sec  = as.numeric(rt$refit_sec %||% NA_real_),
    wall_total_sec  = as.numeric(rt$total_sec %||% r$wall_total_sec %||% NA_real_),
    # CPU timing (CAVEATED: fork workers not reaped -> under-report; use LSF)
    cpu_total_sec   = as.numeric(cpu$total_sec %||% r$cpu_total_sec %||% NA_real_),
    hog_factor      = as.numeric(cpu$hog_factor %||% r$hog_factor  %||% NA_real_),
    scm_workers     = suppressWarnings(as.integer(r$scm_workers %||% NA_integer_)),
    rx_threads      = suppressWarnings(as.integer(r$rx_threads  %||% NA_integer_)),
    message         = as.character(pick("message") %||% NA_character_)
  )
}

.unpack_rse <- function(r, meta) {
  rel <- r$rel_err
  if (is.null(rel) || !nrow(rel)) return(tibble::tibble())
  tibble::tibble(
    sample_N    = r$sample_N    %||% meta$sample_N,
    scenario    = r$scenario_id %||% meta$scenario,
    structure   = r$structure   %||% meta$structure,
    estimator   = r$estimator   %||% meta$estimator,
    outer_opt   = r$outer_opt   %||% meta$outer_opt,
    dataset_id  = r$dataset_id  %||% meta$dataset_id,
    parameter   = as.character(rel$parameter),
    param_class = .classify_param(as.character(rel$parameter)),
    true_value  = as.numeric(rel$true_value),
    estimate    = as.numeric(rel$estimate),
    rel_err     = as.numeric(rel$rel_err),
    rel_err_pct = as.numeric(rel$rel_err_pct %||% (rel$rel_err * 100)),
    sq_rel      = as.numeric(rel$rel_err)^2,
    abs_rel     = abs(as.numeric(rel$rel_err))
  )
}

.unpack_covsel <- function(r, meta) {
  true_set <- .canon_pairs(.true_set_for(r$scenario_id %||% meta$scenario))
  selected <- .canon_pairs(r$scm$selected)
  if (!nrow(true_set) && !nrow(selected)) return(tibble::tibble())
  cmp <- dplyr::full_join(
    dplyr::mutate(true_set, in_true = TRUE),
    dplyr::mutate(selected, in_scm  = TRUE),
    by = c("var", "covar")
  ) |>
    dplyr::mutate(
      in_true = tidyr::replace_na(in_true, FALSE),
      in_scm  = tidyr::replace_na(in_scm,  FALSE),
      verdict = dplyr::case_when(
        in_true &  in_scm ~ "TP",
        in_true & !in_scm ~ "FN",
        !in_true &  in_scm ~ "FP"
      )
    )
  tibble::tibble(
    sample_N   = r$sample_N    %||% meta$sample_N,
    scenario   = r$scenario_id %||% meta$scenario,
    structure  = r$structure   %||% meta$structure,
    estimator  = r$estimator   %||% meta$estimator,
    outer_opt  = r$outer_opt   %||% meta$outer_opt,
    dataset_id = r$dataset_id  %||% meta$dataset_id,
    var        = cmp$var,
    covar      = cmp$covar,
    in_true    = cmp$in_true,
    in_scm     = cmp$in_scm,
    verdict    = cmp$verdict
  )
}

load_one_scm_rds <- function(path, meta) {
  r <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(r)) {
    message(sprintf("[skip] cannot read %s", path))
    return(list(diag = tibble::tibble(), rse = tibble::tibble(),
                covsel = tibble::tibble()))
  }
  true_set <- .true_set_for(r$scenario_id %||% meta$scenario)
  hit <- match_selected_to_truth(r$scm$selected, true_set)
  list(diag   = .unpack_diag(r, meta, hit),
       rse    = .unpack_rse(r, meta),
       covsel = .unpack_covsel(r, meta))
}

# ---- Bulk loader ----------------------------------------------------------
load_all_scm_results <- function(root = "output", sub = "scm_bench",
                                 verbose = TRUE) {
  idx <- discover_scm_files(root, sub = sub)
  if (!nrow(idx)) {
    warning("No SCM bench RDS files found under ", file.path(root, sub))
    return(list(index = idx, diag_long = tibble::tibble(),
                rse_long = tibble::tibble(), covsel_long = tibble::tibble()))
  }
  if (verbose) message(sprintf(
    "Loading %d RDS across %d cells (N x scn x structure x est x opt)%s",
    nrow(idx),
    dplyr::n_distinct(idx[, c("sample_N", "scenario", "structure",
                              "estimator", "outer_opt")]),
    if (any(idx$has_error)) sprintf(" [%d ERROR sidecar(s) present]",
                                    sum(idx$has_error)) else ""))

  loaded <- purrr::map(seq_len(nrow(idx)), function(i) {
    load_one_scm_rds(idx$file_path[i], meta = as.list(idx[i, ]))
  })
  list(
    index       = idx,
    diag_long   = purrr::map_dfr(loaded, "diag"),
    rse_long    = purrr::map_dfr(loaded, "rse"),
    covsel_long = purrr::map_dfr(loaded, "covsel")
  )
}

# ---- Aggregation: convergence + stability + timing per cell ---------------
compute_scm_diag_rates <- function(diag_long) {
  diag_long |>
    dplyr::mutate(
      converged_strict = (converged %in% TRUE) &
                         (cn_below_cutoff %in% TRUE) &
                         !(est_bnd %in% TRUE)
    ) |>
    dplyr::group_by(sample_N, scenario, structure, estimator, outer_opt) |>
    dplyr::summarise(
      n_total             = dplyr::n(),
      n_ok_status         = sum(status == "ok", na.rm = TRUE),
      Converged_pct       = 100 * mean(converged, na.rm = TRUE),
      ConvergedStrict_pct = 100 * mean(converged_strict),
      CovStep_pct         = 100 * mean(cov_ok,    na.rm = TRUE),
      CovDone_pct         = 100 * mean(cov_done,  na.rm = TRUE),
      CNBelowCutoff_pct   = 100 * mean(cn_below_cutoff),
      MinSuc_pct          = 100 * mean(min_suc,   na.rm = TRUE),
      EstBnd_pct          = 100 * mean(est_bnd,   na.rm = TRUE),
      MedCN               = stats::median(cond_num_cor, na.rm = TRUE),
      MaxCN               = suppressWarnings(max(cond_num_cor, na.rm = TRUE)),
      P90CN               = as.numeric(stats::quantile(cond_num_cor, 0.9, na.rm = TRUE)),
      MedObjF             = stats::median(objf, na.rm = TRUE),
      MedNSelected        = stats::median(n_selected, na.rm = TRUE),
      # WALL-clock phase timing (the usable timing numbers)
      MedWallBase_sec     = stats::median(wall_base_sec,  na.rm = TRUE),
      MedWallScm_sec      = stats::median(wall_scm_sec,   na.rm = TRUE),
      MedWallRefit_sec    = stats::median(wall_refit_sec, na.rm = TRUE),
      MedWallTotal_sec    = stats::median(wall_total_sec, na.rm = TRUE),
      MedWallTotal_min    = stats::median(wall_total_sec, na.rm = TRUE) / 60,
      # CPU / hog carried through (CAVEATED: fork CPU not captured; use LSF)
      MedCpuTotal_sec     = stats::median(cpu_total_sec, na.rm = TRUE),
      MedHogFactor        = stats::median(hog_factor,    na.rm = TRUE),
      scm_workers         = dplyr::first(scm_workers),
      rx_threads          = dplyr::first(rx_threads),
      .groups             = "drop"
    ) |>
    dplyr::mutate(dplyr::across(where(is.numeric),
                                \(x) ifelse(is.finite(x), x, NA_real_))) |>
    dplyr::arrange(sample_N, scenario, structure, estimator, outer_opt)
}

# ---- Aggregation: estimation precision per (cell, parameter) --------------
compute_scm_estim <- function(rse_long, diag_long = NULL,
                              mode = c("all", "success", "cond")) {
  mode <- match.arg(mode)
  keys <- c("sample_N", "scenario", "structure", "estimator", "outer_opt",
            "dataset_id")
  df <- rse_long |>
    dplyr::filter(is.finite(true_value), true_value != 0, is.finite(rel_err))

  if (mode == "success") {
    if (is.null(diag_long)) stop("mode='success' requires diag_long")
    keep <- diag_long |>
      dplyr::mutate(
        converged_strict = (converged %in% TRUE) &
                           (cn_below_cutoff %in% TRUE) &
                           !(est_bnd %in% TRUE)
      ) |>
      dplyr::filter(converged_strict) |>
      dplyr::select(dplyr::all_of(keys)) |>
      dplyr::distinct()
    df <- dplyr::semi_join(df, keep, by = keys)
  } else if (mode == "cond") {
    if (is.null(diag_long)) stop("mode='cond' requires diag_long")
    keep <- diag_long |>
      dplyr::filter(exact_match %in% TRUE) |>
      dplyr::select(dplyr::all_of(keys)) |>
      dplyr::distinct()
    df <- dplyr::semi_join(df, keep, by = keys)
  }

  if (!nrow(df))
    return(tibble::tibble(sample_N = integer(), scenario = integer(),
                          structure = character(), estimator = character(),
                          outer_opt = character(), parameter = character(),
                          param_class = character(), true_value = numeric(),
                          n_used = integer(), MedRE_pct = numeric(),
                          MeanRE_pct = numeric(), MARE_pct = numeric(),
                          RMRSE_pct = numeric(), P90AbsRE_pct = numeric()))
  df |>
    dplyr::group_by(sample_N, scenario, structure, estimator, outer_opt,
                    parameter, param_class) |>
    dplyr::summarise(
      true_value   = dplyr::first(true_value),
      n_used       = dplyr::n(),
      MedRE_pct    = 100 * stats::median(rel_err),
      MeanRE_pct   = 100 * mean(rel_err),
      MARE_pct     = 100 * stats::median(abs_rel),
      RMRSE_pct    = 100 * sqrt(mean(sq_rel)),
      P90AbsRE_pct = 100 * as.numeric(stats::quantile(abs_rel, 0.9)),
      .groups      = "drop"
    ) |>
    dplyr::arrange(sample_N, scenario, structure, estimator, outer_opt,
                   param_class, parameter)
}

# ---- Aggregation: Power / PowerCN / PowerMinSuc per cell ------------------
compute_scm_power <- function(diag_long, cn_cor_cut = 1000) {
  diag_long |>
    dplyr::group_by(sample_N, scenario, structure, estimator, outer_opt) |>
    dplyr::summarise(
      N            = dplyr::n(),
      n_exact      = sum(exact_match, na.rm = TRUE),
      n_ok_cn      = sum(cond_num_cor < cn_cor_cut & !is.na(cond_num_cor)),
      n_ok_min     = sum(converged, na.rm = TRUE),
      n_exact_cn   = sum(exact_match &
                           (cond_num_cor < cn_cor_cut & !is.na(cond_num_cor)),
                         na.rm = TRUE),
      n_exact_min  = sum(exact_match & converged, na.rm = TRUE),
      Power        = n_exact / N,
      PowerCN      = ifelse(n_ok_cn  > 0, n_exact_cn  / n_ok_cn,  NA_real_),
      PowerMinSuc  = ifelse(n_ok_min > 0, n_exact_min / n_ok_min, NA_real_),
      .groups      = "drop"
    ) |>
    dplyr::arrange(sample_N, scenario, structure, estimator, outer_opt)
}

# ---- Aggregation: relative power per (cell, k) ----------------------------
compute_scm_relpower <- function(diag_long) {
  diag_long |>
    dplyr::group_by(sample_N, scenario, structure, estimator, outer_opt) |>
    dplyr::group_modify(function(df, key) {
      n_true <- as.integer(stats::median(df$n_true, na.rm = TRUE) %||% 0L)
      if (is.na(n_true) || n_true == 0L)
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
    dplyr::arrange(sample_N, scenario, structure, estimator, outer_opt, k)
}

# ---- Aggregation: per-covariate detection per cell ------------------------
compute_scm_covsel_by_covar <- function(covsel_long) {
  if (!nrow(covsel_long)) return(tibble::tibble())
  covsel_long |>
    dplyr::group_by(sample_N, scenario, structure, estimator, outer_opt,
                    var, covar) |>
    dplyr::summarise(
      n_datasets     = dplyr::n(),
      is_true        = any(in_true %in% TRUE),
      n_detected     = sum(in_scm %in% TRUE),
      detection_rate = mean(in_scm %in% TRUE),
      n_TP           = sum(verdict == "TP", na.rm = TRUE),
      n_FN           = sum(verdict == "FN", na.rm = TRUE),
      n_FP           = sum(verdict == "FP", na.rm = TRUE),
      .groups        = "drop"
    ) |>
    dplyr::arrange(sample_N, scenario, structure, estimator, outer_opt,
                   dplyr::desc(is_true), var, covar)
}

# ---- Top-level driver -----------------------------------------------------
aggregate_scm_bench_run <- function(root          = "output",
                                    sub           = "scm_bench",
                                    out_dir       = file.path(root, "scm_bench_aggregated"),
                                    cn_cor_cut    = 1000,
                                    write_outputs = TRUE,
                                    verbose       = TRUE) {
  loaded <- load_all_scm_results(root, sub = sub, verbose = verbose)
  if (!nrow(loaded$diag_long)) return(invisible(loaded))

  diag_rates    <- compute_scm_diag_rates(loaded$diag_long)
  estim_all     <- compute_scm_estim(loaded$rse_long, mode = "all")
  estim_success <- compute_scm_estim(loaded$rse_long, loaded$diag_long, mode = "success")
  estim_cond    <- compute_scm_estim(loaded$rse_long, loaded$diag_long, mode = "cond")
  power         <- compute_scm_power(loaded$diag_long, cn_cor_cut = cn_cor_cut)
  relpower      <- compute_scm_relpower(loaded$diag_long)
  covsel_by_cov <- compute_scm_covsel_by_covar(loaded$covsel_long)

  if (write_outputs) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(loaded$index,       file.path(out_dir, "scm_file_index.csv"))
    readr::write_csv(loaded$diag_long,   file.path(out_dir, "scm_diag_long.csv"))
    readr::write_csv(loaded$rse_long,    file.path(out_dir, "scm_rse_long.csv"))
    readr::write_csv(loaded$covsel_long, file.path(out_dir, "scm_covsel_long.csv"))
    readr::write_csv(diag_rates,         file.path(out_dir, "scm_diag_rates.csv"))
    readr::write_csv(estim_all,          file.path(out_dir, "scm_estim_all.csv"))
    readr::write_csv(estim_success,      file.path(out_dir, "scm_estim_success.csv"))
    readr::write_csv(estim_cond,         file.path(out_dir, "scm_estim_cond.csv"))
    readr::write_csv(power,              file.path(out_dir, "scm_power.csv"))
    readr::write_csv(relpower,           file.path(out_dir, "scm_relpower.csv"))
    readr::write_csv(covsel_by_cov,      file.path(out_dir, "scm_covsel_by_covar.csv"))
    saveRDS(list(
      index         = loaded$index,
      diag_long     = loaded$diag_long,
      rse_long      = loaded$rse_long,
      covsel_long   = loaded$covsel_long,
      diag_rates    = diag_rates,
      estim_all     = estim_all,
      estim_success = estim_success,
      estim_cond    = estim_cond,
      power         = power,
      relpower      = relpower,
      covsel_by_cov = covsel_by_cov,
      cn_cor_cut    = cn_cor_cut,
      created_at    = Sys.time()
    ), file.path(out_dir, "scm_bench_aggregated.rds"))
    if (verbose)
      message(sprintf("Wrote 11 CSVs + 1 RDS to %s/", out_dir))
  }

  invisible(list(
    index         = loaded$index,
    diag_long     = loaded$diag_long,
    rse_long      = loaded$rse_long,
    covsel_long   = loaded$covsel_long,
    diag_rates    = diag_rates,
    estim_all     = estim_all,
    estim_success = estim_success,
    estim_cond    = estim_cond,
    power         = power,
    relpower      = relpower,
    covsel_by_cov = covsel_by_cov
  ))
}

# ---- CLI dispatcher -------------------------------------------------------
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) > 0L) {
  argv <- commandArgs(trailingOnly = TRUE)
  root <- "output"; sub <- "scm_bench"; out_dir <- NULL
  i <- 1L
  while (i <= length(argv)) {
    if (argv[i] == "--root" && i + 1L <= length(argv)) {
      root <- argv[i + 1L]; i <- i + 2L
    } else if (argv[i] == "--sub" && i + 1L <= length(argv)) {
      sub <- argv[i + 1L]; i <- i + 2L
    } else if (argv[i] == "--out_dir" && i + 1L <= length(argv)) {
      out_dir <- argv[i + 1L]; i <- i + 2L
    } else { i <- i + 1L }
  }
  aggregate_scm_bench_run(
    root    = root,
    sub     = sub,
    out_dir = out_dir %||% file.path(root, "scm_bench_aggregated")
  )
}
