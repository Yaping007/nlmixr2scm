# =============================================================================
# aggregate_psn_scm.R  --  PsN/NONMEM SCM aggregator (schema 2.1 records)
# -----------------------------------------------------------------------------
# The shared aggregator `script/aggregate_scm_estimator2.1.R` was written for the
# nlmixr2 `runSCM` record schema, where per-fit diagnostics live in NESTED
# sub-lists (`r$diag`, `r$diag_t3`) and structural typical values are centred on
# the covariate MEDIAN (so `.backtransform_intercepts()` re-anchors them to
# 70 / 95 post-hoc).
#
# The PsN parser (`parse_psn_scm.R`) writes schema 2.1 records with a DIFFERENT
# shape:
#   * diagnostics are TOP-LEVEL fields (`r$objf`, `r$converged`,
#     `r$cond_num_cor`, `r$min_success`, `r$cov_done`, `r$cn_below_cutoff`);
#   * NONMEM already fits TVCL / TVV2 at the fixed references BW=70 / CrCL=95
#     via the `[code]` section, so NO post-hoc re-centring may be applied
#     (doing so would double-correct the structural error).
#
# Feeding PsN records straight into the runSCM unpacker therefore collapses every
# diagnostic to NA (cond# all-NA -> PowerCN / PowerMinSuc undefined; converged
# NA -> scm_diag_rates.csv full of NA).  This file overrides ONLY the unpack
# layer, reusing every downstream `compute_scm_*` aggregation unchanged, so the
# runSCM aggregator stays untouched for the nlmixr2 workflow.
#
# Usage:
#   source("script/psn_scm/aggregate_psn_scm.R")
#   agg <- aggregate_psn_scm_run(
#     root    = "output/psn_scm_full0727",
#     sub     = "ResforAggregation",
#     out_dir = "output/psn_scm_full0727_aggregated"
#   )
#
# CLI:
#   Rscript script/psn_scm/aggregate_psn_scm.R \
#     --root output/psn_scm_full0727 --sub ResforAggregation \
#     --out_dir output/psn_scm_full0727_aggregated
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(readr)
})

# ---- locate + source the shared aggregator (reuse its compute_* + discovery) -
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L ||
                             (length(a) == 1L && is.na(a))) b else a

.find_shared_aggregator <- function() {
  cands <- c(
    file.path("script", "aggregate_scm_estimator2.1.R"),
    file.path(dirname(getwd()), "script", "aggregate_scm_estimator2.1.R"),
    "aggregate_scm_estimator2.1.R"
  )
  hit <- cands[file.exists(cands)]
  if (!length(hit))
    stop("Cannot locate aggregate_scm_estimator2.1.R; run from the repo root.")
  hit[1]
}
source(.find_shared_aggregator(), local = FALSE)

# =============================================================================
# OVERRIDES  (defined AFTER source() so global-name lookup inside the reused
# `load_one_scm_rds` / `aggregate_scm_bench_run` resolves to THESE versions)
# =============================================================================

# ---- diagnostics from the PsN top-level schema ------------------------------
.unpack_diag <- function(r, meta, hit) {
  num  <- function(x) suppressWarnings(as.numeric(x %||% NA_real_))
  logi <- function(x) suppressWarnings(as.logical(x %||% NA))

  objf_val  <- num(r$objf)
  cn_val    <- num(r$cond_num_cor)
  converged <- logi(r$converged)
  min_succ  <- logi(r$min_success)
  cov_done  <- logi(r$cov_done)
  # `cn_below_cutoff` was computed by the parser at CN_STRICT_CUTOFF; recompute
  # against the aggregation cutoff (1000) for internal consistency, falling back
  # to the stored flag when the cond# is missing.
  cnb <- if (is.finite(cn_val)) cn_val <= 1000 else logi(r$cn_below_cutoff)

  tibble::tibble(
    sample_N        = r$sample_N    %||% meta$sample_N,
    scenario        = r$scenario_id %||% meta$scenario,
    structure       = r$model_type  %||% r$structure %||% meta$structure,
    estimator       = r$estimator   %||% meta$estimator,
    outer_opt       = r$outer_opt   %||% meta$outer_opt,
    dataset_id      = r$dataset_id  %||% meta$dataset_id,
    # PsN records carry no free-text status; derive one so n_ok_status is usable
    status          = dplyr::case_when(
                        isTRUE(converged) ~ "ok",
                        is.na(converged)  ~ NA_character_,
                        TRUE              ~ "not_converged"),
    converged       = converged,
    cn_below_cutoff = cnb,
    # PsN has no separate "cov requested" vs "cov done"; the single cov_done flag
    # (refit covariance step succeeded) serves both CovStep_pct and CovDone_pct.
    cov_ok          = cov_done,
    cov_done        = cov_done,
    objf            = objf_val,
    cond_num_cor    = cn_val,
    cond_num_cor_source = "psn_lst",
    convergence_code    = NA_integer_,
    min_suc         = min_succ,
    cov_step        = cov_done,
    est_bnd         = logi(r$rounding_error),   # rounding == boundary proxy
    phys_bnd        = NA,
    est_family      = "nonmem",
    cov_source      = "psn_cov",
    n_selected      = if (!is.null(r$scm$selected)) nrow(r$scm$selected) else 0L,
    n_true          = nrow(.true_set_for(r$scenario_id %||% meta$scenario)),
    n_true_hit      = hit$n_true_hit,
    n_false_pos     = hit$n_false_pos,
    exact_match     = hit$exact_match,
    # PsN timing is a single total block (base+scm[+refit]); phase split N/A
    wall_base_sec   = NA_real_,
    wall_scm_sec    = NA_real_,
    wall_refit_sec  = NA_real_,
    wall_total_sec  = num(r$runtime$total_sec %||% r$wall_total_sec),
    cpu_total_sec   = num(r$cpu$cpu_sec %||% r$cpu_total_sec),
    hog_factor      = num(r$cpu$hog_factor %||% r$hog_factor),
    scm_workers     = NA_integer_,   # not applicable to PsN/LSF fits
    rx_threads      = NA_integer_,
    message         = as.character(r$message %||% NA_character_)
  )
}

# ---- relative error: PsN rel_err is ALREADY correctly centred (no backtrans) --
.unpack_rse <- function(r, meta) {
  rel <- r$rel_err
  if (is.null(rel) || !nrow(rel)) return(tibble::tibble())
  tibble::tibble(
    sample_N    = r$sample_N    %||% meta$sample_N,
    scenario    = r$scenario_id %||% meta$scenario,
    structure   = r$model_type  %||% r$structure %||% meta$structure,
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
  ) |>
    dplyr::filter(!is.na(parameter))
}

# ---- covsel: reuse the shared canonicalised full-join (via .canon_pairs) -----
# `.unpack_covsel` from the shared file already reads `r$scm$selected` and
# `.true_set_for(r$scenario_id)`, both present in PsN records, and now matches
# case-insensitively (CrCL == CRCL) after the .canon_pairs toupper fix.

# ---- top-level PsN driver ---------------------------------------------------
aggregate_psn_scm_run <- function(root          = "output/psn_scm_full0727",
                                  sub           = "ResforAggregation",
                                  out_dir       = paste0(root, "_aggregated"),
                                  cn_cor_cut    = 1000,
                                  write_outputs = TRUE,
                                  verbose       = TRUE) {
  aggregate_scm_bench_run(
    root          = root,
    sub           = sub,
    out_dir       = out_dir,
    cn_cor_cut    = cn_cor_cut,
    write_outputs = write_outputs,
    verbose       = verbose
  )
}

# ---- CLI --------------------------------------------------------------------
if (identical(environment(), globalenv()) &&
    !interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  getarg <- function(flag, default = NULL) {
    i <- which(args == flag)
    if (length(i) && i < length(args)) args[i + 1L] else default
  }
  root    <- getarg("--root",    "output/psn_scm_full0727")
  sub     <- getarg("--sub",     "ResforAggregation")
  out_dir <- getarg("--out_dir", paste0(root, "_aggregated"))
  cutoff  <- as.numeric(getarg("--cn_cor_cut", "1000"))
  message(sprintf("[aggregate_psn_scm] root=%s sub=%s out=%s cutoff=%g",
                  root, sub, out_dir, cutoff))
  invisible(aggregate_psn_scm_run(root = root, sub = sub, out_dir = out_dir,
                                  cn_cor_cut = cutoff))
}
