# ==============================================================================
# aggregate_vae_covsel.R
# ------------------------------------------------------------------------------
# Aggregate the one-run VAE covariate-selection pilot produced by
# vae_covsel_driver.R (est = "vae", covariateSelection = TRUE).  runSCM is NOT
# involved -- the stepwise search is internal to the VAE training loop.
#
# Directory scanned:
#   output/vae_covsel_pilot/N<N>/scn<SS>_<structure>/covsel/res_ds<DDD>.rds
#
# Grouping cell = (sample_N, scenario, structure).  `estimator` is always "vae"
# and `outer_opt` is NA, so they are carried through but not used as keys.
#
# RECORD SHAPE (schema 2.1, FLAT -- unlike scm_bench which nests under r$test$):
#   r$diag / r$diag_t3  : convergence + stability flags
#   r$rel_err           : tibble(parameter, true_value, estimate, abs_err,
#                                 rel_err, rel_err_pct)
#   r$covsel$selected   : tibble(var, covar, theta_name, estimate)  -- NO shape
#   r$covsel$true_set   : tibble(var, covar, shape)
#   r$fit_runtime_sec, r$status, r$sample_N, r$scenario_id, r$dataset_id, ...
#
# SELECTION SCORING -- on (var, covar) ONLY.
#   VAE has a single fixed covariate form (power for continuous, log-shift for
#   categorical); functional shape is NOT a VAE degree of freedom, and
#   r$covsel$selected carries no `shape` column.  We therefore score TP/FN/FP on
#   (var, covar), mirroring the driver's own verdict.  `true_set` is read from
#   each record (already derived per scenario from the PsN_scenarios indicators).
#
# REFERENCE HANDLING (affects estimation metrics only, not selection).
#   The data-generating model centres continuous covariates at FIXED references
#   (BW/70, CrCL/95).  VAE instead centres at the sample MEAN.  Consequences for
#   the `param_class` tag on rel_err:
#     * covariate_beta      -- the power/cat coefficient is centring-INVARIANT,
#                              so rel_err vs truth is clean and comparable across
#                              VAE / runSCM / truth.  THIS is the headline
#                              estimation metric.
#     * structural_intercept-- TVCL / TVVc are the typical value AT the reference
#                              subject, so they are reference-DEPENDENT.  We now
#                              BACK-TRANSFORM the VAE estimate onto the DGP's
#                              70/95 anchor (see .backtransform_intercepts) using
#                              each fit's own power betas + the per-dataset mean
#                              covariates, so their rel_err is now true bias.
#     * structural_other    -- TVQ / TVVp / TVKA / omegas / ResErr carry no
#                              covariate and are reference-free -> fully clean.
#
# OUTPUTS (to output/vae_covsel_aggregated/):
#   vae_file_index.csv      one row per discovered RDS
#   vae_diag_long.csv       per fit: convergence, stability, timing, selection
#   vae_rse_long.csv        per (fit, parameter): rel_err + param_class
#   vae_covsel_long.csv     per (fit, var, covar): in_true / in_vae / verdict
#   vae_diag_rates.csv      per cell: %converged, CN, timing
#   vae_estim_all.csv       per (cell, parameter): MedRE / MARE / RMRSE (all)
#   vae_estim_success.csv   same, restricted to strictly-converged fits
#   vae_estim_cond.csv      same, restricted to exact-match fits
#   vae_power.csv           per cell: Power / PowerCN / PowerMinSuc
#   vae_relpower.csv        per (cell, k): fraction recovering >= k true covs
#   vae_covsel_by_covar.csv per (cell, var, covar): detection rate / TP-FP-FN
#   vae_covsel_aggregated.rds  bundle of every tibble + params + created_at
#
# Usage:
#   Rscript script/aggregate_vae_covsel.R [--root output] [--out_dir <path>]
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(readr)
})

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

# ---- Parameter classification ---------------------------------------------
# structural_intercept: typical values that carry covariate terms (reference
#   dependent).  structural_other: covariate-free structural / variance / error
#   terms.  covariate_beta: everything else (the promoted covariate coefficients
#   CLBW / CLcrCL / VcBW / VcSEX / CLBMI / VcBMI / VcCrCL / VcRACE / ...).
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

# ---- Backward-compat recovery of $covsel$selected -------------------------
# The updated nlmixr2est (>= 7.0.0) makes `fit$theta` return NULL, so the driver
# stored `covsel$selected = NULL` even though VAE DID promote covariates.  The
# promoted beta_<PARAM>_<COV> coefficients survive in the record's parFixed
# table (rownames = parameter, column "Estimate").  When `covsel$selected` is
# missing/empty we rebuild it from parFixed so old runs re-score correctly --
# no re-fit needed.  Mirrors the driver's own parser.
.recover_selected_from_parfixed <- function(r) {
  pf <- r$parFixed
  if (!is.data.frame(pf) || !nrow(pf) || is.null(rownames(pf))) return(NULL)
  nms  <- rownames(pf)
  hits <- grep("^beta_", nms, value = TRUE)
  if (!length(hits)) return(NULL)
  m     <- regmatches(hits, regexec("^beta_(lTV[[:alnum:]]+)_(.+)$", hits))
  param <- vapply(m, function(x) if (length(x) == 3L) x[2] else NA_character_, character(1))
  covar <- vapply(m, function(x) if (length(x) == 3L) x[3] else NA_character_, character(1))
  var   <- unname(c(lTVCL = "cl", lTVVc = "vc")[param])
  covar <- ifelse(toupper(covar) == "CRCL", "CrCL", covar)
  est   <- if ("Estimate" %in% colnames(pf)) {
    suppressWarnings(as.numeric(pf[hits, "Estimate"]))
  } else rep(NA_real_, length(hits))
  keep <- !is.na(var)
  if (!any(keep)) return(NULL)
  tibble::tibble(var = var[keep], covar = covar[keep],
                 theta_name = hits[keep], estimate = est[keep])
}

# Ensure r$covsel$selected is populated: keep an existing non-empty value,
# otherwise reconstruct from parFixed.  Returns the (possibly patched) record.
.ensure_selected <- function(r) {
  sel <- r$covsel$selected
  if (is.null(sel) || (is.data.frame(sel) && !nrow(sel))) {
    rec <- .recover_selected_from_parfixed(r)
    if (!is.null(rec)) {
      if (is.null(r$covsel)) r$covsel <- list()
      r$covsel$selected <- rec
    }
  }
  r
}

# ---- Backward-compat backfill of covariate-beta rel_err -------------------
# The broken driver (fit$theta -> NULL) also left the covariate_beta rows of
# rel_err with NA estimate/rel_err, because that backfill was keyed off the
# NULL `selected`.  Given a recovered `selected` (var/covar/estimate) we map to
# the FOCEi-style rel_err labels and fill estimate + abs_err + rel_err, exactly
# as the driver does.  Only rows with a non-NA true_value and currently-NA
# estimate are touched; structural rows are never altered.
.vae_true_param <- function(var, covar) {
  key <- paste(var, toupper(covar), sep = "|")
  unname(c(
    "cl|BW"   = "CLBW",
    "cl|CRCL" = "CLcrCL",
    "vc|BW"   = "VcBW",
    "vc|SEX"  = "VcSEX"
  )[key])
}

.ensure_relerr_backfill <- function(r) {
  rel <- r$rel_err
  sel <- r$covsel$selected
  if (is.null(rel) || !nrow(rel) || is.null(sel) || !nrow(sel)) return(r)
  # nothing to do if no covariate_beta rows are missing an estimate
  needs <- rel$parameter %in% c("CLBW", "CLcrCL", "VcBW", "VcSEX") &
           is.na(rel$estimate) & !is.na(rel$true_value)
  if (!any(needs)) return(r)

  sel_named <- sel |>
    dplyr::mutate(parameter = .vae_true_param(var, covar)) |>
    dplyr::filter(!is.na(parameter)) |>
    dplyr::select(parameter, vae_estimate = estimate) |>
    dplyr::distinct(parameter, .keep_all = TRUE)
  if (!nrow(sel_named)) return(r)

  r$rel_err <- rel |>
    dplyr::left_join(sel_named, by = "parameter") |>
    dplyr::mutate(
      estimate    = dplyr::coalesce(estimate, vae_estimate),
      abs_err     = ifelse(is.na(estimate) | is.na(true_value),
                           NA_real_, estimate - true_value),
      rel_err     = ifelse(is.na(estimate) | is.na(true_value) | true_value == 0,
                           NA_real_, (estimate - true_value) / true_value),
      rel_err_pct = rel_err * 100
    ) |>
    dplyr::select(-vae_estimate)
  r
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

# ---- Structural-intercept reference back-transform ------------------------
# The data-generating model defines TVCL / TVVc at FIXED references BW = 70,
# CrCL = 95.  VAE centres continuous covariates at the sample MEAN, so its
# fitted TVCL / TVVc are the typical values at mean(BW) / mean(CrCL), NOT at the
# 70 / 95 anchor -- making their rel_err a REFERENCE ARTEFACT rather than bias.
# We rescale the intercept ESTIMATE back onto the 70 / 95 anchor using the fit's
# OWN power betas (which are themselves reference-invariant):
#   TVCL@70/95 = TVCL_hat * (70/mean(BW))^b_CLBW * (95/mean(CrCL))^b_CLcrCL
#   TVVc@70    = TVVc_hat * (70/mean(BW))^b_VcBW
# (categorical SEX carries no continuous reference, so TVVc uses BW only).
# Covariate betas and covariate-free structural terms are NEVER touched.
# Empirically verified on scn16 / N300: recovers TVCL = 0.594 vs truth 0.6 with
# the MEAN reference (median gives 0.628; uncorrected 0.698).
.COV_REF_BW   <- 70       # DGP BW_REF
.COV_REF_CRCL <- 95       # DGP CRCL_REF
.COV_REF_STAT <- mean     # VAE centres continuous covariates at the sample MEAN

# cache of per-(sample_N, scenario, dataset_id) covariate reference values
.cov_ref_cache <- new.env(parent = emptyenv())

.cov_ref_values <- function(sample_N, scenario, dataset_id,
                            input_root = "Inputdataset") {
  if (is.na(sample_N) || is.na(scenario) || is.na(dataset_id)) return(NULL)
  key <- paste(sample_N, scenario, dataset_id, sep = "|")
  hit <- .cov_ref_cache[[key]]
  if (!is.null(hit)) return(hit$val)          # cached (may be NULL result)
  path <- file.path(input_root, sprintf("sim_obs_N%s", sample_N),
                    sprintf("sim_obs_scenario_%s.rds", scenario))
  out <- NULL
  if (file.exists(path)) {
    d <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(d) && all(c("DATASET", "SUBJECT", "BW", "CrCL") %in% names(d))) {
      d1 <- d[d$DATASET == dataset_id, , drop = FALSE]
      d1 <- d1[!duplicated(d1$SUBJECT), , drop = FALSE]
      if (nrow(d1)) {
        out <- list(BW   = .COV_REF_STAT(d1$BW,   na.rm = TRUE),
                    CrCL = .COV_REF_STAT(d1$CrCL, na.rm = TRUE))
      }
    }
  }
  .cov_ref_cache[[key]] <- list(val = out)     # memoise (including NULL)
  out
}

# Rescale TVCL / TVVc estimate onto the fixed 70 / 95 reference (VAE only).
.backtransform_intercepts <- function(r, meta) {
  estimator <- r$estimator %||% meta$estimator %||% "vae"
  if (!identical(as.character(estimator), "vae")) return(r)
  rel <- r$rel_err
  if (is.null(rel) || !nrow(rel) || is.null(rel$estimate)) return(r)

  gete <- function(p) {
    v <- suppressWarnings(as.numeric(rel$estimate[rel$parameter == p]))
    if (length(v) == 1L) v else NA_real_
  }
  b_clbw <- gete("CLBW"); b_clcr <- gete("CLcrCL"); b_vcbw <- gete("VcBW")

  ref <- .cov_ref_values(r$sample_N    %||% meta$sample_N,
                         r$scenario_id %||% meta$scenario,
                         r$dataset_id  %||% meta$dataset_id)
  if (is.null(ref)) return(r)                  # sim file missing -> leave as-is
  fBW   <- .COV_REF_BW   / ref$BW
  fCRCL <- .COV_REF_CRCL / ref$CrCL

  i_cl <- which(rel$parameter == "TVCL")
  if (length(i_cl) == 1L && is.finite(rel$estimate[i_cl]) &&
      (is.finite(b_clbw) || is.finite(b_clcr))) {
    fac <- 1
    if (is.finite(b_clbw)) fac <- fac * fBW^b_clbw
    if (is.finite(b_clcr)) fac <- fac * fCRCL^b_clcr
    rel$estimate[i_cl] <- rel$estimate[i_cl] * fac
  }
  i_vc <- which(rel$parameter == "TVVc")
  if (length(i_vc) == 1L && is.finite(rel$estimate[i_vc]) && is.finite(b_vcbw)) {
    rel$estimate[i_vc] <- rel$estimate[i_vc] * fBW^b_vcbw
  }

  # recompute error columns on the rescaled intercept estimates
  rel$abs_err     <- ifelse(is.na(rel$estimate) | is.na(rel$true_value),
                            NA_real_, rel$estimate - rel$true_value)
  rel$rel_err     <- ifelse(is.na(rel$estimate) | is.na(rel$true_value) |
                              rel$true_value == 0,
                            NA_real_, (rel$estimate - rel$true_value) / rel$true_value)
  rel$rel_err_pct <- rel$rel_err * 100
  r$rel_err <- rel
  r
}

# ---- File discovery -------------------------------------------------------
discover_vae_files <- function(root = "output",
                               sub  = "vae_covsel_pilot") {
  base <- file.path(root, sub)
  if (!dir.exists(base)) {
    warning("Directory not found: ", base); return(tibble::tibble())
  }
  rds <- list.files(base, pattern = "^res_ds\\d+\\.rds$",
                    recursive = TRUE, full.names = TRUE)
  if (!length(rds)) {
    warning("No res_ds*.rds under ", base); return(tibble::tibble())
  }
  parts <- strsplit(gsub("\\\\", "/", rds), "/", fixed = TRUE)
  covsel_dir <- vapply(parts, function(p) p[length(p) - 1L], character(1)) # "covsel"
  cell_dir   <- vapply(parts, function(p) p[length(p) - 2L], character(1)) # scn<SS>_<struct>
  n_dir      <- vapply(parts, function(p) p[length(p) - 3L], character(1)) # N<N>
  err_txt    <- sub("\\.rds$", "_ERROR.txt", rds)

  tibble::tibble(
    file_path  = rds,
    n_dir      = n_dir,
    cell_dir   = cell_dir,
    covsel_dir = covsel_dir,
    dataset_id = as.integer(sub("^res_ds(\\d+)\\.rds$", "\\1", basename(rds))),
    has_error  = file.exists(err_txt)
  ) |>
    dplyr::mutate(
      sample_N  = suppressWarnings(as.integer(sub("^N", "", n_dir))),
      scenario  = suppressWarnings(as.integer(sub("^scn0*(\\d+)_.*$", "\\1", cell_dir))),
      structure = sub("^scn\\d+_", "", cell_dir),
      estimator = "vae",
      outer_opt = NA_character_
    ) |>
    dplyr::select(sample_N, scenario, structure, estimator, outer_opt,
                  dataset_id, has_error, file_path) |>
    dplyr::arrange(sample_N, scenario, structure, dataset_id)
}

# ---- Unpack one RDS -------------------------------------------------------
.unpack_diag <- function(r, meta, hit) {
  d  <- r$diag    %||% list()
  d3 <- r$diag_t3 %||% list()
  pick <- function(k) d3[[k]] %||% d[[k]]
  cn_val <- suppressWarnings(as.numeric(pick("cond_num_cor") %||% NA_real_))
  tibble::tibble(
    sample_N        = r$sample_N    %||% meta$sample_N,
    scenario        = r$scenario_id %||% meta$scenario,
    structure       = r$structure   %||% meta$structure,
    estimator       = r$estimator   %||% meta$estimator,
    outer_opt       = r$outer_opt   %||% meta$outer_opt,
    dataset_id      = r$dataset_id  %||% meta$dataset_id,
    status          = as.character(r$status %||% NA_character_),
    converged       = as.logical(pick("converged") %||% NA),
    # PMx-standard CN gate (Khandelwal 2019 Fig 7 / NONMEM default = 1000).
    # NA / non-finite -> FALSE so downstream aggregation cannot break.
    cn_below_cutoff = isTRUE(is.finite(cn_val) && cn_val <= 1000),
    cov_ok          = as.logical(pick("cov_ok")  %||% NA),
    objf            = as.numeric(pick("objf")     %||% NA_real_),
    cond_num_cor    = cn_val,
    cond_num_cor_source = as.character(pick("cond_num_cor_source") %||% NA_character_),
    convergence_code    = suppressWarnings(as.integer(pick("convergence_code") %||% NA_integer_)),
    min_suc         = as.logical(d3$min_suc  %||% NA),
    cov_step        = as.logical(d3$cov_step %||% NA),
    est_bnd         = as.logical(d3$est_bnd  %||% NA),
    phys_bnd        = as.logical(d3$phys_bnd %||% NA),
    est_family      = as.character(d3$est_family %||% NA_character_),
    cov_source      = as.character(d3$cov_source %||% NA_character_),
    n_selected      = if (!is.null(r$covsel$selected)) nrow(r$covsel$selected) else 0L,
    n_true          = if (!is.null(r$covsel$true_set)) nrow(r$covsel$true_set) else 0L,
    n_true_hit      = hit$n_true_hit,
    n_false_pos     = hit$n_false_pos,
    exact_match     = hit$exact_match,
    fit_runtime_sec = as.numeric(r$fit_runtime_sec %||% NA_real_),
    message         = as.character(pick("message") %||% NA_character_)
  )
}

.unpack_rse <- function(r, meta) {
  r   <- .backtransform_intercepts(r, meta)
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
  true_set <- .canon_pairs(r$covsel$true_set)
  selected <- .canon_pairs(r$covsel$selected)
  if (!nrow(true_set) && !nrow(selected)) return(tibble::tibble())
  cmp <- dplyr::full_join(
    dplyr::mutate(true_set, in_true = TRUE),
    dplyr::mutate(selected, in_vae  = TRUE),
    by = c("var", "covar")
  ) |>
    dplyr::mutate(
      in_true = tidyr::replace_na(in_true, FALSE),
      in_vae  = tidyr::replace_na(in_vae,  FALSE),
      verdict = dplyr::case_when(
        in_true &  in_vae ~ "TP",
        in_true & !in_vae ~ "FN",
        !in_true &  in_vae ~ "FP"
      )
    )
  tibble::tibble(
    sample_N   = r$sample_N    %||% meta$sample_N,
    scenario   = r$scenario_id %||% meta$scenario,
    structure  = r$structure   %||% meta$structure,
    dataset_id = r$dataset_id  %||% meta$dataset_id,
    var        = cmp$var,
    covar      = cmp$covar,
    in_true    = cmp$in_true,
    in_vae     = cmp$in_vae,
    verdict    = cmp$verdict
  )
}

load_one_vae_rds <- function(path, meta) {
  r <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(r)) {
    message(sprintf("[skip] cannot read %s", path))
    return(list(diag = tibble::tibble(), rse = tibble::tibble(),
                covsel = tibble::tibble()))
  }
  # backfill covsel$selected from parFixed for records written by drivers that
  # hit the fit$theta -> NULL regression (nlmixr2est >= 7.0.0).
  r <- .ensure_selected(r)
  # then backfill the covariate_beta rel_err rows from that recovered selection.
  r <- .ensure_relerr_backfill(r)
  hit <- match_selected_to_truth(r$covsel$selected, r$covsel$true_set)
  list(diag   = .unpack_diag(r, meta, hit),
       rse    = .unpack_rse(r, meta),
       covsel = .unpack_covsel(r, meta))
}

# ---- Bulk loader ----------------------------------------------------------
load_all_vae_results <- function(root = "output", sub = "vae_covsel_pilot",
                                verbose = TRUE) {
  idx <- discover_vae_files(root, sub = sub)
  if (!nrow(idx)) {
    warning("No VAE covsel RDS files found under ", file.path(root, sub))
    return(list(index = idx, diag_long = tibble::tibble(),
                rse_long = tibble::tibble(), covsel_long = tibble::tibble()))
  }
  if (verbose) message(sprintf(
    "Loading %d RDS across %d cells (N x scn x structure)%s",
    nrow(idx),
    dplyr::n_distinct(idx[, c("sample_N", "scenario", "structure")]),
    if (any(idx$has_error)) sprintf(" [%d ERROR sidecar(s) present]",
                                    sum(idx$has_error)) else ""))

  loaded <- purrr::map(seq_len(nrow(idx)), function(i) {
    load_one_vae_rds(idx$file_path[i], meta = as.list(idx[i, ]))
  })
  list(
    index       = idx,
    diag_long   = purrr::map_dfr(loaded, "diag"),
    rse_long    = purrr::map_dfr(loaded, "rse"),
    covsel_long = purrr::map_dfr(loaded, "covsel")
  )
}

# ---- Aggregation: convergence + stability + timing per cell ---------------
compute_vae_diag_rates <- function(diag_long) {
  diag_long |>
    dplyr::mutate(
      converged_strict = (converged %in% TRUE) & (cn_below_cutoff %in% TRUE)
    ) |>
    dplyr::group_by(sample_N, scenario, structure) |>
    dplyr::summarise(
      n_total             = dplyr::n(),
      n_ok_status         = sum(status == "ok", na.rm = TRUE),
      Converged_pct       = 100 * mean(converged, na.rm = TRUE),
      ConvergedStrict_pct = 100 * mean(converged_strict),
      CovStep_pct         = 100 * mean(cov_ok,    na.rm = TRUE),
      CNBelowCutoff_pct   = 100 * mean(cn_below_cutoff),
      MinSuc_pct          = 100 * mean(min_suc,   na.rm = TRUE),
      EstBnd_pct          = 100 * mean(est_bnd,   na.rm = TRUE),
      MedCN               = stats::median(cond_num_cor, na.rm = TRUE),
      MaxCN               = suppressWarnings(max(cond_num_cor, na.rm = TRUE)),
      P90CN               = as.numeric(stats::quantile(cond_num_cor, 0.9, na.rm = TRUE)),
      MedObjF             = stats::median(objf, na.rm = TRUE),
      MedNSelected        = stats::median(n_selected, na.rm = TRUE),
      MedRuntime_sec      = stats::median(fit_runtime_sec, na.rm = TRUE),
      MedRuntime_min      = stats::median(fit_runtime_sec, na.rm = TRUE) / 60,
      .groups             = "drop"
    ) |>
    dplyr::mutate(dplyr::across(where(is.numeric),
                                \(x) ifelse(is.finite(x), x, NA_real_))) |>
    dplyr::arrange(sample_N, scenario, structure)
}

# ---- Aggregation: estimation precision per (cell, parameter) --------------
# rel_err is relative error vs truth (centring-invariant for covariate_beta;
# reference-dependent for structural_intercept -- see header caveat).
compute_vae_estim <- function(rse_long, diag_long = NULL,
                              mode = c("all", "success", "cond")) {
  mode <- match.arg(mode)
  df <- rse_long |>
    dplyr::filter(is.finite(true_value), true_value != 0, is.finite(rel_err))

  if (mode == "success") {
    if (is.null(diag_long)) stop("mode='success' requires diag_long")
    keep <- diag_long |>
      dplyr::mutate(
        converged_strict = (converged %in% TRUE) & (cn_below_cutoff %in% TRUE)
      ) |>
      dplyr::filter(converged_strict) |>
      dplyr::select(sample_N, scenario, structure, dataset_id) |>
      dplyr::distinct()
    df <- dplyr::semi_join(df, keep,
      by = c("sample_N", "scenario", "structure", "dataset_id"))
  } else if (mode == "cond") {
    if (is.null(diag_long)) stop("mode='cond' requires diag_long")
    keep <- diag_long |>
      dplyr::filter(exact_match %in% TRUE) |>
      dplyr::select(sample_N, scenario, structure, dataset_id) |>
      dplyr::distinct()
    df <- dplyr::semi_join(df, keep,
      by = c("sample_N", "scenario", "structure", "dataset_id"))
  }

  if (!nrow(df))
    return(tibble::tibble(sample_N = integer(), scenario = integer(),
                          structure = character(), parameter = character(),
                          param_class = character(), true_value = numeric(),
                          n_used = integer(), MedRE_pct = numeric(),
                          MeanRE_pct = numeric(), MARE_pct = numeric(),
                          RMRSE_pct = numeric(), P90AbsRE_pct = numeric()))
  df |>
    dplyr::group_by(sample_N, scenario, structure, parameter, param_class) |>
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
    dplyr::arrange(sample_N, scenario, structure, param_class, parameter)
}

# ---- Aggregation: Power / PowerCN / PowerMinSuc per cell ------------------
#   Power       = # exact_match          / N
#   PowerCN     = # exact_match & ok_cn  / # ok_cn        (cond_num_cor < cut)
#   PowerMinSuc = # exact_match & ok_min / # ok_min       (converged)
compute_vae_power <- function(diag_long, cn_cor_cut = 1000) {
  diag_long |>
    dplyr::group_by(sample_N, scenario, structure) |>
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
    dplyr::arrange(sample_N, scenario, structure)
}

# ---- Aggregation: relative power per (cell, k) ----------------------------
# For each cell: fraction of datasets recovering >= k of the true covariates,
# for k = 1..n_true.  n_true is read per record (constant within a cell).
compute_vae_relpower <- function(diag_long) {
  diag_long |>
    dplyr::group_by(sample_N, scenario, structure) |>
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
    dplyr::arrange(sample_N, scenario, structure, k)
}

# ---- Aggregation: per-covariate detection per cell ------------------------
# For each (cell, var, covar): whether it is a TRUE covariate, and the fraction
# of datasets in which VAE selected it.  For a true covariate this fraction is
# its per-relationship sensitivity (marginal power); for a distractor it is the
# per-relationship false-positive rate.
compute_vae_covsel_by_covar <- function(covsel_long) {
  if (!nrow(covsel_long)) return(tibble::tibble())
  covsel_long |>
    dplyr::group_by(sample_N, scenario, structure, var, covar) |>
    dplyr::summarise(
      n_datasets     = dplyr::n(),
      is_true        = any(in_true %in% TRUE),
      n_detected     = sum(in_vae %in% TRUE),
      detection_rate = mean(in_vae %in% TRUE),
      n_TP           = sum(verdict == "TP", na.rm = TRUE),
      n_FN           = sum(verdict == "FN", na.rm = TRUE),
      n_FP           = sum(verdict == "FP", na.rm = TRUE),
      .groups        = "drop"
    ) |>
    dplyr::arrange(sample_N, scenario, structure, dplyr::desc(is_true), var, covar)
}

# ---- Top-level driver -----------------------------------------------------
aggregate_vae_covsel_run <- function(root          = "output",
                                     sub           = "vae_covsel_pilot",
                                     out_dir       = file.path(root, "vae_covsel_aggregated"),
                                     cn_cor_cut    = 1000,
                                     write_outputs = TRUE,
                                     verbose       = TRUE) {
  loaded <- load_all_vae_results(root, sub = sub, verbose = verbose)
  if (!nrow(loaded$diag_long)) return(invisible(loaded))

  diag_rates    <- compute_vae_diag_rates(loaded$diag_long)
  estim_all     <- compute_vae_estim(loaded$rse_long, mode = "all")
  estim_success <- compute_vae_estim(loaded$rse_long, loaded$diag_long, mode = "success")
  estim_cond    <- compute_vae_estim(loaded$rse_long, loaded$diag_long, mode = "cond")
  power         <- compute_vae_power(loaded$diag_long, cn_cor_cut = cn_cor_cut)
  relpower      <- compute_vae_relpower(loaded$diag_long)
  covsel_by_cov <- compute_vae_covsel_by_covar(loaded$covsel_long)

  if (write_outputs) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    readr::write_csv(loaded$index,       file.path(out_dir, "vae_file_index.csv"))
    readr::write_csv(loaded$diag_long,   file.path(out_dir, "vae_diag_long.csv"))
    readr::write_csv(loaded$rse_long,    file.path(out_dir, "vae_rse_long.csv"))
    readr::write_csv(loaded$covsel_long, file.path(out_dir, "vae_covsel_long.csv"))
    readr::write_csv(diag_rates,         file.path(out_dir, "vae_diag_rates.csv"))
    readr::write_csv(estim_all,          file.path(out_dir, "vae_estim_all.csv"))
    readr::write_csv(estim_success,      file.path(out_dir, "vae_estim_success.csv"))
    readr::write_csv(estim_cond,         file.path(out_dir, "vae_estim_cond.csv"))
    readr::write_csv(power,              file.path(out_dir, "vae_power.csv"))
    readr::write_csv(relpower,           file.path(out_dir, "vae_relpower.csv"))
    readr::write_csv(covsel_by_cov,      file.path(out_dir, "vae_covsel_by_covar.csv"))
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
    ), file.path(out_dir, "vae_covsel_aggregated.rds"))
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
  root <- "output"; sub <- "vae_covsel_pilot"; out_dir <- NULL
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
  aggregate_vae_covsel_run(
    root    = root,
    sub     = sub,
    out_dir = out_dir %||% file.path(root, "vae_covsel_aggregated")
  )
}
