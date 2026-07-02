library(nlmixr2utils)
library(nlmixr2scm)
library(nlmixr2)
library(rxode2)
library(tidyverse)

true_params <- readRDS("Inputdataset/true_params_long.rds")

## ============================================================================
# Helper1: Convert to NONMEM-format dataset for nlmixr2 fitting / runSCM-------
## ----------------------------------------------------------------------------
##   For each (SCENARIO, DATASET, SUBJECT) we create:
##     - one EVID = 1 dose row at TIME = 0 with AMT = DOSE_MG, CMT = "depot"
##     - 6 EVID = 0 observation rows at sample_times with DV = cp_obs
##   ID is unique within a (SCENARIO, DATASET) fit; for a global ID use
##   paste(SCENARIO, DATASET, SUBJECT, sep = "_").
## ============================================================================
to_nm_dataset <- function(sim_obs) {
  obs_rows <- sim_obs %>%
    dplyr::transmute(
      SCENARIO, DATASET,
      ID   = SUBJECT,
      TIME = time,
      EVID = 0L,
      AMT  = 0,
      CMT  = "central",
      DV   = cp_obs,
      BW, BMI, CrCL, SEX, RACE
    )

  dose_rows <- sim_obs %>%
    dplyr::distinct(SCENARIO, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE) %>%
    dplyr::transmute(
      SCENARIO, DATASET,
      ID   = SUBJECT,
      TIME = 0,
      EVID = 1L,
      AMT  = DOSE_MG,
      CMT  = "depot",
      DV   = NA_real_,
      BW, BMI, CrCL, SEX, RACE
    )

  dplyr::bind_rows(dose_rows, obs_rows) %>%
    dplyr::arrange(SCENARIO, DATASET, ID, TIME, dplyr::desc(EVID))
}

## ============================================================================
## Helper2  Convergence diagnostics
##   converged    : TRUE iff optimizer reported success (fit$convergence == 0)
##                  AND objf is finite. 
##   cond_num_cor : raw condition number of the *correlation* matrix of fixed
##                  effects (lambda_max / lambda_min of cor). Common threshold: < 1000.
##   cov_ok       : TRUE iff fit$cov is a finite-diagonal matrix.  Will be
##                  FALSE when covMethod = "" was used (cov not computed).
##   message      : optimizer exit message, e.g. "Normal exit from bobyqa".#
## ============================================================================
diagnose_fit <- function(fit) {
  if (is.null(fit)) {
    return(list(converged = NA, objf = NA_real_,
                cond_num_cor = NA_real_,
                cov_ok = NA, message = NA_character_))
  }
  conv_code <- if (!is.null(fit$convergence)) fit$convergence else NA_integer_
  cn_cor    <- fit$conditionNumberCor
  cn_cor_val <- if (is.null(cn_cor)) NA_real_ else as.numeric(cn_cor)
  list(
    converged    = isTRUE(conv_code == 0L) &&
                     !is.null(fit$objf) && is.finite(fit$objf),
    objf         = if (!is.null(fit$objf)) fit$objf else NA_real_,
    cond_num_cor = cn_cor_val,
    cov_ok       = isTRUE(!is.null(fit$cov) && all(is.finite(diag(fit$cov)))),
    message      = if (!is.null(fit$message)) as.character(fit$message) else NA_character_
  )
}


## ============================================================================
## Helper3  Extract fixed-effects, random-effects, residual error, and covariate effects into a long tibble
## ============================================================================
extract_params_long <- function(fit, includeCov = TRUE) {
  theta <- fit$theta
  om    <- fit$omega
  ## Fixed-effects: keep the "TV" prefix so names line up with true_params_long()
  ## (TVCL, TVVc, TVQ, TVVp, TVKA).  fit$theta still stores them on the LOG
  ## scale (lTVCL etc.), so we exp() back to natural-scale typical values.
  fe_long <- tibble::tibble(
    parameter = c("TVCL", "TVVc", "TVQ", "TVVp", "TVKA"),
    src_name  = c("lTVCL", "lTVVc", "lTVQ", "lTVVp", "lTVKA")
  ) %>%
    dplyr::mutate(estimate = exp(unname(theta[src_name]))) %>%
    dplyr::select(parameter, estimate)

  rand_long <- tibble::tibble(
    parameter = c("var_CL", "var_Vc", "cov_VcCL"),
    estimate  = c(om["eta.cl", "eta.cl"],
                  om["eta.vc", "eta.vc"],
                  om["eta.vc", "eta.cl"])
  )

  res_long <- tibble::tibble(
    parameter = "ResErr",
    estimate  = unname(theta["prop.err"])
  )
  ## Covariate effects: read from fit$theta if present, else NA.
  ## The regexes below mirror these two naming conventions exactly.
  ##   Naming conventions (verified against runSCM output):
  ##     continuous : cov_<COVAR>_<SHAPE>_<VAR>     e.g. cov_BW_power_cl
  ##     categorical: cov_<COVAR>_<LEVEL>_<VAR>     e.g. cov_SEX_1_vc
  ##   Extra keys (CLBMI, VcBMI, VcCrCL, VcRACE) are included so that
  ##   collinearity-driven false positives at small N (e.g. BMI->CL stealing
  ##   the slot of BW->CL) are surfaced in the wide compare table instead of
  ##   being silently dropped.
  cov_names_in_fit <- names(theta)
  cov_map <- list(
    # --- True relations in the scenario-16 simulator ---------------------
    CLBW   = c("TH_BW_CL",    grep("^cov_(bw|wt)_power_cl$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    CLcrCL = c("TH_CRCL_CL",  grep("^cov_crcl_power_cl$",      cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcBW   = c("TH_BW_VC",    grep("^cov_(bw|wt)_power_vc$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcSEX  = c("TH_SEX_VC",   grep("^cov_sex_[^_]+_vc$",       cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    # --- Distractors / collinear false positives -------------------------
    CLBMI  = c("TH_BMI_CL",   grep("^cov_bmi_power_cl$",       cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcBMI  = c("TH_BMI_VC",   grep("^cov_bmi_power_vc$",       cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcCrCL = c("TH_CRCL_VC",  grep("^cov_crcl_power_vc$",      cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcRACE = c("TH_RACE_VC",  grep("^cov_race_[^_]+_vc$",      cov_names_in_fit, value = TRUE, ignore.case = TRUE))
  )
  cov_long <- tibble::tibble(
    parameter = names(cov_map),
    estimate  = vapply(cov_map, function(nms) {
      hit <- intersect(nms, cov_names_in_fit)
      if (length(hit) == 1L) unname(theta[hit]) else NA_real_
    }, numeric(1))
  )

  out <- dplyr::bind_rows(fe_long, rand_long, res_long, cov_long)
  if (!includeCov) out <- dplyr::filter(out, !parameter %in% names(cov_map))
  out
}

## ============================================================================
## Helper4: Calculate Relative error per parameter (one dataset)
## ============================================================================
rel_err_one <- function(est_long, true_long, scenario_id) {
  true_long %>%
    dplyr::filter(scenario == scenario_id) %>%
    dplyr::select(parameter, true_value) %>%
    dplyr::full_join(est_long, by = "parameter") %>%
    dplyr::mutate(
      abs_err  = estimate - true_value,
      rel_err  = ifelse(is.na(estimate) | is.na(true_value) |
                          true_value == 0, NA_real_,
                        (estimate - true_value) / true_value),
      rel_err_pct = rel_err * 100
    )
}

## ============================================================================
## ---- Helper5: package one runSCM result for downstream comparison --------
##   Returns:
##     - selected:        accepted (var, covar, shape) pairs
##     - step_hist:       full SCM step history (forward + backward summary)
##     - final_est:       extract_params_long() of the final model
##     - rel_err:         relative error vs true scenario parameters
##     - runtime_sec:     wall-clock seconds for the runSCM call
##     - diag:            convergence diagnostics for the final fit
##     - parFixed:        nlmixr2 parFixedDf (Estimate, SE, %RSE, CI) when
##                         the refit succeeded; else NULL
##     - cov_done:        TRUE iff the final_ctrl refit succeeded and the
##                         resulting fit has a populated $cov + parFixedDf.
##                         FALSE when (a) final_ctrl = NULL, or (b) the refit
##                         hit an error.  When FALSE: cond_num_cor stays NA
##                         and PowerCN becomes unavailable for that dataset.
##
##   Covariance computation (the refit path -- restored after the failed
##   getVarCov() experiment).  Why a full refit instead of just cov:
##     * getVarCov() / .setCov() freezes theta + ETAs at the SCM-final point
##       (maxOuterIterations=0L, maxInnerIterations=0L) and only recomputes
##       the Hessian.  But the screening control stops on sigdig=3 with
##       atol=1e-6 / rtol=1e-4, so theta is NOT at a tight-tol stationary
##       point.  The numerical Hessian at a non-stationary theta is ill-
##       conditioned -> R-matrix inversion fails -> fallback to S-only cov
##       -> SEs are 5-50x too small (S inflated by the systematic non-zero
##       gradient).  No cov-time tolerance tightening can fix this -- only
##       moving theta to a tight-tol stationary point does.
##     * The refit (scm_focei_final: sigdig=4, atol=1e-8, rtol=1e-6,
##       covMethod="r,s") re-converges theta + ETAs at tight tols and lands
##       on a stationary point where R is well-conditioned.  Theta shifts
##       are tiny (<3% on covariate effects in pilot ds01) but the cov
##       quality is qualitatively different.
##     * Cost: ~30-60s/ds, negligible against the ~15 min SCM screening.
##
##   final_ctrl:
##     * NULL (default): skip the refit entirely.  selected / step_hist /
##       final_est / rel_err / diag still populate from the screening-tol
##       fit; parFixed = NULL and cov_done = FALSE.
##     * scm_focei_final (or similar): trigger the refit.  Recommended
##       whenever uncertainty / cond_num_cor / PowerCN are needed.
## ============================================================================
package_scm_result <- function(label, scm_res, runtime_sec,
                               true_long = true_params, scenario_id = 9,
                               final_ctrl = NULL) {
  ## Pick the final fit.
  .pickFit <- function(x) {
    if (is.null(x)) return(NULL)
    cand <- if (is.list(x) && length(x) >= 1L) x[[1L]] else x
    if (inherits(cand, "nlmixr2FitCore")) cand else NULL
  }
  final_fit <- .pickFit(scm_res$resBck)
  if (is.null(final_fit)) final_fit <- .pickFit(scm_res$resFwd)

  ## Optional refit with the diagnostic control.  On success the refit
  ## REPLACES final_fit so every downstream extraction (estimates, rel_err,
  ## diag, parFixed) reflects the tight-tol stationary point.
  cov_done <- FALSE
  if (!is.null(final_fit) && !is.null(final_ctrl)) {
    refit <- tryCatch(
      nlmixr2(final_fit$ui, nlme::getData(final_fit),
              est = final_fit$est, control = final_ctrl),
      error = function(e) {
        warning("package_scm_result(): refit failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(refit)) {
      final_fit <- refit
      cov_done  <- !is.null(refit$cov) && all(is.finite(diag(refit$cov)))
    }
  }

  ## --- Final-model covariate relations -
  ##The only reliable source-of-truth for "what's in the final model" is final_fit's
  ## own theta vector, since runSCM names covariate coefficients
  ##     cov_<COVAR>_<SHAPE-or-LEVEL>_<VAR> (e.g. cov_BW_power_cl, cov_SEX_1_vc).
  .parse_cov_theta <- function(nms) {
    hits <- grep("^cov_", nms, value = TRUE)
    if (length(hits) == 0L) {
      return(tibble::tibble(theta_name = character(),
                            covar      = character(),
                            shape      = character(),
                            var        = character()))
    }
    parts <- strsplit(sub("^cov_", "", hits), "_", fixed = TRUE)
    tibble::tibble(
      theta_name = hits,
      covar      = vapply(parts, `[`, character(1), 1L),
      shape      = vapply(parts, function(p) paste(p[-c(1L, length(p))],
                                                   collapse = "_"),
                          character(1)),
      var        = vapply(parts, function(p) p[length(p)], character(1))
    )
  }

  selected <- if (!is.null(final_fit)) {
    .parse_cov_theta(names(final_fit$theta)) %>%
      dplyr::mutate(estimate = unname(final_fit$theta[theta_name])) %>%
      dplyr::select(var, covar, shape, theta_name, estimate)
  } else {
    NULL
  }

  final_est <- if (!is.null(final_fit)) extract_params_long(final_fit) else NULL
  rel_err   <- if (!is.null(final_est)) {
    rel_err_one(final_est, true_long, scenario_id)
  } else NULL
  diag      <- if (!is.null(final_fit)) diagnose_fit(final_fit) else NULL
  ## parFixedDf is only populated when the refit succeeded.
  parFixed  <- if (cov_done) final_fit$parFixedDf else NULL

  ## --- Packaged result -----------------------------------------------------
  ## `step_hist` is the ONLY view of the search trace we keep.
  list(
    label       = label,
    selected    = selected,
    step_hist   = scm_res$summaryTable,
    final_fit   = final_fit,
    final_est   = final_est,
    rel_err     = rel_err,
    diag        = diag,
    parFixed    = parFixed,
    cov_done    = cov_done,
    runtime_sec = runtime_sec
  )
}



## ============================================================================
## ---- Helper 6: runSCM with elapsed-time stash ---------------------------
##   Minimal pass-through to nlmixr2scm::runSCM(): times the call and stashes
##   the elapsed seconds on the returned object as attr(res, "elapsed_s").
##   `label` is accepted for backward compatibility with existing call sites
##   (Part 2 one-shots, the per-ds driver) but is otherwise unused.  
## ============================================================================
runSCM_traced <- function(label, ...) {
  t0  <- Sys.time()
  res <- nlmixr2scm::runSCM(...)
  attr(res, "elapsed_s") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res
}




## ############################################################################
# STAGE 1 -- Scenario 16: all four covariate effects active -------------------
## ############################################################################
##   Truth (Khandelwal 2019, scenario 16):
##     BW   on CL  (power)  TH_BW_CL   = 0.75
##     CrCL on CL  (power)  TH_CRCL_CL = 0.50
##     BW   on Vc  (power)  TH_BW_VC   = 1.00
##     SEX  on Vc  (cat) TH_SEX_VC  = 0.50
##   All four are TRUE positives.  The SCM smoke test must keep these and
##   reject the false-positive distractors (CrCL on Vc, SEX on CL).
##
##   Parallel structure to the scenario-9 block above:
##     Part 1 -- true-model robustness (refexp vs lin)
##     Part 2 -- runSCM feature tests (forward / backward / user / full)
##   All fits use scm_focei_n: bobyqa, sigdig 4, maxOuter/Inner = 2000,
##   covMethod = "" (LRT only needs OFV).  Re-fit with final_focei outside
##   this block if SE / %RSE / parFixedDf are needed for any final model.
## ############################################################################

scm_focei_screen <- nlmixr2est::foceiControl(
  sigdig             = 3,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = FALSE,
  covMethod          = "",            # SKIP: LRT uses OFV only
  stickyRecalcN      = 20,
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl          = rxode2::rxControl(atol = 1e-6, rtol = 1e-4)
)

scm_focei_final <- nlmixr2est::foceiControl(
  sigdig             = 4,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = TRUE,
  covMethod          = "r,s",         # full sandwich for SE / CN
  stickyRecalcN      = 20,
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl          = rxode2::rxControl(atol = 1e-8, rtol = 1e-6)
)

base_2cmt_oral_linCmt <- function() {
  ini({
    lTVCL <- log(0.6)
    lTVQ  <- log(1.8)
    lTVVc <- log(20)
    lTVVp <- log(80)
    lTVKA <- fix(log(0.7))      # KA unidentifiable from this sparse design

    eta.cl + eta.vc ~ c(
      0.1,
      0.02,
      0.1
    )

    prop.err <- 0.1
  })
  model({
    cl <- exp(lTVCL + eta.cl)
    vc <- exp(lTVVc + eta.vc)
    q  <- exp(lTVQ)
    vp <- exp(lTVVp)
    ka <- exp(lTVKA)
    cp <- linCmt()
    cp ~ prop(prop.err)
  })
}

# Part 3: Operating-characteristics pilot (4 datasets, scenario 16, N=80)-----------

##   Wraps the fast-route full-SCM pipeline used above (single dataset) into
##   a per-dataset driver, loops 4 replicates of scenario 16, then computes
##   six OC artefacts:
##     - oc_pilot_power       : Power, PowerCN, PowerMinSuc
##     - oc_pilot_relpower    : relative power at k = 1..Ntrue
##     - oc_pilot_rmrse_uncond: unconditional RMRSE %
##     - oc_pilot_rmrse_cond  : conditional RMRSE % (exact-match runs only)
##     - oc_pilot_timing      : per-dataset + aggregate wall-clock
##     - oc_pilot_per_ds      : audit log (one row per dataset)
##
##   Helpers below are intentionally scenario-agnostic so the same code drops
##   onto the 100 x 16 scale-up.  See the SCALE-UP STUB at the end.
## ============================================================================


# ---- (P3.1) OC helpers (reusable across scenarios) --------------------------
##  Canonical-shape normalizer.  runSCM emits per-LEVEL theta names for
##  categorical covariates (e.g. cov_SEX_1_vc, cov_RACE_2_vc); the parser
##  in package_scm_result() leaves the level digit in `shape`.  For
##  semantic matching against truth (where the relation is one entity)
##  we collapse any all-digit shape to "cat".  This also folds a 3-level
##  cat that selected multiple levels into ONE relation row.
.canon_shape_tbl <- function(rel) {
  rel %>%
    dplyr::mutate(shape = ifelse(grepl("^[0-9]+$", shape), "cat", shape)) %>%
    dplyr::distinct(var, covar, shape)
}

##  Build the truth tibble of "what runSCM should have selected" for a given
##  scenario.  Reads true_long, keeps rows with non-zero true_value AND a
##  parameter name that maps to a (var, covar, shape) triple.  Defaults
##  cover the article scenarios (scn 1-16); pass `shape_map` if you add
##  custom covariate effects later.
extract_true_relations <- function(
    true_long, scenario_id,
    shape_map = list(
      CLBW   = list(var = "cl", covar = "BW",   shape = "power"),
      CLcrCL = list(var = "cl", covar = "CrCL", shape = "power"),
      VcBW   = list(var = "vc", covar = "BW",   shape = "power"),
      VcSEX  = list(var = "vc", covar = "SEX",  shape = "cat")
    )) {
  hits <- true_long %>%
    dplyr::filter(scenario == scenario_id,
                  parameter %in% names(shape_map),
                  !is.na(true_value), true_value != 0) %>%
    dplyr::pull(parameter)
  if (length(hits) == 0L) {
    return(tibble::tibble(var = character(), covar = character(),
                          shape = character()))
  }
  purrr::map_dfr(hits, function(p) {
    m <- shape_map[[p]]
    tibble::tibble(var = m$var, covar = m$covar, shape = m$shape)
  })
}

##  Compare a single run's `selected` tibble against truth.  Returns a
##  list with: n_true_hit, n_false_pos, exact_match (set equality on
##  (var, covar, shape) after canonicalization).
match_selected_to_truth <- function(selected, true_rel) {
  if (is.null(selected) || nrow(selected) == 0L) {
    return(list(n_true_hit  = 0L,
                n_false_pos = 0L,
                exact_match = nrow(true_rel) == 0L))
  }
  sel <- .canon_shape_tbl(selected)
  tru <- .canon_shape_tbl(true_rel)
  true_hit  <- dplyr::inner_join(sel, tru, by = c("var", "covar", "shape"))
  false_pos <- dplyr::anti_join (sel, tru, by = c("var", "covar", "shape"))
  miss      <- dplyr::anti_join (tru, sel, by = c("var", "covar", "shape"))
  list(
    n_true_hit  = nrow(true_hit),
    n_false_pos = nrow(false_pos),
    exact_match = (nrow(false_pos) == 0L) && (nrow(miss) == 0L)
  )
}


# ---- (P3.2) Per-dataset driver ----------------------------------------------
##  Filters the long sim file to one DATASET, fits base with linCmt, runs full
##  SCM via the fast screening control, then attaches a post-hoc covariance to
##  the SCM-final model via nlme::getVarCov() (inside package_scm_result --
##  no re-estimation, ~5-20s/ds vs ~60s/ds for the old refit path).  Wrapped
##  in tryCatch so a single bad dataset returns a stub list (no $test) instead
##  of halting the loop.  `confirm = FALSE` is forced -- no interactive prompts.
##
##  Resilience design:
##    * `true_long` is threaded explicitly so an unloaded `true_params` in
##      the calling env can't cause silent NULL returns from package_scm_result.
##    * SCM result `res_i` is saved to disk IMMEDIATELY after runSCM_traced,
##      before any packaging.  Even if package_scm_result errors later,
##      the (expensive) SCM artefact stays recoverable.
##    * Error handler echoes the exception via warning() AND writes a
##      per-dataset error log to <save_dir>/<ds_tag>_ERROR.txt.
##
##  Resume support (two-tier cache):
##    * Tier 1 -- packaged result.  If `test_<ds_tag>.rds` exists
##      and neither `force_rerun` nor `force_repackage` is TRUE, the driver
##      returns the cached `test_*` object untouched (~ms cost).
##    * Tier 2 -- SCM screening result.  If `res_<ds_tag>.rds`
##      exists (and tier 1 was skipped or bypassed via `force_repackage`),
##      the driver skips base-fit + SCM screening and runs only
##      package_scm_result() (~30-60s for the tight-tol refit).  This is
##      the supported recovery path when package_scm_result() crashed on
##      a previous run, or when `scm_focei_final` settings have changed
##      and only the cov/SE step needs to be redone.
##    * Tier 3 -- full pipeline.  Both caches missing, or `force_rerun =
##      TRUE`.  Runs base fit + SCM screening + packaging end-to-end.
##
##    Flags:
##      - `force_rerun = TRUE`     bypasses BOTH caches -> tier 3.
##      - `force_repackage = TRUE` ignores stale `test_*` but still reuses
##        cached `res_*` -> tier 2 (saves the ~15 min SCM screening).
##        When `res_*` is also missing, falls through to tier 3.
##
##    On tier-3 success the driver unlinks any stale
##    `<ds_tag>_ERROR.txt` from a prior failed run so the log doesn't
##    accumulate misleading sentinels.
##
##  Covariance:
##    * The driver passes `final_ctrl = scm_focei_final` to
##      package_scm_result(), which re-fits the SCM-final model ONCE with
##      tight tols + covMethod="r,s".  Cost: ~30-60s/ds (negligible vs the
##      ~15 min SCM screening).  Yields trustworthy parFixed / cond_num_cor.
##    * Why a full refit rather than getVarCov() on the screening-tol fit:
##      the screening optimum is not a tight-tol stationary point, so the
##      numerical Hessian there is ill-conditioned and the cov calculation
##      falls back to a deflated S-only matrix.  See package_scm_result()
##      docstring for the full explanation.
##
##  saveModels = FALSE inside runSCM:
##    * The package's per-accepted-step audit trail (scm_<ds_tag>/...rds)
##      is disabled to keep the working-set disk small and avoid OneDrive
##      sync overhead on each step.  The only persisted artefact per ds is
##      the packaged `test_<ds_tag>.rds`.
##
##  Parallelism:
##    * `workers` controls INNER parallelism: candidates within one SCM step.
##      For the pilot (sequential outer loop) workers = L is the default.
##      For the scale-up under outer parallelism (future_pmap), set workers
##      = 1L to avoid nested-future thrash.

run_one_dataset_scn16_N80 <- function(ds_id, save_dir,
                                       sim_long     = sim_obs_scn16_N80,
                                       base_fn      = base_2cmt_oral_linCmt,
                                       screen_ctrl  = scm_focei_screen,
                                       vars_vec     = scm16_vars,
                                       covars_vec   = scm16_covars,
                                       catvars_vec  = scm16_catvars,
                                       shapes_vec   = scm16_shapes,
                                       true_long    = true_params,
                                       workers          = 3L,
                                       keep_res         = TRUE,
                                       force_rerun      = FALSE,
                                       force_repackage  = FALSE,
                                       scenario_id      = 16L) {
  ds_tag    <- sprintf("ds%02d", ds_id)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  test_path <- file.path(save_dir, sprintf("test_%s.rds", ds_tag))
  res_path  <- file.path(save_dir, sprintf("res_%s.rds",  ds_tag))
  err_path  <- file.path(save_dir, sprintf("%s_ERROR.txt", ds_tag))

  ## -- Tier 1: packaged result already on disk -> return immediately.
  ##    force_rerun trumps everything; force_repackage skips tier 1 so the
  ##    driver falls through to tier 2 (repackage from cached res_*).
  if (!force_rerun && !force_repackage && file.exists(test_path)) {
    message(sprintf(">>> [%s] tier-1 cached, skipping (%s)",
                    ds_tag, basename(test_path)))
    return(list(
      ds_id       = ds_id,
      ds_tag      = ds_tag,
      test        = readRDS(test_path),
      t_base_sec  = NA_real_,
      t_scm_sec   = NA_real_,
      t_total_sec = NA_real_,
      resumed     = TRUE
    ))
  }

  ## -- Tier 2: SCM screening cached but packaging missing or invalidated.
  ##    Re-run only package_scm_result() on the cached res_*.  This is the
  ##    supported recovery path for a crashed packaging step OR a deliberate
  ##    re-package after scm_focei_final has been tightened/changed.
  ##    Wrapped in its own tryCatch so a refit crash here doesn't lose the
  ##    expensive res_* artefact -- the ERROR.txt is written and the user
  ##    can retry after fixing helpers or final_ctrl.
  if (!force_rerun && file.exists(res_path)) {
    message(sprintf(">>> [%s] tier-2 re-packaging from cached %s",
                    ds_tag, basename(res_path)))
    return(tryCatch({
      res_i     <- readRDS(res_path)
      t_scm_sec <- suppressWarnings(as.numeric(attr(res_i, "elapsed_s")))
      if (!is.finite(t_scm_sec)) {
        cli::cli_warn("[{ds_tag}] cached res_* missing `elapsed_s` attr; t_scm_sec set to NA.")
        t_scm_sec <- NA_real_
      }
      t_pkg <- system.time(
        test_i <- package_scm_result(
          label       = sprintf("scn%02d_%s_N80", scenario_id, ds_tag),
          scm_res     = res_i,
          runtime_sec = t_scm_sec,
          true_long   = true_long,
          scenario_id = scenario_id,
          final_ctrl  = scm_focei_final
        )
      )
      saveRDS(test_i, test_path)
      if (file.exists(err_path)) unlink(err_path)  # stale sentinel from prior fail
      message(sprintf("<<< [%s] tier-2 done; refit %.1fs, cov=%s",
                      ds_tag, as.numeric(t_pkg["elapsed"]),
                      if (isTRUE(test_i$cov_done)) "ok" else "FAIL"))
      list(
        ds_id       = ds_id,
        ds_tag      = ds_tag,
        test        = test_i,
        t_base_sec  = NA_real_,
        t_scm_sec   = t_scm_sec,
        t_total_sec = t_scm_sec,
        resumed     = TRUE
      )
    }, error = function(e) {
      msg <- conditionMessage(e)
      writeLines(c(format(Sys.time()), "tier-2 repackage failed:", msg), err_path)
      warning(sprintf("[%s] tier-2 FAILED: %s (see %s)", ds_tag, msg, err_path),
              call. = FALSE, immediate. = TRUE)
      list(ds_id = ds_id, ds_tag = ds_tag, error = msg, res = NULL,
           resumed = TRUE)
    }))
  }

  ## -- Tier 3: full pipeline (base fit + SCM screening + package).
  message(sprintf("\n>>> [%s] tier-3 starting full pipeline at %s",
                  ds_tag, format(Sys.time(), "%H:%M:%S")))

  ## -- Pre-flight (a): every helper function the driver + package_scm_result()
  ##    transitively call must be in scope.  Promise forcing only catches
  ##    missing *arguments*; missing *functions* fail silently 15 min into
  ##    SCM (cf. the diagnose_fit crash on first dry run).  Checking here
  ##    converts that into an instant hard stop that halts the pilot loop,
  ##    forcing the caller to source the helpers before retrying.
  needed_fns <- c("to_nm_dataset", "runSCM_traced", "package_scm_result",
                  "diagnose_fit", "extract_params_long", "rel_err_one",
                  "match_selected_to_truth")
  missing_fns <- needed_fns[!vapply(needed_fns, exists, logical(1),
                                    mode = "function", inherits = TRUE)]
  if (length(missing_fns)) {
    stop(sprintf("[%s] missing helper function(s): %s -- source PerformanceEvaluation04062026.R first.",
                 ds_tag, paste(missing_fns, collapse = ", ")),
         call. = FALSE)
  }

  ## -- Pre-flight (b): force evaluation of every promise so a missing
  ##    true_params (or any other default argument) fails loudly here,
  ##    not silently 25 min later inside package_scm_result().
  force(true_long); force(sim_long); force(base_fn); force(screen_ctrl)

  res_i <- NULL  # placeholder visible to the error handler

  tryCatch({
    ds_i <- to_nm_dataset(sim_long) %>%
      dplyr::filter(DATASET == ds_id) %>%
      dplyr::select(-SCENARIO, -DATASET) %>%
      dplyr::mutate(ID   = as.integer(ID),
                    SEX  = as.integer(SEX),
                    RACE = as.integer(RACE))

    t_base <- system.time(
      fit_base_i <- nlmixr2(base_fn, ds_i,
                            est = "focei", control = screen_ctrl)
    )
    t_base_sec <- as.numeric(t_base["elapsed"])

    res_i <- runSCM_traced(
      label       = sprintf("scn%02d_%s_N80", scenario_id, ds_tag),
      data        = ds_i,
      fit         = fit_base_i,
      varsVec     = vars_vec,
      covarsVec   = covars_vec,
      catvarsVec  = catvars_vec,
      shapes      = shapes_vec,
      searchType  = "scm",
      control     = screen_ctrl,
      saveModels  = FALSE,       # skip per-step audit trail (OneDrive sync tax)
      workers     = workers,
      print       = 0,           # silence per-iteration FOCEi noise
      maxRetries  = 0L,
      confirm     = FALSE        # no interactive y/n prompt
    )
    t_scm_sec <- as.numeric(attr(res_i, "elapsed_s"))

    ## Save the (expensive) SCM result FIRST -- before any packaging.
    ## Then if package_scm_result fails downstream we can replay it
    ## offline without losing the SCM compute.
    if (keep_res) {
      saveRDS(res_i, file.path(save_dir,
                               sprintf("res_%s.rds", ds_tag)))
    }

    ## package_scm_result() refits the SCM-final model once with scm_focei_final
    ## (tight tols + covMethod="r,s") to obtain trustworthy SE / cond_num_cor.
    ## Returns cov_done = FALSE if the refit hits an error; downstream OC code
    ## handles that by filtering on cond_num_cor < cn_cor_cut.
    test_i <- package_scm_result(
      label         = sprintf("scn%02d_%s_N80", scenario_id, ds_tag),
      scm_res       = res_i,
      runtime_sec   = t_scm_sec,
      true_long     = true_long,         # explicit -- no lazy global lookup
      scenario_id   = scenario_id,
      final_ctrl    = scm_focei_final    # tight-tol refit for cov + SE
    )

    saveRDS(test_i, test_path)
    if (file.exists(err_path)) unlink(err_path)  # clear stale sentinel on recovery

    message(sprintf("<<< [%s] done in %.1f min (base %.1fs + scm %.1fs, cov=%s)",
                    ds_tag, (t_base_sec + t_scm_sec) / 60,
                    t_base_sec, t_scm_sec,
                    if (isTRUE(test_i$cov_done)) "ok" else "FAIL"))

    list(
      ds_id       = ds_id,
      ds_tag      = ds_tag,
      test        = test_i,
      t_base_sec  = t_base_sec,
      t_scm_sec   = t_scm_sec,
      t_total_sec = t_base_sec + t_scm_sec,
      resumed     = FALSE
    )
  }, error = function(e) {
    msg <- conditionMessage(e)
    writeLines(c(format(Sys.time()), msg), err_path)
    warning(sprintf("[%s] FAILED: %s (see %s)", ds_tag, msg, err_path),
            call. = FALSE, immediate. = TRUE)
    ## Even on failure, return enough info to recover.  res_i is the SCM
    ## result if it got that far -- callers can rerun package_scm_result()
    ## offline.  When NULL we know the failure was before/during SCM.
    list(ds_id = ds_id, ds_tag = ds_tag, error = msg, res = res_i,
         resumed = FALSE)
  })
}


DOSE_MG <- 100
scm16_vars       <- c("cl", "vc")
scm16_covars     <- c("BW", "CrCL", "BMI")
scm16_catvars    <- c("SEX", "RACE")
scm16_shapes     <- c("power", "lin")

save_dir <- file.path("output", "output_N80")
sim_obs_scn16_N80 <- readRDS("Inputdataset/sim_obs_N80/sim_obs_scenario_16.rds")


args    <- commandArgs(trailingOnly = TRUE)
ds_id   <- if (length(args) >= 1) as.integer(args[1]) else 1L

message(sprintf("Running dataset ds_id = %d", ds_id))

t_run <- system.time(
  result <- run_one_dataset_scn16_N80(
    ds_id    = ds_id,
    save_dir = save_dir
  )
)

message(sprintf("Dataset %d finished in %.1f min", ds_id, t_run["elapsed"] / 60))