# ============================================================================
# scm_bench_helpers.R
# ----------------------------------------------------------------------------
# Helpers shared by PerformanceEvaluation_scm_bench.R and the aggregator.
# Extracted verbatim (minor comment edits) from
# PerformanceEvaluation_sc_n_ds_HPCE.R so the bench driver can source these
# without triggering the CLI dispatch block at the bottom of that file.
#
# Exports:
#   to_nm_dataset(sim_obs)
#   diagnose_fit(fit)
#   extract_params_long(fit, includeCov = TRUE)
#   rel_err_one(est_long, true_long, scenario_id)
#   runSCM_traced(label, ...)
#   package_scm_result(label, scm_res, runtime_sec, true_long, scenario_id,
#                      final_ctrl = NULL, final_est = "focei")
#
# Base model:
#   base_2cmt_oral_linCmt()
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(nlmixr2)
  library(nlmixr2est)
  library(rxode2)
})

## ---- Base model (2-cmt oral, linCmt) --------------------------------------
base_2cmt_oral_linCmt <- function() {
  ini({
    lTVCL <- log(0.6)
    lTVQ  <- log(1.8)
    lTVVc <- log(20)
    lTVVp <- log(80)
    lTVKA <- fix(log(0.7))

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

## ---- to_nm_dataset --------------------------------------------------------
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

## ---- diagnose_fit ---------------------------------------------------------
diagnose_fit <- function(fit) {
  if (is.null(fit)) {
    return(list(converged = NA, objf = NA_real_,
                cond_num_cor = NA_real_,
                cov_ok = NA, message = NA_character_))
  }
  conv_code  <- if (!is.null(fit$convergence)) fit$convergence else NA_integer_
  cn_cor     <- fit$conditionNumberCor
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

## ---- extract_params_long --------------------------------------------------
extract_params_long <- function(fit, includeCov = TRUE) {
  theta <- fit$theta
  om    <- fit$omega
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
  cov_names_in_fit <- names(theta)
  cov_map <- list(
    CLBW   = c("TH_BW_CL",    grep("^cov_(bw|wt)_power_cl$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    CLcrCL = c("TH_CRCL_CL",  grep("^cov_crcl_power_cl$",      cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcBW   = c("TH_BW_VC",    grep("^cov_(bw|wt)_power_vc$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcSEX  = c("TH_SEX_VC",   grep("^cov_sex_[^_]+_vc$",       cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
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

## ---- rel_err_one ----------------------------------------------------------
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

## ---- runSCM_traced --------------------------------------------------------
runSCM_traced <- function(label, ...) {
  t0  <- Sys.time()
  res <- nlmixr2scm::runSCM(...)
  attr(res, "elapsed_s") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res
}

## ---- package_scm_result ---------------------------------------------------
# BENCH VARIANT: `final_est` selects the estimator used for the tight-tol
# covariance refit.  Default "focei" preserves the legacy behaviour used
# in PerformanceEvaluation_sc_n_ds_HPCE.R.  For the SAEM bench cell the
# driver passes final_est = "foceif" so SAEM's covariance is obtained from
# a fast foceif refit rather than a slow plain-focei refit.
package_scm_result <- function(label, scm_res, runtime_sec,
                               true_long, scenario_id,
                               final_ctrl = NULL,
                               final_est  = "focei") {
  .pickFit <- function(x) {
    if (is.null(x)) return(NULL)
    cand <- if (is.list(x) && length(x) >= 1L) x[[1L]] else x
    if (inherits(cand, "nlmixr2FitCore")) cand else NULL
  }
  final_fit <- .pickFit(scm_res$resBck)
  if (is.null(final_fit)) final_fit <- .pickFit(scm_res$resFwd)

  cov_done      <- FALSE
  refit_runtime <- NA_real_
  if (!is.null(final_fit) && !is.null(final_ctrl)) {
    t0 <- Sys.time()
    refit <- tryCatch(
      nlmixr2(final_fit$ui, nlme::getData(final_fit),
              est = final_est, control = final_ctrl),
      error = function(e) {
        warning("package_scm_result(): refit failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    refit_runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    if (!is.null(refit)) {
      final_fit <- refit
      cov_done  <- !is.null(refit$cov) && all(is.finite(diag(refit$cov)))
    }
  }

  .parse_cov_theta <- function(nms) {
    hits <- grep("^cov_", nms, value = TRUE)
    if (length(hits) == 0L) {
      return(tibble::tibble(theta_name = character(),
                            covar = character(), shape = character(),
                            var = character()))
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
  } else NULL

  final_est_long <- if (!is.null(final_fit)) extract_params_long(final_fit) else NULL
  rel_err        <- if (!is.null(final_est_long))
    rel_err_one(final_est_long, true_long, scenario_id) else NULL
  diag           <- if (!is.null(final_fit)) diagnose_fit(final_fit) else NULL
  parFixed       <- if (cov_done) final_fit$parFixedDf else NULL

  list(
    label             = label,
    selected          = selected,
    step_hist         = scm_res$summaryTable,
    final_fit         = final_fit,
    final_est         = final_est_long,
    rel_err           = rel_err,
    diag              = diag,
    parFixed          = parFixed,
    cov_done          = cov_done,
    runtime_sec       = runtime_sec,
    refit_runtime_sec = refit_runtime,
    refit_estimator   = final_est
  )
}
