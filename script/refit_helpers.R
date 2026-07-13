# ==============================================================================
# refit_helpers.R
# ------------------------------------------------------------------------------
# Self-contained sourceable module of helpers used by the "refit true model"
# workflow (Figure 3 / Table 3 replication).  Every function is a verbatim copy
# from scripts/PerformanceEvaluation04062026.R (see line references below),
# EXCEPT `diagnose_fit_table3()` which is a new extension that returns the seven
# convergence diagnostics from Khandelwal et al. 2019 Table 3 (MinSuc, EstBnd,
# RndErr, ZeroGrad, CovStep, MedCN/StCN via cond_num_cor).
#
# Sourcing order:
#   source("scripts/refit_helpers.R")
#   source("scripts/true_model_factory.R")
#   ...
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

# ---- Constants -------------------------------------------------------------
DOSE_MG <- 100

# ---- PsN scenarios (16 factorial combinations of 4 covariate effects) ------
PsN_scenarios <- data.frame(
  scenario   = 1:16,
  I_BW_CL    = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
  I_CRCL_CL  = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
  I_BW_VC    = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
  I_SEX_VC   = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
)

# ---- Long-format true parameter table --------------------------------------
#   Mirrors scripts/PerformanceEvaluation04062026.R:603 verbatim.
#   Names line up with extract_params_long() column so `rel_err_one()` can
#   left-join without renaming.
true_params_long <- function(scenarios = PsN_scenarios,
                             theta = c(TVCL = 0.6, TVQ = 1.8, TVVc = 20,
                                       TVVp = 80, TVKA = 0.7,
                                       TH_BW_CL = 0.75, TH_CRCL_CL = 0.5,
                                       TH_BW_VC = 1.0, TH_SEX_VC = log(1.5)),
                             omega = c(var_cl = 0.1, var_vc = 0.1,
                                       cov_cl_vc = 0.02),
                             prop_err = 0.1) {
  scenarios |>
    dplyr::rowwise() |>
    dplyr::mutate(
      TVCL     = theta[["TVCL"]],
      TVVc     = theta[["TVVc"]],
      TVQ      = theta[["TVQ"]],
      TVVp     = theta[["TVVp"]],
      TVKA     = theta[["TVKA"]],
      CLBW     = theta[["TH_BW_CL"]]   * I_BW_CL,
      CLcrCL   = theta[["TH_CRCL_CL"]] * I_CRCL_CL,
      VcBW     = theta[["TH_BW_VC"]]   * I_BW_VC,
      VcSEX    = theta[["TH_SEX_VC"]]  * I_SEX_VC,
      var_CL   = omega[["var_cl"]],
      var_Vc   = omega[["var_vc"]],
      cov_VcCL = omega[["cov_cl_vc"]],
      ResErr   = prop_err
    ) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(
      cols      = c(TVCL, TVVc, TVQ, TVVp, TVKA,
                    CLBW, CLcrCL, VcBW, VcSEX,
                    var_CL, var_Vc, cov_VcCL, ResErr),
      names_to  = "parameter",
      values_to = "true_value"
    ) |>
    dplyr::select(scenario, parameter, true_value,
                  I_BW_CL, I_CRCL_CL, I_BW_VC, I_SEX_VC)
}
true_params <- true_params_long()

# ---- Convert one (SCENARIO, DATASET) subset to nlmixr2 dataset -------------
#   Copy of PerformanceEvaluation04062026.R:660 verbatim.
to_nm_dataset <- function(sim_obs) {
  obs_rows <- sim_obs |>
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

  dose_rows <- sim_obs |>
    dplyr::distinct(SCENARIO, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE) |>
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

  dplyr::bind_rows(dose_rows, obs_rows) |>
    dplyr::arrange(SCENARIO, DATASET, ID, TIME, dplyr::desc(EVID))
}

# ---- Load one (SCENARIO, DATASET) slice from a master sim_obs RDS ----------
#   Filters the 7.2M-row master N=300 RDS to the 1800 rows for one fit.
load_scenario_dataset <- function(master_rds, scenario_id, dataset_id) {
  stopifnot(file.exists(master_rds))
  readRDS(master_rds) |>
    dplyr::filter(SCENARIO == scenario_id, DATASET == dataset_id)
}

# ---- Base convergence diagnostics ------------------------------------------
#   Two convergence flags returned:
#     converged  (lenient / numerical): finite OFV + cov_ok + finite cond#
#     cn_below_cutoff:                    cond_num_cor <= CN_STRICT_CUTOFF (1000)
#   The aggregator combines these with est_bnd (from diagnose_fit_table3)
#   to form `converged_strict`, matching NONMEM/Khandelwal Table 3 practice:
#     converged_strict := converged & cn_below_cutoff & !est_bnd
#   The 1000 cutoff matches NONMEM's default "CONDITION NUMBER > 1000" warning
#   and Khandelwal 2019 Fig 7 reporting.
CN_STRICT_CUTOFF <- 1000

## Canonical diagnose_fit (2026-07 rewrite; see scm_bench_helpers.R for notes).
## Kept as duplicate here to avoid a cross-script source dependency in HPCE.
diagnose_fit <- function(fit) {
  if (is.null(fit)) {
    return(list(converged = NA, objf = NA_real_,
                cond_num_cor = NA_real_, cond_num_cor_source = NA_character_,
                cn_below_cutoff = NA,
                cov_ok = NA, convergence_code = NA_integer_,
                message = NA_character_))
  }
  conv_code <- if (!is.null(fit$convergence)) as.integer(fit$convergence) else NA_integer_
  objf_val  <- if (!is.null(fit$objf)) as.numeric(fit$objf) else NA_real_
  cov_ok    <- isTRUE(!is.null(fit$cov) && all(is.finite(diag(fit$cov))))
  cn_native <- fit$conditionNumberCor
  if (!is.null(cn_native) && is.finite(as.numeric(cn_native))) {
    cn_val <- as.numeric(cn_native); cn_src <- "native"
  } else if (cov_ok) {
    cn_val <- tryCatch(kappa(stats::cov2cor(fit$cov), exact = TRUE),
                       error = function(e) NA_real_)
    cn_src <- if (is.finite(cn_val)) "fallback" else NA_character_
  } else {
    cn_val <- NA_real_; cn_src <- NA_character_
  }
  return(list(
    converged           = isTRUE(is.finite(objf_val) && cov_ok && is.finite(cn_val)),
    objf                = objf_val,
    cond_num_cor        = cn_val,
    cond_num_cor_source = cn_src,
    cn_below_cutoff     = isTRUE(is.finite(cn_val) && cn_val <= CN_STRICT_CUTOFF),
    cov_ok              = cov_ok,
    convergence_code    = conv_code,
    message             = if (!is.null(fit$message)) as.character(fit$message) else NA_character_
  ))
}

.diagnose_fit_OLD <- function(fit) {
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

# ---- Extended Table-3 diagnostics ------------------------------------------
#   Maps nlmixr2 (focei) output onto Khandelwal 2019 Table 3 columns:
#     MinSuc    - minimization successful           <- converged
#     EstBnd    - any estimate hit an active bound  <- bounds_spec check
#     RndErr    - rounding / precision warnings     <- grep fit$message
#     ZeroGrad  - zero-gradient exit                <- best-effort probe
#     CovStep   - covariance step successful        <- cov_ok
#     MedCN     - median cond#(cor) across fits     <- aggregated later
#     StCN      - SD of cond#(cor) across fits      <- aggregated later
#
#   `bounds_spec`: named list, one entry per continuous covariate theta
#     that was bounded in this fit's ini block.  Each entry is c(lower, upper)
#     on the natural scale of that theta.  Missing entry = parameter not
#     bounded (config = "none") or categorical/never-bounded parameter.
#     The `boundary_hit_tol` is the absolute distance from the boundary to
#     flag a hit -- 1e-3 matches PsN convention.
diagnose_fit_table3 <- function(fit, bounds_spec = list(),
                                boundary_hit_tol = 1e-3) {
  base <- diagnose_fit(fit)

  if (is.null(fit)) {
    return(c(base, list(
      min_suc       = NA, est_bnd      = NA,
      rnd_err       = NA, zero_grad    = NA,
      cov_step      = NA, phys_bnd     = NA,
      bound_hits    = NA_character_,
      phys_hits     = NA_character_
    )))
  }

  # MinSuc = converged
  min_suc <- base$converged

  # EstBnd: any bounded parameter within `boundary_hit_tol` of its lower/upper?
  bound_hits <- character(0)
  est_bnd    <- FALSE
  if (length(bounds_spec) > 0L) {
    theta <- fit$theta
    for (nm in names(bounds_spec)) {
      lohi <- bounds_spec[[nm]]
      if (!is.numeric(lohi) || length(lohi) != 2L) next
      est <- unname(theta[nm])
      if (is.null(est) || is.na(est)) next
      lo <- lohi[1]; hi <- lohi[2]
      near_lo <- is.finite(lo) && abs(est - lo) < boundary_hit_tol
      near_hi <- is.finite(hi) && abs(est - hi) < boundary_hit_tol
      if (near_lo || near_hi) {
        bound_hits <- c(bound_hits, nm)
        est_bnd    <- TRUE
      }
    }
  } else {
    est_bnd <- NA  # config = "none" has no bounds to hit
  }
  bound_hits_str <- if (length(bound_hits) > 0L)
    paste(bound_hits, collapse = ",") else NA_character_

  # RndErr: search optimizer exit message for precision / rounding cues.
  #   nlmixr2's focei doesn't emit a NONMEM-style "R MATRIX ALGORITHMICALLY
  #   NON-POSITIVE-SEMIDEFINITE" flag, but its message field carries any
  #   numerical-precision warnings from the underlying bobyqa/lbfgsb call.
  rnd_err <- FALSE
  msg <- base$message %||% NA_character_
  if (!is.na(msg) && nzchar(msg)) {
    rnd_err <- grepl("precision|rounding|numerical", msg,
                     perl = TRUE, ignore.case = TRUE)
  }

  # ZeroGrad: nlmixr2 does not expose a stable gradient-norm field across
  # versions.  Best-effort probe: look for common env slots.  Returns NA if
  # unavailable so aggregation can footnote it.
  zero_grad <- NA
  grad_norm <- tryCatch({
    env <- fit$env
    if (!is.null(env)) {
      if (!is.null(env$foceiInfo$grad_norm))            env$foceiInfo$grad_norm
      else if (!is.null(env$grad_norm))                 env$grad_norm
      else if (!is.null(env$.foceiEnv$grad_norm))       env$.foceiEnv$grad_norm
      else NULL
    } else NULL
  }, error = function(e) NULL)
  if (!is.null(grad_norm) && is.finite(grad_norm)) {
    zero_grad <- grad_norm < 1e-6
  }

  # CovStep: already computed as base$cov_ok.
  cov_step <- base$cov_ok

  # PhysBnd (new): flag runs where any covariate theta escapes the
  # physical-plausibility range.  With PsN init 0.001 and no upper bound,
  # nlmixr2's focei can lock onto a degenerate MLE where a categorical
  # covariate theta (e.g. TH_SEX_VC) drifts to 300+ so that predictions for
  # that stratum collapse to zero and the proportional-error sigma^2 -> 0
  # inflates the likelihood to +Inf (log(sigma^2) -> -Inf).  NONMEM's SIGL
  # floor suppresses this pathology; nlmixr2 has no equivalent, so we
  # detect it post-hoc.  Threshold: 10 (matches the "narrow" bound width).
  cov_theta_names <- c("TH_BW_CL", "TH_CRCL_CL", "TH_BW_VC", "TH_SEX_VC")
  cov_theta_names <- intersect(cov_theta_names, names(fit$theta))
  phys_bnd <- FALSE
  phys_hits <- character(0)
  if (length(cov_theta_names) > 0L) {
    for (nm in cov_theta_names) {
      est <- unname(fit$theta[nm])
      if (is.finite(est) && abs(est) > 10) {
        phys_bnd  <- TRUE
        phys_hits <- c(phys_hits, sprintf("%s=%.4g", nm, est))
      }
    }
  }
  phys_hits_str <- if (length(phys_hits) > 0L)
    paste(phys_hits, collapse = ",") else NA_character_

  c(base, list(
    min_suc    = min_suc,
    est_bnd    = est_bnd,
    rnd_err    = rnd_err,
    zero_grad  = zero_grad,
    cov_step   = cov_step,
    phys_bnd   = phys_bnd,
    bound_hits = bound_hits_str,
    phys_hits  = phys_hits_str
  ))
}

# ---- Map fit$theta + fit$omega to article parameter names ------------------
#   Copy of PerformanceEvaluation04062026.R:731 verbatim.
extract_params_long <- function(fit, includeCov = TRUE) {
  theta <- fit$theta
  om    <- fit$omega

  fe_long <- tibble::tibble(
    parameter = c("TVCL", "TVVc", "TVQ", "TVVp", "TVKA"),
    src_name  = c("lTVCL", "lTVVc", "lTVQ", "lTVVp", "lTVKA")
  ) |>
    dplyr::mutate(estimate = exp(unname(theta[src_name]))) |>
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
    CLBW   = c("TH_BW_CL",   grep("^cov_(bw|wt)_power_cl$",  cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    CLcrCL = c("TH_CRCL_CL", grep("^cov_crcl_power_cl$",     cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcBW   = c("TH_BW_VC",   grep("^cov_(bw|wt)_power_vc$",  cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcSEX  = c("TH_SEX_VC",  grep("^cov_sex_[^_]+_vc$",      cov_names_in_fit, value = TRUE, ignore.case = TRUE))
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

# ---- Relative error vs true parameters (one dataset) -----------------------
#   Copy of PerformanceEvaluation04062026.R:793 verbatim.  Sets rel_err = NA
#   when true_value == 0 (the covariate is inactive in this scenario).
rel_err_one <- function(est_long, true_long, scenario_id) {
  true_long |>
    dplyr::filter(scenario == scenario_id) |>
    dplyr::select(parameter, true_value) |>
    dplyr::full_join(est_long, by = "parameter") |>
    dplyr::mutate(
      abs_err     = estimate - true_value,
      rel_err     = ifelse(is.na(estimate) | is.na(true_value) |
                             true_value == 0, NA_real_,
                           (estimate - true_value) / true_value),
      rel_err_pct = rel_err * 100
    )
}

# ---- %||% shim (rlang not required) ----------------------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x
