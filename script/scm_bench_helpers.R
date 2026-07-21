# ==============================================================================
# scm_bench_helpers.R
# ------------------------------------------------------------------------------
# Helper functions for the SCM estimator x optimizer benchmark: base structural
# models, NONMEM-format data conversion, fit diagnostics, parameter extraction,
# runSCM tracing, and the schema-2.1 packaging bridge.
# ==============================================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- base_2cmt_oral_linCmt ------------------------------------------------
## Analytic 2-cmt oral base model (linCmt()). This is the structural base fit
## for the linCmt cells; SCM adds covariate relations on top of cl/vc.
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
# stamp structure so fit_model_type()/assemble_common() resolve model_type
attr(base_2cmt_oral_linCmt, "structure") <- "linCmt"

## ---- base_2cmt_oral_ode ---------------------------------------------------
## Explicit 2-cmt oral ODE surrogate for the linCmt() base model.
##
## Rationale: nlmixr2 structurally forces fast=FALSE (finite-difference outer
## gradient) for any linCmt() model, so the Almquist-2015 analytic gradient
## can never engage.  An explicit d/dt() model is required to unlock
## foceiControl(fast=TRUE) and to give SAEM a well-behaved integrand for
## post-hoc Gaussian-quadrature / IS likelihoods.
##
## Parameterisation is identical to the linCmt() form:
##   cl, vc = clearance / central volume
##   q,  vp = inter-compartmental clearance / peripheral volume
##   ka     = first-order absorption
## Micro-rate form:  k=cl/vc, k12=q/vc, k21=q/vp.
## Compartments named depot / central / periph to match to_nm_dataset()'s
## CMT labels (depot dose, central observation).
base_2cmt_oral_ode <- function() {
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

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot -
                      (cl / vc) * central -
                      (q  / vc) * central +
                      (q  / vp) * periph
    d/dt(periph)  <-  (q / vc) * central -
                      (q / vp) * periph

    cp <- central / vc
    cp ~ prop(prop.err)
  })
}
# stamp structure so fit_model_type()/assemble_common() resolve model_type
attr(base_2cmt_oral_ode, "structure") <- "ode"

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
## Rule change (2026-07): rely on numerical evidence (finite OFV + finite
## covariance + finite condition #), NOT on fit$convergence.  Rationale:
##   - lbfgsb3c returns convergence=8 ("false convergence") for otherwise
##     perfectly good fits, poisoning the old rule.
##   - SAEM has no fit$convergence at all.
## The raw exit code is preserved as `convergence_code` for reference.
##
## Two convergence flags returned:
##   converged        (lenient): finite OFV + cov_ok + finite cond_num_cor
##   cn_below_cutoff:            cond_num_cor <= CN_STRICT_CUTOFF (1000)
## Aggregator combines with est_bnd (from diagnose_fit_table3) to build the
## PMx-standard strict convergence flag.  1000 matches NONMEM's default
## condition-number report threshold and Khandelwal 2019 Fig 7.
##
## cond_num_cor fallback: nlmixr2est's irlsfocei class does not populate
## fit$conditionNumberCor.  When missing, derive from cov via cov2cor+kappa,
## and record via cond_num_cor_source in {"native", "fallback", NA}.
CN_STRICT_CUTOFF <- 1000

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
    cn_val <- tryCatch(
      kappa(stats::cov2cor(fit$cov), exact = TRUE),
      error = function(e) NA_real_
    )
    cn_src <- if (is.finite(cn_val)) "fallback" else NA_character_
  } else {
    cn_val <- NA_real_; cn_src <- NA_character_
  }

  list(
    converged           = isTRUE(is.finite(objf_val) && cov_ok && is.finite(cn_val)),
    objf                = objf_val,
    cond_num_cor        = cn_val,
    cond_num_cor_source = cn_src,
    cn_below_cutoff     = isTRUE(is.finite(cn_val) && cn_val <= CN_STRICT_CUTOFF),
    cov_ok              = cov_ok,
    convergence_code    = conv_code,
    message             = if (!is.null(fit$message)) as.character(fit$message) else NA_character_
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
# Records BOTH wall-clock (Sys.time) and CPU-seconds (proc.time). CPU includes
# the user.child/sys.child columns, where runSCM's forked `workers` credit
# their work once reaped -- so cpu_s captures the parallel candidate fits that
# wall-clock hides.  hog = cpu_s / elapsed_s is the realised parallel speedup.
cpu_secs <- function(pt) {
  unname(sum(pt[c("user.self", "sys.self", "user.child", "sys.child")],
             na.rm = TRUE))
}

runSCM_traced <- function(label, ...) {
  t0  <- Sys.time()
  p0  <- proc.time()
  res <- nlmixr2scm::runSCM(...)
  attr(res, "elapsed_s") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  attr(res, "cpu_s")     <- cpu_secs(proc.time() - p0)
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

## ---- package_scm_schema21 -------------------------------------------------
# SCHEMA-2.1 BRIDGE for the SCM benchmark. Unlike package_scm_result() (which
# emits the ad-hoc 1.x list), this packages a runSCM() result into the SAME
# record shape the single-fit drivers write via assemble_common(), so the
# schema-2.1 aggregators can consume SCM cells without special-casing.
#
# Flow:
#   1. pull the SCM winner  -- forward search leaves the final fit in
#      scm_res$resFwd[[1]] (resBck is NULL for forward-only); fall back to
#      resBck for a backward/bidirectional search.
#   2. tight-tol refit      -- re-fit the winner's ui with refit_ctrl (the
#      CELL'S OWN tier="final" control, e.g. irlsfoceif+lbfgsb3c), dispatched
#      through nlmixr_est_name() so foceif/irlsfoceif hit the right alias.
#      This is what populates cov / SE / cond_num_cor for the record.
#   3. assemble_common()    -- build the estimator-agnostic record core off
#      the refit fit (schema_version, objf, converged, cond_num_cor, rel_err,
#      diag_t3, parFixed, cov, ...).
#   4. runtime override (A) -- assemble_common() stamps a SCALAR
#      fit_runtime_sec; SCM has several phases, so we DROP that scalar and
#      attach a $runtime LIST (base_sec, scm_sec, refit_sec, total_sec).
#   5. $scm block           -- selected covariates (parsed from the winner's
#      cov_<covar>_<shape>_<var> theta names), step_hist = summaryTable, and
#      cov_done. This DUPLICATES nothing on disk: there is no separate scm.rds.
#   6. sidecar write        -- write_fit_sidecar() persists <name>.rds +
#      <name>.fit.rds (the refit winner) + <name>.meta.json.
#
# Args:
#   scm_res         : the runSCM() (or runSCM_traced()) return value.
#   true_mod        : true-model fn (structure/bounds attrs) for assemble_common.
#   scenario_id     : integer scenario for the rel-err join.
#   true_params     : long-format true parameter table.
#   estimator       : grid label of the SEARCH estimator ("focei","irlsfoceif",..)
#   outer_opt       : outer optimizer of the refit ("nlminb","lbfgsb3c",NA).
#   refit_ctrl      : the tier="final" control for the covariance refit. When
#                     NULL, no refit is done and the winner is used as-is.
#   refit_estimator : estimator label for the refit dispatch (defaults to the
#                     cell's own `estimator`, i.e. refined settings of the same
#                     method).
#   runtimes        : list(base_sec=, scm_sec=, base_cpu=, scm_cpu=) timing
#                     inputs. *_sec are wall-clock; *_cpu are CPU-seconds
#                     (self + forked-child). refit_sec/refit_cpu are measured
#                     here; total_sec/total_cpu are the phase sums.
#   identity        : list(sample_N=, dataset_id=, boundary=) extra scalar keys
#                     folded into the record + meta manifest.
#   out_rds         : target .rds path; when non-NULL the three sidecars are
#                     written. When NULL the record is returned without writing.
#
# Returns the schema-2.1 record list (invisibly written to disk when out_rds
# is supplied).
package_scm_schema21 <- function(scm_res, true_mod, scenario_id, true_params,
                                 estimator, outer_opt = NA_character_,
                                 refit_ctrl = NULL, refit_estimator = estimator,
                                 runtimes = list(), identity = list(),
                                 structure = NULL, out_rds = NULL) {

  .pickFit <- function(x) {
    if (is.null(x)) return(NULL)
    cand <- if (is.list(x) && length(x) >= 1L) x[[1L]] else x
    if (inherits(cand, "nlmixr2FitCore")) cand else NULL
  }
  .parse_cov_theta <- function(nms) {
    hits <- grep("^cov_", nms, value = TRUE)
    if (length(hits) == 0L) {
      return(tibble::tibble(theta_name = character(),
                            var = character(), covar = character(),
                            shape = character()))
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

  # 1. SCM winner: for a bidirectional "scm" search the FINAL model is the
  #    backward-eliminated one (resBck) -- forward-selected covariates that
  #    backward drops must NOT appear in `selected`. Prefer resBck; fall back
  #    to resFwd only for a forward-only search (resBck NULL).
  #    2026-07-20: this precedence was previously REVERSED (resFwd first),
  #    which retained backward-dropped covariates and inflated the null-
  #    scenario false-positive rate (Power ~0.9 -> ~0.6). Matches the original
  #    package_scm_result() ordering.
  winner <- .pickFit(scm_res$resBck)
  if (is.null(winner)) winner <- .pickFit(scm_res$resFwd)

  # 2. tight-tol covariance refit with the cell's own final settings.
  refit_sec <- NA_real_
  refit_cpu <- NA_real_
  cov_done  <- FALSE
  if (!is.null(winner) && !is.null(refit_ctrl)) {
    est_dispatch <- nlmixr_est_name(refit_estimator)
    t0 <- Sys.time()
    p0 <- proc.time()
    refit <- tryCatch(
      nlmixr2(winner$ui, nlme::getData(winner),
              est = est_dispatch, control = refit_ctrl),
      error = function(e) {
        warning("package_scm_schema21(): refit failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    refit_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    refit_cpu <- cpu_secs(proc.time() - p0)
    if (!is.null(refit)) {
      winner   <- refit
      cov_done <- !is.null(refit$cov) && all(is.finite(diag(refit$cov)))
    }
  }

  # 3. common schema-2.1 core (built off the refit winner).
  rec <- assemble_common(
    fit         = winner,
    true_mod    = true_mod,
    scenario_id = scenario_id,
    true_params = true_params,
    runtime_sec = NA_real_,             # replaced by $runtime list below
    status      = if (is.null(winner)) "error" else "ok",
    estimator   = estimator
  )

  # 4. runtime override (recommendation A): drop the scalar, attach a list.
  #    Two parallel blocks are kept:
  #      $runtime -> WALL-clock seconds per phase (what the user waits).
  #      $cpu     -> CPU seconds per phase (self + forked-child), i.e. total
  #                  work.  hog_factor = cpu_total / wall_total is the realised
  #                  parallel speedup nlmixr2 delivered (1 = fully serial).
  base_sec <- as.numeric(runtimes$base_sec %||% NA_real_)
  scm_sec  <- as.numeric(runtimes$scm_sec  %||% NA_real_)
  base_cpu <- as.numeric(runtimes$base_cpu %||% NA_real_)
  scm_cpu  <- as.numeric(runtimes$scm_cpu  %||% NA_real_)
  rec$fit_runtime_sec <- NULL
  wall_total <- sum(c(base_sec, scm_sec, refit_sec), na.rm = TRUE)
  cpu_total  <- sum(c(base_cpu, scm_cpu, refit_cpu), na.rm = TRUE)
  rec$runtime <- list(
    base_sec  = base_sec,
    scm_sec   = scm_sec,
    refit_sec = refit_sec,
    total_sec = wall_total
  )
  rec$cpu <- list(
    base_sec   = base_cpu,
    scm_sec    = scm_cpu,
    refit_sec  = refit_cpu,
    total_sec  = cpu_total,
    hog_factor = if (wall_total > 0) cpu_total / wall_total else NA_real_
  )

  # 5. $scm block (subsumes the old standalone scm.rds).
  selected <- if (!is.null(winner)) {
    .parse_cov_theta(names(winner$theta)) %>%
      dplyr::mutate(estimate = unname(winner$theta[theta_name])) %>%
      dplyr::select(var, covar, shape, theta_name, estimate)
  } else NULL
  rec$scm <- list(
    selected  = selected,
    step_hist = scm_res$summaryTable,
    cov_done  = cov_done
  )

  # identity + refit metadata keys (folded in for the record + meta manifest).
  rec$sample_N          <- identity$sample_N  %||% NA_integer_
  rec$scenario_id       <- scenario_id
  rec$dataset_id        <- identity$dataset_id %||% NA_integer_
  rec$estimator         <- estimator
  rec$outer_opt         <- outer_opt
  rec$refit_estimator   <- refit_estimator
  rec$boundary          <- identity$boundary  %||% NA_character_
  rec$structure         <- structure %||% attr(true_mod, "structure") %||% rec$model_type
  rec$refit_runtime_sec <- refit_sec          # populate meta from runtime$refit_sec

  # Parallelism setup actually in effect (recorded, NOT set here): scm_workers
  # is the runSCM fork width; rx_threads is rxode2::getRxThreads() as observed
  # on the compute node.  hog_factor above is the realised speedup.
  rec$scm_workers <- identity$scm_workers %||% NA_integer_
  rec$rx_threads  <- identity$rx_threads  %||% NA_integer_

  # Flat scalar mirrors of the nested $runtime/$cpu blocks so the timing and
  # the parallelism benefit are greppable from the tiny .meta.json manifest
  # (which only carries flat scalars) without opening the .rds.
  rec$wall_total_sec <- wall_total
  rec$cpu_total_sec  <- cpu_total
  rec$hog_factor     <- rec$cpu$hog_factor

  # 6. persist: res_ds*.rds + res_ds*.fit.rds + res_ds*.meta.json.
  if (!is.null(out_rds)) {
    dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
    write_fit_sidecar(rec, out_rds, fit = winner)
  }

  invisible(rec)
}