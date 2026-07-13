# ============================================================================
# estimator_factory.R
# ----------------------------------------------------------------------------
# Central dispatch for the SCM operating-characteristics benchmark.
# Given (est, outer_opt, tier), returns a control object suitable for either
# the SCM screening step (tier="screen") or the final tight-tol refit
# (tier="final"), together with metadata used by the driver.
#
# Method grid (see valid_combos()):
#   focei      x {bobyqa, nlminb, lbfgsb3c}   -- baseline family
#   foceif     x {nlminb, lbfgsb3c}           -- Almquist analytic gradient
#   saem       x  NA                          -- covariance via post-SAEM foceif refit
#
# 2026-07-12: irlsfoceif temporarily removed from the grid; its nlmixr2est
# class does not populate fit$conditionNumberCor and the cov2cor+kappa fallback
# also fails, so it never satisfies the PMx-strict convergence flag. Re-add
# when the upstream cov step is fixed.
#
# VAE was dropped from the benchmark: runSCM() fails to locate cl/vc symbols
# in the vae-returned ui (parameter labels rewritten to lTVCL/lTVVc). Kept in
# nlmixr2est for standalone base fits only.
#
# tier = "screen": tuned sigdig/derivEps/ODE tols per optimizer,
#                  covMethod="", calcTables=FALSE (fast; LRT/OFV only).
# tier = "final" : same tuning + covMethod="r,s", calcTables=TRUE
#                  (SE + cond_num_cor for the SCM winner).
# 2026-07-12: restored two-tier plumbing after Option-A fix -- the loose
# ODE tols were the earlier problem, not the two-tier idea itself. Both tiers
# now share the same tight ODE tols and per-optimizer derivEps/sigdig; they
# differ only in whether the covariance step runs.
# ============================================================================

suppressPackageStartupMessages({
  library(nlmixr2est)
  library(rxode2)
  library(tibble)
})

# ---- Valid (est, outer_opt) cells ------------------------------------------
# Full reference grid.  The submit_all_arrays.sh launcher defaults to a subset
# that EXCLUDES focei (already migrated as focei_bobyqa) to avoid re-running
# it.  Override via ESTIMATORS env var if you want the full 8-cell grid.
valid_combos <- function() {
  tibble::tribble(
    ~estimator,     ~outer_opt,
    "focei",        "bobyqa",
    "focei",        "nlminb",
    "focei",        "lbfgsb3c",
    "foceif",       "nlminb",
    "foceif",       "lbfgsb3c",
    # "irlsfoceif", "lbfgsb3c",  # 2026-07-12: cov step unavailable; excluded
    "saem",         NA_character_
  )
}

is_valid_combo <- function(est, outer_opt) {
  vc <- valid_combos()
  match_na <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
  any(vc$estimator == est & vapply(seq_len(nrow(vc)),
                                   function(i) match_na(vc$outer_opt[i], outer_opt),
                                   logical(1)))
}

# ---- Seed reproducibility --------------------------------------------------
seed_for_dataset <- function(ds_id, base = 1000L) as.integer(base + ds_id)

seed_all <- function(seed, estimator) {
  set.seed(seed)
  invisible(seed)
}

# ---- Control factory -------------------------------------------------------
# Returns list($ctrl, $est, $outer_opt, $tier, $needs_post_refit)
#   needs_post_refit = TRUE for saem final -> use foceif for the covariance
#   step (the driver handles the refit; the estimator_factory just flags it).
make_est_control <- function(est,
                             outer_opt = NA_character_,
                             tier      = c("final", "screen")) {
  tier <- match.arg(tier)

  # Optimizer class: gradient-based outer optimizers need tighter ODE
  # tolerances so that FD gradients are not corrupted by ODE noise, plus a
  # larger FD step (derivEps) to overwhelm residual noise. Derivative-free
  # bobyqa is immune and uses the classic settings.
  is_grad_opt <- !is.na(outer_opt) && outer_opt %in% c("nlminb", "lbfgsb3c")

  sigdig    <- if (is_grad_opt) 5 else 4
  atol      <- if (is_grad_opt) 1e-10 else 1e-8
  rtol      <- if (is_grad_opt) 1e-8  else 1e-6
  rxc       <- rxode2::rxControl(atol = atol, rtol = rtol)

  # Tier-specific: screen skips covariance and table build; final does both.
  covMethod  <- if (tier == "screen") "" else "r,s"
  calcTables <- (tier == "final")

  # Larger FD step for gradient-based outer optimizers (default is ~2e-4,
  # which is inside ODE noise floor at rtol=1e-6). 5e-3 keeps gradient
  # signal well above 1e-8-scale noise for log-scale thetas.
  # 2026-07-13: (B) split derivEps by optimizer class.
  #   - lbfgsb3c: uses a projected-gradient stopping rule that is very
  #     sensitive to FD noise -> keep the larger 5e-3 step.
  #   - nlminb  : uses relative-gradient stopping; 5e-3 combined with
  #     FOCEi's ~0.05-unit inner-ETA ripple caused premature "false
  #     convergence (8)" exits before the LRT could see real slope
  #     signals. 1e-3 gives ~5x sharper gradient without dropping into
  #     the ODE noise floor.
  #   - bobyqa is derivative-free; derivEps ignored.
  derivEps <- if (is.na(outer_opt) || !is_grad_opt) {
    c(2e-4, 2e-4)
  } else if (outer_opt == "nlminb") {
    c(1e-3, 1e-3)
  } else {  # lbfgsb3c
    c(5e-3, 5e-3)
  }

  needs_post_refit <- FALSE

  ctrl <- switch(est,

    focei = nlmixr2est::foceiControl(
      sigdig             = sigdig,
      outerOpt           = outer_opt,
      print              = 0,
      calcTables         = calcTables,
      covMethod          = covMethod,
      stickyRecalcN      = 20,
      maxOuterIterations = 2000,
      maxInnerIterations = 2000,
      derivEps           = derivEps,
      rxControl          = rxc
    ),

    foceif = nlmixr2est::foceiControl(
      interaction        = TRUE,
      sigdig             = sigdig,
      outerOpt           = outer_opt,
      print              = 0,
      calcTables         = calcTables,
      covMethod          = covMethod,
      stickyRecalcN      = 20,
      maxOuterIterations = 2000,
      maxInnerIterations = 2000,
      derivEps           = derivEps,
      rxControl          = rxc
    ),

    irlsfoceif = nlmixr2est::foceiControl(
      interaction        = TRUE,
      sigdig             = sigdig,
      outerOpt           = outer_opt,      # lbfgsb3c only in practice
      print              = 0,
      calcTables         = calcTables,
      covMethod          = covMethod,
      stickyRecalcN      = 20,
      maxOuterIterations = 2000,
      maxInnerIterations = 2000,
      derivEps           = derivEps,
      rxControl          = rxc
    ),

    saem = {
      # SAEM has no outerOpt.
      # 2026-07-13: chain length UNIFIED across screen/final at 500/1000.
      # Rationale: bench_scm at scn16/N80 showed saem_NA selecting the wrong
      # covariate shape (lin vs power) with Power=0. The previous screen
      # chain (300/400) produced OFV noise ~5-10 units, comparable to typical
      # SCM LRT deltas (3.84 for df=1 at alpha=0.05), so shape-vs-shape
      # comparisons became random. 500 burn + 1000 EM is the paper's chain
      # length and yields LRT-grade OFV stability.
      # NOTE: saemControl$covMethod does NOT accept "" (unlike foceiControl);
      # match.arg requires one of {linFim, fim, sa, r,s, r, s}. We use
      # "linFim" for both tiers -- cheap linearized FIM, PMx-standard.
      # Post-SAEM foceif refit still handles r,s cov for the winning model.
      needs_post_refit <- (tier == "final")
      nlmixr2est::saemControl(
        nBurn      = 500L,
        nEm        = 1000L,
        nmc        = 3L,
        print      = 0,
        calcTables = calcTables,
        covMethod  = "linFim",
        rxControl  = rxc
      )
    },

    stop(sprintf("Unknown estimator: '%s'", est))
  )

  list(ctrl = ctrl, est = est, outer_opt = outer_opt, tier = tier,
       needs_post_refit = needs_post_refit)
}

# ---- Build the foceif control used for the post-SAEM covariance refit ----
# 2026-07-12: switched from lbfgsb3c to bobyqa. The bench_refit run at N=80
# scn16 showed lbfgsb3c refits pinned at the boundary (RMRSE ~450%), while
# bobyqa refits recovered the true covariate slopes cleanly (RMRSE ~15-30%).
make_saem_refit_control <- function() {
  make_est_control("foceif", outer_opt = "bobyqa", tier = "final")$ctrl
}
