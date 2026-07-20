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
    "irlsfoceif",   "lbfgsb3c"   # 2026-07-14: re-enabled; dispatches est="ifoceif"
  )
}

is_valid_combo <- function(est, outer_opt) {
  vc <- valid_combos()
  match_na <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
  any(vc$estimator == est & vapply(seq_len(nrow(vc)),
                                   function(i) match_na(vc$outer_opt[i], outer_opt),
                                   logical(1)))
}

# ---- Grid-label -> nlmixr2 est dispatch name -------------------------------
# The benchmark grid uses descriptive labels (foceif, irlsfoceif) that map onto
# nlmixr2's internal est= aliases. With the console now loading nlmixr2est 6.2.0
# (scratch lib, guaranteed by the project .Rprofile), these aliases dispatch
# correctly:
#   getValidNlmixrControl(., "ifoceif") -> class ifoceiControl  (verified live)
#
# We dispatch them DIRECTLY (not via est="focei") because the alias sets
# fast=TRUE, which turns on the genuine Almquist ANALYTIC outer gradient. The
# earlier est="focei" + muModel="irls" recipe did NOT set fast=TRUE, so it
# silently fell back to a FINITE-DIFFERENCE gradient (header showed only
# "mu: irls", no "grad: analytic") and, as a side effect, left covMethod empty
# so SEs never propagated. Verified scn16/ds1/N80:
#   est="ifoceif"          -> grad: analytic, covMethod=r,s (native SE),  307 s
#   est="focei"+irls recipe-> grad: FD,       covMethod="" (cov2se SE),   205 s
# The analytic gradient is the whole point of this estimator, so we take the
# direct alias despite the ~50% longer runtime (the extra time is the symbolic
# sensitivity-model setup, ~106 s, not optExpression).
#
# optExpression=FALSE (set in make_est_control) is still REQUIRED: it prevents
# the parallel-CSE daemon deadlock on the ODE analytic-gradient build. It does
# NOT disable the analytic gradient (proven: the direct ifoceif run had
# optExpression=FALSE yet still reported grad: analytic).
#
# saem/focei pass through unchanged. Requires nlmixr2est >= 6.2.0.
nlmixr_est_name <- function(estimator) {
  switch(estimator,
    foceif     = "foceif",    # FOCEi + interaction, analytic gradient (fast=TRUE)
    irlsfoceif = "ifoceif",   # IRLS-FOCEi + interaction, analytic gradient
    estimator                 # focei, saem pass through unchanged
  )
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
  # 2026-07-13: reverted (B) per-optimizer split; nlminb=1e-3 did not
  # help SCM Power and made bench_refit worse. Return to uniform 5e-3
  # for all gradient optimizers.
  derivEps <- if (is_grad_opt) c(5e-3, 5e-3) else c(2e-4, 2e-4)

  needs_post_refit <- FALSE

  ctrl <- switch(est,

    focei = nlmixr2est::foceiControl(
      sigdig             = sigdig,
      outerOpt           = outer_opt,
      print              = 0,
      printNcol          = 10L,   # explicit: iterPrintControl derives ncol from
                                  # getOption("width") when NULL, which can be <1
                                  # under a narrow console and errors the build.
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
      printNcol          = 10L,   # explicit; see focei branch note.
      calcTables         = calcTables,
      covMethod          = covMethod,
      stickyRecalcN      = 20,
      maxOuterIterations = 2000,
      maxInnerIterations = 2000,
      derivEps           = derivEps,
      # 2026-07-14: optExpression=TRUE (common-subexpression elimination ON).
      # An earlier build DEADLOCKED here ("optimizing duplicate expressions in
      # FOCEi outer gradient model (N chunks, M daemons)") because the CSE pass
      # spawned parallel daemons. nlmixr2est 6.2.0 rewrote that pass to be
      # SERIAL + chunked (.foceiAnalyticAugModelDirs: a plain vapply over
      # 40-line chunks, no daemons), which fixes the deadlock. Direct probe at
      # scn16/ds1/N80 (est="ifoceif"): optExpression=TRUE = 278.7 s vs FALSE =
      # 307.0 s (9% faster), identical objf (-1591.4887), gradient still
      # analytic. So CSE is back ON. fallbackFD=TRUE remains the safety net if
      # the sensitivity equations fail to solve.
      optExpression      = TRUE,
      fallbackFD         = TRUE,
      rxControl          = rxc
    ),

    irlsfoceif = nlmixr2est::foceiControl(
      interaction        = TRUE,
      muModel            = "irls",         # IRLS mu-referenced regression
      sigdig             = sigdig,
      outerOpt           = outer_opt,      # lbfgsb3c only in practice
      print              = 0,
      printNcol          = 10L,   # explicit; see focei branch note.
      calcTables         = calcTables,
      covMethod          = covMethod,
      stickyRecalcN      = 20,
      maxOuterIterations = 2000,
      maxInnerIterations = 2000,
      derivEps           = derivEps,
      # See foceif branch: optExpression=TRUE (CSE on) is 9% faster than FALSE
      # on 6.2.0's serial-chunked gradient build and no longer deadlocks;
      # fallbackFD=TRUE is the safety net. muModel="irls" + interaction=TRUE +
      # est="ifoceif" (via nlmixr_est_name) gives the genuine IRLS-FOCEi
      # analytic-gradient algorithm with native r,s covariance.
      optExpression      = TRUE,
      fallbackFD         = TRUE,
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
      # length and yields LRT-grade OFV stability. NEVER cut iterations to
      # save time -- a 300/500 tuning shifted TH_BW_VC 0.11->0.39 and moved
      # OFV by 38 units (a convergence failure, not a speedup).
      #
      # 2026-07-14: GQ (Gaussian quadrature -2LL) is DELIBERATELY OFF.
      # logLik=FALSE => fit$objf is the FOCEi-approximation OBJF, which is
      # what runSCM's LRT consumes. It is always present, consistent across
      # candidates, and cross-estimator comparable. The approximation error
      # cancels in the nested delta, so GQ adds nothing to covariate
      # selection and only costs time. Final results table reports the
      # FOCEi-approx -2LL/AIC/BIC (decision (a)).
      #
      # 2026-07-14: ODE tolerances kept TIGHT (atol=1e-8, rtol=1e-6). A loose
      # tol experiment (atol=1e-6/rtol=1e-4) with iterations held fixed at
      # the baseline 300/500 moved the OFV by 44 units (-1511.9 -> -1467.9).
      # SAEM's stochastic approximation amplifies ODE-solve error, so even a
      # modest tol loosening corrupts the fit enough to flip covariate
      # decisions. There is NO safe tol speedup here -- accuracy wins.
      #
      # 2026-07-14: cores = 1L (single-threaded), DELIBERATELY. On the HPCE
      # bench array each task is given a fixed thread budget (scm=3 SCM
      # workers x rxThreads=2). Letting SAEM grab extra cores would (a)
      # oversubscribe the node against the array's thread accounting and
      # (b) make SAEM wall-time incomparable to the single-threaded FOCEi
      # family -- an unfair speed comparison. Multi-core rxode2 solving can
      # also reorder the stochastic-approximation across threads, which is
      # not a like-for-like algorithm. Keep SAEM single-threaded.
      # NOTE: saemControl$covMethod does NOT accept "" (unlike foceiControl);
      # match.arg requires one of {linFim, fim, sa, r,s, r, s}. We use
      # "linFim" for both tiers -- cheap (~0.25 s) linearized FIM giving the
      # condition number. Post-SAEM focei refit (bobyqa) handles the r,s cov
      # for the winning model.
      needs_post_refit <- (tier == "final")
      saem_rxc   <- rxode2::rxControl(atol = 1e-8, rtol = 1e-6, cores = 1L)
      nlmixr2est::saemControl(
        nBurn      = 500L,
        nEm        = 1000L,
        nmc        = 3L,
        print      = 0,
        logLik     = FALSE,      # GQ OFF -- FOCEi-approx OBJF only
        calcTables = calcTables,
        covMethod  = "linFim",
        rxControl  = saem_rxc
      )
    },

    stop(sprintf("Unknown estimator: '%s'", est))
  )

  list(ctrl = ctrl, est = est, outer_opt = outer_opt, tier = tier,
       needs_post_refit = needs_post_refit)
}

# ---- Build the FOCEi control used for the post-SAEM covariance refit -------
# 2026-07-12: switched from lbfgsb3c to bobyqa. The bench_refit run at N=80
# scn16 showed lbfgsb3c refits pinned at the boundary (RMRSE ~450%), while
# bobyqa refits recovered the true covariate slopes cleanly (RMRSE ~15-30%).
# 2026-07-15: est is plain "focei" (NOT "foceif"). The "f" in foceif = fast
# (analytic gradient), and both focei/foceif already use interaction=TRUE. Since
# this refit uses bobyqa (derivative-free), the analytic gradient is downgraded
# to fast=FALSE and is never used -- so foceif+bobyqa was byte-identical to
# focei+bobyqa. Using "focei" makes the refit honest: FOCEi-with-interaction,
# FD, bobyqa. It still computes the r,s covariance for the SAEM point estimate.
make_saem_refit_control <- function() {
  make_est_control("focei", outer_opt = "bobyqa", tier = "final")$ctrl
}
