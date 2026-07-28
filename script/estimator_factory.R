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
#   irlsfocei  x {bobyqa}                      -- IRLS mu-referenced FOCEi (non-fast, est=ifocei)
#   irlsfoceif x {lbfgsb3c}                    -- IRLS mu-referenced FOCEi (fast, est=ifoceif)
#   saem       x  NA                          -- covariance via post-SAEM foceif refit
#
# ACTIVE COMPARISON (2026-07-23, manager-agreed): focei x bobyqa vs
#   irlsfocei x bobyqa -- same derivative-free outer optimiser, differing ONLY
#   in whether the mu-referenced IRLS profiling of the population/covariate
#   thetas is on. This isolates the IRLS speed-up. `irlsfocei` is the NON-fast
#   IRLS grid label (dispatches est="ifocei"): bobyqa is derivative-free so the
#   fast `irlsfoceif` (est="ifoceif") variant's analytic sensitivity-model build
#   (~106 s) would be pure overhead. All other (est, opt) cells remain valid
#   combos but are PARKED (not in the default sweep).
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
    "irlsfocei",    "bobyqa",    # 2026-07-23 ACTIVE: IRLS-vs-FOCEi speed test
                                 # under the SAME derivative-free outer optimiser
                                 # (non-fast, dispatches est="ifocei").
    "irlsfoceif",   "lbfgsb3c"   # PARKED: fast IRLS grid label (dispatches
                                 # est="ifoceif"), for the analytic-gradient outer.
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
# The benchmark grid uses descriptive labels (foceif, irlsfocei, irlsfoceif)
# that map onto nlmixr2's internal est= aliases (verified live on nlmixr2est
# 7.0.1):
#   getValidNlmixrCtl.foceif -> foceiControl   (foceif,  fast=TRUE)
#   getValidNlmixrCtl.ifocei -> ifoceiControl  (ifocei,  fast=FALSE, IRLS mu-ref)
#   getValidNlmixrCtl.ifoceif-> ifoceiControl  (ifoceif, fast=TRUE,  IRLS mu-ref)
#
# foceif dispatches DIRECTLY (not via est="focei") because the alias sets
# fast=TRUE, turning on the Almquist ANALYTIC outer gradient consumed by its
# gradient-based optimisers (nlminb/lbfgsb3c).
#
# irlsfocei (ACTIVE, x bobyqa) is the NON-fast IRLS grid label -> est="ifocei".
# bobyqa is derivative-free -- it never consumes an outer gradient -- so the
# fast irlsfoceif (est="ifoceif") variant would BUILD the symbolic sensitivity
# model (~106 s setup) only to discard it. ifocei gives the identical IRLS
# mu-referenced estimation and native r,s covariance with none of that overhead.
#
# irlsfoceif (PARKED, x lbfgsb3c) is the FAST IRLS grid label -> est="ifoceif";
# kept for a possible analytic-gradient-outer re-activation.
#
# saem/focei pass through unchanged. Requires nlmixr2est >= 7.0.1.
#
# NOTE: nlmixr2est 7.0.0+ rejects the RAW grid labels "irlsfocei"/"irlsfoceif"
# as est strings (no such method), so every nlmixr2()/base-fit call site MUST
# route through nlmixr_est_name(). "foceif" slips through only because it happens
# to also be a valid est string.
nlmixr_est_name <- function(estimator) {
  switch(estimator,
    foceif     = "foceif",    # FOCEi + interaction, analytic gradient (fast=TRUE)
    irlsfocei  = "ifocei",    # 2026-07-23 ACTIVE: plain (non-fast) IRLS-FOCEi.
    irlsfoceif = "ifoceif",   # PARKED: fast IRLS-FOCEi + interaction.
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
                             outer_opt     = NA_character_,
                             tier          = c("final", "screen"),
                             screen_sigdig = NA_real_,
                             screen_atol   = NA_real_,
                             screen_rtol   = NA_real_,
                             warm          = "calc",
                             fix1          = TRUE) {
  tier <- match.arg(tier)
  # fix1 toggle (2026-07-22): when TRUE, the analytic-gradient branches
  # (foceif/irlsfoceif) use sensitivity-matched ODE tolerances (rxc_sens);
  # when FALSE they use the default rxControl (rxc), whose sensitivity tols are
  # loosened 10x. Exposed as an argument so a diagnostic can A/B the effect of
  # Fix 1 on a single (est, opt) cell. Has no effect on focei/saem.
  # Inner-Hessian seeding for the n1qn1 inner problem. nlmixr2est 6.2.0 changed
  # the DEFAULT from the classic self-initialized Hessian ("save", used by the
  # ORIGINAL SCM run) to a recomputed Hessian ("calc"). "calc" perturbs each
  # candidate fit's OFV at the ~1e-2 level, which flips near-threshold SCM
  # selections and inflates the null-scenario false-positive rate. warm="save"
  # reproduces the original behaviour. Applied to BOTH tiers so the whole
  # search is consistent.
  warm <- match.arg(warm, c("calc", "save"))

  # Optimizer class: gradient-based outer optimizers need tighter ODE
  # tolerances so that FD gradients are not corrupted by ODE noise, plus a
  # larger FD step (derivEps) to overwhelm residual noise. Derivative-free
  # bobyqa is immune and uses the classic settings.
  is_grad_opt <- !is.na(outer_opt) && outer_opt %in% c("nlminb", "lbfgsb3c")

  sigdig    <- if (is_grad_opt) 5 else 4
  atol      <- if (is_grad_opt) 1e-10 else 1e-8
  rtol      <- if (is_grad_opt) 1e-8  else 1e-6

  # A/B knob: on the SCREEN tier ONLY, optionally coarsen the candidate-LRT
  # optimization precision to match the ORIGINAL pipeline (sigdig=3, atol=1e-6,
  # rtol=1e-4). Tighter screening (the current default sigdig=4) realizes more
  # of each spurious covariate's dOFV past the chi^2_1=3.84 forward-LRT
  # threshold, inflating the null-scenario false-positive rate. The FINAL tier
  # (covariance refit) is deliberately untouched so it stays tight.
  if (tier == "screen") {
    if (!is.na(screen_sigdig)) sigdig <- screen_sigdig
    if (!is.na(screen_atol))   atol   <- screen_atol
    if (!is.na(screen_rtol))   rtol   <- screen_rtol
  }
  rxc       <- rxode2::rxControl(atol = atol, rtol = rtol)

  # FIX 1 (2026-07-22): sensitivity-matched rxControl for the ANALYTIC-gradient
  # estimators (foceif, irlsfoceif; est aliases set fast=TRUE).
  #
  # By default rxControl loosens the forward-sensitivity tolerances 10x relative
  # to the state solve (maxAtolRtolFactor=0.1): with atol=1e-10/rtol=1e-8 the
  # sensitivity ODEs are actually solved at atolSens=1e-9/rtolSens=1e-7. That
  # looser sensitivity solve is what trips the Almquist analytic OUTER gradient
  # into "could not be solved at this point" -> fallbackFD swaps in an FD
  # gradient for the affected iterations. A quasi-Newton outer optimizer then
  # sees a gradient that flips between the analytic and FD scales, poisoning its
  # inverse-Hessian (Hessian resets, "last objf not at minimum", near-threshold
  # SCM selection flips).
  #
  # Matching atolSens/rtolSens (and their steady-state counterparts) to the
  # tight STATE tols keeps the sensitivity solve accurate enough that the
  # analytic gradient stays solvable, so fallbackFD stops firing. Verified live
  # (scn16 ODE, ds001, cold-start candidate): identical OBJF to the loose-sens
  # run and ~5x faster (17.6 s vs 87.6 s), with NO analytic->FD switch.
  # Applies only to the analytic-gradient branches; focei (fast=FALSE) and saem
  # are unaffected and keep the default rxc.
  rxc_sens  <- rxode2::rxControl(atol = atol, rtol = rtol,
                                 atolSens = atol, rtolSens = rtol,
                                 ssAtolSens = atol, ssRtolSens = rtol)

  # fix1 selects which rxControl the analytic-gradient branches receive.
  rxc_analytic <- if (isTRUE(fix1)) rxc_sens else rxc

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
      warm               = warm,
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
      warm               = warm,
      rxControl          = rxc_analytic   # FIX 1 (toggled): sens-matched tols
    ),

    # irlsfocei (ACTIVE, non-fast) and irlsfoceif (PARKED, fast) share one
    # control body: the fast-ness is carried by the est= alias (ifocei vs
    # ifoceif) that nlmixr_est_name() picks. We build the NATIVE ifoceiControl
    # (nlmixr2est 7.0.1) rather than a plain foceiControl(muModel="irls"): the
    # latter is only accepted via an auto-conversion that emits
    # `.minfo("converting foceiControl to ifoceiControl")` every fit. ifoceiControl
    # forwards ... to foceiControl and HARD-CODES muModel="irls" (so we omit it),
    # keeping interaction=TRUE. The parked fast path (getValidNlmixrCtl.ifoceif =
    # .foceiFastCtl(control, ifoceiControl)) is built from this SAME control class
    # and just flips fast=TRUE, so one body serves both est aliases cleanly.
    irlsfocei  = ,
    irlsfoceif = nlmixr2est::ifoceiControl(
      # 2026-07-23: the ACTIVE cell is irlsfocei x bobyqa (est="ifocei"). bobyqa
      # is derivative-free so no outer gradient is consumed; the mu-referenced
      # IRLS profiling of the population/covariate thetas (the whole point of the
      # i* family) is active regardless of fast=, so irlsfocei x bobyqa is the
      # correct IRLS counterpart to the plain focei x bobyqa baseline. The parked
      # irlsfoceif x lbfgsb3c cell reuses this body but dispatches est="ifoceif".
      interaction        = TRUE,           # ifoceiControl forces muModel="irls"
      sigdig             = sigdig,
      outerOpt           = outer_opt,      # bobyqa (active) or lbfgsb3c (parked)
      print              = 0,
      printNcol          = 10L,   # explicit; see focei branch note.
      calcTables         = calcTables,
      covMethod          = covMethod,
      stickyRecalcN      = 20,
      maxOuterIterations = 2000,
      maxInnerIterations = 2000,
      derivEps           = derivEps,
      # optExpression/fallbackFD are inert on the non-fast ifocei path (no
      # analytic outer gradient is built under fast=FALSE), but are kept TRUE as
      # harmless safe defaults so the parked fast irlsfoceif cell needs no
      # control change. Together with est="ifocei"/"ifoceif" (via
      # nlmixr_est_name) this yields the IRLS-FOCEi algorithm with native r,s
      # covariance; bobyqa supplies the derivative-free outer search.
      optExpression      = TRUE,
      fallbackFD         = TRUE,
      warm               = warm,
      rxControl          = rxc_analytic   # FIX 1 (toggled): sens-matched tols
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
