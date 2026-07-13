# ==============================================================================
# true_model_factory.R
# ------------------------------------------------------------------------------
# Scenario- and boundary-aware factory for the "true" 2-cmt oral model used in
# the Figure 3 / Table 3 refit robustness study.
#
# Convention (from Khandelwal et al. 2019, Fig 7 and Eqs 1-3):
#   * Continuous covariates enter as power terms on natural scale, but are
#     re-parameterised inside the log-additive nlmixr2 form:
#       CL = exp(lTVCL + TH_BW_CL * log(BW/70) + TH_CRCL_CL * log(CrCL/95) + eta.cl)
#     so TH_BW_CL, TH_CRCL_CL, TH_BW_VC are dimensionless exponents in
#     the natural scale of the paper's Eq 1-2.  Boundaries `c(lower, upper)`
#     are applied to these thetas directly, matching PsN's $THETA syntax.
#   * Categorical covariate SEX enters as a linear additive term on the log
#     scale:  Vc = ... * exp(TH_SEX_VC * SEX).  Per user spec, TH_SEX_VC is
#     ALWAYS unbounded regardless of the boundary config.
#   * lTVCL / lTVQ / lTVVc / lTVVp are log-scale typical values (unbounded);
#     lTVKA is FIXED at log(0.7) because the 6-point sampling schedule does
#     not identify absorption.
#
# Boundary configurations (applied only to continuous covariate thetas):
#   * "none"   - bare initial estimate, no bounds
#   * "wide"   - c(-1e5, init, 1e5)   (PsN default from paper Fig 7)
#   * "narrow" - c(-10,  init, 10)    (paper's narrow tier)
#   * "tight"  - c(-5,   init, 5)     (sensitivity tier added 2026-07-09)
#
# Initial estimates for continuous covariate thetas: 0.001 (PsN convention,
# paper Fig 7).  This is deliberately far from the true value so that the
# boundary can actively influence the optimizer trajectory.  Categorical
# TH_SEX_VC uses the same 0.001 magnitude.
#
# Model body uses `linCmt()` for the analytic 2-cmt oral solution (~5-15x
# faster than the ODE form, verified in scripts/PerformanceEvaluation04062026.R).
# ==============================================================================

PSN_INIT_CONT   <- 0.5                            # continuous cov init (2026-07-12: was 0.1; still plateau-trapped for nlminb/lbfgsb3c)
PSN_INIT_CAT    <- log(1.5)                       # categorical cov init ~= 0.405 (matches truth for TH_SEX_VC scenarios)
BOUNDARY_WIDE   <- c(-1e5, 1e5)                   # PsN default
BOUNDARY_NARROW <- c(-10, 10)                     # paper's narrow tier
BOUNDARY_TIGHT  <- c(-5, 5)                       # added 2026-07-09: tighter
                                                  # sensitivity tier below
                                                  # BOUNDARY_NARROW

# ---- Build one scenario x boundary UI function -----------------------------
#   Returns:
#     * an anonymous nlmixr2 UI function (has both ini() and model() blocks)
#     * with attributes:
#         scenario_id           <- int
#         boundary              <- "none" | "wide" | "narrow"
#         bounds_spec           <- named list of c(lower, upper) for bounded
#                                  continuous cov thetas (empty for "none"
#                                  or scenarios with no continuous covs)
#         estimated_cont_params <- character vector of continuous cov names
#         categorical_params    <- character vector of cat cov names
#         fn_text               <- assembled source text (for debugging)
make_true_model <- function(scenario_id, boundary = c("none", "wide", "narrow", "tight"),
                            scenarios = PsN_scenarios) {
  boundary <- match.arg(boundary)
  scn <- scenarios[scenarios$scenario == scenario_id, , drop = FALSE]
  if (nrow(scn) != 1L) {
    stop("scenario_id = ", scenario_id, " not found in PsN_scenarios.")
  }

  # ---- ini() lines: base thetas (always present, unbounded) ----------------
  ini_lines <- c(
    "lTVCL <- log(0.6)",
    "lTVQ  <- log(1.8)",
    "lTVVc <- log(20)",
    "lTVVp <- log(80)",
    "lTVKA <- fix(log(0.7))"
  )

  # ---- Helper: emit one bounded / unbounded theta assignment ---------------
  emit_bounded <- function(name, init, boundary) {
    switch(boundary,
      "none"   = sprintf("%s <- %.6g", name, init),
      "wide"   = sprintf("%s <- c(%.6g, %.6g, %.6g)",
                         name, BOUNDARY_WIDE[1], init, BOUNDARY_WIDE[2]),
      "narrow" = sprintf("%s <- c(%.6g, %.6g, %.6g)",
                         name, BOUNDARY_NARROW[1], init, BOUNDARY_NARROW[2]),
      "tight"  = sprintf("%s <- c(%.6g, %.6g, %.6g)",
                         name, BOUNDARY_TIGHT[1],  init, BOUNDARY_TIGHT[2])
    )
  }
  bounds_of <- function(boundary) {
    switch(boundary,
      "none"   = c(-Inf, Inf),
      "wide"   = BOUNDARY_WIDE,
      "narrow" = BOUNDARY_NARROW,
      "tight"  = BOUNDARY_TIGHT
    )
  }

  # ---- Emit continuous covariate thetas conditionally ----------------------
  bounds_spec       <- list()
  estimated_cont    <- character(0)
  cov_body_terms_cl <- character(0)
  cov_body_terms_vc <- character(0)

  if (scn$I_BW_CL == 1L) {
    ini_lines <- c(ini_lines, emit_bounded("TH_BW_CL", PSN_INIT_CONT, boundary))
    cov_body_terms_cl <- c(cov_body_terms_cl, "TH_BW_CL * log(BW / 70)")
    if (boundary != "none") bounds_spec$TH_BW_CL <- bounds_of(boundary)
    estimated_cont <- c(estimated_cont, "TH_BW_CL")
  }
  if (scn$I_CRCL_CL == 1L) {
    ini_lines <- c(ini_lines,
                   emit_bounded("TH_CRCL_CL", PSN_INIT_CONT, boundary))
    cov_body_terms_cl <- c(cov_body_terms_cl, "TH_CRCL_CL * log(CrCL / 95)")
    if (boundary != "none") bounds_spec$TH_CRCL_CL <- bounds_of(boundary)
    estimated_cont <- c(estimated_cont, "TH_CRCL_CL")
  }
  if (scn$I_BW_VC == 1L) {
    ini_lines <- c(ini_lines, emit_bounded("TH_BW_VC", PSN_INIT_CONT, boundary))
    cov_body_terms_vc <- c(cov_body_terms_vc, "TH_BW_VC * log(BW / 70)")
    if (boundary != "none") bounds_spec$TH_BW_VC <- bounds_of(boundary)
    estimated_cont <- c(estimated_cont, "TH_BW_VC")
  }

  # ---- Categorical covariate theta ----------------------------------------
  # 2026-07-12: previously always unbounded; that let lbfgsb3c/saem diverge to
  # 1e13 on VcSEX. Now bounded by same tier as continuous covs.
  categorical <- character(0)
  if (scn$I_SEX_VC == 1L) {
    ini_lines <- c(ini_lines, emit_bounded("TH_SEX_VC", PSN_INIT_CAT, boundary))
    cov_body_terms_vc <- c(cov_body_terms_vc, "TH_SEX_VC * SEX")
    if (boundary != "none") bounds_spec$TH_SEX_VC <- bounds_of(boundary)
    categorical <- "TH_SEX_VC"
  }

  # ---- Omega + residual error (always present) ------------------------------
  ini_lines <- c(
    ini_lines,
    "eta.cl + eta.vc ~ c(0.1, 0.02, 0.1)",
    "prop.err <- 0.1"
  )

  # ---- model() lines: log-additive form + linCmt() ------------------------
  cl_rhs <- paste(c("lTVCL", cov_body_terms_cl), collapse = " + ")
  vc_rhs <- paste(c("lTVVc", cov_body_terms_vc), collapse = " + ")
  model_lines <- c(
    sprintf("lTVCL_typ <- %s", cl_rhs),
    sprintf("lTVVc_typ <- %s", vc_rhs),
    "cl        <- exp(lTVCL_typ + eta.cl)",
    "vc        <- exp(lTVVc_typ + eta.vc)",
    "q         <- exp(lTVQ)",
    "vp        <- exp(lTVVp)",
    "ka        <- exp(lTVKA)",
    "cp        <- linCmt()",
    "cp ~ prop(prop.err)"
  )

  # ---- Assemble function source text and eval() to build a real closure ---
  fn_text <- paste0(
    "function() {\n",
    "  ini({\n    ", paste(ini_lines,   collapse = "\n    "), "\n  })\n",
    "  model({\n    ", paste(model_lines, collapse = "\n    "), "\n  })\n",
    "}"
  )
  fn <- eval(parse(text = fn_text), envir = globalenv())

  attr(fn, "scenario_id")           <- scenario_id
  attr(fn, "boundary")              <- boundary
  attr(fn, "bounds_spec")           <- bounds_spec
  attr(fn, "estimated_cont_params") <- estimated_cont
  attr(fn, "categorical_params")    <- categorical
  attr(fn, "fn_text")               <- fn_text
  fn
}
