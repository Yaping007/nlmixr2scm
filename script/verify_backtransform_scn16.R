# ==============================================================================
# verify_backtransform_scn16.R
# ------------------------------------------------------------------------------
# PURPOSE
#   Empirically verify the covariate-effect back-transform correction used to
#   report structural typical values at the DGP reference centers (BW@70,
#   CrCL@95) from a model that was centered somewhere else (e.g. the subject
#   median, as runSCM does by default, or the sample mean, as the VAE does).
#
# STRATEGY (Point 3 = Option B: DIRECT injection, no SCM selection)
#   For each of the first 5 scenario-16 datasets (N=300) we fit the *known*
#   true covariate structure twice, differing ONLY in the centering constant:
#
#     (A) FORCED center  : log(BW/70), log(CrCL/95)   -> exp(lTVCL) is TVCL@(70,95)
#                                                        exp(lTVVc) is TVVc@70   BY CONSTRUCTION
#     (B) MEDIAN center  : log(BW/medBW), log(CrCL/medCrCL)  (runSCM default)
#
#   The back-transform applied to run (B) must reproduce run (A)'s intercepts:
#
#       TV@ref = exp(l_c) * prod_k (ref_k / c_k)^theta_k        (continuous only)
#
#   Categorical SEX (center 0) is excluded from the product.
#
# CHECKS (Point 1 = Option B: also verify slope invariance)
#   1. Intercept identity     : exp(lTVCL)_A  ~=  backtransform( exp(lTVCL)_B )
#                               exp(lTVVc)_A  ~=  backtransform( exp(lTVVc)_B )
#   2. Slope invariance       : theta_hat_A  ~=  theta_hat_B      (centering-invariant)
#   3. Ground-truth recovery  : run (A) intercepts ~= DGP truth (TVCL=0.6, TVVc=20)
#
# scn16 truth: TVCL=0.6, TVVc=20; TH_BW_CL=0.75, TH_CRCL_CL=0.5,
#              TH_BW_VC=1.0, TH_SEX_VC=log(1.5). DGP centers BW@70, CrCL@95.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(nlmixr2est)
})

sd_ <- "script"
source(file.path(sd_, "scm_bench_helpers.R"), chdir = FALSE)  # base models, to_nm_dataset, extract_params_long
source(file.path(sd_, "estimator_factory.R"), chdir = FALSE)  # make_est_control, nlmixr_est_name

# ---- constants ---------------------------------------------------------------
if (!exists("DOSE_MG"))    DOSE_MG    <- 100
REF_BW    <- 70
REF_CRCL  <- 95
N_DATASETS <- 5
SIM_RDS   <- "Inputdataset/sim_obs_N300/sim_obs_scenario_16.rds"

TRUE_TVCL     <- 0.6
TRUE_TVVc     <- 20
TRUE_TH_BW_CL <- 0.75
TRUE_TH_CR_CL <- 0.5
TRUE_TH_BW_VC <- 1.0
TRUE_TH_SEX   <- log(1.5)

# ---- covariate UI builder (direct injection of the scn16 structure) ----------
# Uses the internal single-pass rebuilder so all four covariates are injected
# onto the clean base without the tainted-UI problem of a per-covariate loop.
# theta names: cov_BW_power_cl, cov_CrCL_power_cl, cov_BW_power_vc, cov_SEX_cat_vc
build_scn16_ui <- function(base_fn, cBW, cCRCL) {
  base_ui <- rxode2::assertRxUi(base_fn)
  pairs_df <- data.frame(
    var     = c("cl", "cl", "vc", "vc"),
    covar   = c("BW_power", "CrCL_power", "BW_power", "SEX_cat"),
    covExpr = c(sprintf("log(BW/%.10g)",   cBW),
                sprintf("log(CrCL/%.10g)", cCRCL),
                sprintf("log(BW/%.10g)",   cBW),
                "SEX"),
    init    = c(0.75, 0.5, 1.0, 0.4),   # near-truth starts; centering-independent
    lower   = rep(-5, 4),
    upper   = rep( 5, 4),
    stringsAsFactors = FALSE
  )
  nlmixr2scm:::.rebuildUiFromPairs(base_ui, pairs_df)
}

# Pull the four covariate slopes by exact theta name (Point 1 = B: no reliance
# on extract_params_long's narrow 'power' regex; read them directly).
SLOPE_NAMES <- c(TH_BW_CL = "cov_BW_power_cl",
                 TH_CR_CL = "cov_CrCL_power_cl",
                 TH_BW_VC = "cov_BW_power_vc",
                 TH_SEX_VC = "cov_SEX_cat_vc")

grab <- function(fit) {
  th <- fit$theta
  list(
    TVCL_raw = exp(unname(th["lTVCL"])),   # intercept at the fit's own center
    TVVc_raw = exp(unname(th["lTVVc"])),
    th_bw_cl = unname(th[SLOPE_NAMES["TH_BW_CL"]]),
    th_cr_cl = unname(th[SLOPE_NAMES["TH_CR_CL"]]),
    th_bw_vc = unname(th[SLOPE_NAMES["TH_BW_VC"]]),
    th_sex   = unname(th[SLOPE_NAMES["TH_SEX_VC"]]),
    objf     = as.numeric(fit$objf)
  )
}

# ---- driver ------------------------------------------------------------------
sim_all <- readRDS(SIM_RDS)
ctrl    <- make_est_control("focei", "bobyqa", tier = "final")$ctrl

rows <- list()
for (dsid in seq_len(N_DATASETS)) {
  cat(sprintf("\n===== dataset %d =====\n", dsid))

  ds <- to_nm_dataset(sim_all) |>
    dplyr::filter(DATASET == dsid) |>
    dplyr::select(-SCENARIO, -DATASET) |>
    dplyr::mutate(ID = as.integer(ID),
                  SEX = as.integer(SEX),
                  RACE = as.integer(RACE)) |>
    dplyr::filter(!(EVID == 0L & TIME == 0))

  subj  <- ds[!duplicated(ds$ID), ]
  medBW <- median(subj$BW)
  medCR <- median(subj$CrCL)
  cat(sprintf("  median BW = %.4f   median CrCL = %.4f\n", medBW, medCR))

  ui_forced <- build_scn16_ui(base_2cmt_oral_linCmt, REF_BW, REF_CRCL)
  ui_median <- build_scn16_ui(base_2cmt_oral_linCmt, medBW, medCR)

  fitA <- tryCatch(nlmixr2(ui_forced, ds, est = "focei", control = ctrl),
                   error = function(e) {cat("  [A forced] FAIL:", conditionMessage(e), "\n"); NULL})
  fitB <- tryCatch(nlmixr2(ui_median, ds, est = "focei", control = ctrl),
                   error = function(e) {cat("  [B median] FAIL:", conditionMessage(e), "\n"); NULL})
  if (is.null(fitA) || is.null(fitB)) next

  A <- grab(fitA); B <- grab(fitB)

  # Back-transform run B's intercepts to the DGP reference centers.
  #   CL depends on BW (@70) and CrCL (@95); Vc depends on BW (@70) only.
  TVCL_bt <- B$TVCL_raw *
    (REF_BW  / medBW)^B$th_bw_cl *
    (REF_CRCL / medCR)^B$th_cr_cl
  TVVc_bt <- B$TVVc_raw *
    (REF_BW / medBW)^B$th_bw_vc

  rows[[length(rows) + 1]] <- tibble::tibble(
    dataset   = dsid,
    medBW     = medBW,   medCrCL = medCR,
    # intercepts
    TVCL_forced = A$TVCL_raw,  TVCL_median_raw = B$TVCL_raw,  TVCL_backtrans = TVCL_bt,
    TVVc_forced = A$TVVc_raw,  TVVc_median_raw = B$TVVc_raw,  TVVc_backtrans = TVVc_bt,
    # slope invariance (A vs B)
    th_bw_cl_A = A$th_bw_cl,   th_bw_cl_B = B$th_bw_cl,
    th_cr_cl_A = A$th_cr_cl,   th_cr_cl_B = B$th_cr_cl,
    th_bw_vc_A = A$th_bw_vc,   th_bw_vc_B = B$th_bw_vc,
    th_sex_A   = A$th_sex,     th_sex_B   = B$th_sex,
    objf_A = A$objf,           objf_B = B$objf
  )
}

res <- dplyr::bind_rows(rows)

# ---- summaries ---------------------------------------------------------------
cat("\n\n================ INTERCEPT RECOVERY (target: TVCL=0.6, TVVc=20) ================\n")
res |>
  dplyr::transmute(
    dataset,
    TVCL_forced, TVCL_backtrans,
    TVCL_bt_err_pct = 100 * (TVCL_backtrans - TVCL_forced) / TVCL_forced,
    TVVc_forced, TVVc_backtrans,
    TVVc_bt_err_pct = 100 * (TVVc_backtrans - TVVc_forced) / TVVc_forced
  ) |>
  as.data.frame() |>
  print(digits = 5)

cat("\n---- naive (uncorrected) median intercept vs forced, % high ----\n")
res |>
  dplyr::transmute(
    dataset,
    TVCL_naive_err_pct = 100 * (TVCL_median_raw - TVCL_forced) / TVCL_forced,
    TVVc_naive_err_pct = 100 * (TVVc_median_raw - TVVc_forced) / TVVc_forced
  ) |>
  as.data.frame() |>
  print(digits = 4)

cat("\n================ SLOPE INVARIANCE (A forced vs B median) ================\n")
res |>
  dplyr::transmute(
    dataset,
    th_bw_cl_A, th_bw_cl_B,
    th_cr_cl_A, th_cr_cl_B,
    th_bw_vc_A, th_bw_vc_B,
    th_sex_A,   th_sex_B,
    dObjf = objf_B - objf_A
  ) |>
  as.data.frame() |>
  print(digits = 5)

cat("\nDone. `res` holds the full per-dataset comparison.\n")
