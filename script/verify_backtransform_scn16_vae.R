# ==============================================================================
# verify_backtransform_scn16_vae.R
# ------------------------------------------------------------------------------
# PURPOSE
#   VAE analog of verify_backtransform_scn16.R.  The VAE encodes a continuous
#   covariate as log(v / mean(v)) -- i.e. it CENTERS ON THE SAMPLE MEAN, not the
#   median (runSCM) and not the DGP clinical reference (70 / 95).  Source:
#   nlmixr2est/R/vaeData.R, .vaeCovariateSearch():
#       .covPop[j] <- mean(v); .covMat[, j] <- log(v / .covPop[j])
#   Categorical indicators (0/1) stay raw (center 0), same as SEX here.
#
#   This script verifies the SAME back-transform identity for the VAE's mean
#   center:
#       TV@ref = exp(l_c) * prod_k (ref_k / c_k)^theta_k        (continuous only)
#   with c_k = mean(cov_k).  To isolate the centering question from the VAE's
#   torch estimator, we hold the ESTIMATOR fixed (FOCEi/bobyqa) and only swap
#   the centering constant.  Any estimator-specific behaviour of the actual VAE
#   is orthogonal to whether the mean-center back-transform is algebraically
#   correct -- which is exactly what is under test.
#
# THREE centering variants, identical model + estimator otherwise:
#     (A) FORCED : log(BW/70),      log(CrCL/95)        -> TV at DGP reference
#     (B) MEDIAN : log(BW/medBW),   log(CrCL/medCrCL)   -> runSCM default
#     (C) MEAN   : log(BW/meanBW),  log(CrCL/meanCrCL)  -> VAE convention
#
#   Back-transform (B) and (C) to the reference and check both reproduce (A)
#   and the DGP truth (TVCL=0.6, TVVc=20).  Slopes must be invariant across all
#   three centers.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(nlmixr2est)
})

sd_ <- "script"
source(file.path(sd_, "scm_bench_helpers.R"), chdir = FALSE)  # base models, to_nm_dataset
source(file.path(sd_, "estimator_factory.R"), chdir = FALSE)  # make_est_control

if (!exists("DOSE_MG")) DOSE_MG <- 100
REF_BW     <- 70
REF_CRCL   <- 95
N_DATASETS <- 5
SIM_RDS    <- "Inputdataset/sim_obs_N300/sim_obs_scenario_16.rds"

TRUE_TVCL <- 0.6
TRUE_TVVc <- 20

# ---- covariate UI builder (direct covExpr injection; shape-agnostic) ---------
build_scn16_ui <- function(base_fn, cBW, cCRCL) {
  base_ui <- rxode2::assertRxUi(base_fn)
  pairs_df <- data.frame(
    var     = c("cl", "cl", "vc", "vc"),
    covar   = c("BW_power", "CrCL_power", "BW_power", "SEX_cat"),
    covExpr = c(sprintf("log(BW/%.10g)",   cBW),
                sprintf("log(CrCL/%.10g)", cCRCL),
                sprintf("log(BW/%.10g)",   cBW),
                "SEX"),
    init    = c(0.75, 0.5, 1.0, 0.4),
    lower   = rep(-5, 4),
    upper   = rep( 5, 4),
    stringsAsFactors = FALSE
  )
  nlmixr2scm:::.rebuildUiFromPairs(base_ui, pairs_df)
}

SLOPE_NAMES <- c(TH_BW_CL = "cov_BW_power_cl",
                 TH_CR_CL = "cov_CrCL_power_cl",
                 TH_BW_VC = "cov_BW_power_vc",
                 TH_SEX_VC = "cov_SEX_cat_vc")

grab <- function(fit) {
  th <- fit$theta
  list(
    TVCL_raw = exp(unname(th["lTVCL"])),
    TVVc_raw = exp(unname(th["lTVVc"])),
    th_bw_cl = unname(th[SLOPE_NAMES["TH_BW_CL"]]),
    th_cr_cl = unname(th[SLOPE_NAMES["TH_CR_CL"]]),
    th_bw_vc = unname(th[SLOPE_NAMES["TH_BW_VC"]]),
    th_sex   = unname(th[SLOPE_NAMES["TH_SEX_VC"]]),
    objf     = as.numeric(fit$objf)
  )
}

# Back-transform CL and Vc intercepts from an arbitrary center to the reference.
bt <- function(g, cBW, cCR) {
  list(
    TVCL = g$TVCL_raw * (REF_BW / cBW)^g$th_bw_cl * (REF_CRCL / cCR)^g$th_cr_cl,
    TVVc = g$TVVc_raw * (REF_BW / cBW)^g$th_bw_vc
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

  subj    <- ds[!duplicated(ds$ID), ]
  medBW   <- median(subj$BW);  medCR   <- median(subj$CrCL)
  meanBW  <- mean(subj$BW);    meanCR  <- mean(subj$CrCL)
  cat(sprintf("  BW  median=%.4f mean=%.4f | CrCL median=%.4f mean=%.4f\n",
              medBW, meanBW, medCR, meanCR))

  fitA <- tryCatch(nlmixr2(build_scn16_ui(base_2cmt_oral_linCmt, REF_BW,  REF_CRCL),  ds, "focei", control = ctrl),
                   error = function(e) {cat("  [A forced] FAIL:", conditionMessage(e), "\n"); NULL})
  fitB <- tryCatch(nlmixr2(build_scn16_ui(base_2cmt_oral_linCmt, medBW,   medCR),     ds, "focei", control = ctrl),
                   error = function(e) {cat("  [B median] FAIL:", conditionMessage(e), "\n"); NULL})
  fitC <- tryCatch(nlmixr2(build_scn16_ui(base_2cmt_oral_linCmt, meanBW,  meanCR),    ds, "focei", control = ctrl),
                   error = function(e) {cat("  [C mean]   FAIL:", conditionMessage(e), "\n"); NULL})
  if (is.null(fitA) || is.null(fitB) || is.null(fitC)) next

  A <- grab(fitA); B <- grab(fitB); C <- grab(fitC)
  Bbt <- bt(B, medBW,  medCR)
  Cbt <- bt(C, meanBW, meanCR)

  rows[[length(rows) + 1L]] <- tibble::tibble(
    dataset = dsid,
    medBW = medBW, meanBW = meanBW, medCrCL = medCR, meanCrCL = meanCR,
    # intercepts at own center
    TVCL_forced = A$TVCL_raw, TVCL_median_raw = B$TVCL_raw, TVCL_mean_raw = C$TVCL_raw,
    TVVc_forced = A$TVVc_raw, TVVc_median_raw = B$TVVc_raw, TVVc_mean_raw = C$TVVc_raw,
    # back-transformed to reference
    TVCL_median_bt = Bbt$TVCL, TVCL_mean_bt = Cbt$TVCL,
    TVVc_median_bt = Bbt$TVVc, TVVc_mean_bt = Cbt$TVVc,
    # slopes across the three centers (VAE-relevant = mean, C)
    th_bw_cl_A = A$th_bw_cl, th_bw_cl_C = C$th_bw_cl,
    th_cr_cl_A = A$th_cr_cl, th_cr_cl_C = C$th_cr_cl,
    th_bw_vc_A = A$th_bw_vc, th_bw_vc_C = C$th_bw_vc,
    th_sex_A   = A$th_sex,   th_sex_C   = C$th_sex,
    objf_A = A$objf, objf_C = C$objf
  )
}

res_vae <- dplyr::bind_rows(rows)

# ---- summaries ---------------------------------------------------------------
cat("\n\n===== VAE (mean-center) BACK-TRANSFORM RECOVERY vs forced ref =====\n")
res_vae |>
  dplyr::transmute(
    dataset,
    TVCL_forced, TVCL_mean_bt,
    TVCL_bt_err_pct = 100 * (TVCL_mean_bt - TVCL_forced) / TVCL_forced,
    TVVc_forced, TVVc_mean_bt,
    TVVc_bt_err_pct = 100 * (TVVc_mean_bt - TVVc_forced) / TVVc_forced
  ) |> as.data.frame() |> print(digits = 5)

cat("\n---- naive UNCORRECTED mean-center intercept vs forced, % high ----\n")
res_vae |>
  dplyr::transmute(
    dataset,
    TVCL_naive_err_pct = 100 * (TVCL_mean_raw - TVCL_forced) / TVCL_forced,
    TVVc_naive_err_pct = 100 * (TVVc_mean_raw - TVVc_forced) / TVVc_forced
  ) |> as.data.frame() |> print(digits = 4)

cat("\n---- median vs mean back-transform (both should hit the reference) ----\n")
res_vae |>
  dplyr::transmute(
    dataset,
    TVCL_median_bt, TVCL_mean_bt,
    TVVc_median_bt, TVVc_mean_bt
  ) |> as.data.frame() |> print(digits = 5)

cat("\n===== SLOPE INVARIANCE (forced A vs mean C) =====\n")
res_vae |>
  dplyr::transmute(
    dataset,
    th_bw_cl_A, th_bw_cl_C,
    th_cr_cl_A, th_cr_cl_C,
    th_bw_vc_A, th_bw_vc_C,
    th_sex_A,   th_sex_C,
    dObjf = objf_C - objf_A
  ) |> as.data.frame() |> print(digits = 5)

cat("\nDone. `res_vae` holds the full per-dataset comparison.\n")
