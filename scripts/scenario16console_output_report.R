# R 4.5.3 started.
# 
# R version 4.5.3 (2026-03-11 ucrt) -- "Reassured Reassurer"
# Copyright (C) 2026 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64
# 
# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.
# 
#   Natural language support but running in an English locale
# 
# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.
# 
# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.
# 
library(testthat)
library(nlmixr2scm)
devtools::load_all("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
# ℹ Loading nlmixr2scm
library(nlmixr2utils)
rxode2::rxClean()
library(nlmixr2)
# ── Attaching packages ─────────────────────────────────────────────────────────────────────────── nlmixr2 5.0.0 ──
# ✔ lotri        1.0.4     ✔ nlmixr2extra 5.1.0
# ✔ nlmixr2data  2.0.9     ✔ nlmixr2plot  5.0.1
# ✔ nlmixr2est   6.0.2     ✔ rxode2       5.1.3
# ── Optional Packages Loaded/Ignored ───────────────────────────────────────────────────────────── nlmixr2 5.0.0 ──
# ✖ babelmixr2     ✖ nonmem2rx
# ✖ ggPMX     ✖ posologyr
# ✖ monolix2rx     ✖ shinyMixR
# ✖ nlmixr2lib     ✖ xpose.nlmixr2
# ✖ nlmixr2rpt     
# ── Conflicts ─────────────────────────────────────────────────────────────────────────────── nlmixr2conflicts() ──
# ✖ rxode2::boxCox()     masks nlmixr2est::boxCox()
# ✖ rxode2::yeoJohnson() masks nlmixr2est::yeoJohnson()
library(rxode2)
library(tidyverse)
# ── Attaching core tidyverse packages ────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
# ✔ dplyr     1.2.1     ✔ readr     2.2.0
# ✔ forcats   1.0.1     ✔ stringr   1.6.0
# ✔ ggplot2   4.0.3     ✔ tibble    3.3.1
# ✔ lubridate 1.9.5     ✔ tidyr     1.3.2
# ✔ purrr     1.2.2     
# ── Conflicts ──────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
# ✖ readr::edition_get()   masks testthat::edition_get()
# ✖ dplyr::filter()        masks stats::filter()
# ✖ dplyr::lag()           masks stats::lag()
# ✖ readr::local_edition() masks testthat::local_edition()
# ℹ Use the conflicted package ([object Object])  to force all conflicts to become errors

scm_focei_n <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = "r,s",         # SCM doesn't need cov matrix for LRT
  maxOuterIterations = 2000,
  maxInnerIterations = 2000
)
candidate_pairs_test_corContinusous_cat <- list(
  list(var = "cl", covar = "BW",   shapes = "power"),
  list(var = "cl", covar = "CrCL", shapes = "power"),
  list(var = "vc", covar = "BW",   shapes = "power"),
  list(var = "vc", covar = "BMI",  shapes = "power"),
  list(var = "vc", covar = "CrCL", shapes = "power"),
  list(var = "vc", covar = "SEX",  shapes = "cat"),
  list(var = "vc", covar = "RACE", shapes = "cat")
)
scm16_vars       <- c("cl", "vc")
scm16_covars     <- c("BW", "CrCL", "BMI")
scm16_catvars    <- c("SEX", "RACE")
scm16_shapes     <- c("power", "lin")


#### Forward selection (cartesian product)-------------
res16_N80_fwd_auto <- runSCM_traced(
  label       = "scn16_N80_forward_auto",
  data        = ds16_01_N80,
  fit         = fit_base16_N80,
  varsVec     = scm16_vars,
  covarsVec   = scm16_covars,
  catvarsVec  = scm16_catvars,
  shapes      = scm16_shapes,
  searchType  = "forward",
  control     = scm_focei_n,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
# 
# === [scn16_N80_forward_auto] 09:37:07 starting forward  | NA candidate(s), 3 worker(s) ===
# ── SCM Summary ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -4395.917
# ℹ Base model params : 9
# ℹ Search type : forward
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (cl, vc)
# ℹ Covariates : 8
# (BW_power, BW_lin, CrCL_power, CrCL_lin, BMI_power, BMI_lin, SEX_1, RACE_0)
# ℹ Total candidates : 16
# ── Categorical covariates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ SEX: reference = '0'
# Indicators (1): SEX_1
# ℹ RACE: reference = '1'
# Indicators (1): RACE_0
# ── Relationships to test ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ cl [power]
# 2. BW_lin ~ cl [lin]
# 3. BW_power ~ vc [power]
# 4. BW_lin ~ vc [lin]
# 5. CrCL_power ~ cl [power]
# 6. CrCL_lin ~ cl [lin]
# 7. CrCL_power ~ vc [power]
# 8. CrCL_lin ~ vc [lin]
# 9. BMI_power ~ cl [power]
# 10. BMI_lin ~ cl [lin]
# 11. BMI_power ~ vc [power]
# 12. BMI_lin ~ vc [lin]
# 13. SEX_1 ~ cl [cat]
# 14. SEX_1 ~ vc [cat]
# 15. RACE_0 ~ cl [cat]
# 16. RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Proceed with SCM? [y/N]: 
# 
# ── SCM Step Summary ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.917   -4419.326    23.409     0.0000  Added
# Forward    2     CrCL_power~cl [power]   -4419.326   -4449.376    30.050     0.0000  Added
# Forward    3     BMI_power~cl [power]    -4449.376   -4460.772    11.397     0.0007  Added
# Forward    4     SEX_1~vc [cat]          -4460.772   -4472.988    12.216     0.0005  Added
# Forward    5     BW_power~cl [power]     -4472.988   -4474.801     1.813     0.1782  Not selected
# ── SCM All Candidates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.917   -4419.326    23.409     0.0000  Added
# Forward    1     CrCL_power~cl [power]   -4395.917   -4417.850    21.933     0.0000  Not selected
# Forward    1     SEX_1~vc [cat]          -4395.917   -4417.620    21.703     0.0000  Not selected
# Forward    1     BMI_power~cl [power]    -4395.917   -4410.297    14.380     0.0001  Not selected
# Forward    1     BW_power~cl [power]     -4395.917   -4409.977    14.060     0.0002  Not selected
# Forward    1     BW_lin~vc [lin]         -4395.917   -4407.907    11.991     0.0005  Not selected
# Forward    1     BMI_lin~cl [lin]        -4395.917   -4407.220    11.304     0.0008  Not selected
# Forward    1     BMI_power~vc [power]    -4395.917   -4402.681     6.764     0.0093  Not selected
# Forward    1     CrCL_lin~cl [lin]       -4395.917   -4402.052     6.135     0.0133  Not selected
# Forward    1     BMI_lin~vc [lin]        -4395.917   -4400.480     4.564     0.0327  Not selected
# Forward    1     BW_lin~cl [lin]         -4395.917   -4399.456     3.539     0.0599  Not selected
# Forward    1     RACE_0~vc [cat]         -4395.917   -4397.600     1.684     0.1945  Not selected
# Forward    1     CrCL_power~vc [power]   -4395.917   -4397.298     1.381     0.2399  Not selected
# Forward    1     SEX_1~cl [cat]          -4395.917   -4397.186     1.269     0.2599  Not selected
# Forward    1     CrCL_lin~vc [lin]       -4395.917   -4397.092     1.175     0.2784  Not selected
# Forward    1     RACE_0~cl [cat]         -4395.917   -4395.958     0.041     0.8395  Not selected
# 
# Forward    2     CrCL_power~cl [power]   -4419.326   -4449.376    30.050     0.0000  Added
# Forward    2     BW_power~cl [power]     -4419.326   -4446.817    27.491     0.0000  Not selected
# Forward    2     BMI_power~cl [power]    -4419.326   -4445.125    25.799     0.0000  Not selected
# Forward    2     BMI_lin~cl [lin]        -4419.326   -4430.887    11.562     0.0007  Not selected
# Forward    2     SEX_1~vc [cat]          -4419.326   -4430.806    11.481     0.0007  Not selected
# Forward    2     BMI_power~vc [power]    -4419.326   -4423.801     4.476     0.0344  Not selected
# Forward    2     BW_lin~cl [lin]         -4419.326   -4412.542    -6.784     1.0000  Not selected
# Forward    2     CrCL_lin~cl [lin]       -4419.326   -4414.328    -4.998     1.0000  Not selected
# Forward    2     CrCL_power~vc [power]   -4419.326   -4419.067    -0.258     1.0000  Not selected
# Forward    2     CrCL_lin~vc [lin]       -4419.326   -4409.292   -10.034     1.0000  Not selected
# Forward    2     BMI_lin~vc [lin]        -4419.326   -4155.830  -263.496     1.0000  Not selected
# Forward    2     SEX_1~cl [cat]          -4419.326   -4417.222    -2.104     1.0000  Not selected
# Forward    2     RACE_0~cl [cat]         -4419.326   -4417.892    -1.434     1.0000  Not selected
# Forward    2     RACE_0~vc [cat]         -4419.326   -4418.708    -0.618     1.0000  Not selected
# 
# Forward    3     BMI_power~cl [power]    -4449.376   -4460.772    11.397     0.0007  Added
# Forward    3     SEX_1~vc [cat]          -4449.376   -4460.175    10.800     0.0010  Not selected
# Forward    3     BW_power~cl [power]     -4449.376   -4459.550    10.175     0.0014  Not selected
# Forward    3     BMI_power~vc [power]    -4449.376   -4453.928     4.553     0.0329  Not selected
# Forward    3     RACE_0~vc [cat]         -4449.376   -4449.490     0.115     0.7348  Not selected
# Forward    3     CrCL_power~vc [power]   -4449.376   -4449.463     0.088     0.7674  Not selected
# Forward    3     BW_lin~cl [lin]         -4449.376   -4435.910   -13.465     1.0000  Not selected
# Forward    3     CrCL_lin~vc [lin]       -4449.376   -4432.652   -16.724     1.0000  Not selected
# Forward    3     BMI_lin~cl [lin]        -4449.376   -4443.717    -5.658     1.0000  Not selected
# Forward    3     BMI_lin~vc [lin]        -4449.376   -4436.301   -13.075     1.0000  Not selected
# Forward    3     SEX_1~cl [cat]          -4449.376   -4447.443    -1.933     1.0000  Not selected
# Forward    3     RACE_0~cl [cat]         -4449.376   -4449.176    -0.199     1.0000  Not selected
# 
# Forward    4     SEX_1~vc [cat]          -4460.772   -4472.988    12.216     0.0005  Added
# Forward    4     BMI_power~vc [power]    -4460.772   -4465.991     5.219     0.0223  Not selected
# Forward    4     BW_power~cl [power]     -4460.772   -4464.523     3.751     0.0528  Not selected
# Forward    4     CrCL_power~vc [power]   -4460.772   -4463.363     2.591     0.1075  Not selected
# Forward    4     RACE_0~vc [cat]         -4460.772   -4463.064     2.291     0.1301  Not selected
# Forward    4     SEX_1~cl [cat]          -4460.772   -4461.641     0.869     0.3513  Not selected
# Forward    4     RACE_0~cl [cat]         -4460.772   -4461.551     0.778     0.3776  Not selected
# Forward    4     BW_lin~cl [lin]         -4460.772   -4446.277   -14.495     1.0000  Not selected
# Forward    4     CrCL_lin~vc [lin]       -4460.772   -4445.974   -14.798     1.0000  Not selected
# Forward    4     BMI_lin~vc [lin]        -4460.772   -4446.690   -14.082     1.0000  Not selected
# 
# Forward    5     BW_power~cl [power]     -4472.988   -4474.801     1.813     0.1782  Not selected
# Forward    5     CrCL_power~vc [power]   -4472.988   -4473.779     0.790     0.3740  Not selected
# Forward    5     RACE_0~cl [cat]         -4472.988   -4473.429     0.441     0.5066  Not selected
# Forward    5     BMI_power~vc [power]    -4472.988   -4473.322     0.334     0.5634  Not selected
# Forward    5     RACE_0~vc [cat]         -4472.988   -4473.197     0.209     0.6474  Not selected
# Forward    5     BW_lin~cl [lin]         -4472.988   -4443.287   -29.701     1.0000  Not selected
# Forward    5     CrCL_lin~vc [lin]       -4472.988   -4442.698   -30.291     1.0000  Not selected
# Forward    5     BMI_lin~vc [lin]        -4472.988   -4447.857   -25.131     1.0000  Not selected
# Forward    5     SEX_1~cl [cat]          -4472.988   -4472.527    -0.461     1.0000  Not selected
# ── Final model ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Retained:
# BW_power~vc [power]
# CrCL_power~cl [power]
# BMI_power~cl [power]
# SEX_1~vc [cat]
# ℹ Final model OFV: -4472.988
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [scn16_N80_forward_auto] 10:11:42 DONE | elapsed 34.6 min (2075 s) ===
# 
# There were 21 warnings (use warnings() to see them)
warnings()
# Warning messages:
# 1: ! BW_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4412.542 > -4419.326). Accepting best available result (attempt 1/1, dOFV = -6.784).
# 2: ! CrCL_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4414.328 > -4419.326). Accepting best available result (attempt 1/1, dOFV = -4.998).
# 3: ! CrCL_power ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4419.067 > -4419.326). Accepting best available result (attempt 1/1, dOFV = -0.258).
# 4: ! CrCL_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4409.292 > -4419.326). Accepting best available result (attempt 1/1, dOFV = -10.034).
# 5: ! BMI_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4155.83 > -4419.326). Accepting best available result (attempt 1/1, dOFV = -263.496).
# 6: ! SEX_1 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4417.222 > -4419.326). Accepting best available result (attempt 1/1, dOFV = -2.104).
# 7: ! RACE_0 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4417.892 > -4419.326). Accepting best available result (attempt 1/1, dOFV = -1.434).
# 8: ! RACE_0 ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4418.708 > -4419.326). Accepting best available result (attempt 1/1, dOFV = -0.618).
# 9: ! BW_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4435.91 > -4449.376). Accepting best available result (attempt 1/1, dOFV = -13.465).
# 10: ! CrCL_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4432.652 > -4449.376). Accepting best available result (attempt 1/1, dOFV = -16.724).
# 11: ! BMI_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4443.717 > -4449.376). Accepting best available result (attempt 1/1, dOFV = -5.658).
# 12: ! BMI_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4436.301 > -4449.376). Accepting best available result (attempt 1/1, dOFV = -13.075).
# 13: ! SEX_1 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4447.443 > -4449.376). Accepting best available result (attempt 1/1, dOFV = -1.933).
# 14: ! RACE_0 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4449.176 > -4449.376). Accepting best available result (attempt 1/1, dOFV = -0.199).
# 15: ! BW_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4446.277 > -4460.772). Accepting best available result (attempt 1/1, dOFV = -14.495).
# 16: ! CrCL_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4445.974 > -4460.772). Accepting best available result (attempt 1/1, dOFV = -14.798).
# 17: ! BMI_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4446.69 > -4460.772). Accepting best available result (attempt 1/1, dOFV = -14.082).
# 18: ! BW_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4443.287 > -4472.988). Accepting best available result (attempt 1/1, dOFV = -29.701).
# 19: ! CrCL_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4442.698 > -4472.988). Accepting best available result (attempt 1/1, dOFV = -30.291).
# 20: ! BMI_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4447.857 > -4472.988). Accepting best available result (attempt 1/1, dOFV = -25.131).
# 21: ! SEX_1 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4472.527 > -4472.988). Accepting best available result (attempt 1/1, dOFV = -0.461).
t16_N80_fwd_auto    <- attr(res16_N80_fwd_auto, "elapsed_s")
saveRDS(res16_N80_fwd_auto, file.path(stage1_dir16_N80, "res_fwd_auto.rds"))
test16_N80_fwd_auto <- package_scm_result(
  "scn16_N80_forward_auto", res16_N80_fwd_auto, t16_N80_fwd_auto,
  scenario_id = 16
)
saveRDS(test16_N80_fwd_auto, file.path(stage1_dir16_N80, "test_fwd_auto.rds"))
test16_N80_fwd_auto
# $label
# [1] "scn16_N80_forward_auto"
# 
# $selected
#          step      covar var shape      objf deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames covarEffect bsvReduction
# eta.vc      1   BW_power  vc power -4419.326 23.40874 -3519.145 -3481.580        10 3.841459 1.309820e-06      yes    forward   cov_BW_power_vc   1.0990235    43.701685
# eta.cl21    2 CrCL_power  cl power -4449.376 30.05001 -3547.195 -3505.457        11 3.841459 4.210461e-08      yes    forward cov_CrCL_power_cl   0.6142960    33.315191
# eta.cl22    3  BMI_power  cl power -4460.772 11.39668 -3556.591 -3510.680        12 3.841459 7.357546e-04      yes    forward  cov_BMI_power_cl   0.7450663     6.095665
# eta.vc43    4      SEX_1  vc   cat -4472.988 12.21602 -3566.807 -3516.722        13 3.841459 4.738083e-04      yes    forward      cov_SEX_1_vc   0.2795625    15.222192
# 
# $step_hist
#          step      covar var shape      objf      deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames   covarEffect  bsvReduction
# eta.cl      1   BW_power  cl power -4409.977   14.06031796 -3509.796 -3472.232        10 3.841459 1.770399e-04       no    forward   cov_BW_power_cl  0.6900469793  1.920616e+01
# eta.cl1     1     BW_lin  cl   lin -4399.456    3.53887630 -3499.275 -3461.711        10 3.841459 5.994606e-02       no    forward     cov_BW_lin_cl  0.0013054830 -4.631446e-14
# eta.vc      1   BW_power  vc power -4419.326   23.40874072 -3519.145 -3481.580        10 3.841459 1.309820e-06      yes    forward   cov_BW_power_vc  1.0990234884  4.370168e+01
# eta.vc1     1     BW_lin  vc   lin -4407.907   11.99050609 -3507.726 -3470.162        10 3.841459 5.347227e-04       no    forward     cov_BW_lin_vc  0.0082510107  1.515788e-04
# eta.cl2     1 CrCL_power  cl power -4417.850   21.93274247 -3517.669 -3480.104        10 3.841459 2.823747e-06       no    forward cov_CrCL_power_cl  0.5089612002  2.869984e+01
# eta.cl3     1   CrCL_lin  cl   lin -4402.052    6.13482109 -3501.871 -3464.307        10 3.841459 1.325450e-02       no    forward   cov_CrCL_lin_cl  0.0009985787 -4.631446e-14
# eta.vc2     1 CrCL_power  vc power -4397.298    1.38098727 -3497.117 -3459.553        10 3.841459 2.399331e-01       no    forward cov_CrCL_power_vc  0.1524307040  8.544358e+00
# eta.vc3     1   CrCL_lin  vc   lin -4397.092    1.17501651 -3496.911 -3459.347        10 3.841459 2.783729e-01       no    forward   cov_CrCL_lin_vc  0.0013340691  5.015468e-05
# eta.cl4     1  BMI_power  cl power -4410.297   14.38020067 -3510.116 -3472.552        10 3.841459 1.493646e-04       no    forward  cov_BMI_power_cl  0.7133239611  2.477168e+01
# eta.cl5     1    BMI_lin  cl   lin -4407.220   11.30351364 -3507.039 -3469.475        10 3.841459 7.736055e-04       no    forward    cov_BMI_lin_cl  0.0198558940 -1.477100e-03
# eta.vc4     1  BMI_power  vc power -4402.681    6.76379867 -3502.500 -3464.936        10 3.841459 9.302553e-03       no    forward  cov_BMI_power_vc  0.6115931974  1.938594e+01
# eta.vc5     1    BMI_lin  vc   lin -4400.480    4.56352038 -3500.299 -3462.735        10 3.841459 3.265987e-02       no    forward    cov_BMI_lin_vc  0.0137283161  1.454456e-03
# eta.cl6     1      SEX_1  cl   cat -4397.186    1.26913478 -3497.005 -3459.441        10 3.841459 2.599291e-01       no    forward      cov_SEX_1_cl -0.1110616793  5.920366e-01
# eta.vc6     1      SEX_1  vc   cat -4417.620   21.70348080 -3517.439 -3479.875        10 3.841459 3.182134e-06       no    forward      cov_SEX_1_vc  0.4315517983  2.390614e+01
# eta.cl7     1     RACE_0  cl   cat -4395.958    0.04100203 -3495.777 -3458.213        10 3.841459 8.395340e-01       no    forward     cov_RACE_0_cl -0.0381982480 -3.611577e-02
# eta.vc7     1     RACE_0  vc   cat -4397.600    1.68354054 -3497.419 -3459.855        10 3.841459 1.944547e-01       no    forward     cov_RACE_0_vc -0.2195091970  3.314031e+00
# eta.cl8     2   BW_power  cl power -4446.817   27.49133112 -3544.636 -3502.898        11 3.841459 1.578002e-07       no    forward   cov_BW_power_cl  0.9889863276  2.669405e+01
# eta.cl11    2     BW_lin  cl   lin -4412.542   -6.78368045 -3510.361 -3468.623        11 3.841459 1.000000e+00       no    forward     cov_BW_lin_cl  0.0013054830  1.999157e+00
# eta.cl21    2 CrCL_power  cl power -4449.376   30.05001098 -3547.195 -3505.457        11 3.841459 4.210461e-08      yes    forward cov_CrCL_power_cl  0.6142959938  3.331519e+01
# eta.cl31    2   CrCL_lin  cl   lin -4414.328   -4.99792091 -3512.147 -3470.409        11 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_cl  0.0009985787  1.998089e+00
# eta.vc8     2 CrCL_power  vc power -4419.067   -0.25826434 -3516.886 -3475.148        11 3.841459 1.000000e+00       no    forward cov_CrCL_power_vc -0.0550059725 -1.096379e+00
# eta.vc11    2   CrCL_lin  vc   lin -4409.292  -10.03350800 -3507.111 -3465.373        11 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc -0.0027860703 -7.762521e+01
# eta.cl41    2  BMI_power  cl power -4445.125   25.79926168 -3542.944 -3501.206        11 3.841459 3.788334e-07       no    forward  cov_BMI_power_cl  1.0326307193  2.823060e+01
# eta.cl51    2    BMI_lin  cl   lin -4430.887   11.56166369 -3528.706 -3486.968        11 3.841459 6.732559e-04       no    forward    cov_BMI_lin_cl  0.0208207236  2.014886e+00
# eta.vc21    2  BMI_power  vc power -4423.801    4.47569660 -3521.620 -3479.882        11 3.841459 3.438018e-02       no    forward  cov_BMI_power_vc -0.7598683210  7.775842e+00
# eta.vc31    2    BMI_lin  vc   lin -4155.830 -263.49564775 -3253.649 -3211.911        11 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0162853895 -7.753985e+01
# eta.cl61    2      SEX_1  cl   cat -4417.222   -2.10377016 -3515.041 -3473.303        11 3.841459 1.000000e+00       no    forward      cov_SEX_1_cl -0.0676786318 -5.394655e+00
# eta.vc41    2      SEX_1  vc   cat -4430.806   11.48076162 -3528.625 -3486.887        11 3.841459 7.032031e-04       no    forward      cov_SEX_1_vc  0.3113991290  1.462851e+01
# eta.cl71    2     RACE_0  cl   cat -4417.892   -1.43368279 -3515.711 -3473.973        11 3.841459 1.000000e+00       no    forward     cov_RACE_0_cl -0.1291309973 -3.442525e+00
# eta.vc51    2     RACE_0  vc   cat -4418.708   -0.61801033 -3516.527 -3474.789        11 3.841459 1.000000e+00       no    forward     cov_RACE_0_vc -0.0012446642 -4.059950e+00
# eta.cl9     3   BW_power  cl power -4459.550   10.17474926 -3555.369 -3509.458        12 3.841459 1.423771e-03       no    forward   cov_BW_power_cl  0.6875057315  1.162559e+01
# eta.cl12    3     BW_lin  cl   lin -4435.910  -13.46520579 -3531.729 -3485.818        12 3.841459 1.000000e+00       no    forward     cov_BW_lin_cl  0.0013054830 -4.696127e+01
# eta.vc9     3 CrCL_power  vc power -4449.463    0.08751069 -3545.282 -3499.370        12 3.841459 7.673660e-01       no    forward cov_CrCL_power_vc  0.0390217591  2.804494e-01
# eta.vc12    3   CrCL_lin  vc   lin -4432.652  -16.72398895 -3528.471 -3482.559        12 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc  0.0001377488 -7.770817e+01
# eta.cl22    3  BMI_power  cl power -4460.772   11.39668060 -3556.591 -3510.680        12 3.841459 7.357546e-04      yes    forward  cov_BMI_power_cl  0.7450662510  6.095665e+00
# eta.cl32    3    BMI_lin  cl   lin -4443.717   -5.65844713 -3539.536 -3493.624        12 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0176752579 -4.694159e+01
# eta.vc22    3  BMI_power  vc power -4453.928    4.55259802 -3549.747 -3503.835        12 3.841459 3.286884e-02       no    forward  cov_BMI_power_vc -0.7895331256  7.243455e+00
# eta.vc32    3    BMI_lin  vc   lin -4436.301  -13.07483033 -3532.120 -3486.208        12 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0122819281 -7.770456e+01
# eta.cl42    3      SEX_1  cl   cat -4447.443   -1.93262057 -3543.262 -3497.350        12 3.841459 1.000000e+00       no    forward      cov_SEX_1_cl -0.1017437904 -4.840909e+00
# eta.vc42    3      SEX_1  vc   cat -4460.175   10.79955128 -3555.994 -3510.082        12 3.841459 1.015247e-03       no    forward      cov_SEX_1_vc  0.3047727151  1.493336e+01
# eta.cl52    3     RACE_0  cl   cat -4449.176   -0.19916433 -3544.995 -3499.084        12 3.841459 1.000000e+00       no    forward     cov_RACE_0_cl -0.0952731188 -3.753656e+00
# eta.vc52    3     RACE_0  vc   cat -4449.490    0.11474706 -3545.309 -3499.398        12 3.841459 7.348031e-01       no    forward     cov_RACE_0_vc  0.0276511974  1.298605e+00
# eta.cl10    4   BW_power  cl power -4464.523    3.75096729 -3558.342 -3508.257        13 3.841459 5.277696e-02       no    forward   cov_BW_power_cl  0.3414203235  1.142124e+01
# eta.cl13    4     BW_lin  cl   lin -4446.277  -14.49544396 -3540.096 -3490.010        13 3.841459 1.000000e+00       no    forward     cov_BW_lin_cl  0.0013054830 -5.650276e+01
# eta.vc10    4 CrCL_power  vc power -4463.363    2.59096852 -3557.182 -3507.097        13 3.841459 1.074746e-01       no    forward cov_CrCL_power_vc  0.0182662448 -2.380539e+00
# eta.vc13    4   CrCL_lin  vc   lin -4445.974  -14.79839898 -3539.793 -3489.707        13 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc -0.0004384248 -8.346818e+01
# eta.vc23    4  BMI_power  vc power -4465.991    5.21851116 -3559.810 -3509.724        13 3.841459 2.234768e-02       no    forward  cov_BMI_power_vc -0.7278845145  5.638787e+00
# eta.vc33    4    BMI_lin  vc   lin -4446.690  -14.08180631 -3540.509 -3490.424        13 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0045055498 -8.329374e+01
# eta.cl23    4      SEX_1  cl   cat -4461.641    0.86859471 -3555.460 -3505.374        13 3.841459 3.513446e-01       no    forward      cov_SEX_1_cl -0.0574827874  7.264655e-01
# eta.vc43    4      SEX_1  vc   cat -4472.988   12.21602257 -3566.807 -3516.722        13 3.841459 4.738083e-04      yes    forward      cov_SEX_1_vc  0.2795625388  1.522219e+01
# eta.cl33    4     RACE_0  cl   cat -4461.551    0.77839644 -3555.370 -3505.284        13 3.841459 3.776320e-01       no    forward     cov_RACE_0_cl  0.0425240688  5.176501e+00
# eta.vc53    4     RACE_0  vc   cat -4463.064    2.29128810 -3556.883 -3506.797        13 3.841459 1.301019e-01       no    forward     cov_RACE_0_vc  0.0075359762 -4.958800e+00
# eta.cl14    5   BW_power  cl power -4474.801    1.81259775 -3566.620 -3512.361        14 3.841459 1.781969e-01       no    forward   cov_BW_power_cl  0.5087533243  3.885274e+00
# eta.cl15    5     BW_lin  cl   lin -4443.287  -29.70101362 -3535.106 -3480.847        14 3.841459 1.000000e+00       no    forward     cov_BW_lin_cl  0.0013054830 -6.353133e+01
# eta.vc14    5 CrCL_power  vc power -4473.779    0.79043488 -3565.598 -3511.338        14 3.841459 3.739687e-01       no    forward cov_CrCL_power_vc  0.0654826294 -3.024820e+00
# eta.vc15    5   CrCL_lin  vc   lin -4442.698  -30.29058484 -3534.517 -3480.257        14 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc  0.0002280385 -1.164373e+02
# eta.vc24    5  BMI_power  vc power -4473.322    0.33388057 -3565.141 -3510.882        14 3.841459 5.633830e-01       no    forward  cov_BMI_power_vc  0.3309706853 -5.668959e+00
# eta.vc34    5    BMI_lin  vc   lin -4447.857  -25.13096168 -3539.676 -3485.417        14 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0025959731 -1.163146e+02
# eta.cl24    5      SEX_1  cl   cat -4472.527   -0.46112420 -3564.346 -3510.087        14 3.841459 1.000000e+00       no    forward      cov_SEX_1_cl -0.0032807318 -1.019441e+00
# eta.cl34    5     RACE_0  cl   cat -4473.429    0.44114203 -3565.248 -3510.989        14 3.841459 5.065718e-01       no    forward     cov_RACE_0_cl  0.0405829507  4.106049e+00
# eta.vc44    5     RACE_0  vc   cat -4473.197    0.20916377 -3565.016 -3510.757        14 3.841459 6.474236e-01       no    forward     cov_RACE_0_vc -0.0202027868 -3.522031e-01
# 
# $final_fit
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#            OBJF       AIC       BIC Log-likelihood Condition#(Cov) Condition#(Cor)
# FOCEi -4472.988 -3566.807 -3516.722       1795.404        369.2086        7.793718
# 
# ── Time (sec $time): ──
# 
#             setup optimize covariance preprocess postprocess    other
# elapsed 0.0669905 94.57312   35.48047       0.17        0.06 15.90942
# 
# ── Population Parameters ($parFixed or $parFixedDf): ──
# 
#                      Est.      SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# lTVCL             -0.3821 0.03583  9.378 0.6824 (0.6361, 0.7321)     34.1      2.81% 
# lTVQ               0.6209 0.01857  2.991    1.861 (1.794, 1.929)                     
# lTVVc               3.158 0.05164  1.635    23.53 (21.26, 26.04)     30.1      13.1% 
# lTVVp               4.397 0.01153 0.2622    81.22 (79.41, 83.08)                     
# lTVKA             -0.3567   FIXED  FIXED                 -0.3567                     
# prop.err           0.1063                                 0.1063                     
# cov_BW_power_vc    0.9935  0.1654  16.65  0.9935 (0.6693, 1.318)                     
# cov_CrCL_power_cl  0.4421 0.09766  22.09 0.4421 (0.2506, 0.6335)                     
# cov_BMI_power_cl   0.6544  0.1671  25.53 0.6544 (0.3269, 0.9819)                     
# cov_SEX_1_vc       0.2796 0.07535  26.95 0.2796 (0.1319, 0.4272)                     
#  
#   Covariance Type ($covMethod): r,s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.283   
#  
# 
#   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found ($runInfo):
#    • gradient problems with covariance; see $scaleInfo 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring ($censInformation): No censoring
#   Minimization message ($message):  
#     Normal exit from bobyqa 
# 
# $final_est
# # A tibble: 17 × 2
#    parameter estimate
#    <chr>        <dbl>
#  1 TVCL        0.682 
#  2 TVVc       23.5   
#  3 TVQ         1.86  
#  4 TVVp       81.2   
#  5 TVKA        0.7   
#  6 var_CL      0.110 
#  7 var_Vc      0.0869
#  8 cov_VcCL    0.0277
#  9 ResErr      0.106 
# 10 CLBW       NA     
# 11 CLcrCL      0.442 
# 12 VcBW        0.993 
# 13 VcSEX       0.280 
# 14 CLBMI       0.654 
# 15 VcBMI      NA     
# 16 VcCrCL     NA     
# 17 VcRACE     NA     
# 
# $rel_err
# # A tibble: 17 × 6
#    parameter true_value estimate   abs_err   rel_err rel_err_pct
#    <chr>          <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
#  1 TVCL           0.6     0.682   8.24e- 2  1.37e- 1    1.37e+ 1
#  2 TVVc          20      23.5     3.53e+ 0  1.76e- 1    1.76e+ 1
#  3 TVQ            1.8     1.86    6.05e- 2  3.36e- 2    3.36e+ 0
#  4 TVVp          80      81.2     1.22e+ 0  1.53e- 2    1.53e+ 0
#  5 TVKA           0.7     0.7     3.33e-16  4.76e-16    4.76e-14
#  6 CLBW           0.75   NA      NA        NA          NA       
#  7 CLcrCL         0.5     0.442  -5.79e- 2 -1.16e- 1   -1.16e+ 1
#  8 VcBW           1       0.993  -6.54e- 3 -6.54e- 3   -6.54e- 1
#  9 VcSEX          0.405   0.280  -1.26e- 1 -3.11e- 1   -3.11e+ 1
# 10 var_CL         0.1     0.110   9.94e- 3  9.94e- 2    9.94e+ 0
# 11 var_Vc         0.1     0.0869 -1.31e- 2 -1.31e- 1   -1.31e+ 1
# 12 cov_VcCL       0.02    0.0277  7.68e- 3  3.84e- 1    3.84e+ 1
# 13 ResErr         0.1     0.106   6.26e- 3  6.26e- 2    6.26e+ 0
# 14 CLBMI         NA       0.654  NA        NA          NA       
# 15 VcBMI         NA      NA      NA        NA          NA       
# 16 VcCrCL        NA      NA      NA        NA          NA       
# 17 VcRACE        NA      NA      NA        NA          NA       
# 
# $diag
# $diag$converged
# [1] TRUE
# 
# $diag$objf
# [1] -4472.988
# 
# $diag$cond_num
# [1] 369.2086
# 
# $diag$cond_num_sqrt
# [1] 19.2148
# 
# $diag$cov_ok
# [1] TRUE
# 
# $diag$message
# [1] "Normal exit from bobyqa"
# 
# 
# $parFixed
# NULL
# 
# $final_fit_refit
# NULL
# 
# $runtime_sec
# [1] 2074.707


### Backward elimination (includedrelations)-----------
res16_N80_bck <- runSCM_traced(
  label             = "scn16_N80_backward",
  data              = ds16_01_N80,
  fit               = fit_base16_N80,
  pairsVec          = candidate_pairs_test_corContinusous_cat,
  catvarsVec        = scm16_catvars,
  searchType        = "backward",
  includedRelations = candidate_pairs_test_corContinusous_cat,
  control           = scm_focei_n,
  saveModels        = FALSE,
  workers           = 3L,
  print             = 100,
  maxRetries        = 0L
)
# 
# === [scn16_N80_backward] 10:38:17 starting backward | 7 candidate(s), 3 worker(s) ===
# ── SCM Summary ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -4395.917
# ℹ Base model params : 9
# ℹ Search type : backward
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (cl, vc)
# ℹ Covariates : 5
# (BW_power, CrCL_power, BMI_power, SEX_1, RACE_0)
# ℹ Total candidates : 7
# ── Categorical covariates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ SEX: reference = '0'
# Indicators (1): SEX_1
# ℹ RACE: reference = '1'
# Indicators (1): RACE_0
# ── Relationships to test ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ cl [power]
# 2. CrCL_power ~ cl [power]
# 3. BW_power ~ vc [power]
# 4. BMI_power ~ vc [power]
# 5. CrCL_power ~ vc [power]
# 6. SEX_1 ~ vc [cat]
# 7. RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Proceed with SCM? [y/N]: 

# ── SCM Step Summary ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Backward   1     BMI_power~vc [power]    -4474.442   -4474.449    -0.007     1.0000  Removed
# Backward   2     CrCL_power~vc [power]   -4474.449   -4474.088     0.361     0.5482  Removed
# Backward   3     RACE_0~vc [cat]         -4474.088   -4474.661    -0.572     1.0000  Removed
# Backward   4     SEX_1~vc [cat]          -4474.661   -4463.204    11.457     0.0007  Retained
# ── SCM All Candidates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Backward   1     BW_power~cl [power]     -4474.442   -4460.494    13.948     0.0002  Retained
# Backward   1     CrCL_power~cl [power]   -4474.442   -4461.948    12.494     0.0004  Retained
# Backward   1     SEX_1~vc [cat]          -4474.442   -4466.433     8.009     0.0047  Retained
# Backward   1     BW_power~vc [power]     -4474.442   -4470.988     3.454     0.0631  Retained
# Backward   1     CrCL_power~vc [power]   -4474.442   -4474.226     0.215     0.6425  Retained
# Backward   1     RACE_0~vc [cat]         -4474.442   -4474.306     0.135     0.7128  Retained
# Backward   1     BMI_power~vc [power]    -4474.442   -4474.449    -0.007     1.0000  Removed
# 
# Backward   2     BW_power~vc [power]     -4474.449   -4457.794    16.655     0.0000  Retained
# Backward   2     SEX_1~vc [cat]          -4474.449   -4460.167    14.282     0.0002  Retained
# Backward   2     BW_power~cl [power]     -4474.449   -4461.110    13.339     0.0003  Retained
# Backward   2     CrCL_power~cl [power]   -4474.449   -4462.590    11.859     0.0006  Retained
# Backward   2     RACE_0~vc [cat]         -4474.449   -4473.716     0.733     0.3918  Retained
# Backward   2     CrCL_power~vc [power]   -4474.449   -4474.088     0.361     0.5482  Removed
# 
# Backward   3     SEX_1~vc [cat]          -4474.088   -4206.019   268.069     0.0000  Retained
# Backward   3     BW_power~vc [power]     -4474.088   -4448.899    25.189     0.0000  Retained
# Backward   3     BW_power~cl [power]     -4474.088   -4460.125    13.963     0.0002  Retained
# Backward   3     CrCL_power~cl [power]   -4474.088   -4462.193    11.895     0.0006  Retained
# Backward   3     RACE_0~vc [cat]         -4474.088   -4474.661    -0.572     1.0000  Removed
# 
# Backward   4     BW_power~vc [power]     -4474.661   -4447.761    26.900     0.0000  Retained
# Backward   4     BW_power~cl [power]     -4474.661   -4460.632    14.028     0.0002  Retained
# Backward   4     CrCL_power~cl [power]   -4474.661   -4462.439    12.222     0.0005  Retained
# Backward   4     SEX_1~vc [cat]          -4474.661   -4463.204    11.457     0.0007  Retained
# ── Final model ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Removed:
# BMI_power~vc [power]
# CrCL_power~vc [power]
# RACE_0~vc [cat]
# ℹ Final model OFV: -4474.661
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [scn16_N80_backward] 11:00:33 DONE | elapsed 22.3 min (1336 s) ===
# 
t16_N80_bck    <- attr(res16_N80_bck, "elapsed_s")
saveRDS(res16_N80_bck, file.path(stage1_dir16_N80, "res_bck.rds"))
test16_N80_bck <- package_scm_result(
  "scn16_N80_backward_only", res16_N80_bck, t16_N80_bck,
  scenario_id = 16
)
saveRDS(test16_N80_bck, file.path(stage1_dir16_N80, "test_bck.rds"))
test16_N80_bck
# $label
# [1] "scn16_N80_backward_only"
# 
# $selected
#          step      covar var shape      objf    deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames covarEffect bsvReduction
# eta.cl      1   BW_power  cl power -4460.494  13.9482116 -3550.313 -3491.880        15 6.634897 1.879163e-04 retained   backward   cov_BW_power_cl  0.66980238   18.0616584
# eta.cl1     1 CrCL_power  cl power -4461.948  12.4940909 -3551.767 -3493.334        15 6.634897 4.082412e-04 retained   backward cov_CrCL_power_cl  0.40261792   19.2908783
# eta.vc      1   BW_power  vc power -4470.988   3.4536562 -3560.807 -3502.374        15 6.634897 6.311202e-02 retained   backward   cov_BW_power_vc  0.80301404   14.1187136
# eta.vc2     1 CrCL_power  vc power -4474.226   0.2154624 -3564.045 -3505.612        15 6.634897 6.425193e-01 retained   backward cov_CrCL_power_vc  0.06604757   -2.5624983
# eta.vc3     1      SEX_1  vc   cat -4466.433   8.0086334 -3556.252 -3497.819        15 6.634897 4.655486e-03 retained   backward      cov_SEX_1_vc  0.34218598    5.1602269
# eta.vc4     1     RACE_0  vc   cat -4474.306   0.1354741 -3564.125 -3505.692        15 6.634897 7.128225e-01 retained   backward     cov_RACE_0_vc -0.03766625   -4.5091165
# eta.cl2     2   BW_power  cl power -4461.110  13.3387276 -3552.929 -3498.670        14 6.634897 2.599807e-04 retained   backward   cov_BW_power_cl  0.67070010   14.9521291
# eta.cl11    2 CrCL_power  cl power -4462.590  11.8593100 -3554.409 -3500.149        14 6.634897 5.737380e-04 retained   backward cov_CrCL_power_cl  0.40002069   12.7276946
# eta.vc5     2   BW_power  vc power -4457.794  16.6551665 -3549.613 -3495.354        14 6.634897 4.482804e-05 retained   backward   cov_BW_power_vc  0.93032147   29.2729456
# eta.vc21    2      SEX_1  vc   cat -4460.167  14.2822476 -3551.986 -3497.726        14 6.634897 1.573418e-04 retained   backward      cov_SEX_1_vc  0.31073711   10.3839030
# eta.vc31    2     RACE_0  vc   cat -4473.716   0.7334241 -3565.535 -3511.275        14 6.634897 3.917755e-01 retained   backward     cov_RACE_0_vc -0.03373208    0.9364553
# eta.cl3     3   BW_power  cl power -4460.125  13.9634618 -3553.944 -3503.859        13 6.634897 1.863981e-04 retained   backward   cov_BW_power_cl  0.67527796   16.3481872
# eta.cl12    3 CrCL_power  cl power -4462.193  11.8949121 -3556.013 -3505.927        13 6.634897 5.628748e-04 retained   backward cov_CrCL_power_cl  0.38580907   14.7485948
# eta.vc6     3   BW_power  vc power -4448.899  25.1891421 -3542.718 -3492.633        13 6.634897 5.197399e-07 retained   backward   cov_BW_power_vc  0.90403520   37.3507884
# eta.vc12    3      SEX_1  vc   cat -4206.019 268.0693487 -3299.838 -3249.753        13 6.634897 0.000000e+00 retained   backward      cov_SEX_1_vc  0.32578199   14.5688477
# eta.cl4     4   BW_power  cl power -4460.632  14.0284226 -3556.451 -3510.540        12 6.634897 1.800681e-04 retained   backward   cov_BW_power_cl  0.67695055   16.0630662
# eta.cl13    4 CrCL_power  cl power -4462.439  12.2215529 -3558.258 -3512.346        12 6.634897 4.724059e-04 retained   backward cov_CrCL_power_cl  0.38836920   14.0800763
# eta.vc7     4   BW_power  vc power -4447.761  26.8997054 -3543.580 -3497.668        12 6.634897 2.142915e-07 retained   backward   cov_BW_power_vc  1.00482271   34.1831131
# eta.vc13    4      SEX_1  vc   cat -4463.204  11.4565135 -3559.023 -3513.112        12 6.634897 7.124384e-04 retained   backward      cov_SEX_1_vc  0.30473293   18.1051691
# 
# $step_hist
#          step      covar var shape      objf      deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames covarEffect bsvReduction
# eta.cl      1   BW_power  cl power -4460.494  13.948211552 -3550.313 -3491.880        15 6.634897 1.879163e-04 retained   backward   cov_BW_power_cl  0.66980238   18.0616584
# eta.cl1     1 CrCL_power  cl power -4461.948  12.494090898 -3551.767 -3493.334        15 6.634897 4.082412e-04 retained   backward cov_CrCL_power_cl  0.40261792   19.2908783
# eta.vc      1   BW_power  vc power -4470.988   3.453656162 -3560.807 -3502.374        15 6.634897 6.311202e-02 retained   backward   cov_BW_power_vc  0.80301404   14.1187136
# eta.vc1     1  BMI_power  vc power -4474.449  -0.007083477 -3564.268 -3505.835        15 6.634897 1.000000e+00  dropped   backward  cov_BMI_power_vc  0.12591270   -4.9324214
# eta.vc2     1 CrCL_power  vc power -4474.226   0.215462437 -3564.045 -3505.612        15 6.634897 6.425193e-01 retained   backward cov_CrCL_power_vc  0.06604757   -2.5624983
# eta.vc3     1      SEX_1  vc   cat -4466.433   8.008633358 -3556.252 -3497.819        15 6.634897 4.655486e-03 retained   backward      cov_SEX_1_vc  0.34218598    5.1602269
# eta.vc4     1     RACE_0  vc   cat -4474.306   0.135474117 -3564.125 -3505.692        15 6.634897 7.128225e-01 retained   backward     cov_RACE_0_vc -0.03766625   -4.5091165
# eta.cl2     2   BW_power  cl power -4461.110  13.338727648 -3552.929 -3498.670        14 6.634897 2.599807e-04 retained   backward   cov_BW_power_cl  0.67070010   14.9521291
# eta.cl11    2 CrCL_power  cl power -4462.590  11.859309952 -3554.409 -3500.149        14 6.634897 5.737380e-04 retained   backward cov_CrCL_power_cl  0.40002069   12.7276946
# eta.vc5     2   BW_power  vc power -4457.794  16.655166509 -3549.613 -3495.354        14 6.634897 4.482804e-05 retained   backward   cov_BW_power_vc  0.93032147   29.2729456
# eta.vc11    2 CrCL_power  vc power -4474.088   0.360545729 -3565.907 -3511.648        14 6.634897 5.482033e-01  dropped   backward cov_CrCL_power_vc  0.05910443    5.4952047
# eta.vc21    2      SEX_1  vc   cat -4460.167  14.282247561 -3551.986 -3497.726        14 6.634897 1.573418e-04 retained   backward      cov_SEX_1_vc  0.31073711   10.3839030
# eta.vc31    2     RACE_0  vc   cat -4473.716   0.733424134 -3565.535 -3511.275        14 6.634897 3.917755e-01 retained   backward     cov_RACE_0_vc -0.03373208    0.9364553
# eta.cl3     3   BW_power  cl power -4460.125  13.963461757 -3553.944 -3503.859        13 6.634897 1.863981e-04 retained   backward   cov_BW_power_cl  0.67527796   16.3481872
# eta.cl12    3 CrCL_power  cl power -4462.193  11.894912125 -3556.013 -3505.927        13 6.634897 5.628748e-04 retained   backward cov_CrCL_power_cl  0.38580907   14.7485948
# eta.vc6     3   BW_power  vc power -4448.899  25.189142068 -3542.718 -3492.633        13 6.634897 5.197399e-07 retained   backward   cov_BW_power_vc  0.90403520   37.3507884
# eta.vc12    3      SEX_1  vc   cat -4206.019 268.069348730 -3299.838 -3249.753        13 6.634897 0.000000e+00 retained   backward      cov_SEX_1_vc  0.32578199   14.5688477
# eta.vc22    3     RACE_0  vc   cat -4474.661  -0.572257057 -3568.480 -3518.394        13 6.634897 1.000000e+00  dropped   backward     cov_RACE_0_vc -0.04371761   -4.3706273
# eta.cl4     4   BW_power  cl power -4460.632  14.028422586 -3556.451 -3510.540        12 6.634897 1.800681e-04 retained   backward   cov_BW_power_cl  0.67695055   16.0630662
# eta.cl13    4 CrCL_power  cl power -4462.439  12.221552926 -3558.258 -3512.346        12 6.634897 4.724059e-04 retained   backward cov_CrCL_power_cl  0.38836920   14.0800763
# eta.vc7     4   BW_power  vc power -4447.761  26.899705361 -3543.580 -3497.668        12 6.634897 2.142915e-07 retained   backward   cov_BW_power_vc  1.00482271   34.1831131
# eta.vc13    4      SEX_1  vc   cat -4463.204  11.456513521 -3559.023 -3513.112        12 6.634897 7.124384e-04 retained   backward      cov_SEX_1_vc  0.30473293   18.1051691
# 
# $final_fit
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#            OBJF      AIC       BIC Log-likelihood Condition#(Cov) Condition#(Cor)
# FOCEi -4474.661 -3568.48 -3518.394        1796.24        391.6634        8.294813
# 
# ── Time (sec $time): ──
# 
#             setup optimize covariance preprocess postprocess    other
# elapsed 0.0657362  37.5714   20.38189       0.19        0.05 15.41097
# 
# ── Population Parameters ($parFixed or $parFixedDf): ──
# 
#                      Est.      SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# lTVCL             -0.3793  0.0357  9.414  0.6844 (0.6381, 0.734)     33.3      1.42% 
# lTVQ               0.6201 0.01986  3.203    1.859 (1.788, 1.933)                     
# lTVVc               3.146 0.05485  1.743    23.25 (20.88, 25.89)     30.3      12.7% 
# lTVVp               4.397 0.01223 0.2781    81.19 (79.27, 83.16)                     
# lTVKA             -0.3567   FIXED  FIXED                 -0.3567                     
# prop.err           0.1054                                 0.1054                     
# cov_BW_power_cl     0.677  0.1741  25.72   0.677 (0.3357, 1.018)                     
# cov_CrCL_power_cl  0.3884  0.1014  26.12 0.3884 (0.1896, 0.5872)                     
# cov_BW_power_vc     1.005  0.1752  17.43   1.005 (0.6615, 1.348)                     
# cov_SEX_1_vc       0.3047 0.08162  26.78 0.3047 (0.1448, 0.4647)                     
#  
#   Covariance Type ($covMethod): r,s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.244   
#  
# 
#   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found ($runInfo):
#    • gradient problems with covariance; see $scaleInfo 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring ($censInformation): No censoring
#   Minimization message ($message):  
#     Normal exit from bobyqa 
# 
# $final_est
# # A tibble: 17 × 2
#    parameter estimate
#    <chr>        <dbl>
#  1 TVCL        0.684 
#  2 TVVc       23.2   
#  3 TVQ         1.86  
#  4 TVVp       81.2   
#  5 TVKA        0.7   
#  6 var_CL      0.105 
#  7 var_Vc      0.0880
#  8 cov_VcCL    0.0235
#  9 ResErr      0.105 
# 10 CLBW        0.677 
# 11 CLcrCL      0.388 
# 12 VcBW        1.00  
# 13 VcSEX       0.305 
# 14 CLBMI      NA     
# 15 VcBMI      NA     
# 16 VcCrCL     NA     
# 17 VcRACE     NA     
# 
# $rel_err
# # A tibble: 17 × 6
#    parameter true_value estimate   abs_err   rel_err rel_err_pct
#    <chr>          <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
#  1 TVCL           0.6     0.684   8.44e- 2  1.41e- 1    1.41e+ 1
#  2 TVVc          20      23.2     3.25e+ 0  1.62e- 1    1.62e+ 1
#  3 TVQ            1.8     1.86    5.92e- 2  3.29e- 2    3.29e+ 0
#  4 TVVp          80      81.2     1.19e+ 0  1.49e- 2    1.49e+ 0
#  5 TVKA           0.7     0.7     3.33e-16  4.76e-16    4.76e-14
#  6 CLBW           0.75    0.677  -7.30e- 2 -9.74e- 2   -9.74e+ 0
#  7 CLcrCL         0.5     0.388  -1.12e- 1 -2.23e- 1   -2.23e+ 1
#  8 VcBW           1       1.00    4.82e- 3  4.82e- 3    4.82e- 1
#  9 VcSEX          0.405   0.305  -1.01e- 1 -2.48e- 1   -2.48e+ 1
# 10 var_CL         0.1     0.105   5.41e- 3  5.41e- 2    5.41e+ 0
# 11 var_Vc         0.1     0.0880 -1.20e- 2 -1.20e- 1   -1.20e+ 1
# 12 cov_VcCL       0.02    0.0235  3.53e- 3  1.77e- 1    1.77e+ 1
# 13 ResErr         0.1     0.105   5.45e- 3  5.45e- 2    5.45e+ 0
# 14 CLBMI         NA      NA      NA        NA          NA       
# 15 VcBMI         NA      NA      NA        NA          NA       
# 16 VcCrCL        NA      NA      NA        NA          NA       
# 17 VcRACE        NA      NA      NA        NA          NA       
# 
# $diag
# $diag$converged
# [1] TRUE
# 
# $diag$objf
# [1] -4474.661
# 
# $diag$cond_num
# [1] 391.6634
# 
# $diag$cond_num_sqrt
# [1] 19.79049
# 
# $diag$cov_ok
# [1] TRUE
# 
# $diag$message
# [1] "Normal exit from bobyqa"
# 
# 
# $parFixed
# NULL
# 
# $final_fit_refit
# NULL
# 
# $runtime_sec
# [1] 1335.51



### 3. FullSCM (cartesian product) ---------

scm_focei_production <- nlmixr2est::foceiControl(
  sigdig             = 3,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = FALSE,
  covMethod          = "r,s",
  stickyRecalcN      = 20,                      
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl          = rxode2::rxControl(atol = 1e-8, rtol = 1e-6)
)
res16_N80_full <- runSCM_traced(
  label       = "scn16_N80_full",
  data        = ds16_01_N80,
  fit         = fit_base16_N80,
  varsVec     = scm16_vars,
  covarsVec   = scm16_covars,
  catvarsVec  = scm16_catvars,
  shapes      = scm16_shapes,
  searchType  = "scm",
  control     = scm_focei_production,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
# 
# === [scn16_N80_full] 12:22:09 starting scm      | NA candidate(s), 3 worker(s) ===
# ── SCM Summary ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -4395.917
# ℹ Base model params : 9
# ℹ Search type : scm
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (cl, vc)
# ℹ Covariates : 8
# (BW_power, BW_lin, CrCL_power, CrCL_lin, BMI_power, BMI_lin, SEX_1, RACE_0)
# ℹ Total candidates : 16
# ── Categorical covariates ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ SEX: reference = '0'
# Indicators (1): SEX_1
# ℹ RACE: reference = '1'
# Indicators (1): RACE_0
# ── Relationships to test ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ cl [power]
# 2. BW_lin ~ cl [lin]
# 3. BW_power ~ vc [power]
# 4. BW_lin ~ vc [lin]
# 5. CrCL_power ~ cl [power]
# 6. CrCL_lin ~ cl [lin]
# 7. CrCL_power ~ vc [power]
# 8. CrCL_lin ~ vc [lin]
# 9. BMI_power ~ cl [power]
# 10. BMI_lin ~ cl [lin]
# 11. BMI_power ~ vc [power]
# 12. BMI_lin ~ vc [lin]
# 13. SEX_1 ~ cl [cat]
# 14. SEX_1 ~ vc [cat]
# 15. RACE_0 ~ cl [cat]
# 16. RACE_0 ~ vc [cat]
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Proceed with SCM? [y/N]: 

# ── SCM Step Summary ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     SEX_1~vc [cat]          -4395.917   -4418.035    22.118     0.0000  Added
# Forward    2     CrCL_power~cl [power]   -4418.035   -4441.016    22.980     0.0000  Added
# Forward    3     BW_power~vc [power]     -4441.016   -4460.205    19.189     0.0000  Added
# Forward    4     BW_power~cl [power]     -4460.205   -4473.735    13.530     0.0002  Added
# Forward    5     SEX_1~cl [cat]          -4473.735   -4476.139     2.404     0.1210  Not selected
# Backward   1     SEX_1~vc [cat]          -4473.735   -4463.201    10.534     0.0012  Retained
# ── SCM All Candidates ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     SEX_1~vc [cat]          -4395.917   -4418.035    22.118     0.0000  Added
# Forward    1     CrCL_power~cl [power]   -4395.917   -4417.702    21.785     0.0000  Not selected
# Forward    1     BW_power~vc [power]     -4395.917   -4416.723    20.806     0.0000  Not selected
# Forward    1     BW_power~cl [power]     -4395.917   -4409.959    14.042     0.0002  Not selected
# Forward    1     BMI_power~cl [power]    -4395.917   -4408.143    12.226     0.0005  Not selected
# Forward    1     BW_lin~vc [lin]         -4395.917   -4407.905    11.988     0.0005  Not selected
# Forward    1     BMI_power~vc [power]    -4395.917   -4402.605     6.688     0.0097  Not selected
# Forward    1     CrCL_lin~cl [lin]       -4395.917   -4402.052     6.135     0.0133  Not selected
# Forward    1     BMI_lin~vc [lin]        -4395.917   -4400.387     4.471     0.0345  Not selected
# Forward    1     BMI_lin~cl [lin]        -4395.917   -4399.718     3.801     0.0512  Not selected
# Forward    1     BW_lin~cl [lin]         -4395.917   -4399.455     3.538     0.0600  Not selected
# Forward    1     RACE_0~vc [cat]         -4395.917   -4397.603     1.686     0.1941  Not selected
# Forward    1     SEX_1~cl [cat]          -4395.917   -4397.339     1.423     0.2330  Not selected
# Forward    1     CrCL_power~vc [power]   -4395.917   -4397.320     1.403     0.2362  Not selected
# Forward    1     CrCL_lin~vc [lin]       -4395.917   -4397.100     1.183     0.2767  Not selected
# Forward    1     RACE_0~cl [cat]         -4395.917   -4395.946     0.030     0.8632  Not selected
# 
# Forward    2     CrCL_power~cl [power]   -4418.035   -4441.016    22.980     0.0000  Added
# Forward    2     BW_power~cl [power]     -4418.035   -4435.842    17.806     0.0000  Not selected
# Forward    2     BMI_power~cl [power]    -4418.035   -4431.082    13.047     0.0003  Not selected
# Forward    2     BW_power~vc [power]     -4418.035   -4431.012    12.977     0.0003  Not selected
# Forward    2     BMI_power~vc [power]    -4418.035   -4428.918    10.883     0.0010  Not selected
# Forward    2     RACE_0~vc [cat]         -4418.035   -4419.900     1.864     0.1721  Not selected
# Forward    2     CrCL_power~vc [power]   -4418.035   -4418.611     0.575     0.4481  Not selected
# Forward    2     BW_lin~cl [lin]         -4418.035   -4402.478   -15.557     1.0000  Not selected
# Forward    2     BW_lin~vc [lin]         -4418.035   -4398.874   -19.161     1.0000  Not selected
# Forward    2     CrCL_lin~cl [lin]       -4418.035   -4404.340   -13.695     1.0000  Not selected
# Forward    2     CrCL_lin~vc [lin]       -4418.035   -4397.909   -20.127     1.0000  Not selected
# Forward    2     BMI_lin~cl [lin]        -4418.035   -4402.045   -15.991     1.0000  Not selected
# Forward    2     BMI_lin~vc [lin]        -4418.035   -4399.289   -18.747     1.0000  Not selected
# Forward    2     SEX_1~cl [cat]          -4418.035   -4416.235    -1.800     1.0000  Not selected
# Forward    2     RACE_0~cl [cat]         -4418.035   -4417.522    -0.513     1.0000  Not selected
# 
# Forward    3     BW_power~vc [power]     -4441.016   -4460.205    19.189     0.0000  Added
# Forward    3     BMI_power~vc [power]    -4441.016   -4457.261    16.245     0.0001  Not selected
# Forward    3     CrCL_power~vc [power]   -4441.016   -4449.057     8.042     0.0046  Not selected
# Forward    3     BW_power~cl [power]     -4441.016   -4448.471     7.455     0.0063  Not selected
# Forward    3     BMI_power~cl [power]    -4441.016   -4444.902     3.886     0.0487  Not selected
# Forward    3     RACE_0~vc [cat]         -4441.016   -4443.025     2.010     0.1563  Not selected
# Forward    3     RACE_0~cl [cat]         -4441.016   -4441.066     0.050     0.8226  Not selected
# Forward    3     BW_lin~cl [lin]         -4441.016   -4414.523   -26.493     1.0000  Not selected
# Forward    3     BW_lin~vc [lin]         -4441.016   -4415.640   -25.376     1.0000  Not selected
# Forward    3     CrCL_lin~vc [lin]       -4441.016   -4417.267   -23.749     1.0000  Not selected
# Forward    3     BMI_lin~cl [lin]        -4441.016   -4414.580   -26.436     1.0000  Not selected
# Forward    3     BMI_lin~vc [lin]        -4441.016   -4416.119   -24.896     1.0000  Not selected
# Forward    3     SEX_1~cl [cat]          -4441.016   -4441.005    -0.011     1.0000  Not selected
# 
# Forward    4     BW_power~cl [power]     -4460.205   -4473.735    13.530     0.0002  Added
# Forward    4     BMI_power~cl [power]    -4460.205   -4473.149    12.944     0.0003  Not selected
# Forward    4     CrCL_power~vc [power]   -4460.205   -4460.806     0.601     0.4382  Not selected
# Forward    4     RACE_0~vc [cat]         -4460.205   -4460.326     0.121     0.7279  Not selected
# Forward    4     RACE_0~cl [cat]         -4460.205   -4460.228     0.023     0.8808  Not selected
# Forward    4     BW_lin~cl [lin]         -4460.205   -4432.742   -27.463     1.0000  Not selected
# Forward    4     CrCL_lin~vc [lin]       -4460.205   -4429.444   -30.761     1.0000  Not selected
# Forward    4     BMI_lin~cl [lin]        -4460.205   -4432.851   -27.354     1.0000  Not selected
# Forward    4     BMI_power~vc [power]    -4460.205   -4459.852    -0.353     1.0000  Not selected
# Forward    4     BMI_lin~vc [lin]        -4460.205   -4433.864   -26.341     1.0000  Not selected
# Forward    4     SEX_1~cl [cat]          -4460.205   -4460.158    -0.047     1.0000  Not selected
# 
# Forward    5     SEX_1~cl [cat]          -4473.735   -4476.139     2.404     0.1210  Not selected
# Forward    5     RACE_0~cl [cat]         -4473.735   -4474.666     0.931     0.3346  Not selected
# Forward    5     RACE_0~vc [cat]         -4473.735   -4474.639     0.905     0.3416  Not selected
# Forward    5     BMI_power~cl [power]    -4473.735   -4474.576     0.841     0.3590  Not selected
# Forward    5     BMI_power~vc [power]    -4473.735   -4474.380     0.646     0.4217  Not selected
# Forward    5     CrCL_power~vc [power]   -4473.735   -4474.307     0.572     0.4495  Not selected
# Forward    5     CrCL_lin~vc [lin]       -4473.735   -4441.218   -32.517     1.0000  Not selected
# Forward    5     BMI_lin~cl [lin]        -4473.735   -4441.935   -31.799     1.0000  Not selected
# Forward    5     BMI_lin~vc [lin]        -4473.735   -4448.242   -25.493     1.0000  Not selected
# 
# Backward   1     BW_power~vc [power]     -4473.735   -4448.197    25.538     0.0000  Retained
# Backward   1     BW_power~cl [power]     -4473.735   -4460.814    12.921     0.0003  Retained
# Backward   1     CrCL_power~cl [power]   -4473.735   -4461.305    12.430     0.0004  Retained
# Backward   1     SEX_1~vc [cat]          -4473.735   -4463.201    10.534     0.0012  Retained
# ── Final model ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Retained:
# SEX_1~vc [cat]
# CrCL_power~cl [power]
# BW_power~vc [power]
# BW_power~cl [power]
# ℹ Final model OFV: -4473.735
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [scn16_N80_full] 12:53:36 DONE | elapsed 31.4 min (1887 s) ===
# 
# There were 23 warnings (use warnings() to see them)
 warnings()
# Warning messages:
# 1: ! BW_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4402.478 > -4418.035). Accepting best available result (attempt 1/1, dOFV = -15.557).
# 2: ! BW_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4398.874 > -4418.035). Accepting best available result (attempt 1/1, dOFV = -19.161).
# 3: ! CrCL_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4404.34 > -4418.035). Accepting best available result (attempt 1/1, dOFV = -13.695).
# 4: ! CrCL_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4397.909 > -4418.035). Accepting best available result (attempt 1/1, dOFV = -20.127).
# 5: ! BMI_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4402.045 > -4418.035). Accepting best available result (attempt 1/1, dOFV = -15.991).
# 6: ! BMI_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4399.289 > -4418.035). Accepting best available result (attempt 1/1, dOFV = -18.747).
# 7: ! SEX_1 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4416.235 > -4418.035). Accepting best available result (attempt 1/1, dOFV = -1.8).
# 8: ! RACE_0 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4417.522 > -4418.035). Accepting best available result (attempt 1/1, dOFV = -0.513).
# 9: ! BW_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4414.523 > -4441.016). Accepting best available result (attempt 1/1, dOFV = -26.493).
# 10: ! BW_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4415.64 > -4441.016). Accepting best available result (attempt 1/1, dOFV = -25.376).
# 11: ! CrCL_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4417.267 > -4441.016). Accepting best available result (attempt 1/1, dOFV = -23.749).
# 12: ! BMI_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4414.58 > -4441.016). Accepting best available result (attempt 1/1, dOFV = -26.436).
# 13: ! BMI_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4416.119 > -4441.016). Accepting best available result (attempt 1/1, dOFV = -24.896).
# 14: ! SEX_1 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4441.005 > -4441.016). Accepting best available result (attempt 1/1, dOFV = -0.011).
# 15: ! BW_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4432.742 > -4460.205). Accepting best available result (attempt 1/1, dOFV = -27.463).
# 16: ! CrCL_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4429.444 > -4460.205). Accepting best available result (attempt 1/1, dOFV = -30.761).
# 17: ! BMI_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4432.851 > -4460.205). Accepting best available result (attempt 1/1, dOFV = -27.354).
# 18: ! BMI_power ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4459.852 > -4460.205). Accepting best available result (attempt 1/1, dOFV = -0.353).
# 19: ! BMI_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4433.864 > -4460.205). Accepting best available result (attempt 1/1, dOFV = -26.341).
# 20: ! SEX_1 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4460.158 > -4460.205). Accepting best available result (attempt 1/1, dOFV = -0.047).
# 21: ! CrCL_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4441.218 > -4473.735). Accepting best available result (attempt 1/1, dOFV = -32.517).
# 22: ! BMI_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4441.935 > -4473.735). Accepting best available result (attempt 1/1, dOFV = -31.799).
# 23: ! BMI_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4448.242 > -4473.735). Accepting best available result (attempt 1/1, dOFV = -25.493).
t16_N80_full    <- attr(res16_N80_full, "elapsed_s")
saveRDS(res16_N80_full, file.path(stage1_dir16_N80, "res_full.rds"))
test16_N80_full <- package_scm_result(
  "scn16_N80_full_scm", res16_N80_full, t16_N80_full,
  scenario_id = 16
)
test16_N80_full
# $label
# [1] "scn16_N80_full_scm"
# 
# $selected
#          step      covar var shape      objf deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames covarEffect bsvReduction
# eta.vc6     1      SEX_1  vc   cat -4418.035 22.11846 -3517.854 -3480.290        10 3.841459 2.563327e-06      yes    forward      cov_SEX_1_vc   0.4280355     22.02129
# eta.cl21    2 CrCL_power  cl power -4441.016 22.98030 -3538.835 -3497.097        11 3.841459 1.636697e-06      yes    forward cov_CrCL_power_cl   0.5258416     31.42834
# eta.vc9     3   BW_power  vc power -4460.205 19.18945 -3556.024 -3510.112        12 3.841459 1.183656e-05      yes    forward   cov_BW_power_vc   0.8636865     36.42384
# eta.cl10    4   BW_power  cl power -4473.735 13.52989 -3567.554 -3517.468        13 3.841459 2.347935e-04      yes    forward   cov_BW_power_cl   0.6957325     15.90497
# eta.cl16    1   BW_power  cl power -4460.814 12.92091 -3556.633 -3510.721        12 6.634897 3.249322e-04 retained   backward   cov_BW_power_cl   0.6957325     13.21296
# eta.vc16    1   BW_power  vc power -4448.197 25.53795 -3544.016 -3498.104        12 6.634897 4.337657e-07 retained   backward   cov_BW_power_vc   1.0023455     35.29656
# eta.cl17    1 CrCL_power  cl power -4461.305 12.42960 -3557.124 -3511.213        12 6.634897 4.225831e-04 retained   backward cov_CrCL_power_cl   0.3811412     15.93696
# eta.vc17    1      SEX_1  vc   cat -4463.201 10.53378 -3559.020 -3513.108        12 6.634897 1.172120e-03 retained   backward      cov_SEX_1_vc   0.3103543     16.68256
# 
# $step_hist
#          step      covar var shape      objf     deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames   covarEffect  bsvReduction
# eta.cl      1   BW_power  cl power -4409.959  14.04182376 -3509.778 -3472.214        10 3.841459 1.787895e-04       no    forward   cov_BW_power_cl  0.6995574106  1.905230e+01
# eta.cl1     1     BW_lin  cl   lin -4399.455   3.53843291 -3499.274 -3461.710        10 3.841459 5.996209e-02       no    forward     cov_BW_lin_cl  0.0013054830 -4.631446e-14
# eta.vc      1   BW_power  vc power -4416.723  20.80614291 -3516.542 -3478.978        10 3.841459 5.081970e-06       no    forward   cov_BW_power_vc  0.9315814727  4.280486e+01
# eta.vc1     1     BW_lin  vc   lin -4407.905  11.98826875 -3507.724 -3470.160        10 3.841459 5.353650e-04       no    forward     cov_BW_lin_vc  0.0080191779  4.131633e-09
# eta.cl2     1 CrCL_power  cl power -4417.702  21.78508557 -3517.521 -3479.957        10 3.841459 3.049616e-06       no    forward cov_CrCL_power_cl  0.4975848794  2.833045e+01
# eta.cl3     1   CrCL_lin  cl   lin -4402.052   6.13535997 -3501.871 -3464.307        10 3.841459 1.325046e-02       no    forward   cov_CrCL_lin_cl  0.0009985787  6.373240e-04
# eta.vc2     1 CrCL_power  vc power -4397.320   1.40302124 -3497.139 -3459.575        10 3.841459 2.362184e-01       no    forward cov_CrCL_power_vc  0.1558478084  8.222043e+00
# eta.vc3     1   CrCL_lin  vc   lin -4397.100   1.18306144 -3496.919 -3459.355        10 3.841459 2.767336e-01       no    forward   cov_CrCL_lin_vc  0.0015063985  3.743230e-03
# eta.cl4     1  BMI_power  cl power -4408.143  12.22646403 -3507.962 -3470.398        10 3.841459 4.711641e-04       no    forward  cov_BMI_power_cl  0.6779588384  7.871260e+00
# eta.cl5     1    BMI_lin  cl   lin -4399.718   3.80137326 -3499.537 -3461.973        10 3.841459 5.121057e-02       no    forward    cov_BMI_lin_cl  0.0036900369 -4.631446e-14
# eta.vc4     1  BMI_power  vc power -4402.605   6.68834979 -3502.424 -3464.860        10 3.841459 9.704499e-03       no    forward  cov_BMI_power_vc  0.5997085820  2.046187e+01
# eta.vc5     1    BMI_lin  vc   lin -4400.387   4.47056843 -3500.206 -3462.642        10 3.841459 3.448352e-02       no    forward    cov_BMI_lin_vc  0.0158661073  2.858648e-07
# eta.cl6     1      SEX_1  cl   cat -4397.339   1.42268807 -3497.158 -3459.594        10 3.841459 2.329615e-01       no    forward      cov_SEX_1_cl -0.1171779182  1.984261e+00
# eta.vc6     1      SEX_1  vc   cat -4418.035  22.11846142 -3517.854 -3480.290        10 3.841459 2.563327e-06      yes    forward      cov_SEX_1_vc  0.4280354916  2.202129e+01
# eta.cl7     1     RACE_0  cl   cat -4395.946   0.02967167 -3495.765 -3458.201        10 3.841459 8.632373e-01       no    forward     cov_RACE_0_cl -0.0326594042 -8.227126e-01
# eta.vc7     1     RACE_0  vc   cat -4397.603   1.68585968 -3497.422 -3459.858        10 3.841459 1.941477e-01       no    forward     cov_RACE_0_vc -0.2190620351  3.905938e+00
# eta.cl8     2   BW_power  cl power -4435.842  17.80647842 -3533.661 -3491.923        11 3.841459 2.445492e-05       no    forward   cov_BW_power_cl  0.7743234917  2.895541e+01
# eta.cl11    2     BW_lin  cl   lin -4402.478 -15.55710416 -3500.297 -3458.559        11 3.841459 1.000000e+00       no    forward     cov_BW_lin_cl  0.0013054830  2.005322e-01
# eta.vc8     2   BW_power  vc power -4431.012  12.97654490 -3528.831 -3487.093        11 3.841459 3.154175e-04       no    forward   cov_BW_power_vc  0.7652090032  3.689476e+01
# eta.vc11    2     BW_lin  vc   lin -4398.874 -19.16107171 -3496.693 -3454.955        11 3.841459 1.000000e+00       no    forward     cov_BW_lin_vc  0.0023035747 -2.823997e+01
# eta.cl21    2 CrCL_power  cl power -4441.016  22.98030344 -3538.835 -3497.097        11 3.841459 1.636697e-06      yes    forward cov_CrCL_power_cl  0.5258415576  3.142834e+01
# eta.cl31    2   CrCL_lin  cl   lin -4404.340 -13.69543171 -3502.159 -3460.421        11 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_cl  0.0009985787  2.005322e-01
# eta.vc21    2 CrCL_power  vc power -4418.611   0.57540977 -3516.430 -3474.692        11 3.841459 4.481172e-01       no    forward cov_CrCL_power_vc  0.1361521972  1.163869e+01
# eta.vc31    2   CrCL_lin  vc   lin -4397.909 -20.12663177 -3495.728 -3453.990        11 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc  0.0004723927 -2.823987e+01
# eta.cl41    2  BMI_power  cl power -4431.082  13.04679922 -3528.901 -3487.163        11 3.841459 3.038032e-04       no    forward  cov_BMI_power_cl  0.6886322326  2.423673e+01
# eta.cl51    2    BMI_lin  cl   lin -4402.045 -15.99056272 -3499.864 -3458.126        11 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369  2.005322e-01
# eta.vc41    2  BMI_power  vc power -4428.918  10.88304889 -3526.737 -3484.999        11 3.841459 9.704836e-04       no    forward  cov_BMI_power_vc  0.6899692943  3.049462e+01
# eta.vc51    2    BMI_lin  vc   lin -4399.289 -18.74654095 -3497.108 -3455.370        11 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc  0.0094069779 -2.824096e+01
# eta.cl61    2      SEX_1  cl   cat -4416.235  -1.79980840 -3514.054 -3472.317        11 3.841459 1.000000e+00       no    forward      cov_SEX_1_cl  0.0324978091  2.663351e+00
# eta.cl71    2     RACE_0  cl   cat -4417.522  -0.51326058 -3515.341 -3473.603        11 3.841459 1.000000e+00       no    forward     cov_RACE_0_cl -0.0563050433  3.774089e+00
# eta.vc61    2     RACE_0  vc   cat -4419.900   1.86442300 -3517.719 -3475.981        11 3.841459 1.721152e-01       no    forward     cov_RACE_0_vc -0.1984795796  7.059647e+00
# eta.cl9     3   BW_power  cl power -4448.471   7.45540268 -3544.290 -3498.378        12 3.841459 6.324633e-03       no    forward   cov_BW_power_cl  0.4960898366  1.375284e+01
# eta.cl12    3     BW_lin  cl   lin -4414.523 -26.49299530 -3510.342 -3464.430        12 3.841459 1.000000e+00       no    forward     cov_BW_lin_cl  0.0013054830 -4.554039e+01
# eta.vc9     3   BW_power  vc power -4460.205  19.18945310 -3556.024 -3510.112        12 3.841459 1.183656e-05      yes    forward   cov_BW_power_vc  0.8636864970  3.642384e+01
# eta.vc12    3     BW_lin  vc   lin -4415.640 -25.37583940 -3511.459 -3465.547        12 3.841459 1.000000e+00       no    forward     cov_BW_lin_vc  0.0040567315 -3.376716e+01
# eta.vc22    3 CrCL_power  vc power -4449.057   8.04164199 -3544.876 -3498.965        12 3.841459 4.571408e-03       no    forward cov_CrCL_power_vc  0.3426113487  1.205548e+01
# eta.vc32    3   CrCL_lin  vc   lin -4417.267 -23.74889794 -3513.086 -3467.174        12 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc  0.0029581115 -3.376716e+01
# eta.cl22    3  BMI_power  cl power -4444.902   3.88628990 -3540.721 -3494.809        12 3.841459 4.868185e-02       no    forward  cov_BMI_power_cl  0.4021207027 -9.750290e+00
# eta.cl32    3    BMI_lin  cl   lin -4414.580 -26.43602459 -3510.399 -3464.487        12 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369 -4.554039e+01
# eta.vc42    3  BMI_power  vc power -4457.261  16.24527840 -3553.080 -3507.168        12 3.841459 5.564810e-05       no    forward  cov_BMI_power_vc  0.8014753478  3.116005e+01
# eta.vc52    3    BMI_lin  vc   lin -4416.119 -24.89609252 -3511.938 -3466.027        12 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc  0.0129753362 -3.376267e+01
# eta.cl42    3      SEX_1  cl   cat -4441.005  -0.01083529 -3536.824 -3490.912        12 3.841459 1.000000e+00       no    forward      cov_SEX_1_cl -0.0084251904  4.335413e-01
# eta.cl52    3     RACE_0  cl   cat -4441.066   0.05026516 -3536.885 -3490.973        12 3.841459 8.226025e-01       no    forward     cov_RACE_0_cl -0.0256432973 -3.779853e-01
# eta.vc62    3     RACE_0  vc   cat -4443.025   2.00981503 -3538.844 -3492.933        12 3.841459 1.562844e-01       no    forward     cov_RACE_0_vc -0.2156726489  6.227989e+00
# eta.cl10    4   BW_power  cl power -4473.735  13.52989247 -3567.554 -3517.468        13 3.841459 2.347935e-04      yes    forward   cov_BW_power_cl  0.6957325184  1.590497e+01
# eta.cl13    4     BW_lin  cl   lin -4432.742 -27.46337746 -3526.561 -3476.475        13 3.841459 1.000000e+00       no    forward     cov_BW_lin_cl  0.0013054830 -4.202642e+01
# eta.vc10    4 CrCL_power  vc power -4460.806   0.60099005 -3554.625 -3504.540        13 3.841459 4.382005e-01       no    forward cov_CrCL_power_vc  0.1083106787 -9.752339e-01
# eta.vc13    4   CrCL_lin  vc   lin -4429.444 -30.76083265 -3523.263 -3473.178        13 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc  0.0005423428 -1.104034e+02
# eta.cl23    4  BMI_power  cl power -4473.149  12.94404879 -3566.968 -3516.883        13 3.841459 3.209401e-04       no    forward  cov_BMI_power_cl  0.6590758392  1.581745e+01
# eta.cl33    4    BMI_lin  cl   lin -4432.851 -27.35366724 -3526.670 -3476.585        13 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369 -4.202642e+01
# eta.vc23    4  BMI_power  vc power -4459.852  -0.35266891 -3553.671 -3503.586        13 3.841459 1.000000e+00       no    forward  cov_BMI_power_vc  0.0326707957 -2.862777e-01
# eta.vc33    4    BMI_lin  vc   lin -4433.864 -26.34081861 -3527.683 -3477.598        13 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0081663821 -1.095506e+02
# eta.cl43    4      SEX_1  cl   cat -4460.158  -0.04716810 -3553.977 -3503.891        13 3.841459 1.000000e+00       no    forward      cov_SEX_1_cl -0.0157455485 -8.210606e-02
# eta.cl53    4     RACE_0  cl   cat -4460.228   0.02250538 -3554.047 -3503.961        13 3.841459 8.807505e-01       no    forward     cov_RACE_0_cl -0.0861632304  1.056077e+00
# eta.vc43    4     RACE_0  vc   cat -4460.326   0.12107565 -3554.145 -3504.060        13 3.841459 7.278708e-01       no    forward     cov_RACE_0_vc -0.0239746921  1.697130e-01
# eta.vc14    5 CrCL_power  vc power -4474.307   0.57191291 -3566.126 -3511.867        14 3.841459 4.494998e-01       no    forward cov_CrCL_power_vc  0.0519406290  1.468331e+00
# eta.vc15    5   CrCL_lin  vc   lin -4441.218 -32.51708759 -3533.037 -3478.778        14 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc  0.0001446700 -1.108123e+02
# eta.cl14    5  BMI_power  cl power -4474.576   0.84137897 -3566.395 -3512.136        14 3.841459 3.590027e-01       no    forward  cov_BMI_power_cl  0.2403941440  1.266162e+00
# eta.cl15    5    BMI_lin  cl   lin -4441.935 -31.79940339 -3533.755 -3479.495        14 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369 -6.888799e+01
# eta.vc24    5  BMI_power  vc power -4474.380   0.64551533 -3566.199 -3511.940        14 3.841459 4.217206e-01       no    forward  cov_BMI_power_vc -0.0469589379 -1.408461e+00
# eta.vc34    5    BMI_lin  vc   lin -4448.242 -25.49264868 -3540.061 -3485.802        14 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0041237348 -1.079906e+02
# eta.cl24    5      SEX_1  cl   cat -4476.139   2.40424545 -3567.958 -3513.699        14 3.841459 1.210065e-01       no    forward      cov_SEX_1_cl -0.1300668457 -7.090384e-02
#  [ reached 'max' / getOption("max.print") -- omitted 6 rows ]
# 
# $final_fit
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#            OBJF       AIC       BIC Log-likelihood Condition#(Cov) Condition#(Cor)
# FOCEi -4473.735 -3567.554 -3517.468       1795.777        413.0033        8.595527
# 
# ── Time (sec $time): ──
# 
#           setup optimize covariance preprocess postprocess    other
# elapsed 0.06373 81.00773   36.39802       0.18        0.11 22.69051
# 
# ── Population Parameters ($parFixed or $parFixedDf): ──
# 
#                     Est.     SE  %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# lTVCL             -0.379 0.0368  9.71    0.685 (0.637, 0.736)     33.5      1.76% 
# lTVQ               0.621 0.0193  3.11       1.86 (1.79, 1.93)                     
# lTVVc               3.14 0.0567   1.8       23.2 (20.8, 25.9)     30.5      13.6% 
# lTVVp                4.4 0.0125 0.285       81.2 (79.2, 83.2)                     
# lTVKA             -0.357  FIXED FIXED                  -0.357                     
# prop.err           0.106                                0.106                     
# cov_SEX_1_vc        0.31 0.0846  27.3     0.31 (0.144, 0.476)                     
# cov_CrCL_power_cl  0.381  0.109  28.7    0.381 (0.167, 0.595)                     
# cov_BW_power_vc        1  0.173  17.2         1 (0.664, 1.34)                     
# cov_BW_power_cl    0.696  0.188    27     0.696 (0.328, 1.06)                     
#  
#   Covariance Type ($covMethod): r,s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.316  
#  
# 
#   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found ($runInfo):
#    • gradient problems with covariance; see $scaleInfo 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring ($censInformation): No censoring
#   Minimization message ($message):  
#     Normal exit from bobyqa 
# 
# $final_est
# # A tibble: 17 × 2
#    parameter estimate
#    <chr>        <dbl>
#  1 TVCL        0.685 
#  2 TVVc       23.2   
#  3 TVQ         1.86  
#  4 TVVp       81.2   
#  5 TVKA        0.7   
#  6 var_CL      0.106 
#  7 var_Vc      0.0892
#  8 cov_VcCL    0.0308
#  9 ResErr      0.106 
# 10 CLBW        0.696 
# 11 CLcrCL      0.381 
# 12 VcBW        1.00  
# 13 VcSEX       0.310 
# 14 CLBMI      NA     
# 15 VcBMI      NA     
# 16 VcCrCL     NA     
# 17 VcRACE     NA     
# 
# $rel_err
# # A tibble: 17 × 6
#    parameter true_value estimate   abs_err   rel_err rel_err_pct
#    <chr>          <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
#  1 TVCL           0.6     0.685   8.45e- 2  1.41e- 1    1.41e+ 1
#  2 TVVc          20      23.2     3.22e+ 0  1.61e- 1    1.61e+ 1
#  3 TVQ            1.8     1.86    6.14e- 2  3.41e- 2    3.41e+ 0
#  4 TVVp          80      81.2     1.22e+ 0  1.53e- 2    1.53e+ 0
#  5 TVKA           0.7     0.7     3.33e-16  4.76e-16    4.76e-14
#  6 CLBW           0.75    0.696  -5.43e- 2 -7.24e- 2   -7.24e+ 0
#  7 CLcrCL         0.5     0.381  -1.19e- 1 -2.38e- 1   -2.38e+ 1
#  8 VcBW           1       1.00    2.35e- 3  2.35e- 3    2.35e- 1
#  9 VcSEX          0.405   0.310  -9.51e- 2 -2.35e- 1   -2.35e+ 1
# 10 var_CL         0.1     0.106   6.45e- 3  6.45e- 2    6.45e+ 0
# 11 var_Vc         0.1     0.0892 -1.08e- 2 -1.08e- 1   -1.08e+ 1
# 12 cov_VcCL       0.02    0.0308  1.08e- 2  5.40e- 1    5.40e+ 1
# 13 ResErr         0.1     0.106   5.82e- 3  5.82e- 2    5.82e+ 0
# 14 CLBMI         NA      NA      NA        NA          NA       
# 15 VcBMI         NA      NA      NA        NA          NA       
# 16 VcCrCL        NA      NA      NA        NA          NA       
# 17 VcRACE        NA      NA      NA        NA          NA       
# 
# $diag
# $diag$converged
# [1] TRUE
# 
# $diag$objf
# [1] -4473.735
# 
# $diag$cond_num
# [1] 413.0033
# 
# $diag$cond_num_sqrt
# [1] 20.32248
# 
# $diag$cov_ok
# [1] TRUE
# 
# $diag$message
# [1] "Normal exit from bobyqa"
# 
# 
# $parFixed
# NULL
# 
# $final_fit_refit
# NULL
# 
# $runtime_sec
# [1] 1886.927
saveRDS(test16_N80_full, file.path(stage1_dir16_N80, "test_full.rds"))


### 4.1 userUser specified relations ------------
scm_focei_production4 <- nlmixr2est::foceiControl(
  sigdig             = 4,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = FALSE,
  covMethod          = "r,s",
  stickyRecalcN      = 20,                      
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl          = rxode2::rxControl(atol = 1e-8, rtol = 1e-6) #Even with this setting, dofv still can be negative
)

res16_N80_user <- runSCM_traced(
  label       = "scn16_N80_user",
  data        = ds16_01_N80,
  fit         = fit_base16_N80,
  pairsVec    = list(list(var = "vc", covar = "BW",   shapes = "power"),
                     list(var = "cl", covar = "BW",   shapes = "power"),
                     list(var = "cl", covar = "CrCL", shapes = "power"),
                     list(var = "vc", covar = "SEX",  shapes = "cat")),
  catvarsVec  = "SEX",     # required: declares SEX so runSCM builds the SEX_<level> dummies
  searchType  = "scm",
  control     = scm_focei_production4,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
# ── SCM Step Summary ───────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.917   -4419.286    23.369     0.0000  Added
# Forward    2     BW_power~cl [power]     -4419.286   -4450.821    31.536     0.0000  Added
# Forward    3     CrCL_power~cl [power]   -4450.821   -4462.977    12.156     0.0005  Added
# Forward    4     SEX_1~vc [cat]          -4462.977   -4474.174    11.197     0.0008  Added
# Backward   1     SEX_1~vc [cat]          -4474.174   -4463.198    10.977     0.0009  Retained
# ── SCM All Candidates ─────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.917   -4419.286    23.369     0.0000  Added
# Forward    1     SEX_1~vc [cat]          -4395.917   -4418.114    22.197     0.0000  Not selected
# Forward    1     CrCL_power~cl [power]   -4395.917   -4417.807    21.891     0.0000  Not selected
# Forward    1     BW_power~cl [power]     -4395.917   -4409.699    13.782     0.0002  Not selected
# 
# Forward    2     BW_power~cl [power]     -4419.286   -4450.821    31.536     0.0000  Added
# Forward    2     CrCL_power~cl [power]   -4419.286   -4445.905    26.619     0.0000  Not selected
# Forward    2     SEX_1~vc [cat]          -4419.286   -4431.047    11.761     0.0006  Not selected
# 
# Forward    3     CrCL_power~cl [power]   -4450.821   -4462.977    12.156     0.0005  Added
# Forward    3     SEX_1~vc [cat]          -4450.821   -4461.750    10.929     0.0009  Not selected
# 
# Forward    4     SEX_1~vc [cat]          -4462.977   -4474.174    11.197     0.0008  Added
# 
# Backward   1     BW_power~vc [power]     -4474.174   -4448.458    25.716     0.0000  Retained
# Backward   1     BW_power~cl [power]     -4474.174   -4460.474    13.700     0.0002  Retained
# Backward   1     CrCL_power~cl [power]   -4474.174   -4462.724    11.451     0.0007  Retained
# Backward   1     SEX_1~vc [cat]          -4474.174   -4463.198    10.977     0.0009  Retained
# ── Final model ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Retained:
# BW_power~vc [power]
# BW_power~cl [power]
# CrCL_power~cl [power]
# SEX_1~vc [cat]
# ℹ Final model OFV: -4474.174
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [scn16_N80_user] 15:19:59 DONE | elapsed 9.6 min (576 s) ===
# 
t16_N80_user    <- attr(res16_N80_user, "elapsed_s")
saveRDS(res16_N80_user, file.path(stage1_dir16_N80, "res_user.rds"))
test16_N80_user <- package_scm_result(
  "scn16_N80_user_BWonVc", res16_N80_user, t16_N80_user,
  scenario_id = 16
)
saveRDS(test16_N80_user, file.path(stage1_dir16_N80, "test_user.rds"))
test16_N80_user
# $label
# [1] "scn16_N80_user_BWonVc"
# 
# $selected
#          step      covar var shape      objf deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames covarEffect bsvReduction
# eta.vc      1   BW_power  vc power -4419.286 23.36881 -3519.105 -3481.541        10 3.841459 1.337292e-06      yes    forward   cov_BW_power_vc   1.0876951     42.80153
# eta.cl2     2   BW_power  cl power -4450.821 31.53557 -3548.640 -3506.902        11 3.841459 1.958205e-08      yes    forward   cov_BW_power_cl   1.0144435     33.57192
# eta.cl3     3 CrCL_power  cl power -4462.977 12.15594 -3558.796 -3512.884        12 3.841459 4.893170e-04      yes    forward cov_CrCL_power_cl   0.3951763     14.50952
# eta.vc4     4      SEX_1  vc   cat -4474.174 11.19711 -3567.993 -3517.908        13 3.841459 8.192491e-04      yes    forward      cov_SEX_1_vc   0.3127793     15.70659
# eta.vc5     1   BW_power  vc power -4448.458 25.71637 -3544.277 -3498.365        12 6.634897 3.954587e-07 retained   backward   cov_BW_power_vc   0.9921970     36.99777
# eta.cl4     1   BW_power  cl power -4460.474 13.69999 -3556.293 -3510.382        12 6.634897 2.144553e-04 retained   backward   cov_BW_power_cl   0.6862169     14.82973
# eta.cl12    1 CrCL_power  cl power -4462.724 11.45057 -3558.543 -3512.631        12 6.634897 7.147220e-04 retained   backward cov_CrCL_power_cl   0.3858241     13.78382
# eta.vc11    1      SEX_1  vc   cat -4463.198 10.97661 -3559.017 -3513.105        12 6.634897 9.226901e-04 retained   backward      cov_SEX_1_vc   0.3127793     18.59475
# 
# $step_hist
#          step      covar var shape      objf deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames covarEffect bsvReduction
# eta.vc      1   BW_power  vc power -4419.286 23.36881 -3519.105 -3481.541        10 3.841459 1.337292e-06      yes    forward   cov_BW_power_vc   1.0876951    42.801525
# eta.cl      1   BW_power  cl power -4409.699 13.78229 -3509.518 -3471.954        10 3.841459 2.052617e-04       no    forward   cov_BW_power_cl   0.6781998    17.068237
# eta.cl1     1 CrCL_power  cl power -4417.807 21.89056 -3517.626 -3480.062        10 3.841459 2.886499e-06       no    forward cov_CrCL_power_cl   0.5080737    30.035367
# eta.vc1     1      SEX_1  vc   cat -4418.114 22.19693 -3517.933 -3480.369        10 3.841459 2.460675e-06       no    forward      cov_SEX_1_vc   0.4290643    24.752773
# eta.cl2     2   BW_power  cl power -4450.821 31.53557 -3548.640 -3506.902        11 3.841459 1.958205e-08      yes    forward   cov_BW_power_cl   1.0144435    33.571917
# eta.cl11    2 CrCL_power  cl power -4445.905 26.61944 -3543.724 -3501.986        11 3.841459 2.477381e-07       no    forward cov_CrCL_power_cl   0.5975543    28.941427
# eta.vc2     2      SEX_1  vc   cat -4431.047 11.76105 -3528.866 -3487.128        11 3.841459 6.048298e-04       no    forward      cov_SEX_1_vc   0.2941354    16.031311
# eta.cl3     3 CrCL_power  cl power -4462.977 12.15594 -3558.796 -3512.884        12 3.841459 4.893170e-04      yes    forward cov_CrCL_power_cl   0.3951763    14.509516
# eta.vc3     3      SEX_1  vc   cat -4461.750 10.92896 -3557.569 -3511.657        12 3.841459 9.467270e-04       no    forward      cov_SEX_1_vc   0.3121157     8.810766
# eta.vc4     4      SEX_1  vc   cat -4474.174 11.19711 -3567.993 -3517.908        13 3.841459 8.192491e-04      yes    forward      cov_SEX_1_vc   0.3127793    15.706594
# eta.vc5     1   BW_power  vc power -4448.458 25.71637 -3544.277 -3498.365        12 6.634897 3.954587e-07 retained   backward   cov_BW_power_vc   0.9921970    36.997769
# eta.cl4     1   BW_power  cl power -4460.474 13.69999 -3556.293 -3510.382        12 6.634897 2.144553e-04 retained   backward   cov_BW_power_cl   0.6862169    14.829730
# eta.cl12    1 CrCL_power  cl power -4462.724 11.45057 -3558.543 -3512.631        12 6.634897 7.147220e-04 retained   backward cov_CrCL_power_cl   0.3858241    13.783823
# eta.vc11    1      SEX_1  vc   cat -4463.198 10.97661 -3559.017 -3513.105        12 6.634897 9.226901e-04 retained   backward      cov_SEX_1_vc   0.3127793    18.594750
# 
# $final_fit
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#            OBJF       AIC       BIC Log-likelihood Condition#(Cov) Condition#(Cor)
# FOCEi -4474.174 -3567.993 -3517.908       1795.997        415.5014        8.634761
# 
# ── Time (sec $time): ──
# 
#             setup optimize covariance preprocess postprocess    other
# elapsed 0.0173598 39.23995   19.75486       0.09        0.03 15.75783
# 
# ── Population Parameters ($parFixed or $parFixedDf): ──
# 
#                      Est.      SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# lTVCL             -0.3801 0.03675  9.668 0.6838 (0.6363, 0.7349)     33.5      1.88% 
# lTVQ               0.6218  0.0195  3.136    1.862 (1.792, 1.935)                     
# lTVVc               3.142 0.05838  1.858    23.15 (20.65, 25.96)     30.3      13.1% 
# lTVVp               4.397 0.01251 0.2845    81.21 (79.24, 83.22)                     
# lTVKA             -0.3567   FIXED  FIXED                 -0.3567                     
# prop.err           0.1058                                 0.1058                     
# cov_BW_power_vc    0.9922   0.174  17.53  0.9922 (0.6512, 1.333)                     
# cov_BW_power_cl    0.6862  0.1888  27.52  0.6862 (0.3161, 1.056)                     
# cov_CrCL_power_cl  0.3858  0.1108  28.72  0.3858 (0.1687, 0.603)                     
# cov_SEX_1_vc       0.3128 0.08333  26.64 0.3128 (0.1495, 0.4761)                     
#  
#   Covariance Type ($covMethod): r,s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.287   
#  
# 
#   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found ($runInfo):
#    • gradient problems with covariance; see $scaleInfo 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring ($censInformation): No censoring
#   Minimization message ($message):  
#     Normal exit from bobyqa 
# 
# $final_est
# # A tibble: 17 × 2
#    parameter estimate
#    <chr>        <dbl>
#  1 TVCL        0.684 
#  2 TVVc       23.1   
#  3 TVQ         1.86  
#  4 TVVp       81.2   
#  5 TVKA        0.7   
#  6 var_CL      0.107 
#  7 var_Vc      0.0877
#  8 cov_VcCL    0.0277
#  9 ResErr      0.106 
# 10 CLBW        0.686 
# 11 CLcrCL      0.386 
# 12 VcBW        0.992 
# 13 VcSEX       0.313 
# 14 CLBMI      NA     
# 15 VcBMI      NA     
# 16 VcCrCL     NA     
# 17 VcRACE     NA     
# 
# $rel_err
# # A tibble: 17 × 6
#    parameter true_value estimate   abs_err   rel_err rel_err_pct
#    <chr>          <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
#  1 TVCL           0.6     0.684   8.38e- 2  1.40e- 1    1.40e+ 1
#  2 TVVc          20      23.1     3.15e+ 0  1.57e- 1    1.57e+ 1
#  3 TVQ            1.8     1.86    6.23e- 2  3.46e- 2    3.46e+ 0
#  4 TVVp          80      81.2     1.21e+ 0  1.51e- 2    1.51e+ 0
#  5 TVKA           0.7     0.7     3.33e-16  4.76e-16    4.76e-14
#  6 CLBW           0.75    0.686  -6.38e- 2 -8.50e- 2   -8.50e+ 0
#  7 CLcrCL         0.5     0.386  -1.14e- 1 -2.28e- 1   -2.28e+ 1
#  8 VcBW           1       0.992  -7.80e- 3 -7.80e- 3   -7.80e- 1
#  9 VcSEX          0.405   0.313  -9.27e- 2 -2.29e- 1   -2.29e+ 1
# 10 var_CL         0.1     0.107   6.55e- 3  6.55e- 2    6.55e+ 0
# 11 var_Vc         0.1     0.0877 -1.23e- 2 -1.23e- 1   -1.23e+ 1
# 12 cov_VcCL       0.02    0.0277  7.73e- 3  3.86e- 1    3.86e+ 1
# 13 ResErr         0.1     0.106   5.85e- 3  5.85e- 2    5.85e+ 0
# 14 CLBMI         NA      NA      NA        NA          NA       
# 15 VcBMI         NA      NA      NA        NA          NA       
# 16 VcCrCL        NA      NA      NA        NA          NA       
# 17 VcRACE        NA      NA      NA        NA          NA       
# 
# $diag
# $diag$converged
# [1] TRUE
# 
# $diag$objf
# [1] -4474.174
# 
# $diag$cond_num
# [1] 415.5014
# 
# $diag$cond_num_sqrt
# [1] 20.38385
# 
# $diag$cov_ok
# [1] TRUE
# 
# $diag$message
# [1] "Normal exit from bobyqa"
# 
# 
# $parFixed
# NULL
# 
# $final_fit_refit
# NULL
# 
# $runtime_sec
# [1] 576.022


### 4.2 user specified relations with wrong pairs-----------
### Tried tighter focei sigdig and explicit rxcontrol - for wrong pairs, still fails to converge. 
scm_focei_production4 <- nlmixr2est::foceiControl(
  sigdig             = 4,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = FALSE,
  covMethod          = "r,s",
  stickyRecalcN      = 20,                      
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl          = rxode2::rxControl(atol = 1e-8, rtol = 1e-6) #Even with this setting, dofv still can be negative
)
res16_N80_user_wr <- runSCM_traced(
  label       = "scn16_N80_user",
  data        = ds16_01_N80,
  fit         = fit_base16_N80,
  pairsVec    = list(list(var = "vc", covar = "BW",   shapes = "power"),
                     list(var = "cl", covar = "BMI",   shapes = "lin"),
                     list(var = "cl", covar = "CrCL", shapes = "power"),
                     list(var = "vc", covar = "RACE",  shapes = "cat")),
  catvarsVec  = "RACE",     # required: declares SEX so runSCM builds the SEX_<level> dummies
  searchType  = "scm",
  control     = scm_focei_production4,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
# 
# === [scn16_N80_user] 15:47:43 starting scm      | 4 candidate(s), 3 worker(s) ===
# ── SCM Summary ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -4395.917
# ℹ Base model params : 9
# ℹ Search type : scm
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (vc, cl)
# ℹ Covariates : 4
# (BW_power, BMI_lin, CrCL_power, RACE_0)
# ℹ Total candidates : 4
# ── Categorical covariates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ RACE: reference = '1'
# Indicators (1): RACE_0
# ── Relationships to test ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ vc [power]
# 2. BMI_lin ~ cl [lin]
# 3. CrCL_power ~ cl [power]
# 4. RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Proceed with SCM? [y/N]: 
# 
# ── SCM Step Summary ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.917   -4419.286    23.369     0.0000  Added
# Forward    2     CrCL_power~cl [power]   -4419.286   -4445.905    26.619     0.0000  Added
# Forward    3     RACE_0~vc [cat]         -4445.905   -4449.382     3.477     0.0622  Not selected
# Backward   1     CrCL_power~cl [power]   -4445.905   -4419.276    26.629     0.0000  Retained
# ── SCM All Candidates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.917   -4419.286    23.369     0.0000  Added
# Forward    1     CrCL_power~cl [power]   -4395.917   -4417.807    21.891     0.0000  Not selected
# Forward    1     BMI_lin~cl [lin]        -4395.917   -4405.725     9.808     0.0017  Not selected
# Forward    1     RACE_0~vc [cat]         -4395.917   -4397.612     1.695     0.1930  Not selected
# 
# Forward    2     CrCL_power~cl [power]   -4419.286   -4445.905    26.619     0.0000  Added
# Forward    2     BMI_lin~cl [lin]        -4419.286   -4429.930    10.644     0.0011  Not selected
# Forward    2     RACE_0~vc [cat]         -4419.286   -4419.122    -0.163     1.0000  Not selected
# 
# Forward    3     RACE_0~vc [cat]         -4445.905   -4449.382     3.477     0.0622  Not selected
# Forward    3     BMI_lin~cl [lin]        -4445.905   -4443.817    -2.088     1.0000  Not selected
# 
# Backward   1     BW_power~vc [power]     -4445.905   -4418.152    27.753     0.0000  Retained
# Backward   1     CrCL_power~cl [power]   -4445.905   -4419.276    26.629     0.0000  Retained
# ── Final model ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Retained:
# BW_power~vc [power]
# CrCL_power~cl [power]
# ℹ Final model OFV: -4445.905
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [scn16_N80_user] 15:53:11 DONE | elapsed 5.5 min (327 s) ===
# 
# Warning messages:
# 1: ! RACE_0 ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4419.122 > -4419.286). Accepting best available result (attempt 1/1, dOFV = -0.163). 
# 2: ! BMI_lin ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-4443.817 > -4445.905). Accepting best available result (attempt 1/1, dOFV = -2.088). 




## 5 Speed up runscm -------------
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

fit_base16_N80_linCmt <- readRDS(file.path(stage1_dir16_N80, "fit_base_N80_linCmt.rds"))
## 5.1 test16_N80_fwd_auto_fast---------
res16_N80_fwd_auto_fast <- runSCM_traced(
  label       = "scn16_N80_forward_auto_fast",
  data        = ds16_01_N80,
  fit         = fit_base16_N80_linCmt,
  varsVec     = scm16_vars,
  covarsVec   = scm16_covars,
  catvarsVec  = scm16_catvars,
  shapes      = scm16_shapes,
  searchType  = "forward",
  control     = scm_focei_screen,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
# ── SCM Step Summary ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.942   -4419.355    23.413     0.0000  Added
# Forward    2     BW_power~cl [power]     -4419.355   -4451.133    31.778     0.0000  Added
# Forward    3     CrCL_power~cl [power]   -4451.133   -4463.213    12.080     0.0005  Added
# Forward    4     SEX_1~vc [cat]          -4463.213   -4473.985    10.772     0.0010  Added
# Forward    5     SEX_1~cl [cat]          -4473.985   -4477.111     3.126     0.0771  Not selected
# ── SCM All Candidates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.942   -4419.355    23.413     0.0000  Added
# Forward    1     CrCL_power~cl [power]   -4395.942   -4418.167    22.225     0.0000  Not selected
# Forward    1     SEX_1~vc [cat]          -4395.942   -4417.810    21.869     0.0000  Not selected
# Forward    1     BW_power~cl [power]     -4395.942   -4413.574    17.632     0.0000  Not selected
# Forward    1     BMI_power~cl [power]    -4395.942   -4411.689    15.748     0.0001  Not selected
# Forward    1     BW_lin~vc [lin]         -4395.942   -4407.744    11.802     0.0006  Not selected
# Forward    1     BMI_power~vc [power]    -4395.942   -4402.714     6.773     0.0093  Not selected
# Forward    1     CrCL_lin~cl [lin]       -4395.942   -4402.045     6.104     0.0135  Not selected
# Forward    1     BMI_lin~vc [lin]        -4395.942   -4400.265     4.323     0.0376  Not selected
# Forward    1     BMI_lin~cl [lin]        -4395.942   -4399.699     3.757     0.0526  Not selected
# Forward    1     BW_lin~cl [lin]         -4395.942   -4399.413     3.471     0.0625  Not selected
# Forward    1     SEX_1~cl [cat]          -4395.942   -4398.338     2.397     0.1216  Not selected
# Forward    1     RACE_0~vc [cat]         -4395.942   -4397.632     1.690     0.1936  Not selected
# Forward    1     CrCL_power~vc [power]   -4395.942   -4397.422     1.480     0.2237  Not selected
# Forward    1     CrCL_lin~vc [lin]       -4395.942   -4397.022     1.081     0.2985  Not selected
# Forward    1     RACE_0~cl [cat]         -4395.942   -4396.003     0.061     0.8044  Not selected
# 
# Forward    2     BW_power~cl [power]     -4419.355   -4451.133    31.778     0.0000  Added
# Forward    2     CrCL_power~cl [power]   -4419.355   -4449.473    30.118     0.0000  Not selected
# Forward    2     BMI_power~cl [power]    -4419.355   -4445.847    26.491     0.0000  Not selected
# Forward    2     SEX_1~vc [cat]          -4419.355   -4431.041    11.686     0.0006  Not selected
# Forward    2     BMI_power~vc [power]    -4419.355   -4424.067     4.712     0.0300  Not selected
# Forward    2     RACE_0~cl [cat]         -4419.355   -4420.092     0.737     0.3906  Not selected
# Forward    2     CrCL_power~vc [power]   -4419.355   -4419.392     0.036     0.8486  Not selected
# Forward    2     RACE_0~vc [cat]         -4419.355   -4419.365     0.010     0.9189  Not selected
# Forward    2     SEX_1~cl [cat]          -4419.355   -4419.360     0.004     0.9466  Not selected
# Forward    2     BW_lin~cl [lin]         -4419.355   -4411.990    -7.365     1.0000  Not selected
# Forward    2     CrCL_lin~cl [lin]       -4419.355   -4413.772    -5.583     1.0000  Not selected
# Forward    2     CrCL_lin~vc [lin]       -4419.355   -4409.013   -10.342     1.0000  Not selected
# Forward    2     BMI_lin~cl [lin]        -4419.355   -4411.920    -7.435     1.0000  Not selected
# Forward    2     BMI_lin~vc [lin]        -4419.355   -4412.245    -7.110     1.0000  Not selected
# 
# Forward    3     CrCL_power~cl [power]   -4451.133   -4463.213    12.080     0.0005  Added
# Forward    3     SEX_1~vc [cat]          -4451.133   -4462.964    11.831     0.0006  Not selected
# Forward    3     BMI_power~vc [power]    -4451.133   -4455.826     4.693     0.0303  Not selected
# Forward    3     SEX_1~cl [cat]          -4451.133   -4455.119     3.986     0.0459  Not selected
# Forward    3     BMI_power~cl [power]    -4451.133   -4451.929     0.796     0.3724  Not selected
# Forward    3     RACE_0~cl [cat]         -4451.133   -4451.740     0.607     0.4359  Not selected
# Forward    3     CrCL_power~vc [power]   -4451.133   -4451.174     0.041     0.8397  Not selected
# Forward    3     RACE_0~vc [cat]         -4451.133   -4451.155     0.022     0.8826  Not selected
# Forward    3     CrCL_lin~cl [lin]       -4451.133   -4438.660   -12.473     1.0000  Not selected
# Forward    3     CrCL_lin~vc [lin]       -4451.133   -4435.959   -15.174     1.0000  Not selected
# Forward    3     BMI_lin~cl [lin]        -4451.133   -4435.805   -15.328     1.0000  Not selected
# Forward    3     BMI_lin~vc [lin]        -4451.133   -4436.340   -14.793     1.0000  Not selected
# 
# Forward    4     SEX_1~vc [cat]          -4463.213   -4473.985    10.772     0.0010  Added
# Forward    4     BMI_power~vc [power]    -4463.213   -4468.048     4.835     0.0279  Not selected
# Forward    4     SEX_1~cl [cat]          -4463.213   -4466.968     3.754     0.0527  Not selected
# Forward    4     BMI_power~cl [power]    -4463.213   -4464.579     1.366     0.2426  Not selected
# Forward    4     RACE_0~cl [cat]         -4463.213   -4463.407     0.194     0.6598  Not selected
# Forward    4     RACE_0~vc [cat]         -4463.213   -4463.244     0.031     0.8605  Not selected
# Forward    4     CrCL_power~vc [power]   -4463.213   -4463.222     0.008     0.9267  Not selected
# Forward    4     CrCL_lin~vc [lin]       -4463.213   -4443.547   -19.666     1.0000  Not selected
# Forward    4     BMI_lin~cl [lin]        -4463.213   -4444.580   -18.634     1.0000  Not selected
# Forward    4     BMI_lin~vc [lin]        -4463.213   -4445.040   -18.173     1.0000  Not selected
# 
# Forward    5     SEX_1~cl [cat]          -4473.985   -4477.111     3.126     0.0771  Not selected
# Forward    5     BMI_power~cl [power]    -4473.985   -4475.795     1.810     0.1785  Not selected
# Forward    5     RACE_0~cl [cat]         -4473.985   -4475.217     1.232     0.2670  Not selected
# Forward    5     CrCL_power~vc [power]   -4473.985   -4475.164     1.179     0.2775  Not selected
# Forward    5     RACE_0~vc [cat]         -4473.985   -4475.005     1.020     0.3124  Not selected
# Forward    5     BMI_power~vc [power]    -4473.985   -4474.996     1.012     0.3145  Not selected
# Forward    5     CrCL_lin~vc [lin]       -4473.985   -4440.855   -33.130     1.0000  Not selected
# Forward    5     BMI_lin~cl [lin]        -4473.985   -4441.595   -32.390     1.0000  Not selected
# Forward    5     BMI_lin~vc [lin]        -4473.985   -4442.137   -31.848     1.0000  Not selected
# ── Final model ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Retained:
# BW_power~vc [power]
# BW_power~cl [power]
# CrCL_power~cl [power]
# SEX_1~vc [cat]
# ℹ Final model OFV: -4473.985
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [scn16_N80_forward_auto_fast] 19:49:39 DONE | elapsed 19.2 min (1151 s) ===
# 
# There were 15 warnings (use warnings() to see them)
saveRDS(res16_N80_fwd_auto_fast,
        file.path(stage1_dir16_N80, "res_fwd_auto_fast.rds"))
test16_N80_fwd_auto_fast <- package_scm_result(
  "scn16_N80_forward_auto_fast",
  res16_N80_fwd_auto_fast, t16_N80_fwd_auto_fast,
  scenario_id   = 16,
  refit_control = scm_focei_final
)

test16_N80_fwd_auto_fast
# $label
# [1] "scn16_N80_forward_auto_fast"
# 
# $selected
# # A tibble: 4 × 5
#   var   covar shape theta_name        estimate
#   <chr> <chr> <chr> <chr>                <dbl>
# 1 vc    BW    power cov_BW_power_vc      1.00 
# 2 cl    BW    power cov_BW_power_cl      0.684
# 3 cl    CrCL  power cov_CrCL_power_cl    0.389
# 4 vc    SEX   1     cov_SEX_1_vc         0.301
# 
# $step_hist
#          step      covar var shape      objf      deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames   covarEffect
# eta.cl      1   BW_power  cl power -4413.574  17.631967332 -3513.393 -3475.829        10 3.841459 2.680447e-05       no    forward   cov_BW_power_cl  0.9084817332
# eta.cl1     1     BW_lin  cl   lin -4399.413   3.471055472 -3499.232 -3461.668        10 3.841459 6.245145e-02       no    forward     cov_BW_lin_cl  0.0013054830
# eta.vc      1   BW_power  vc power -4419.355  23.413461400 -3519.174 -3481.610        10 3.841459 1.306610e-06      yes    forward   cov_BW_power_vc  1.1120120021
# eta.vc1     1     BW_lin  vc   lin -4407.744  11.802031905 -3507.563 -3469.999        10 3.841459 5.916611e-04       no    forward     cov_BW_lin_vc  0.0078896306
# eta.cl2     1 CrCL_power  cl power -4418.167  22.225292875 -3517.986 -3480.422        10 3.841459 2.424587e-06       no    forward cov_CrCL_power_cl  0.5240732977
# eta.cl3     1   CrCL_lin  cl   lin -4402.045   6.103574150 -3501.864 -3464.300        10 3.841459 1.349088e-02       no    forward   cov_CrCL_lin_cl  0.0009985787
# eta.vc2     1 CrCL_power  vc power -4397.422   1.480429634 -3497.241 -3459.677        10 3.841459 2.237073e-01       no    forward cov_CrCL_power_vc  0.1747862892
# eta.vc3     1   CrCL_lin  vc   lin -4397.022   1.080830179 -3496.841 -3459.277        10 3.841459 2.985119e-01       no    forward   cov_CrCL_lin_vc  0.0016314046
# eta.cl4     1  BMI_power  cl power -4411.689  15.747608463 -3511.508 -3473.944        10 3.841459 7.237938e-05       no    forward  cov_BMI_power_cl  0.8019273708
# eta.cl5     1    BMI_lin  cl   lin -4399.699   3.756972789 -3499.518 -3461.954        10 3.841459 5.258771e-02       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc4     1  BMI_power  vc power -4402.714   6.772613843 -3502.533 -3464.969        10 3.841459 9.256717e-03       no    forward  cov_BMI_power_vc  0.6286474645
# eta.vc5     1    BMI_lin  vc   lin -4400.265   4.323282992 -3500.084 -3462.520        10 3.841459 3.759432e-02       no    forward    cov_BMI_lin_vc  0.0155213049
# eta.cl6     1      SEX_1  cl   cat -4398.338   2.396771020 -3498.157 -3460.593        10 3.841459 1.215860e-01       no    forward      cov_SEX_1_cl -0.1348785294
# eta.vc6     1      SEX_1  vc   cat -4417.810  21.868641154 -3517.629 -3480.065        10 3.841459 2.919664e-06       no    forward      cov_SEX_1_vc  0.4301377843
# eta.cl7     1     RACE_0  cl   cat -4396.003   0.061364394 -3495.822 -3458.258        10 3.841459 8.043523e-01       no    forward     cov_RACE_0_cl -0.0380652094
# eta.vc7     1     RACE_0  vc   cat -4397.632   1.690371210 -3497.451 -3459.887        10 3.841459 1.935520e-01       no    forward     cov_RACE_0_vc -0.2180457935
# eta.cl8     2   BW_power  cl power -4451.133  31.777972138 -3548.952 -3507.214        11 3.841459 1.728419e-08      yes    forward   cov_BW_power_cl  1.0090481488
# eta.cl11    2     BW_lin  cl   lin -4411.990  -7.364878408 -3509.809 -3468.071        11 3.841459 1.000000e+00       no    forward     cov_BW_lin_cl  0.0013054830
# eta.cl21    2 CrCL_power  cl power -4449.473  30.118256639 -3547.292 -3505.555        11 3.841459 4.064865e-08       no    forward cov_CrCL_power_cl  0.6153718887
# eta.cl31    2   CrCL_lin  cl   lin -4413.772  -5.583387416 -3511.591 -3469.853        11 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_cl  0.0009985787
# eta.vc8     2 CrCL_power  vc power -4419.392   0.036422416 -3517.211 -3475.473        11 3.841459 8.486458e-01       no    forward cov_CrCL_power_vc -0.0295864121
# eta.vc11    2   CrCL_lin  vc   lin -4409.013 -10.342497829 -3506.832 -3465.094        11 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc -0.0028816572
# eta.cl41    2  BMI_power  cl power -4445.847  26.491420092 -3543.666 -3501.928        11 3.841459 2.647112e-07       no    forward  cov_BMI_power_cl  1.0443050096
# eta.cl51    2    BMI_lin  cl   lin -4411.920  -7.435376841 -3509.739 -3468.001        11 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc21    2  BMI_power  vc power -4424.067   4.711567226 -3521.886 -3480.148        11 3.841459 2.996033e-02       no    forward  cov_BMI_power_vc -0.7893514874
# eta.vc31    2    BMI_lin  vc   lin -4412.245  -7.110264252 -3510.064 -3468.326        11 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0161404715
# eta.cl61    2      SEX_1  cl   cat -4419.360   0.004479203 -3517.179 -3475.441        11 3.841459 9.466399e-01       no    forward      cov_SEX_1_cl  0.0062176780
# eta.vc41    2      SEX_1  vc   cat -4431.041  11.685864939 -3528.860 -3487.122        11 3.841459 6.297670e-04       no    forward      cov_SEX_1_vc  0.3091304608
# eta.cl71    2     RACE_0  cl   cat -4420.092   0.736938966 -3517.911 -3476.173        11 3.841459 3.906432e-01       no    forward     cov_RACE_0_cl -0.1450503430
# eta.vc51    2     RACE_0  vc   cat -4419.365   0.010371563 -3517.184 -3475.447        11 3.841459 9.188830e-01       no    forward     cov_RACE_0_vc  0.0188040797
# eta.cl9     3 CrCL_power  cl power -4463.213  12.080331332 -3559.032 -3513.121        12 3.841459 5.095654e-04      yes    forward cov_CrCL_power_cl  0.3949712524
# eta.cl12    3   CrCL_lin  cl   lin -4438.660 -12.473497115 -3534.479 -3488.567        12 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_cl  0.0009985787
# eta.vc9     3 CrCL_power  vc power -4451.174   0.040935381 -3546.993 -3501.081        12 3.841459 8.396627e-01       no    forward cov_CrCL_power_vc -0.0262543920
# eta.vc12    3   CrCL_lin  vc   lin -4435.959 -15.174122371 -3531.778 -3485.866        12 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc -0.0016914004
# eta.cl22    3  BMI_power  cl power -4451.929   0.795666291 -3547.748 -3501.836        12 3.841459 3.723922e-01       no    forward  cov_BMI_power_cl  0.3091569738
# eta.cl32    3    BMI_lin  cl   lin -4435.805 -15.328351920 -3531.624 -3485.712        12 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc22    3  BMI_power  vc power -4455.826   4.692735177 -3551.645 -3505.733        12 3.841459 3.029039e-02       no    forward  cov_BMI_power_vc -0.7952841184
# eta.vc32    3    BMI_lin  vc   lin -4436.340 -14.792962526 -3532.159 -3486.247        12 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0073954280
# eta.cl42    3      SEX_1  cl   cat -4455.119   3.986077188 -3550.938 -3505.027        12 3.841459 4.587776e-02       no    forward      cov_SEX_1_cl -0.2314962091
# eta.vc42    3      SEX_1  vc   cat -4462.964  11.831031480 -3558.783 -3512.871        12 3.841459 5.825172e-04       no    forward      cov_SEX_1_vc  0.3073305474
# eta.cl52    3     RACE_0  cl   cat -4451.740   0.607007795 -3547.559 -3501.647        12 3.841459 4.359167e-01       no    forward     cov_RACE_0_cl  0.1093257544
# eta.vc52    3     RACE_0  vc   cat -4451.155   0.021794720 -3546.974 -3501.062        12 3.841459 8.826345e-01       no    forward     cov_RACE_0_vc  0.0217399098
# eta.vc10    4 CrCL_power  vc power -4463.222   0.008459491 -3557.041 -3506.955        13 3.841459 9.267175e-01       no    forward cov_CrCL_power_vc  0.0141391079
# eta.vc13    4   CrCL_lin  vc   lin -4443.547 -19.665952598 -3537.366 -3487.281        13 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc -0.0004684074
# eta.cl10    4  BMI_power  cl power -4464.579   1.365697884 -3558.398 -3508.313        13 3.841459 2.425524e-01       no    forward  cov_BMI_power_cl  0.3892663061
# eta.cl13    4    BMI_lin  cl   lin -4444.580 -18.633643027 -3538.399 -3488.313        13 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc23    4  BMI_power  vc power -4468.048   4.835006793 -3561.867 -3511.782        13 3.841459 2.788753e-02       no    forward  cov_BMI_power_vc -0.7887056815
# eta.vc33    4    BMI_lin  vc   lin -4445.040 -18.173391160 -3538.859 -3488.774        13 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0078427518
# eta.cl23    4      SEX_1  cl   cat -4466.968   3.754156805 -3560.787 -3510.701        13 3.841459 5.267636e-02       no    forward      cov_SEX_1_cl -0.2020882634
# eta.vc43    4      SEX_1  vc   cat -4473.985  10.771500113 -3567.804 -3517.718        13 3.841459 1.030749e-03      yes    forward      cov_SEX_1_vc  0.3081906737
# eta.cl33    4     RACE_0  cl   cat -4463.407   0.193720354 -3557.226 -3507.141        13 3.841459 6.598381e-01       no    forward     cov_RACE_0_cl  0.0597881204
# eta.vc53    4     RACE_0  vc   cat -4463.244   0.030867276 -3557.063 -3506.978        13 3.841459 8.605368e-01       no    forward     cov_RACE_0_vc  0.0263377498
# eta.vc14    5 CrCL_power  vc power -4475.164   1.179234401 -3566.983 -3512.724        14 3.841459 2.775119e-01       no    forward cov_CrCL_power_vc  0.0524172300
# eta.vc15    5   CrCL_lin  vc   lin -4440.855 -33.129717629 -3532.674 -3478.415        14 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc  0.0001307784
# eta.cl14    5  BMI_power  cl power -4475.795   1.809707362 -3567.614 -3513.354        14 3.841459 1.785433e-01       no    forward  cov_BMI_power_cl  0.2937043767
# eta.cl15    5    BMI_lin  cl   lin -4441.595 -32.389587586 -3533.414 -3479.155        14 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc24    5  BMI_power  vc power -4474.996   1.011506230 -3566.815 -3512.556        14 3.841459 3.145423e-01       no    forward  cov_BMI_power_vc  0.0298406592
# eta.vc34    5    BMI_lin  vc   lin -4442.137 -31.848255087 -3533.956 -3479.696        14 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0067708754
# eta.cl24    5      SEX_1  cl   cat -4477.111   3.125840013 -3568.930 -3514.671        14 3.841459 7.706015e-02       no    forward      cov_SEX_1_cl -0.1325983626
# eta.cl34    5     RACE_0  cl   cat -4475.217   1.232027424 -3567.036 -3512.777        14 3.841459 2.670131e-01       no    forward     cov_RACE_0_cl  0.0635971984
# eta.vc44    5     RACE_0  vc   cat -4475.005   1.020446988 -3566.824 -3512.565        14 3.841459 3.124130e-01       no    forward     cov_RACE_0_vc -0.0167253220
#           bsvReduction
# eta.cl    3.301946e+01
# eta.cl1  -1.997907e-13
# eta.vc    4.331752e+01
# eta.vc1   5.250347e-09
# eta.cl2   3.104095e+01
# eta.cl3  -1.997907e-13
# eta.vc2   9.174201e+00
# eta.vc3   2.330973e-03
# eta.cl4   2.735838e+01
# eta.cl5   4.155184e-02
# eta.vc4   2.018488e+01
# eta.vc5   2.420442e-07
# eta.cl6  -1.299902e+00
# eta.vc6   2.522587e+01
# eta.cl7   4.570263e-01
# eta.vc7   4.171365e+00
# eta.cl8   3.334699e+01
# eta.cl11  1.532540e-01
# eta.cl21  3.142949e+01
# eta.cl31  1.521120e-01
# eta.vc8  -2.659845e-01
# eta.vc11 -7.642134e+01
# eta.cl41  2.543172e+01
# eta.cl51  1.612334e-01
# eta.vc21  6.528204e+00
# eta.vc31 -7.628704e+01
# eta.cl61  7.775346e-02
# eta.vc41  1.321021e+01
# eta.cl71  7.270450e-01
# eta.vc51 -2.371172e-01
# eta.cl9   1.405358e+01
# eta.cl12 -4.980151e+01
# eta.vc9   1.971837e-01
# eta.vc12 -7.647146e+01
# eta.cl22  2.262342e-01
# eta.cl32 -4.978882e+01
# eta.vc22  7.303346e+00
# eta.vc32 -7.639440e+01
# eta.cl42  2.186661e-02
# eta.vc42  1.496117e+01
# eta.cl52  7.280446e-01
# eta.vc52  7.685914e-02
# eta.vc10 -1.135475e-01
# eta.vc13 -7.639553e+01
# eta.cl10  1.181416e+00
# eta.cl13 -7.428164e+01
# eta.vc23  7.499391e+00
# eta.vc33 -7.631242e+01
# eta.cl23  1.522263e-01
# eta.vc43  1.926044e+01
# eta.cl33  7.715297e-02
# eta.vc53 -3.592586e-02
# eta.vc14 -5.130787e+00
# eta.vc15 -1.184738e+02
# eta.cl14  3.793556e+00
# eta.cl15 -6.869606e+01
# eta.vc24 -5.566974e+00
# eta.vc34 -1.183825e+02
# eta.cl24  4.367539e+00
# eta.cl34  3.523347e+00
# eta.vc44 -5.125540e+00
# 
# $final_fit
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#        OBJF       AIC       BIC Log-likelihood Condition#(Cov) Condition#(Cor)
# FOCEi -4475 -3568.819 -3518.734        1796.41        412.9822        9.408444
# 
# ── Time (sec $time): ──
# 
#             setup optimize covariance preprocess postprocess table     other
# elapsed 0.0261655 48.34642   11.26855       0.18        0.01  0.08 0.1988633
# 
# ── Population Parameters ($parFixed or $parFixedDf): ──
# 
#                                                   Parameter    Est.      SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# lTVCL                                                       -0.3802 0.03607  9.488 0.6837 (0.6371, 0.7338)     33.0     0.608% 
# lTVQ                                                         0.6192  0.0187   3.02    1.857 (1.791, 1.927)                     
# lTVVc                                                         3.148 0.05791   1.84    23.28 (20.78, 26.08)     30.7      12.9% 
# lTVVp                                                         4.396 0.01193 0.2713    81.16 (79.28, 83.08)                     
# lTVKA             KA unidentifiable from this sparse design -0.3567   FIXED  FIXED                 -0.3567                     
# prop.err                                                     0.1054                                 0.1054                     
# cov_BW_power_vc                                               1.002  0.1723   17.2     1.002 (0.664, 1.34)                     
# cov_BW_power_cl                                              0.6842  0.1809  26.44  0.6842 (0.3296, 1.039)                     
# cov_CrCL_power_cl                                            0.3892  0.1072  27.55  0.3892 (0.179, 0.5994)                     
# cov_SEX_1_vc                                                 0.3005 0.08561  28.49 0.3005 (0.1327, 0.4683)                     
#  
#   Covariance Type ($covMethod): r,s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.172   
#  
# 
#   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found ($runInfo):
#    • gradient problems with covariance; see $scaleInfo 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring ($censInformation): No censoring
#   Minimization message ($message):  
#     Normal exit from bobyqa 
# 
# ── Fit Data (object is a modified tibble): ──
# # A tibble: 480 × 32
#   ID     TIME    DV  PRED      RES     WRES IPRED  IRES IWRES CPRED     CRES    CWRES eta.cl eta.vc     depot central peripheral1 BW_power_cl CrCL_power_cl cov_cl    cl
#   <fct> <dbl> <dbl> <dbl>    <dbl>    <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl>  <dbl>     <dbl>   <dbl>       <dbl>       <dbl>         <dbl>  <dbl> <dbl>
# 1 1      0     0     0    0        0         0    0     0      0     0        0        0.156  0.144 100           0           0        -0.199      -0.00868 -0.207 0.649
# 2 1      7.06  2.77  2.77 0.000163 0.000377  2.57 0.202 0.748  2.78 -0.00426 -0.00979  0.156  0.144   0.712      51.7        34.6      -0.199      -0.00868 -0.207 0.649
# 3 1     14.1   1.58  1.36 0.224    0.997     1.33 0.255 1.82   1.37  0.216    1.01     0.156  0.144   0.00507    26.7        51.7      -0.199      -0.00868 -0.207 0.649
# # ℹ 477 more rows
# # ℹ 11 more variables: BW_power_vc <dbl>, SEX_1_vc <dbl>, cov_vc <dbl>, vc <dbl>, q <dbl>, vp <dbl>, tad <dbl>, dosenum <dbl>, SEX_1 <dbl>, CrCL <dbl>, BW <dbl>
# # ℹ Use `print(n = ...)` to see more rows
# 
# $final_est
# # A tibble: 17 × 2
#    parameter estimate
#    <chr>        <dbl>
#  1 TVCL        0.684 
#  2 TVVc       23.3   
#  3 TVQ         1.86  
#  4 TVVp       81.2   
#  5 TVKA        0.7   
#  6 var_CL      0.103 
#  7 var_Vc      0.0903
#  8 cov_VcCL    0.0166
#  9 ResErr      0.105 
# 10 CLBW        0.684 
# 11 CLcrCL      0.389 
# 12 VcBW        1.00  
# 13 VcSEX       0.301 
# 14 CLBMI      NA     
# 15 VcBMI      NA     
# 16 VcCrCL     NA     
# 17 VcRACE     NA     
# 
# $rel_err
# # A tibble: 17 × 6
#    parameter true_value estimate   abs_err   rel_err rel_err_pct
#    <chr>          <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
#  1 TVCL           0.6     0.684   8.37e- 2  1.40e- 1    1.40e+ 1
#  2 TVVc          20      23.3     3.28e+ 0  1.64e- 1    1.64e+ 1
#  3 TVQ            1.8     1.86    5.75e- 2  3.19e- 2    3.19e+ 0
#  4 TVVp          80      81.2     1.16e+ 0  1.45e- 2    1.45e+ 0
#  5 TVKA           0.7     0.7     3.33e-16  4.76e-16    4.76e-14
#  6 CLBW           0.75    0.684  -6.58e- 2 -8.77e- 2   -8.77e+ 0
#  7 CLcrCL         0.5     0.389  -1.11e- 1 -2.22e- 1   -2.22e+ 1
#  8 VcBW           1       1.00    1.77e- 3  1.77e- 3    1.77e- 1
#  9 VcSEX          0.405   0.301  -1.05e- 1 -2.59e- 1   -2.59e+ 1
# 10 var_CL         0.1     0.103   3.32e- 3  3.32e- 2    3.32e+ 0
# 11 var_Vc         0.1     0.0903 -9.71e- 3 -9.71e- 2   -9.71e+ 0
# 12 cov_VcCL       0.02    0.0166 -3.39e- 3 -1.69e- 1   -1.69e+ 1
# 13 ResErr         0.1     0.105   5.40e- 3  5.40e- 2    5.40e+ 0
# 14 CLBMI         NA      NA      NA        NA          NA       
# 15 VcBMI         NA      NA      NA        NA          NA       
# 16 VcCrCL        NA      NA      NA        NA          NA       
# 17 VcRACE        NA      NA      NA        NA          NA       
# 
# $diag
# $diag$converged
# [1] TRUE
# 
# $diag$objf
# [1] -4475
# 
# $diag$cond_num
# [1] 412.9822
# 
# $diag$cond_num_sqrt
# [1] 20.32196
# 
# $diag$cov_ok
# [1] TRUE
# 
# $diag$message
# [1] "Normal exit from bobyqa"
# 
# 
# $parFixed
#                                                   Parameter   Estimate         SE       %RSE Back-transformed   CI Lower   CI Upper BSV(CV%) Shrink(SD)%
# lTVCL                                                       -0.3801686 0.03607010  9.4879230        0.6837461  0.6370770  0.7338340 32.99140   0.6079082
# lTVQ                                                         0.6192196 0.01870190  3.0202367        1.8574779  1.7906247  1.9268271       NA          NA
# lTVVc                                                        3.1476690 0.05790706  1.8396808       23.2817310 20.7837873 26.0798955 30.73936  12.9195978
# lTVVp                                                        4.3964033 0.01192604  0.2712681       81.1584413 79.2833939 83.0778334       NA          NA
# lTVKA             KA unidentifiable from this sparse design -0.3566749         NA         NA       -0.3566749         NA         NA       NA          NA
# prop.err                                                     0.1053955         NA         NA        0.1053955         NA         NA       NA          NA
# cov_BW_power_vc                                              1.0017715 0.17233943 17.2034673        1.0017715  0.6639924  1.3395505       NA          NA
# cov_BW_power_cl                                              0.6841879 0.18090041 26.4401660        0.6841879  0.3296296  1.0387461       NA          NA
# cov_CrCL_power_cl                                            0.3892281 0.10724026 27.5520302        0.3892281  0.1790411  0.5994152       NA          NA
# cov_SEX_1_vc                                                 0.3005307 0.08561311 28.4873142        0.3005307  0.1327320  0.4683293       NA          NA
# 
# $refit_done
# [1] TRUE
# 
# $runtime_sec
# [1] 1150.83
saveRDS(test16_N80_fwd_auto_fast,
        file.path(stage1_dir16_N80, "test_fwd_auto_fast.rds"))


#5.2 res16_N80_bck_fast--------------
res16_N80_bck_fast <- runSCM_traced(
  label             = "scn16_N80_backward_fast",
  data              = ds16_01_N80,
  fit               = fit_base16_N80_linCmt,
  pairsVec          = candidate_pairs_test_corContinusous_cat,
  catvarsVec        = scm16_catvars,
  searchType        = "backward",
  includedRelations = candidate_pairs_test_corContinusous_cat,
  control           = scm_focei_screen,
  saveModels        = FALSE,
  workers           = 3L,
  print             = 100,
  maxRetries        = 0L
)
# 
# === [scn16_N80_backward_fast] 20:22:28 starting backward | 7 candidate(s), 3 worker(s) ===
# ── SCM Summary ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -4395.942
# ℹ Base model params : 9
# ℹ Search type : backward
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (cl, vc)
# ℹ Covariates : 5
# (BW_power, CrCL_power, BMI_power, SEX_1, RACE_0)
# ℹ Total candidates : 7
# ── Categorical covariates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ SEX: reference = '0'
# Indicators (1): SEX_1
# ℹ RACE: reference = '1'
# Indicators (1): RACE_0
# ── Relationships to test ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ cl [power]
# 2. CrCL_power ~ cl [power]
# 3. BW_power ~ vc [power]
# 4. BMI_power ~ vc [power]
# 5. CrCL_power ~ vc [power]
# 6. SEX_1 ~ vc [cat]
# 7. RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Proceed with SCM? [y/N]: 

# ── SCM Step Summary ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Backward   1     BMI_power~vc [power]    -4475.199   -4475.178     0.021     0.8861  Removed
# Backward   2     RACE_0~vc [cat]         -4475.178   -4475.059     0.120     0.7292  Removed
# Backward   3     CrCL_power~vc [power]   -4475.059   -4474.998     0.060     0.8062  Removed
# Backward   4     SEX_1~vc [cat]          -4474.998   -4463.216    11.783     0.0006  Retained
# ── SCM All Candidates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Backward   1     BW_power~cl [power]     -4475.199   -4461.415    13.784     0.0002  Retained
# Backward   1     CrCL_power~cl [power]   -4475.199   -4462.317    12.882     0.0003  Retained
# Backward   1     SEX_1~vc [cat]          -4475.199   -4468.003     7.196     0.0073  Retained
# Backward   1     BW_power~vc [power]     -4475.199   -4472.112     3.087     0.0789  Retained
# Backward   1     CrCL_power~vc [power]   -4475.199   -4474.621     0.578     0.4470  Retained
# Backward   1     RACE_0~vc [cat]         -4475.199   -4475.065     0.134     0.7146  Retained
# Backward   1     BMI_power~vc [power]    -4475.199   -4475.178     0.021     0.8861  Removed
# 
# Backward   2     BW_power~vc [power]     -4475.178   -4447.101    28.077     0.0000  Retained
# Backward   2     BW_power~cl [power]     -4475.178   -4461.043    14.135     0.0002  Retained
# Backward   2     CrCL_power~cl [power]   -4475.178   -4462.378    12.800     0.0003  Retained
# Backward   2     SEX_1~vc [cat]          -4475.178   -4463.229    11.950     0.0005  Retained
# Backward   2     CrCL_power~vc [power]   -4475.178   -4473.628     1.551     0.2130  Retained
# Backward   2     RACE_0~vc [cat]         -4475.178   -4475.059     0.120     0.7292  Removed
# 
# Backward   3     BW_power~vc [power]     -4475.059   -4449.735    25.323     0.0000  Retained
# Backward   3     BW_power~cl [power]     -4475.059   -4461.397    13.661     0.0002  Retained
# Backward   3     CrCL_power~cl [power]   -4475.059   -4462.829    12.230     0.0005  Retained
# Backward   3     SEX_1~vc [cat]          -4475.059   -4463.178    11.881     0.0006  Retained
# Backward   3     CrCL_power~vc [power]   -4475.059   -4474.998     0.060     0.8062  Removed
# 
# Backward   4     BW_power~vc [power]     -4474.998   -4448.692    26.306     0.0000  Retained
# Backward   4     BW_power~cl [power]     -4474.998   -4460.999    14.000     0.0002  Retained
# Backward   4     CrCL_power~cl [power]   -4474.998   -4463.025    11.973     0.0005  Retained
# Backward   4     SEX_1~vc [cat]          -4474.998   -4463.216    11.783     0.0006  Retained
# ── Final model ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Removed:
# BMI_power~vc [power]
# RACE_0~vc [cat]
# CrCL_power~vc [power]
# ℹ Final model OFV: -4474.998
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [scn16_N80_backward_fast] 20:35:46 DONE | elapsed 13.3 min (798 s) ===
# 
t16_N80_bck_fast <- attr(res16_N80_bck_fast, "elapsed_s")
saveRDS(res16_N80_bck_fast,
        file.path(stage1_dir16_N80, "res_bck_fast.rds"))
test16_N80_bck_fast <- package_scm_result(
  "scn16_N80_backward_only_fast",
  res16_N80_bck_fast, t16_N80_bck_fast,
  scenario_id   = 16,
  refit_control = scm_focei_final
)
# calculating covariance matrix
# [====|====|====|====|====|====|====|====|====|====] 100%; 0:00:11 done
# → Calculating residuals/tables
# ✔ done
saveRDS(test16_N80_bck_fast,
        file.path(stage1_dir16_N80, "test_bck_fast.rds"))
test16_N80_bck_fast
# $label
# [1] "scn16_N80_backward_only_fast"
# 
# $selected
# # A tibble: 4 × 5
#   var   covar shape theta_name        estimate
#   <chr> <chr> <chr> <chr>                <dbl>
# 1 cl    BW    power cov_BW_power_cl      0.684
# 2 cl    CrCL  power cov_CrCL_power_cl    0.389
# 3 vc    BW    power cov_BW_power_vc      1.00 
# 4 vc    SEX   1     cov_SEX_1_vc         0.300
# 
# $step_hist
#          step      covar var shape      objf    deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames covarEffect bsvReduction
# eta.cl      1   BW_power  cl power -4461.415 13.78351197 -3551.234 -3492.801        15 6.634897 2.051286e-04 retained   backward   cov_BW_power_cl  0.67677124   16.8425962
# eta.cl1     1 CrCL_power  cl power -4462.317 12.88235168 -3552.136 -3493.703        15 6.634897 3.316952e-04 retained   backward cov_CrCL_power_cl  0.39690016   16.3249635
# eta.vc      1   BW_power  vc power -4472.112  3.08660833 -3561.931 -3503.498        15 6.634897 7.893919e-02 retained   backward   cov_BW_power_vc  0.87728060    5.7627668
# eta.vc1     1  BMI_power  vc power -4475.178  0.02052460 -3564.997 -3506.564        15 6.634897 8.860816e-01  dropped   backward  cov_BMI_power_vc  0.06763612    1.5135321
# eta.vc2     1 CrCL_power  vc power -4474.621  0.57836114 -3564.440 -3506.007        15 6.634897 4.469554e-01 retained   backward cov_CrCL_power_vc  0.05700645   -0.7143426
# eta.vc3     1      SEX_1  vc   cat -4468.003  7.19583319 -3557.822 -3499.389        15 6.634897 7.307305e-03 retained   backward      cov_SEX_1_vc  0.31720340    8.4725594
# eta.vc4     1     RACE_0  vc   cat -4475.065  0.13368230 -3564.884 -3506.451        15 6.634897 7.146442e-01 retained   backward     cov_RACE_0_vc -0.02661441    0.4824961
# eta.cl2     2   BW_power  cl power -4461.043 14.13527301 -3552.862 -3498.603        14 6.634897 1.701237e-04 retained   backward   cov_BW_power_cl  0.67713178   17.8340257
# eta.cl11    2 CrCL_power  cl power -4462.378 12.80008164 -3554.197 -3499.938        14 6.634897 3.466042e-04 retained   backward cov_CrCL_power_cl  0.39668698   15.3940451
# eta.vc5     2   BW_power  vc power -4447.101 28.07732840 -3538.920 -3484.661        14 6.634897 1.165634e-07 retained   backward   cov_BW_power_vc  0.94296285   47.8336332
# eta.vc11    2 CrCL_power  vc power -4473.628  1.55094434 -3565.447 -3511.187        14 6.634897 2.129961e-01 retained   backward cov_CrCL_power_vc  0.05424424   -4.8942599
# eta.vc21    2      SEX_1  vc   cat -4463.229 11.94991823 -3555.048 -3500.788        14 6.634897 5.464978e-04 retained   backward      cov_SEX_1_vc  0.30735210   14.8390502
# eta.vc31    2     RACE_0  vc   cat -4475.059  0.11984568 -3566.878 -3512.618        14 6.634897 7.292019e-01  dropped   backward     cov_RACE_0_vc -0.02577667   -1.0786679
# eta.cl3     3   BW_power  cl power -4461.397 13.66121049 -3555.216 -3505.131        13 6.634897 2.189303e-04 retained   backward   cov_BW_power_cl  0.67395624   15.6384041
# eta.cl12    3 CrCL_power  cl power -4462.829 12.22977136 -3556.648 -3506.562        13 6.634897 4.703296e-04 retained   backward cov_CrCL_power_cl  0.39601461   13.8570038
# eta.vc6     3   BW_power  vc power -4449.735 25.32320233 -3543.554 -3493.469        13 6.634897 4.848414e-07 retained   backward   cov_BW_power_vc  0.97463710   41.8040402
# eta.vc12    3 CrCL_power  vc power -4474.998  0.06018588 -3568.817 -3518.732        13 6.634897 8.062024e-01  dropped   backward cov_CrCL_power_vc  0.04554222    0.4668776
# eta.vc22    3      SEX_1  vc   cat -4463.178 11.88094231 -3556.997 -3506.911        13 6.634897 5.671124e-04 retained   backward      cov_SEX_1_vc  0.30546642   14.5128782
# eta.cl4     4   BW_power  cl power -4460.999 13.99956111 -3556.818 -3510.906        12 6.634897 1.828533e-04 retained   backward   cov_BW_power_cl  0.68507441   16.0308254
# eta.cl13    4 CrCL_power  cl power -4463.025 11.97297591 -3558.844 -3512.933        12 6.634897 5.397767e-04 retained   backward cov_CrCL_power_cl  0.38873508   14.3425839
# eta.vc7     4   BW_power  vc power -4448.692 26.30625061 -3544.511 -3498.600        12 6.634897 2.913442e-07 retained   backward   cov_BW_power_vc  1.00253154   34.4302317
# eta.vc13    4      SEX_1  vc   cat -4463.216 11.78259690 -3559.035 -3513.123        12 6.634897 5.978702e-04 retained   backward      cov_SEX_1_vc  0.30090536   14.9131803
# 
# $final_fit
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#        OBJF       AIC       BIC Log-likelihood Condition#(Cov) Condition#(Cor)
# FOCEi -4475 -3568.819 -3518.734        1796.41        412.9404         9.40415
# 
# ── Time (sec $time): ──
# 
#             setup optimize covariance preprocess postprocess table     other
# elapsed 0.0368028 13.72226   11.24431        0.2        0.03  0.05 0.2566217
# 
# ── Population Parameters ($parFixed or $parFixedDf): ──
# 
#                                                   Parameter    Est.      SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# lTVCL                                                       -0.3802 0.03607  9.488 0.6837 (0.6371, 0.7338)     33.0     0.622% 
# lTVQ                                                         0.6192  0.0187   3.02    1.857 (1.791, 1.927)                     
# lTVVc                                                         3.148  0.0579  1.839    23.28 (20.78, 26.08)     30.7      12.9% 
# lTVVp                                                         4.396 0.01193 0.2713    81.16 (79.28, 83.08)                     
# lTVKA             KA unidentifiable from this sparse design -0.3567   FIXED  FIXED                 -0.3567                     
# prop.err                                                     0.1054                                 0.1054                     
# cov_BW_power_cl                                              0.6841  0.1809  26.44  0.6841 (0.3296, 1.039)                     
# cov_CrCL_power_cl                                            0.3893  0.1072  27.55 0.3893 (0.1791, 0.5995)                     
# cov_BW_power_vc                                               1.002  0.1723  17.19    1.002 (0.6646, 1.34)                     
# cov_SEX_1_vc                                                 0.3005  0.0856  28.49 0.3005 (0.1327, 0.4682)                     
#  
#   Covariance Type ($covMethod): r,s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.172   
#  
# 
#   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found ($runInfo):
#    • gradient problems with covariance; see $scaleInfo 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring ($censInformation): No censoring
#   Minimization message ($message):  
#     Normal exit from bobyqa 
# 
# ── Fit Data (object is a modified tibble): ──
# # A tibble: 480 × 32
#   ID     TIME    DV  PRED      RES    WRES IPRED  IRES IWRES CPRED     CRES   CWRES eta.cl eta.vc   depot central peripheral1 BW_power_cl CrCL_power_cl cov_cl    cl BW_power_vc
#   <fct> <dbl> <dbl> <dbl>    <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>   <dbl>  <dbl>  <dbl>   <dbl>   <dbl>       <dbl>       <dbl>         <dbl>  <dbl> <dbl>       <dbl>
# 1 1      0     0     0     0       0        0    0     0      0     0        0       0.156  0.145 1   e+2     0           0        -0.199      -0.00868 -0.207 0.649      -0.291
# 2 1      7.06  2.77  2.77  4.39e-5 1.01e-4  2.57 0.202 0.747  2.78 -0.00439 -0.0101  0.156  0.145 7.12e-1    51.7        34.6      -0.199      -0.00868 -0.207 0.649      -0.291
# 3 1     14.1   1.58  1.36  2.25e-1 9.97e-1  1.33 0.255 1.82   1.37  0.216    1.01    0.156  0.145 5.07e-3    26.7        51.7      -0.199      -0.00868 -0.207 0.649      -0.291
# # ℹ 477 more rows
# # ℹ 10 more variables: SEX_1_vc <dbl>, cov_vc <dbl>, vc <dbl>, q <dbl>, vp <dbl>, tad <dbl>, dosenum <dbl>, SEX_1 <dbl>, CrCL <dbl>, BW <dbl>
# # ℹ Use `print(n = ...)` to see more rows
# 
# $final_est
# # A tibble: 17 × 2
#    parameter estimate
#    <chr>        <dbl>
#  1 TVCL        0.684 
#  2 TVVc       23.3   
#  3 TVQ         1.86  
#  4 TVVp       81.2   
#  5 TVKA        0.7   
#  6 var_CL      0.103 
#  7 var_Vc      0.0903
#  8 cov_VcCL    0.0166
#  9 ResErr      0.105 
# 10 CLBW        0.684 
# 11 CLcrCL      0.389 
# 12 VcBW        1.00  
# 13 VcSEX       0.300 
# 14 CLBMI      NA     
# 15 VcBMI      NA     
# 16 VcCrCL     NA     
# 17 VcRACE     NA     
# 
# $rel_err
# # A tibble: 17 × 6
#    parameter true_value estimate   abs_err   rel_err rel_err_pct
#    <chr>          <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
#  1 TVCL           0.6     0.684   8.37e- 2  1.40e- 1    1.40e+ 1
#  2 TVVc          20      23.3     3.28e+ 0  1.64e- 1    1.64e+ 1
#  3 TVQ            1.8     1.86    5.75e- 2  3.19e- 2    3.19e+ 0
#  4 TVVp          80      81.2     1.16e+ 0  1.45e- 2    1.45e+ 0
#  5 TVKA           0.7     0.7     3.33e-16  4.76e-16    4.76e-14
#  6 CLBW           0.75    0.684  -6.59e- 2 -8.79e- 2   -8.79e+ 0
#  7 CLcrCL         0.5     0.389  -1.11e- 1 -2.21e- 1   -2.21e+ 1
#  8 VcBW           1       1.00    2.30e- 3  2.30e- 3    2.30e- 1
#  9 VcSEX          0.405   0.300  -1.05e- 1 -2.59e- 1   -2.59e+ 1
# 10 var_CL         0.1     0.103   3.35e- 3  3.35e- 2    3.35e+ 0
# 11 var_Vc         0.1     0.0903 -9.72e- 3 -9.72e- 2   -9.72e+ 0
# 12 cov_VcCL       0.02    0.0166 -3.39e- 3 -1.70e- 1   -1.70e+ 1
# 13 ResErr         0.1     0.105   5.39e- 3  5.39e- 2    5.39e+ 0
# 14 CLBMI         NA      NA      NA        NA          NA       
# 15 VcBMI         NA      NA      NA        NA          NA       
# 16 VcCrCL        NA      NA      NA        NA          NA       
# 17 VcRACE        NA      NA      NA        NA          NA       
# 
# $diag
# $diag$converged
# [1] TRUE
# 
# $diag$objf
# [1] -4475
# 
# $diag$cond_num
# [1] 412.9404
# 
# $diag$cond_num_sqrt
# [1] 20.32094
# 
# $diag$cov_ok
# [1] TRUE
# 
# $diag$message
# [1] "Normal exit from bobyqa"
# 
# 
# $parFixed
#                                                   Parameter   Estimate         SE      %RSE Back-transformed   CI Lower   CI Upper BSV(CV%) Shrink(SD)%
# lTVCL                                                       -0.3801644 0.03606887  9.487703        0.6837490  0.6370812  0.7338353  32.9967   0.6223343
# lTVQ                                                         0.6192247 0.01870076  3.020028        1.8574873  1.7906378  1.9268325       NA          NA
# lTVVc                                                        3.1476581 0.05789671  1.839358       23.2814773 20.7839823 26.0790822  30.7371  12.9167750
# lTVVp                                                        4.3964081 0.01192543  0.271254       81.1588309 79.2838688 83.0781334       NA          NA
# lTVKA             KA unidentifiable from this sparse design -0.3566749         NA        NA       -0.3566749         NA         NA       NA          NA
# prop.err                                                     0.1053940         NA        NA        0.1053940         NA         NA       NA          NA
# cov_BW_power_cl                                              0.6841017 0.18088949 26.441902        0.6841017  0.3295648  1.0386386       NA          NA
# cov_CrCL_power_cl                                            0.3893021 0.10723628 27.545774        0.3893021  0.1791229  0.5994814       NA          NA
# cov_BW_power_vc                                              1.0022993 0.17231930 17.192400        1.0022993  0.6645596  1.3400389       NA          NA
# cov_SEX_1_vc                                                 0.3004529 0.08560353 28.491497        0.3004529  0.1326731  0.4682328       NA          NA
# 
# $refit_done
# [1] TRUE
# 
# $runtime_sec
# [1] 797.9679

#5.3 res16_N80_full_fast----------------
res16_N80_full_fast <- runSCM_traced(
  label       = "scn16_N80_full_fast",
  data        = ds16_01_N80,
  fit         = fit_base16_N80_linCmt,
  varsVec     = scm16_vars,
  covarsVec   = scm16_covars,
  catvarsVec  = scm16_catvars,
  shapes      = scm16_shapes,
  searchType  = "scm",
  control     = scm_focei_screen,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
# 
# === [scn16_N80_full_fast] 20:37:12 starting scm      | NA candidate(s), 3 worker(s) ===
# ── SCM Summary ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -4395.942
# ℹ Base model params : 9
# ℹ Search type : scm
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (cl, vc)
# ℹ Covariates : 8
# (BW_power, BW_lin, CrCL_power, CrCL_lin, BMI_power, BMI_lin, SEX_1, RACE_0)
# ℹ Total candidates : 16
# ── Categorical covariates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ SEX: reference = '0'
# Indicators (1): SEX_1
# ℹ RACE: reference = '1'
# Indicators (1): RACE_0
# ── Relationships to test ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ cl [power]
# 2. BW_lin ~ cl [lin]
# 3. BW_power ~ vc [power]
# 4. BW_lin ~ vc [lin]
# 5. CrCL_power ~ cl [power]
# 6. CrCL_lin ~ cl [lin]
# 7. CrCL_power ~ vc [power]
# 8. CrCL_lin ~ vc [lin]
# 9. BMI_power ~ cl [power]
# 10. BMI_lin ~ cl [lin]
# 11. BMI_power ~ vc [power]
# 12. BMI_lin ~ vc [lin]
# 13. SEX_1 ~ cl [cat]
# 14. SEX_1 ~ vc [cat]
# 15. RACE_0 ~ cl [cat]
# 16. RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Proceed with SCM? [y/N]: 


# ── SCM Step Summary ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.942   -4419.355    23.413     0.0000  Added
# Forward    2     BW_power~cl [power]     -4419.355   -4451.133    31.778     0.0000  Added
# Forward    3     CrCL_power~cl [power]   -4451.133   -4463.213    12.080     0.0005  Added
# Forward    4     SEX_1~vc [cat]          -4463.213   -4473.985    10.772     0.0010  Added
# Forward    5     SEX_1~cl [cat]          -4473.985   -4477.111     3.126     0.0771  Not selected
# Backward   1     SEX_1~vc [cat]          -4473.985   -4463.214    10.770     0.0010  Retained
# ── SCM All Candidates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.942   -4419.355    23.413     0.0000  Added
# Forward    1     CrCL_power~cl [power]   -4395.942   -4418.167    22.225     0.0000  Not selected
# Forward    1     SEX_1~vc [cat]          -4395.942   -4417.810    21.869     0.0000  Not selected
# Forward    1     BW_power~cl [power]     -4395.942   -4413.574    17.632     0.0000  Not selected
# Forward    1     BMI_power~cl [power]    -4395.942   -4411.689    15.748     0.0001  Not selected
# Forward    1     BW_lin~vc [lin]         -4395.942   -4407.744    11.802     0.0006  Not selected
# Forward    1     BMI_power~vc [power]    -4395.942   -4402.714     6.773     0.0093  Not selected
# Forward    1     CrCL_lin~cl [lin]       -4395.942   -4402.045     6.104     0.0135  Not selected
# Forward    1     BMI_lin~vc [lin]        -4395.942   -4400.265     4.323     0.0376  Not selected
# Forward    1     BMI_lin~cl [lin]        -4395.942   -4399.699     3.757     0.0526  Not selected
# Forward    1     BW_lin~cl [lin]         -4395.942   -4399.413     3.471     0.0625  Not selected
# Forward    1     SEX_1~cl [cat]          -4395.942   -4398.338     2.397     0.1216  Not selected
# Forward    1     RACE_0~vc [cat]         -4395.942   -4397.632     1.690     0.1936  Not selected
# Forward    1     CrCL_power~vc [power]   -4395.942   -4397.422     1.480     0.2237  Not selected
# Forward    1     CrCL_lin~vc [lin]       -4395.942   -4397.022     1.081     0.2985  Not selected
# Forward    1     RACE_0~cl [cat]         -4395.942   -4396.003     0.061     0.8044  Not selected
# 
# Forward    2     BW_power~cl [power]     -4419.355   -4451.133    31.778     0.0000  Added
# Forward    2     CrCL_power~cl [power]   -4419.355   -4449.473    30.118     0.0000  Not selected
# Forward    2     BMI_power~cl [power]    -4419.355   -4445.847    26.491     0.0000  Not selected
# Forward    2     SEX_1~vc [cat]          -4419.355   -4431.041    11.686     0.0006  Not selected
# Forward    2     BMI_power~vc [power]    -4419.355   -4424.067     4.712     0.0300  Not selected
# Forward    2     RACE_0~cl [cat]         -4419.355   -4420.092     0.737     0.3906  Not selected
# Forward    2     CrCL_power~vc [power]   -4419.355   -4419.392     0.036     0.8486  Not selected
# Forward    2     RACE_0~vc [cat]         -4419.355   -4419.365     0.010     0.9189  Not selected
# Forward    2     SEX_1~cl [cat]          -4419.355   -4419.360     0.004     0.9466  Not selected
# Forward    2     BW_lin~cl [lin]         -4419.355   -4411.990    -7.365     1.0000  Not selected
# Forward    2     CrCL_lin~cl [lin]       -4419.355   -4413.772    -5.583     1.0000  Not selected
# Forward    2     CrCL_lin~vc [lin]       -4419.355   -4409.013   -10.342     1.0000  Not selected
# Forward    2     BMI_lin~cl [lin]        -4419.355   -4411.920    -7.435     1.0000  Not selected
# Forward    2     BMI_lin~vc [lin]        -4419.355   -4412.245    -7.110     1.0000  Not selected
# 
# Forward    3     CrCL_power~cl [power]   -4451.133   -4463.213    12.080     0.0005  Added
# Forward    3     SEX_1~vc [cat]          -4451.133   -4462.964    11.831     0.0006  Not selected
# Forward    3     BMI_power~vc [power]    -4451.133   -4455.826     4.693     0.0303  Not selected
# Forward    3     SEX_1~cl [cat]          -4451.133   -4455.119     3.986     0.0459  Not selected
# Forward    3     BMI_power~cl [power]    -4451.133   -4451.929     0.796     0.3724  Not selected
# Forward    3     RACE_0~cl [cat]         -4451.133   -4451.740     0.607     0.4359  Not selected
# Forward    3     CrCL_power~vc [power]   -4451.133   -4451.174     0.041     0.8397  Not selected
# Forward    3     RACE_0~vc [cat]         -4451.133   -4451.155     0.022     0.8826  Not selected
# Forward    3     CrCL_lin~cl [lin]       -4451.133   -4438.660   -12.473     1.0000  Not selected
# Forward    3     CrCL_lin~vc [lin]       -4451.133   -4435.959   -15.174     1.0000  Not selected
# Forward    3     BMI_lin~cl [lin]        -4451.133   -4435.805   -15.328     1.0000  Not selected
# Forward    3     BMI_lin~vc [lin]        -4451.133   -4436.340   -14.793     1.0000  Not selected
# 
# Forward    4     SEX_1~vc [cat]          -4463.213   -4473.985    10.772     0.0010  Added
# Forward    4     BMI_power~vc [power]    -4463.213   -4468.048     4.835     0.0279  Not selected
# Forward    4     SEX_1~cl [cat]          -4463.213   -4466.968     3.754     0.0527  Not selected
# Forward    4     BMI_power~cl [power]    -4463.213   -4464.579     1.366     0.2426  Not selected
# Forward    4     RACE_0~cl [cat]         -4463.213   -4463.407     0.194     0.6598  Not selected
# Forward    4     RACE_0~vc [cat]         -4463.213   -4463.244     0.031     0.8605  Not selected
# Forward    4     CrCL_power~vc [power]   -4463.213   -4463.222     0.008     0.9267  Not selected
# Forward    4     CrCL_lin~vc [lin]       -4463.213   -4443.547   -19.666     1.0000  Not selected
# Forward    4     BMI_lin~cl [lin]        -4463.213   -4444.580   -18.634     1.0000  Not selected
# Forward    4     BMI_lin~vc [lin]        -4463.213   -4445.040   -18.173     1.0000  Not selected
# 
# Forward    5     SEX_1~cl [cat]          -4473.985   -4477.111     3.126     0.0771  Not selected
# Forward    5     BMI_power~cl [power]    -4473.985   -4475.795     1.810     0.1785  Not selected
# Forward    5     RACE_0~cl [cat]         -4473.985   -4475.217     1.232     0.2670  Not selected
# Forward    5     CrCL_power~vc [power]   -4473.985   -4475.164     1.179     0.2775  Not selected
# Forward    5     RACE_0~vc [cat]         -4473.985   -4475.005     1.020     0.3124  Not selected
# Forward    5     BMI_power~vc [power]    -4473.985   -4474.996     1.012     0.3145  Not selected
# Forward    5     CrCL_lin~vc [lin]       -4473.985   -4440.855   -33.130     1.0000  Not selected
# Forward    5     BMI_lin~cl [lin]        -4473.985   -4441.595   -32.390     1.0000  Not selected
# Forward    5     BMI_lin~vc [lin]        -4473.985   -4442.137   -31.848     1.0000  Not selected
# 
# Backward   1     BW_power~vc [power]     -4473.985   -4448.694    25.291     0.0000  Retained
# Backward   1     BW_power~cl [power]     -4473.985   -4461.011    12.974     0.0003  Retained
# Backward   1     CrCL_power~cl [power]   -4473.985   -4463.034    10.951     0.0009  Retained
# Backward   1     SEX_1~vc [cat]          -4473.985   -4463.214    10.770     0.0010  Retained
# ── Final model ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Retained:
# BW_power~vc [power]
# BW_power~cl [power]
# CrCL_power~cl [power]
# SEX_1~vc [cat]
# ℹ Final model OFV: -4473.985
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [scn16_N80_full_fast] 20:51:44 DONE | elapsed 14.5 min (871 s) ===
# 
# There were 15 warnings (use warnings() to see them)
t16_N80_full_fast <- attr(res16_N80_full_fast, "elapsed_s")
saveRDS(res16_N80_full_fast,
        file.path(stage1_dir16_N80, "res_full_fast.rds"))
test16_N80_full_fast <- package_scm_result(
  "scn16_N80_full_scm_fast",
  res16_N80_full_fast, t16_N80_full_fast,
  scenario_id   = 16,
  refit_control = scm_focei_final
)
# calculating covariance matrix
# done
# → Calculating residuals/tables
# ✔ done
saveRDS(test16_N80_full_fast,
        file.path(stage1_dir16_N80, "test_full_fast.rds"))

test16_N80_full_fast
# $label
# [1] "scn16_N80_full_scm_fast"
# 
# $selected
# # A tibble: 4 × 5
#   var   covar shape theta_name        estimate
#   <chr> <chr> <chr> <chr>                <dbl>
# 1 vc    BW    power cov_BW_power_vc      1.00 
# 2 cl    BW    power cov_BW_power_cl      0.684
# 3 cl    CrCL  power cov_CrCL_power_cl    0.389
# 4 vc    SEX   1     cov_SEX_1_vc         0.301
# 
# $step_hist
#          step      covar var shape      objf      deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames   covarEffect
# eta.cl      1   BW_power  cl power -4413.574  17.631967332 -3513.393 -3475.829        10 3.841459 2.680447e-05       no    forward   cov_BW_power_cl  0.9084817332
# eta.cl1     1     BW_lin  cl   lin -4399.413   3.471055472 -3499.232 -3461.668        10 3.841459 6.245145e-02       no    forward     cov_BW_lin_cl  0.0013054830
# eta.vc      1   BW_power  vc power -4419.355  23.413461400 -3519.174 -3481.610        10 3.841459 1.306610e-06      yes    forward   cov_BW_power_vc  1.1120120021
# eta.vc1     1     BW_lin  vc   lin -4407.744  11.802031905 -3507.563 -3469.999        10 3.841459 5.916611e-04       no    forward     cov_BW_lin_vc  0.0078896306
# eta.cl2     1 CrCL_power  cl power -4418.167  22.225292875 -3517.986 -3480.422        10 3.841459 2.424587e-06       no    forward cov_CrCL_power_cl  0.5240732977
# eta.cl3     1   CrCL_lin  cl   lin -4402.045   6.103574150 -3501.864 -3464.300        10 3.841459 1.349088e-02       no    forward   cov_CrCL_lin_cl  0.0009985787
# eta.vc2     1 CrCL_power  vc power -4397.422   1.480429634 -3497.241 -3459.677        10 3.841459 2.237073e-01       no    forward cov_CrCL_power_vc  0.1747862892
# eta.vc3     1   CrCL_lin  vc   lin -4397.022   1.080830179 -3496.841 -3459.277        10 3.841459 2.985119e-01       no    forward   cov_CrCL_lin_vc  0.0016314046
# eta.cl4     1  BMI_power  cl power -4411.689  15.747608463 -3511.508 -3473.944        10 3.841459 7.237938e-05       no    forward  cov_BMI_power_cl  0.8019273708
# eta.cl5     1    BMI_lin  cl   lin -4399.699   3.756972789 -3499.518 -3461.954        10 3.841459 5.258771e-02       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc4     1  BMI_power  vc power -4402.714   6.772613843 -3502.533 -3464.969        10 3.841459 9.256717e-03       no    forward  cov_BMI_power_vc  0.6286474645
# eta.vc5     1    BMI_lin  vc   lin -4400.265   4.323282992 -3500.084 -3462.520        10 3.841459 3.759432e-02       no    forward    cov_BMI_lin_vc  0.0155213049
# eta.cl6     1      SEX_1  cl   cat -4398.338   2.396771020 -3498.157 -3460.593        10 3.841459 1.215860e-01       no    forward      cov_SEX_1_cl -0.1348785294
# eta.vc6     1      SEX_1  vc   cat -4417.810  21.868641154 -3517.629 -3480.065        10 3.841459 2.919664e-06       no    forward      cov_SEX_1_vc  0.4301377843
# eta.cl7     1     RACE_0  cl   cat -4396.003   0.061364394 -3495.822 -3458.258        10 3.841459 8.043523e-01       no    forward     cov_RACE_0_cl -0.0380652094
# eta.vc7     1     RACE_0  vc   cat -4397.632   1.690371210 -3497.451 -3459.887        10 3.841459 1.935520e-01       no    forward     cov_RACE_0_vc -0.2180457935
# eta.cl8     2   BW_power  cl power -4451.133  31.777972138 -3548.952 -3507.214        11 3.841459 1.728419e-08      yes    forward   cov_BW_power_cl  1.0090481488
# eta.cl11    2     BW_lin  cl   lin -4411.990  -7.364878408 -3509.809 -3468.071        11 3.841459 1.000000e+00       no    forward     cov_BW_lin_cl  0.0013054830
# eta.cl21    2 CrCL_power  cl power -4449.473  30.118256639 -3547.292 -3505.555        11 3.841459 4.064865e-08       no    forward cov_CrCL_power_cl  0.6153718887
# eta.cl31    2   CrCL_lin  cl   lin -4413.772  -5.583387416 -3511.591 -3469.853        11 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_cl  0.0009985787
# eta.vc8     2 CrCL_power  vc power -4419.392   0.036422416 -3517.211 -3475.473        11 3.841459 8.486458e-01       no    forward cov_CrCL_power_vc -0.0295864121
# eta.vc11    2   CrCL_lin  vc   lin -4409.013 -10.342497829 -3506.832 -3465.094        11 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc -0.0028816572
# eta.cl41    2  BMI_power  cl power -4445.847  26.491420092 -3543.666 -3501.928        11 3.841459 2.647112e-07       no    forward  cov_BMI_power_cl  1.0443050096
# eta.cl51    2    BMI_lin  cl   lin -4411.920  -7.435376841 -3509.739 -3468.001        11 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc21    2  BMI_power  vc power -4424.067   4.711567226 -3521.886 -3480.148        11 3.841459 2.996033e-02       no    forward  cov_BMI_power_vc -0.7893514874
# eta.vc31    2    BMI_lin  vc   lin -4412.245  -7.110264252 -3510.064 -3468.326        11 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0161404715
# eta.cl61    2      SEX_1  cl   cat -4419.360   0.004479203 -3517.179 -3475.441        11 3.841459 9.466399e-01       no    forward      cov_SEX_1_cl  0.0062176780
# eta.vc41    2      SEX_1  vc   cat -4431.041  11.685864939 -3528.860 -3487.122        11 3.841459 6.297670e-04       no    forward      cov_SEX_1_vc  0.3091304608
# eta.cl71    2     RACE_0  cl   cat -4420.092   0.736938966 -3517.911 -3476.173        11 3.841459 3.906432e-01       no    forward     cov_RACE_0_cl -0.1450503430
# eta.vc51    2     RACE_0  vc   cat -4419.365   0.010371563 -3517.184 -3475.447        11 3.841459 9.188830e-01       no    forward     cov_RACE_0_vc  0.0188040797
# eta.cl9     3 CrCL_power  cl power -4463.213  12.080331332 -3559.032 -3513.121        12 3.841459 5.095654e-04      yes    forward cov_CrCL_power_cl  0.3949712524
# eta.cl12    3   CrCL_lin  cl   lin -4438.660 -12.473497115 -3534.479 -3488.567        12 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_cl  0.0009985787
# eta.vc9     3 CrCL_power  vc power -4451.174   0.040935381 -3546.993 -3501.081        12 3.841459 8.396627e-01       no    forward cov_CrCL_power_vc -0.0262543920
# eta.vc12    3   CrCL_lin  vc   lin -4435.959 -15.174122371 -3531.778 -3485.866        12 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc -0.0016914004
# eta.cl22    3  BMI_power  cl power -4451.929   0.795666291 -3547.748 -3501.836        12 3.841459 3.723922e-01       no    forward  cov_BMI_power_cl  0.3091569738
# eta.cl32    3    BMI_lin  cl   lin -4435.805 -15.328351920 -3531.624 -3485.712        12 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc22    3  BMI_power  vc power -4455.826   4.692735177 -3551.645 -3505.733        12 3.841459 3.029039e-02       no    forward  cov_BMI_power_vc -0.7952841184
# eta.vc32    3    BMI_lin  vc   lin -4436.340 -14.792962526 -3532.159 -3486.247        12 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0073954280
# eta.cl42    3      SEX_1  cl   cat -4455.119   3.986077188 -3550.938 -3505.027        12 3.841459 4.587776e-02       no    forward      cov_SEX_1_cl -0.2314962091
# eta.vc42    3      SEX_1  vc   cat -4462.964  11.831031480 -3558.783 -3512.871        12 3.841459 5.825172e-04       no    forward      cov_SEX_1_vc  0.3073305474
# eta.cl52    3     RACE_0  cl   cat -4451.740   0.607007795 -3547.559 -3501.647        12 3.841459 4.359167e-01       no    forward     cov_RACE_0_cl  0.1093257544
# eta.vc52    3     RACE_0  vc   cat -4451.155   0.021794720 -3546.974 -3501.062        12 3.841459 8.826345e-01       no    forward     cov_RACE_0_vc  0.0217399098
# eta.vc10    4 CrCL_power  vc power -4463.222   0.008459491 -3557.041 -3506.955        13 3.841459 9.267175e-01       no    forward cov_CrCL_power_vc  0.0141391079
# eta.vc13    4   CrCL_lin  vc   lin -4443.547 -19.665952598 -3537.366 -3487.281        13 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc -0.0004684074
# eta.cl10    4  BMI_power  cl power -4464.579   1.365697884 -3558.398 -3508.313        13 3.841459 2.425524e-01       no    forward  cov_BMI_power_cl  0.3892663061
# eta.cl13    4    BMI_lin  cl   lin -4444.580 -18.633643027 -3538.399 -3488.313        13 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc23    4  BMI_power  vc power -4468.048   4.835006793 -3561.867 -3511.782        13 3.841459 2.788753e-02       no    forward  cov_BMI_power_vc -0.7887056815
# eta.vc33    4    BMI_lin  vc   lin -4445.040 -18.173391160 -3538.859 -3488.774        13 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0078427518
# eta.cl23    4      SEX_1  cl   cat -4466.968   3.754156805 -3560.787 -3510.701        13 3.841459 5.267636e-02       no    forward      cov_SEX_1_cl -0.2020882634
# eta.vc43    4      SEX_1  vc   cat -4473.985  10.771500113 -3567.804 -3517.718        13 3.841459 1.030749e-03      yes    forward      cov_SEX_1_vc  0.3081906737
# eta.cl33    4     RACE_0  cl   cat -4463.407   0.193720354 -3557.226 -3507.141        13 3.841459 6.598381e-01       no    forward     cov_RACE_0_cl  0.0597881204
# eta.vc53    4     RACE_0  vc   cat -4463.244   0.030867276 -3557.063 -3506.978        13 3.841459 8.605368e-01       no    forward     cov_RACE_0_vc  0.0263377498
# eta.vc14    5 CrCL_power  vc power -4475.164   1.179234401 -3566.983 -3512.724        14 3.841459 2.775119e-01       no    forward cov_CrCL_power_vc  0.0524172300
# eta.vc15    5   CrCL_lin  vc   lin -4440.855 -33.129717629 -3532.674 -3478.415        14 3.841459 1.000000e+00       no    forward   cov_CrCL_lin_vc  0.0001307784
# eta.cl14    5  BMI_power  cl power -4475.795   1.809707362 -3567.614 -3513.354        14 3.841459 1.785433e-01       no    forward  cov_BMI_power_cl  0.2937043767
# eta.cl15    5    BMI_lin  cl   lin -4441.595 -32.389587586 -3533.414 -3479.155        14 3.841459 1.000000e+00       no    forward    cov_BMI_lin_cl  0.0036900369
# eta.vc24    5  BMI_power  vc power -4474.996   1.011506230 -3566.815 -3512.556        14 3.841459 3.145423e-01       no    forward  cov_BMI_power_vc  0.0298406592
# eta.vc34    5    BMI_lin  vc   lin -4442.137 -31.848255087 -3533.956 -3479.696        14 3.841459 1.000000e+00       no    forward    cov_BMI_lin_vc -0.0067708754
# eta.cl24    5      SEX_1  cl   cat -4477.111   3.125840013 -3568.930 -3514.671        14 3.841459 7.706015e-02       no    forward      cov_SEX_1_cl -0.1325983626
# eta.cl34    5     RACE_0  cl   cat -4475.217   1.232027424 -3567.036 -3512.777        14 3.841459 2.670131e-01       no    forward     cov_RACE_0_cl  0.0635971984
# eta.vc44    5     RACE_0  vc   cat -4475.005   1.020446988 -3566.824 -3512.565        14 3.841459 3.124130e-01       no    forward     cov_RACE_0_vc -0.0167253220
# eta.cl16    1   BW_power  cl power -4461.011  12.974384411 -3556.830 -3510.918        12 6.634897 3.157817e-04 retained   backward   cov_BW_power_cl  0.6954515011
#           bsvReduction
# eta.cl    3.301946e+01
# eta.cl1  -1.997907e-13
# eta.vc    4.331752e+01
# eta.vc1   5.250347e-09
# eta.cl2   3.104095e+01
# eta.cl3  -1.997907e-13
# eta.vc2   9.174201e+00
# eta.vc3   2.330973e-03
# eta.cl4   2.735838e+01
# eta.cl5   4.155184e-02
# eta.vc4   2.018488e+01
# eta.vc5   2.420442e-07
# eta.cl6  -1.299902e+00
# eta.vc6   2.522587e+01
# eta.cl7   4.570263e-01
# eta.vc7   4.171365e+00
# eta.cl8   3.334699e+01
# eta.cl11  1.532540e-01
# eta.cl21  3.142949e+01
# eta.cl31  1.521120e-01
# eta.vc8  -2.659845e-01
# eta.vc11 -7.642134e+01
# eta.cl41  2.543172e+01
# eta.cl51  1.612334e-01
# eta.vc21  6.528204e+00
# eta.vc31 -7.628704e+01
# eta.cl61  7.775346e-02
# eta.vc41  1.321021e+01
# eta.cl71  7.270450e-01
# eta.vc51 -2.371172e-01
# eta.cl9   1.405358e+01
# eta.cl12 -4.980151e+01
# eta.vc9   1.971837e-01
# eta.vc12 -7.647146e+01
# eta.cl22  2.262342e-01
# eta.cl32 -4.978882e+01
# eta.vc22  7.303346e+00
# eta.vc32 -7.639440e+01
# eta.cl42  2.186661e-02
# eta.vc42  1.496117e+01
# eta.cl52  7.280446e-01
# eta.vc52  7.685914e-02
# eta.vc10 -1.135475e-01
# eta.vc13 -7.639553e+01
# eta.cl10  1.181416e+00
# eta.cl13 -7.428164e+01
# eta.vc23  7.499391e+00
# eta.vc33 -7.631242e+01
# eta.cl23  1.522263e-01
# eta.vc43  1.926044e+01
# eta.cl33  7.715297e-02
# eta.vc53 -3.592586e-02
# eta.vc14 -5.130787e+00
# eta.vc15 -1.184738e+02
# eta.cl14  3.793556e+00
# eta.cl15 -6.869606e+01
# eta.vc24 -5.566974e+00
# eta.vc34 -1.183825e+02
# eta.cl24  4.367539e+00
# eta.cl34  3.523347e+00
# eta.vc44 -5.125540e+00
# eta.cl16  1.322813e+01
#  [ reached 'max' / getOption("max.print") -- omitted 3 rows ]
# 
# $final_fit
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#        OBJF       AIC       BIC Log-likelihood Condition#(Cov) Condition#(Cor)
# FOCEi -4475 -3568.819 -3518.734        1796.41        412.9822        9.408444
# 
# ── Time (sec $time): ──
# 
#             setup optimize covariance preprocess postprocess table     other
# elapsed 0.0068215 49.56159   11.03013       0.16        0.01  0.07 0.2414577
# 
# ── Population Parameters ($parFixed or $parFixedDf): ──
# 
#                                                   Parameter    Est.      SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# lTVCL                                                       -0.3802 0.03607  9.488 0.6837 (0.6371, 0.7338)     33.0     0.608% 
# lTVQ                                                         0.6192  0.0187   3.02    1.857 (1.791, 1.927)                     
# lTVVc                                                         3.148 0.05791   1.84    23.28 (20.78, 26.08)     30.7      12.9% 
# lTVVp                                                         4.396 0.01193 0.2713    81.16 (79.28, 83.08)                     
# lTVKA             KA unidentifiable from this sparse design -0.3567   FIXED  FIXED                 -0.3567                     
# prop.err                                                     0.1054                                 0.1054                     
# cov_BW_power_vc                                               1.002  0.1723   17.2     1.002 (0.664, 1.34)                     
# cov_BW_power_cl                                              0.6842  0.1809  26.44  0.6842 (0.3296, 1.039)                     
# cov_CrCL_power_cl                                            0.3892  0.1072  27.55  0.3892 (0.179, 0.5994)                     
# cov_SEX_1_vc                                                 0.3005 0.08561  28.49 0.3005 (0.1327, 0.4683)                     
#  
#   Covariance Type ($covMethod): r,s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.172   
#  
# 
#   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found ($runInfo):
#    • gradient problems with covariance; see $scaleInfo 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring ($censInformation): No censoring
#   Minimization message ($message):  
#     Normal exit from bobyqa 
# 
# ── Fit Data (object is a modified tibble): ──
# # A tibble: 480 × 32
#   ID     TIME    DV  PRED     RES    WRES IPRED  IRES IWRES CPRED     CRES    CWRES eta.cl eta.vc   depot central peripheral1 BW_power_cl CrCL_power_cl cov_cl    cl BW_power_vc
#   <fct> <dbl> <dbl> <dbl>   <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl>  <dbl>   <dbl>   <dbl>       <dbl>       <dbl>         <dbl>  <dbl> <dbl>       <dbl>
# 1 1      0     0     0    0       0        0    0     0      0     0        0        0.156  0.144 1   e+2     0           0        -0.199      -0.00868 -0.207 0.649      -0.291
# 2 1      7.06  2.77  2.77 1.63e-4 3.77e-4  2.57 0.202 0.748  2.78 -0.00426 -0.00979  0.156  0.144 7.12e-1    51.7        34.6      -0.199      -0.00868 -0.207 0.649      -0.291
# 3 1     14.1   1.58  1.36 2.24e-1 9.97e-1  1.33 0.255 1.82   1.37  0.216    1.01     0.156  0.144 5.07e-3    26.7        51.7      -0.199      -0.00868 -0.207 0.649      -0.291
# # ℹ 477 more rows
# # ℹ 10 more variables: SEX_1_vc <dbl>, cov_vc <dbl>, vc <dbl>, q <dbl>, vp <dbl>, tad <dbl>, dosenum <dbl>, SEX_1 <dbl>, CrCL <dbl>, BW <dbl>
# # ℹ Use `print(n = ...)` to see more rows
# 
# $final_est
# # A tibble: 17 × 2
#    parameter estimate
#    <chr>        <dbl>
#  1 TVCL        0.684 
#  2 TVVc       23.3   
#  3 TVQ         1.86  
#  4 TVVp       81.2   
#  5 TVKA        0.7   
#  6 var_CL      0.103 
#  7 var_Vc      0.0903
#  8 cov_VcCL    0.0166
#  9 ResErr      0.105 
# 10 CLBW        0.684 
# 11 CLcrCL      0.389 
# 12 VcBW        1.00  
# 13 VcSEX       0.301 
# 14 CLBMI      NA     
# 15 VcBMI      NA     
# 16 VcCrCL     NA     
# 17 VcRACE     NA     
# 
# $rel_err
# # A tibble: 17 × 6
#    parameter true_value estimate   abs_err   rel_err rel_err_pct
#    <chr>          <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
#  1 TVCL           0.6     0.684   8.37e- 2  1.40e- 1    1.40e+ 1
#  2 TVVc          20      23.3     3.28e+ 0  1.64e- 1    1.64e+ 1
#  3 TVQ            1.8     1.86    5.75e- 2  3.19e- 2    3.19e+ 0
#  4 TVVp          80      81.2     1.16e+ 0  1.45e- 2    1.45e+ 0
#  5 TVKA           0.7     0.7     3.33e-16  4.76e-16    4.76e-14
#  6 CLBW           0.75    0.684  -6.58e- 2 -8.77e- 2   -8.77e+ 0
#  7 CLcrCL         0.5     0.389  -1.11e- 1 -2.22e- 1   -2.22e+ 1
#  8 VcBW           1       1.00    1.77e- 3  1.77e- 3    1.77e- 1
#  9 VcSEX          0.405   0.301  -1.05e- 1 -2.59e- 1   -2.59e+ 1
# 10 var_CL         0.1     0.103   3.32e- 3  3.32e- 2    3.32e+ 0
# 11 var_Vc         0.1     0.0903 -9.71e- 3 -9.71e- 2   -9.71e+ 0
# 12 cov_VcCL       0.02    0.0166 -3.39e- 3 -1.69e- 1   -1.69e+ 1
# 13 ResErr         0.1     0.105   5.40e- 3  5.40e- 2    5.40e+ 0
# 14 CLBMI         NA      NA      NA        NA          NA       
# 15 VcBMI         NA      NA      NA        NA          NA       
# 16 VcCrCL        NA      NA      NA        NA          NA       
# 17 VcRACE        NA      NA      NA        NA          NA       
# 
# $diag
# $diag$converged
# [1] TRUE
# 
# $diag$objf
# [1] -4475
# 
# $diag$cond_num
# [1] 412.9822
# 
# $diag$cond_num_sqrt
# [1] 20.32196
# 
# $diag$cov_ok
# [1] TRUE
# 
# $diag$message
# [1] "Normal exit from bobyqa"
# 
# 
# $parFixed
#                                                   Parameter   Estimate         SE       %RSE Back-transformed   CI Lower   CI Upper BSV(CV%) Shrink(SD)%
# lTVCL                                                       -0.3801686 0.03607010  9.4879230        0.6837461  0.6370770  0.7338340 32.99140   0.6079082
# lTVQ                                                         0.6192196 0.01870190  3.0202367        1.8574779  1.7906247  1.9268271       NA          NA
# lTVVc                                                        3.1476690 0.05790706  1.8396808       23.2817310 20.7837873 26.0798955 30.73936  12.9195978
# lTVVp                                                        4.3964033 0.01192604  0.2712681       81.1584413 79.2833939 83.0778334       NA          NA
# lTVKA             KA unidentifiable from this sparse design -0.3566749         NA         NA       -0.3566749         NA         NA       NA          NA
# prop.err                                                     0.1053955         NA         NA        0.1053955         NA         NA       NA          NA
# cov_BW_power_vc                                              1.0017715 0.17233943 17.2034673        1.0017715  0.6639924  1.3395505       NA          NA
# cov_BW_power_cl                                              0.6841879 0.18090041 26.4401660        0.6841879  0.3296296  1.0387461       NA          NA
# cov_CrCL_power_cl                                            0.3892281 0.10724026 27.5520302        0.3892281  0.1790411  0.5994152       NA          NA
# cov_SEX_1_vc                                                 0.3005307 0.08561311 28.4873142        0.3005307  0.1327320  0.4683293       NA          NA
# 
# $refit_done
# [1] TRUE
# 
# $runtime_sec
# [1] 871.1197

# 5.4 res16_N80_user_4true fast--------
res16_N80_user_fast <- runSCM_traced(
  label       = "scn16_N80_user_fast",
  data        = ds16_01_N80,
  fit         = fit_base16_N80_linCmt,
  pairsVec    = list(list(var = "vc", covar = "BW",   shapes = "power"),
                     list(var = "cl", covar = "BW",   shapes = "power"),
                     list(var = "cl", covar = "CrCL", shapes = "power"),
                     list(var = "vc", covar = "SEX",  shapes = "cat")),
  catvarsVec  = "SEX",
  searchType  = "scm",
  control     = scm_focei_screen,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
# 
# === [scn16_N80_user_fast] 21:16:13 starting scm      | 4 candidate(s), 3 worker(s) ===
# ── SCM Summary ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -4395.942
# ℹ Base model params : 9
# ℹ Search type : scm
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (vc, cl)
# ℹ Covariates : 3
# (BW_power, CrCL_power, SEX_1)
# ℹ Total candidates : 4
# ── Categorical covariates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ SEX: reference = '0'
# Indicators (1): SEX_1
# ── Relationships to test ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ vc [power]
# 2. BW_power ~ cl [power]
# 3. CrCL_power ~ cl [power]
# 4. SEX_1 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Proceed with SCM? [y/N]: 
# ── SCM Step Summary ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.942   -4419.355    23.413     0.0000  Added
# Forward    2     BW_power~cl [power]     -4419.355   -4451.133    31.778     0.0000  Added
# Forward    3     CrCL_power~cl [power]   -4451.133   -4463.213    12.080     0.0005  Added
# Forward    4     SEX_1~vc [cat]          -4463.213   -4473.985    10.772     0.0010  Added
# Backward   1     SEX_1~vc [cat]          -4473.985   -4463.214    10.770     0.0010  Retained
# ── SCM All Candidates ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                  Ref OFV         OFV      dOFV    p-value  Decision 
# --------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~vc [power]     -4395.942   -4419.355    23.413     0.0000  Added
# Forward    1     CrCL_power~cl [power]   -4395.942   -4418.167    22.225     0.0000  Not selected
# Forward    1     SEX_1~vc [cat]          -4395.942   -4417.810    21.869     0.0000  Not selected
# Forward    1     BW_power~cl [power]     -4395.942   -4413.574    17.632     0.0000  Not selected
# 
# Forward    2     BW_power~cl [power]     -4419.355   -4451.133    31.778     0.0000  Added
# Forward    2     CrCL_power~cl [power]   -4419.355   -4449.473    30.118     0.0000  Not selected
# Forward    2     SEX_1~vc [cat]          -4419.355   -4431.041    11.686     0.0006  Not selected
# 
# Forward    3     CrCL_power~cl [power]   -4451.133   -4463.213    12.080     0.0005  Added
# Forward    3     SEX_1~vc [cat]          -4451.133   -4462.964    11.831     0.0006  Not selected
# 
# Forward    4     SEX_1~vc [cat]          -4463.213   -4473.985    10.772     0.0010  Added
# 
# Backward   1     BW_power~vc [power]     -4473.985   -4448.694    25.291     0.0000  Retained
# Backward   1     BW_power~cl [power]     -4473.985   -4461.011    12.974     0.0003  Retained
# Backward   1     CrCL_power~cl [power]   -4473.985   -4463.034    10.951     0.0009  Retained
# Backward   1     SEX_1~vc [cat]          -4473.985   -4463.214    10.770     0.0010  Retained
# ── Final model ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Retained:
# BW_power~vc [power]
# BW_power~cl [power]
# CrCL_power~cl [power]
# SEX_1~vc [cat]
# ℹ Final model OFV: -4473.985
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [scn16_N80_user_fast] 21:20:58 DONE | elapsed 4.7 min (285 s) ===
# 
t16_N80_user_fast <- attr(res16_N80_user_fast, "elapsed_s")
saveRDS(res16_N80_user_fast,
        file.path(stage1_dir16_N80, "res_user_fast.rds"))
test16_N80_user_fast <- package_scm_result(
  "scn16_N80_user_4tr_fast",
  res16_N80_user_fast, t16_N80_user_fast,
  scenario_id   = 16,
  refit_control = scm_focei_final
)
# calculating covariance matrix
# done
# → Calculating residuals/tables
# ✔ done
saveRDS(test16_N80_user_fast,
        file.path(stage1_dir16_N80, "test_user_fast.rds"))

test16_N80_user_fast
# $label
# [1] "scn16_N80_user_4tr_fast"
# 
# $selected
# # A tibble: 4 × 5
#   var   covar shape theta_name        estimate
#   <chr> <chr> <chr> <chr>                <dbl>
# 1 vc    BW    power cov_BW_power_vc      1.00 
# 2 cl    BW    power cov_BW_power_cl      0.684
# 3 cl    CrCL  power cov_CrCL_power_cl    0.389
# 4 vc    SEX   1     cov_SEX_1_vc         0.301
# 
# $step_hist
#          step      covar var shape      objf deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames covarEffect bsvReduction
# eta.vc      1   BW_power  vc power -4419.355 23.41346 -3519.174 -3481.610        10 3.841459 1.306610e-06      yes    forward   cov_BW_power_vc   1.1120120     43.31752
# eta.cl      1   BW_power  cl power -4413.574 17.63197 -3513.393 -3475.829        10 3.841459 2.680447e-05       no    forward   cov_BW_power_cl   0.9084817     33.01946
# eta.cl1     1 CrCL_power  cl power -4418.167 22.22529 -3517.986 -3480.422        10 3.841459 2.424587e-06       no    forward cov_CrCL_power_cl   0.5240733     31.04095
# eta.vc1     1      SEX_1  vc   cat -4417.810 21.86864 -3517.629 -3480.065        10 3.841459 2.919664e-06       no    forward      cov_SEX_1_vc   0.4301378     25.22587
# eta.cl2     2   BW_power  cl power -4451.133 31.77797 -3548.952 -3507.214        11 3.841459 1.728419e-08      yes    forward   cov_BW_power_cl   1.0090481     33.34699
# eta.cl11    2 CrCL_power  cl power -4449.473 30.11826 -3547.292 -3505.555        11 3.841459 4.064865e-08       no    forward cov_CrCL_power_cl   0.6153719     31.42949
# eta.vc2     2      SEX_1  vc   cat -4431.041 11.68586 -3528.860 -3487.122        11 3.841459 6.297670e-04       no    forward      cov_SEX_1_vc   0.3091305     13.21021
# eta.cl3     3 CrCL_power  cl power -4463.213 12.08033 -3559.032 -3513.121        12 3.841459 5.095654e-04      yes    forward cov_CrCL_power_cl   0.3949713     14.05358
# eta.vc3     3      SEX_1  vc   cat -4462.964 11.83103 -3558.783 -3512.871        12 3.841459 5.825172e-04       no    forward      cov_SEX_1_vc   0.3073305     14.96117
# eta.vc4     4      SEX_1  vc   cat -4473.985 10.77150 -3567.804 -3517.718        13 3.841459 1.030749e-03      yes    forward      cov_SEX_1_vc   0.3081907     19.26044
# eta.vc5     1   BW_power  vc power -4448.694 25.29123 -3544.513 -3498.601        12 6.634897 4.929446e-07 retained   backward   cov_BW_power_vc   1.0068294     37.93463
# eta.cl4     1   BW_power  cl power -4461.011 12.97438 -3556.830 -3510.918        12 6.634897 3.157817e-04 retained   backward   cov_BW_power_cl   0.6954515     13.22813
# eta.cl12    1 CrCL_power  cl power -4463.034 10.95089 -3558.853 -3512.941        12 6.634897 9.355864e-04 retained   backward cov_CrCL_power_cl   0.3814720     10.85347
# eta.vc11    1      SEX_1  vc   cat -4463.214 10.77042 -3559.034 -3513.122        12 6.634897 1.031350e-03 retained   backward      cov_SEX_1_vc   0.3081907     19.25957
# 
# $final_fit
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#        OBJF       AIC       BIC Log-likelihood Condition#(Cov) Condition#(Cor)
# FOCEi -4475 -3568.819 -3518.734        1796.41        412.9822        9.408444
# 
# ── Time (sec $time): ──
# 
#             setup optimize covariance preprocess postprocess table     other
# elapsed 0.0269752 48.55235   11.56133       0.22        0.03  0.08 0.2593399
# 
# ── Population Parameters ($parFixed or $parFixedDf): ──
# 
#                                                   Parameter    Est.      SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# lTVCL                                                       -0.3802 0.03607  9.488 0.6837 (0.6371, 0.7338)     33.0     0.608% 
# lTVQ                                                         0.6192  0.0187   3.02    1.857 (1.791, 1.927)                     
# lTVVc                                                         3.148 0.05791   1.84    23.28 (20.78, 26.08)     30.7      12.9% 
# lTVVp                                                         4.396 0.01193 0.2713    81.16 (79.28, 83.08)                     
# lTVKA             KA unidentifiable from this sparse design -0.3567   FIXED  FIXED                 -0.3567                     
# prop.err                                                     0.1054                                 0.1054                     
# cov_BW_power_vc                                               1.002  0.1723   17.2     1.002 (0.664, 1.34)                     
# cov_BW_power_cl                                              0.6842  0.1809  26.44  0.6842 (0.3296, 1.039)                     
# cov_CrCL_power_cl                                            0.3892  0.1072  27.55  0.3892 (0.179, 0.5994)                     
# cov_SEX_1_vc                                                 0.3005 0.08561  28.49 0.3005 (0.1327, 0.4683)                     
#  
#   Covariance Type ($covMethod): r,s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.172   
#  
# 
#   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found ($runInfo):
#    • gradient problems with covariance; see $scaleInfo 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring ($censInformation): No censoring
#   Minimization message ($message):  
#     Normal exit from bobyqa 
# 
# ── Fit Data (object is a modified tibble): ──
# # A tibble: 480 × 32
#   ID     TIME    DV  PRED     RES    WRES IPRED  IRES IWRES CPRED     CRES    CWRES eta.cl eta.vc   depot central peripheral1 BW_power_cl CrCL_power_cl cov_cl    cl BW_power_vc
#   <fct> <dbl> <dbl> <dbl>   <dbl>   <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl>    <dbl>  <dbl>  <dbl>   <dbl>   <dbl>       <dbl>       <dbl>         <dbl>  <dbl> <dbl>       <dbl>
# 1 1      0     0     0    0       0        0    0     0      0     0        0        0.156  0.144 1   e+2     0           0        -0.199      -0.00868 -0.207 0.649      -0.291
# 2 1      7.06  2.77  2.77 1.63e-4 3.77e-4  2.57 0.202 0.748  2.78 -0.00426 -0.00979  0.156  0.144 7.12e-1    51.7        34.6      -0.199      -0.00868 -0.207 0.649      -0.291
# 3 1     14.1   1.58  1.36 2.24e-1 9.97e-1  1.33 0.255 1.82   1.37  0.216    1.01     0.156  0.144 5.07e-3    26.7        51.7      -0.199      -0.00868 -0.207 0.649      -0.291
# # ℹ 477 more rows
# # ℹ 10 more variables: SEX_1_vc <dbl>, cov_vc <dbl>, vc <dbl>, q <dbl>, vp <dbl>, tad <dbl>, dosenum <dbl>, SEX_1 <dbl>, CrCL <dbl>, BW <dbl>
# # ℹ Use `print(n = ...)` to see more rows
# 
# $final_est
# # A tibble: 17 × 2
#    parameter estimate
#    <chr>        <dbl>
#  1 TVCL        0.684 
#  2 TVVc       23.3   
#  3 TVQ         1.86  
#  4 TVVp       81.2   
#  5 TVKA        0.7   
#  6 var_CL      0.103 
#  7 var_Vc      0.0903
#  8 cov_VcCL    0.0166
#  9 ResErr      0.105 
# 10 CLBW        0.684 
# 11 CLcrCL      0.389 
# 12 VcBW        1.00  
# 13 VcSEX       0.301 
# 14 CLBMI      NA     
# 15 VcBMI      NA     
# 16 VcCrCL     NA     
# 17 VcRACE     NA     
# 
# $rel_err
# # A tibble: 17 × 6
#    parameter true_value estimate   abs_err   rel_err rel_err_pct
#    <chr>          <dbl>    <dbl>     <dbl>     <dbl>       <dbl>
#  1 TVCL           0.6     0.684   8.37e- 2  1.40e- 1    1.40e+ 1
#  2 TVVc          20      23.3     3.28e+ 0  1.64e- 1    1.64e+ 1
#  3 TVQ            1.8     1.86    5.75e- 2  3.19e- 2    3.19e+ 0
#  4 TVVp          80      81.2     1.16e+ 0  1.45e- 2    1.45e+ 0
#  5 TVKA           0.7     0.7     3.33e-16  4.76e-16    4.76e-14
#  6 CLBW           0.75    0.684  -6.58e- 2 -8.77e- 2   -8.77e+ 0
#  7 CLcrCL         0.5     0.389  -1.11e- 1 -2.22e- 1   -2.22e+ 1
#  8 VcBW           1       1.00    1.77e- 3  1.77e- 3    1.77e- 1
#  9 VcSEX          0.405   0.301  -1.05e- 1 -2.59e- 1   -2.59e+ 1
# 10 var_CL         0.1     0.103   3.32e- 3  3.32e- 2    3.32e+ 0
# 11 var_Vc         0.1     0.0903 -9.71e- 3 -9.71e- 2   -9.71e+ 0
# 12 cov_VcCL       0.02    0.0166 -3.39e- 3 -1.69e- 1   -1.69e+ 1
# 13 ResErr         0.1     0.105   5.40e- 3  5.40e- 2    5.40e+ 0
# 14 CLBMI         NA      NA      NA        NA          NA       
# 15 VcBMI         NA      NA      NA        NA          NA       
# 16 VcCrCL        NA      NA      NA        NA          NA       
# 17 VcRACE        NA      NA      NA        NA          NA       
# 
# $diag
# $diag$converged
# [1] TRUE
# 
# $diag$objf
# [1] -4475
# 
# $diag$cond_num
# [1] 412.9822
# 
# $diag$cond_num_sqrt
# [1] 20.32196
# 
# $diag$cov_ok
# [1] TRUE
# 
# $diag$message
# [1] "Normal exit from bobyqa"
# 
# 
# $parFixed
#                                                   Parameter   Estimate         SE       %RSE Back-transformed   CI Lower   CI Upper BSV(CV%) Shrink(SD)%
# lTVCL                                                       -0.3801686 0.03607010  9.4879230        0.6837461  0.6370770  0.7338340 32.99140   0.6079082
# lTVQ                                                         0.6192196 0.01870190  3.0202367        1.8574779  1.7906247  1.9268271       NA          NA
# lTVVc                                                        3.1476690 0.05790706  1.8396808       23.2817310 20.7837873 26.0798955 30.73936  12.9195978
# lTVVp                                                        4.3964033 0.01192604  0.2712681       81.1584413 79.2833939 83.0778334       NA          NA
# lTVKA             KA unidentifiable from this sparse design -0.3566749         NA         NA       -0.3566749         NA         NA       NA          NA
# prop.err                                                     0.1053955         NA         NA        0.1053955         NA         NA       NA          NA
# cov_BW_power_vc                                              1.0017715 0.17233943 17.2034673        1.0017715  0.6639924  1.3395505       NA          NA
# cov_BW_power_cl                                              0.6841879 0.18090041 26.4401660        0.6841879  0.3296296  1.0387461       NA          NA
# cov_CrCL_power_cl                                            0.3892281 0.10724026 27.5520302        0.3892281  0.1790411  0.5994152       NA          NA
# cov_SEX_1_vc                                                 0.3005307 0.08561311 28.4873142        0.3005307  0.1327320  0.4683293       NA          NA
# 
# $refit_done
# [1] TRUE
# 
# $runtime_sec
# [1] 284.9108



# 6.1 Pilot Tset- smoke run on 1 ds01 ---------
parallel::detectCores(logical = FALSE)
# [1] 4
.canon_shape_tbl <- function(rel) {
  rel %>%
    dplyr::mutate(shape = ifelse(grepl("^[0-9]+$", shape), "cat", shape)) %>%
    dplyr::distinct(var, covar, shape)
}
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
compute_power_block <- function(per_ds, n_true, cn_cor_cut = 1000) {
  N      <- nrow(per_ds)
  ok_cn  <- per_ds$cond_num_cor < cn_cor_cut & !is.na(per_ds$cond_num_cor)
  ok_min <- per_ds$converged
  ok_min[is.na(ok_min)] <- FALSE

  power_main <- tibble::tibble(
    metric = c("Power", "PowerCN", "PowerMinSuc"),
    num    = c(
      sum(per_ds$exact_match, na.rm = TRUE),
      sum(per_ds$exact_match & ok_cn,  na.rm = TRUE),
      sum(per_ds$exact_match & ok_min, na.rm = TRUE)
    ),
    denom  = c(N, sum(ok_cn), sum(ok_min))
  ) %>%
    dplyr::mutate(value = num / denom)

  rel_power <- tibble::tibble(
    k                   = seq_len(n_true),
    n_at_least_k        = vapply(seq_len(n_true), function(k)
      sum(per_ds$n_true_hit >= k, na.rm = TRUE), integer(1)),
    fraction_at_least_k = vapply(seq_len(n_true), function(k)
      mean(per_ds$n_true_hit >= k, na.rm = TRUE), numeric(1))
  )

  list(power = power_main, rel_power = rel_power)
}
compute_rmrse_block <- function(per_ds, true_long, scenario_id,
                                pop_params, cov_params) {
  truth <- true_long %>%
    dplyr::filter(scenario == scenario_id,
                  parameter %in% c(pop_params, cov_params)) %>%
    dplyr::select(parameter, true_value)

  long <- purrr::imap_dfr(per_ds$final_est, function(est, i) {
    if (is.null(est)) return(NULL)
    dplyr::mutate(est,
                  dataset_id  = per_ds$dataset_id[i],
                  exact_match = per_ds$exact_match[i])
  }) %>%
    dplyr::inner_join(truth, by = "parameter") %>%
    dplyr::mutate(
      rel_err = (estimate - true_value) / true_value,
      sq_rel  = rel_err^2,
      abs_rel = abs(rel_err)
    )

  rmrse_one <- function(df) {
    df <- dplyr::filter(df, is.finite(sq_rel))
    if (nrow(df) == 0L) {
      tibble::tibble(RMRSE_pct = NA_real_, MARE_pct = NA_real_, n_used = 0L)
    } else {
      tibble::tibble(
        RMRSE_pct = 100 * sqrt(mean(df$sq_rel)),
        MARE_pct  = 100 * stats::median(df$abs_rel),
        n_used    = nrow(df)
      )
    }
  }

  uncond <- long %>%
    dplyr::group_by(parameter) %>%
    dplyr::group_modify(~ rmrse_one(.x)) %>%
    dplyr::ungroup()

  cond <- long %>%
    dplyr::filter(exact_match) %>%
    dplyr::group_by(parameter) %>%
    dplyr::group_modify(~ rmrse_one(.x)) %>%
    dplyr::ungroup()

  list(rmrse_uncond = uncond, rmrse_cond = cond)
}
run_one_dataset_scn16_N80 <- function(ds_id, save_dir,
                                       sim_long     = sim_obs_scn16_N80,
                                       base_fn      = base_2cmt_oral_linCmt,
                                       screen_ctrl  = scm_focei_screen,
                                       vars_vec     = scm16_vars,
                                       covars_vec   = scm16_covars,
                                       catvars_vec  = scm16_catvars,
                                       shapes_vec   = scm16_shapes,
                                       true_long    = true_params,
                                       workers      = 3L,
                                       keep_res     = TRUE,
                                       force_rerun  = FALSE,
                                       scenario_id  = 16L) {
  ds_tag <- sprintf("ds%02d", ds_id)

  ## -- Resume check FIRST: before any expensive work or promise forcing
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  test_path <- file.path(save_dir, sprintf("test_full_fast_%s.rds", ds_tag))
  if (!force_rerun && file.exists(test_path)) {
    message(sprintf(">>> [%s] cached, skipping (delete %s or pass force_rerun=TRUE to redo)",
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

  message(sprintf("\n>>> [%s] starting dataset %d at %s",
                  ds_tag, ds_id, format(Sys.time(), "%H:%M:%S")))

  ## -- Pre-flight: surface missing-dependency errors BEFORE expensive SCM.
  ##    Force evaluation of every promise so a missing true_params (or any
  ##    other default argument) fails loudly here, not silently 25 min later.
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
      label       = sprintf("scn%02d_%s_N80_full_fast", scenario_id, ds_tag),
      data        = ds_i,
      fit         = fit_base_i,
      varsVec     = vars_vec,
      covarsVec   = covars_vec,
      catvarsVec  = catvars_vec,
      shapes      = shapes_vec,
      searchType  = "scm",
      control     = screen_ctrl,
      saveModels  = FALSE,
      workers     = workers,
      print       = 100,
      maxRetries  = 0L,
      confirm     = FALSE        # no interactive y/n prompt
    )
    t_scm_sec <- as.numeric(attr(res_i, "elapsed_s"))

    ## Save the (expensive) SCM result FIRST -- before any packaging.
    ## Then if package_scm_result fails downstream we can replay it
    ## offline without losing the SCM compute.
    if (keep_res) {
      saveRDS(res_i, file.path(save_dir,
                               sprintf("res_full_fast_%s.rds", ds_tag)))
    }

    test_i <- package_scm_result(
      label         = sprintf("scn%02d_%s_N80_full_fast", scenario_id, ds_tag),
      scm_res       = res_i,
      runtime_sec   = t_scm_sec,
      true_long     = true_long,         # explicit -- no lazy global lookup
      scenario_id   = scenario_id
    )

    saveRDS(test_i, test_path)

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
    err_path <- file.path(save_dir, sprintf("%s_ERROR.txt", ds_tag))
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
out_dir <- "simulated_virtual_dataset"
out_dir_v2_N80   <- "simulated_virtual_dataset_eta_filtered_N80"
stage1_dir16_N80 <- file.path(out_dir_v2_N80 , "stage1_smoke_scn16_ds01_N80")
sim_obs_scn16_N80 <- readRDS(file.path(out_dir_v2_N80, "sim_obs_scenario_16.rds"))
stage1_pilot_scn16_N80 <- file.path(out_dir_v2_N80,
                                     "stage1_pilot_scn16_N80")
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
scm16_vars       <- c("cl", "vc")
scm16_covars     <- c("BW", "CrCL", "BMI")
scm16_catvars    <- c("SEX", "RACE")
scm16_shapes     <- c("power", "lin")
## Verify all required symbols and that ds01 needs running (no cached test file)
need <- c("run_one_dataset_scn16_N80", "sim_obs_scn16_N80",
          "base_2cmt_oral_linCmt", "scm_focei_screen", "scm_focei_final",
          "scm16_vars", "scm16_covars", "scm16_catvars", "scm16_shapes",
          "to_nm_dataset", "true_params", "stage1_pilot_scn16_N80",
          "package_scm_result", "runSCM_traced")
have <- vapply(need, exists, logical(1))
missing <- need[!have]

cat("All required symbols present:", all(have), "\n")
if (length(missing)) cat("MISSING:", paste(missing, collapse = ", "), "\n")

## Check if ds01 already has a cached test file (would trigger resume, not actual run)
test_path <- file.path(stage1_pilot_scn16_N80, "test_full_fast_ds01.rds")
cat("\nstage1_pilot_scn16_N80 =", stage1_pilot_scn16_N80, "\n")
cat("ds01 test file exists?     ", file.exists(test_path), "\n")
if (file.exists(test_path)) {
  cat("File size:", round(file.size(test_path) / 1024, 1), "KB\n")
  cat("Modified:", format(file.mtime(test_path)), "\n")
}

## Confirm function signature includes the new args
cat("\nDriver formals:\n")
print(names(formals(run_one_dataset_scn16_N80)))
# All required symbols present: TRUE 
# 
# stage1_pilot_scn16_N80 = simulated_virtual_dataset_eta_filtered_N80/stage1_pilot_scn16_N80 
# ds01 test file exists?      FALSE 
# 
# Driver formals:
#  [1] "ds_id"        "save_dir"     "sim_long"     "base_fn"      "screen_ctrl"  "final_ctrl"   "vars_vec"     "covars_vec"   "catvars_vec"  "shapes_vec"   "true_long"   
# [12] "refit_policy" "workers"      "keep_res"     "force_rerun"  "scenario_id" 
## Smoke test: run ds01 with all defaults.
##   refit_policy = "full" (covMethod = "r,s" + tables)
##   workers      = 3L     (inner parallelism for the 4-core box)
##   keep_res     = TRUE   (keep res_full_fast_ds01.rds for debugging)
##   force_rerun  = FALSE  (would skip if cache existed; cache does not)
t_smoke_ds01 <- system.time(
  smoke_ds01 <- run_one_dataset_scn16_N80(
    ds_id    = 1L,
    save_dir = stage1_pilot_scn16_N80
  )
)

cat("\n\n=== Smoke test wall-clock:", round(t_smoke_ds01["elapsed"] / 60, 2), "min ===\n")
cat("=== Return structure: ===\n")
str(smoke_ds01, max.level = 1)
# 
# >>> [ds01] starting dataset 1 at 09:23:25
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# === [scn16_ds01_N80_full_fast] 09:23:37 starting scm      | NA candidate(s), 3 worker(s) ===
# ── SCM Summary ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -4395.942
# ℹ Base model params : 9
# ℹ Search type : scm
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (cl, vc)
# ℹ Covariates : 8
# (BW_power, BW_lin, CrCL_power, CrCL_lin, BMI_power, BMI_lin, SEX_1, RACE_0)
# ℹ Total candidates : 16
# ── Categorical covariates ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ SEX: reference = '0'
# Indicators (1): SEX_1
# ℹ RACE: reference = '1'
# Indicators (1): RACE_0
# ── Relationships to test ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ cl [power]
# 2. BW_lin ~ cl [lin]
# 3. BW_power ~ vc [power]
# 4. BW_lin ~ vc [lin]
# 5. CrCL_power ~ cl [power]
# 6. CrCL_lin ~ cl [lin]
# 7. CrCL_power ~ vc [power]
# 8. CrCL_lin ~ vc [lin]
# 9. BMI_power ~ cl [power]
# 10. BMI_lin ~ cl [lin]
# 11. BMI_power ~ vc [power]
# 12. BMI_lin ~ vc [lin]
# 13. SEX_1 ~ cl [cat]
# 14. SEX_1 ~ vc [cat]
# 15. RACE_0 ~ cl [cat]
# 16. RACE_0 ~ vc [cat]


# 6.2 Pilot Test-continued for ds02-04 -----------
run_one_dataset_scn16_N80 <- function(ds_id, save_dir,
                                       sim_long     = sim_obs_scn16_N80,
                                       base_fn      = base_2cmt_oral_linCmt,
                                       screen_ctrl  = scm_focei_screen,
                                       vars_vec     = scm16_vars,
                                       covars_vec   = scm16_covars,
                                       catvars_vec  = scm16_catvars,
                                       shapes_vec   = scm16_shapes,
                                       true_long    = true_params,
                                       workers      = 3L,
                                       keep_res     = TRUE,
                                       force_rerun  = FALSE,
                                       scenario_id  = 16L) {
  ds_tag <- sprintf("ds%02d", ds_id)

  ## -- Resume check FIRST: before any expensive work or promise forcing
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  test_path <- file.path(save_dir, sprintf("test_full_fast_%s.rds", ds_tag))
  if (!force_rerun && file.exists(test_path)) {
    message(sprintf(">>> [%s] cached, skipping (delete %s or pass force_rerun=TRUE to redo)",
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

  message(sprintf("\n>>> [%s] starting dataset %d at %s",
                  ds_tag, ds_id, format(Sys.time(), "%H:%M:%S")))

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
      label       = sprintf("scn%02d_%s_N80_full_fast", scenario_id, ds_tag),
      data        = ds_i,
      fit         = fit_base_i,
      varsVec     = vars_vec,
      covarsVec   = covars_vec,
      catvarsVec  = catvars_vec,
      shapes      = shapes_vec,
      searchType  = "scm",
      control     = screen_ctrl,
      saveModels  = FALSE,
      workers     = workers,
      print       = 100,
      maxRetries  = 0L,
      confirm     = FALSE        # no interactive y/n prompt
    )
    t_scm_sec <- as.numeric(attr(res_i, "elapsed_s"))

    ## Save the (expensive) SCM result FIRST -- before any packaging.
    ## Then if package_scm_result fails downstream we can replay it
    ## offline without losing the SCM compute.
    if (keep_res) {
      saveRDS(res_i, file.path(save_dir,
                               sprintf("res_full_fast_%s.rds", ds_tag)))
    }

    test_i <- package_scm_result(
      label         = sprintf("scn%02d_%s_N80_full_fast", scenario_id, ds_tag),
      scm_res       = res_i,
      runtime_sec   = t_scm_sec,
      true_long     = true_long,         # explicit -- no lazy global lookup
      scenario_id   = scenario_id
    )

    saveRDS(test_i, test_path)

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
    err_path <- file.path(save_dir, sprintf("%s_ERROR.txt", ds_tag))
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
t_pilot_total <- system.time(
  scm_pilot_scn16_N80 <- purrr::map(
    1:4,
    run_one_dataset_scn16_N80,
    save_dir = stage1_pilot_scn16_N80
  ) %>%
    rlang::set_names(sprintf("ds%02d", 1:4))
)
# >>> [ds01] cached, skipping (delete test_full_fast_ds01.rds or pass force_rerun=TRUE to redo)
# 
# >>> [ds02] starting dataset 2 at 09:53:39
# ℹ parameter labels from comments will be replaced by 'label()'
# done