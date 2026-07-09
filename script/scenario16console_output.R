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
diagnose_fit <- function(fit) {
  if (is.null(fit)) {
    return(list(converged = NA, objf = NA_real_,
                cond_num = NA_real_, cond_num_sqrt = NA_real_,
                cov_ok = NA, message = NA_character_))
  }
  conv_code <- if (!is.null(fit$convergence)) fit$convergence else NA_integer_
  cn <- fit$conditionNumber
  if (is.null(cn)) cn <- fit$conditionNumberCov
  if (is.null(cn)) cn <- fit$conditionNumberTheta
  cn_val <- if (is.null(cn)) NA_real_ else as.numeric(cn)
  list(
    converged     = isTRUE(conv_code == 0L) &&
                      !is.null(fit$objf) && is.finite(fit$objf),
    objf          = if (!is.null(fit$objf)) fit$objf else NA_real_,
    cond_num      = cn_val,
    cond_num_sqrt = if (is.na(cn_val)) NA_real_ else sqrt(cn_val),
    cov_ok        = isTRUE(!is.null(fit$cov) && all(is.finite(diag(fit$cov)))),
    message       = if (!is.null(fit$message)) as.character(fit$message) else NA_character_
  )
}
extract_params_long <- function(fit, includeCov = TRUE) {
  theta <- fit$theta
  om    <- fit$omega

  fe_long <- tibble::tibble(
    parameter = c("CL", "Vc", "Q", "Vp", "KA"),
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

  ## Covariate effects: read from fit$theta if present, else NA.
  ##   The theta names emitted by runSCM() preserve the original case of the
  ##   data column (e.g. "cov_BW_power_cl" when the data has column "BW"),
  ##   so the regex must match case-insensitively to catch both "BW" and
  ##   "wt" / "bw" parameterisations.
  cov_names_in_fit <- names(theta)
  cov_map <- list(
    CLBW   = c("TH_BW_CL",   grep("^cov_(bw|wt)_power_cl$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    CLcrCL = c("TH_CRCL_CL", grep("^cov_crcl_power_cl$",      cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcBW   = c("TH_BW_VC",   grep("^cov_(bw|wt)_power_vc$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcSEX  = c("TH_SEX_VC",  grep("^cov_sex(_male)?_cat_vc$", cov_names_in_fit, value = TRUE, ignore.case = TRUE))
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
fit_base <- readRDS(file.path(stage1_dir, "fit_base.rds"))
# Error:
# ! object 'stage1_dir' not found
#     ▆
#  1. ├─base::readRDS(file.path(stage1_dir, "fit_base.rds"))
#  2. └─base::file.path(stage1_dir, "fit_base.rds")
stage1_dir <- file.path(out_dir_v2, "stage1_smoke_scn09_ds01")
# Error:
# ! object 'out_dir_v2' not found
#     ▆
#  1. └─base::file.path(out_dir_v2, "stage1_smoke_scn09_ds01")
out_dir <- "simulated_virtual_dataset"
true_params <- readRDS(file.path(out_dir, "true_params_long.rds"))
out_dir_v2 <- "simulated_virtual_dataset_eta_filtered"
stage1_dir <- file.path(out_dir_v2, "stage1_smoke_scn09_ds01")
fit_base <- readRDS(file.path(stage1_dir, "fit_base.rds"))
fit_base_cov <- readRDS(file.path(stage1_dir, "fit_base_cov.rds"))
fit_base
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#            OBJF       AIC       BIC Log-likelihood
# FOCEi -16722.83 -13398.65 -13354.69       6707.326
# 
# ── Time (sec fit_base$time): ──
# 
#            setup optimize preprocess postprocess   other
# elapsed 0.042359 39.83114       0.11        0.02 13.7365
# 
# ── Population Parameters (fit_base$parFixed or fit_base$parFixedDf): ──
# 
#             Est. Back-transformed BSV(CV%) Shrink(SD)%
# lTVCL    -0.4121           0.6623     37.5      1.28% 
# lTVQ      0.6165            1.852                     
# lTVVc      3.058            21.29     31.5      12.9% 
# lTVVp      4.382            80.01                     
# lTVKA    -0.3567          -0.3567                     
# prop.err  0.1008           0.1008                     
#  
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.196   
#  
# 
#   Full BSV covariance (fit_base$omega) or correlation (fit_base$omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in fit_base$shrink 
#   Information about run found (fit_base$runInfo):
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring (fit_base$censInformation): No censoring
#   Minimization message (fit_base$message):  
#     Normal exit from bobyqa 
fit_base_cov
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#            OBJF       AIC       BIC Log-likelihood Condition#(Cov) Condition#(Cor)
# FOCEi -16722.83 -13398.65 -13354.69       6707.326        21.07669        2.116988
# 
# ── Time (sec fit_base_cov$time): ──
# 
#             setup optimize covariance preprocess postprocess table     other
# elapsed 0.0382303 38.70701   19.54924       0.09        0.02  0.15 0.3155251
# 
# ── Population Parameters (fit_base_cov$parFixed or fit_base_cov$parFixedDf): ──
# 
#             Est.       SE   %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
# lTVCL    -0.4121  0.02107  5.112 0.6623 (0.6355, 0.6902)     37.5      1.28% 
# lTVQ      0.6165 0.007683  1.246    1.852 (1.825, 1.881)                     
# lTVVc      3.058  0.02197 0.7183    21.29 (20.39, 22.22)     31.5      12.9% 
# lTVVp      4.382  0.00561  0.128     80.01 (79.14, 80.9)                     
# lTVKA    -0.3567    FIXED  FIXED                 -0.3567                     
# prop.err  0.1008                                  0.1008                     
#  
#   Covariance Type (fit_base_cov$covMethod): r,s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.196   
#  
# 
#   Full BSV covariance (fit_base_cov$omega) or correlation (fit_base_cov$omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in fit_base_cov$shrink 
#   Information about run found (fit_base_cov$runInfo):
#    • gradient problems with covariance; see $scaleInfo 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring (fit_base_cov$censInformation): No censoring
#   Minimization message (fit_base_cov$message):  
#     Normal exit from bobyqa 
# 
# ── Fit Data (object fit_base_cov is a modified tibble): ──
# # A tibble: 1,800 × 26
#   ID     TIME    DV  PRED     RES   WRES IPRED    IRES  IWRES CPRED   CRES CWRES eta.cl eta.vc     depot central
#   <fct> <dbl> <dbl> <dbl>   <dbl>  <dbl> <dbl>   <dbl>  <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>     <dbl>   <dbl>
# 1 1      0     0     0     0       0      0     0       0      0     0      0    -0.532  0.675 100           0  
# 2 1      7.06  1.71  2.51 -0.801  -1.79   1.79 -0.0789 -0.438  2.61 -0.898 -2.03 -0.532  0.675   0.712      74.7
# 3 1     14.1   1.28  1.34 -0.0565 -0.260  1.33 -0.0501 -0.373  1.64 -0.360 -1.52 -0.532  0.675   0.00507    55.6
# # ℹ 1,797 more rows
# # ℹ 10 more variables: peripheral <dbl>, cl <dbl>, vc <dbl>, q <dbl>, vp <dbl>, k10 <dbl>, k12 <dbl>, k21 <dbl>,
# #   tad <dbl>, dosenum <dbl>
# # ℹ Use `print(n = ...)` to see more rows
scm_focei <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = ""         # SCM doesn't need cov matrix for LRT
)
scm_focei_maxiteration <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  covMethod  = ""         # SCM doesn't need cov matrix for LRT
)
final_focei <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = TRUE,
  covMethod  = "r,s"
)
ds01 <- readRDS(file.path(stage1_dir, "nm_scn09_ds01.rds"))
candidate_pairs_test <- list(
  list(var = "cl", covar = "BW",   shapes = "power"),
  list(var = "cl", covar = "CrCL", shapes = "power"),
 list(var = "vc", covar = "BW",   shapes = "power"),
  list(var = "vc", covar = "CrCL", shapes = "power")
)
runSCM_traced <- function(label, ...) {
  dots      <- list(...)
  n_pairs   <- if (!is.null(dots$pairsVec))   length(dots$pairsVec) else NA
  n_workers <- if (!is.null(dots$workers))    dots$workers          else 1L
  search    <- if (!is.null(dots$searchType)) dots$searchType       else "scm"

  cat(sprintf(
    "\n=== [%s] %s starting %-8s | %d candidate(s), %d worker(s) ===\n",
    label, format(Sys.time(), "%H:%M:%S"), search, n_pairs, n_workers
  ), file = stderr())

  t0  <- Sys.time()
  res <- withCallingHandlers(
    do.call(nlmixr2scm::runSCM, dots),
    message = function(m) {
      txt <- conditionMessage(m)
      ## Tag only step/candidate boundary messages; let cli rules pass through.
      if (grepl("step\\s+\\d+,\\s*candidate|forward search|backward search",
                txt, ignore.case = TRUE)) {
        elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("[%s | +%5.1f min] ",
                    format(Sys.time(), "%H:%M:%S"), elapsed_min),
            file = stderr())
      }
      ## Don't muffle -- let the original cli message print as usual.
    }
  )
  elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  cat(sprintf(
    "=== [%s] %s DONE | elapsed %.1f min (%.0f s) ===\n\n",
    label, format(Sys.time(), "%H:%M:%S"),
    elapsed_s / 60, elapsed_s
  ), file = stderr())

  attr(res, "elapsed_s") <- elapsed_s
  res
}
out_dir_v2   <- "simulated_virtual_dataset_eta_filtered"
stage1_dir16 <- file.path(out_dir_v2, "stage1_smoke_scn16_ds01")
if (!dir.exists(stage1_dir16)) dir.create(stage1_dir16, recursive = TRUE)
scm_focei_n <- nlmixr2est::foceiControl(
  sigdig             = 4,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = FALSE,    # SCM doesn't need IPRED / CWRES tables
  covMethod          = "",       # SCM doesn't need cov matrix for LRT
  maxOuterIterations = 2000,
  maxInnerIterations = 2000
)
sim_obs_scn16 <- readRDS(file.path(out_dir_v2, "sim_obs_scenario_16.rds"))
ds16_01 <- to_nm_dataset(sim_obs_scn16) %>%
  dplyr::filter(DATASET == 1) %>%
  dplyr::select(-SCENARIO, -DATASET) %>%
  dplyr::mutate(
    ID   = as.integer(ID),
    SEX  = as.integer(SEX),
    RACE = as.integer(RACE)
  )
# Error in `to_nm_dataset()`:
# ! could not find function "to_nm_dataset"
#     ▆
#  1. ├─... %>% ...
#  2. ├─dplyr::mutate(...)
#  3. ├─dplyr::select(., -SCENARIO, -DATASET)
#  4. └─dplyr::filter(., DATASET == 1)
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
ds16_01 <- to_nm_dataset(sim_obs_scn16) %>%
  dplyr::filter(DATASET == 1) %>%
  dplyr::select(-SCENARIO, -DATASET) %>%
  dplyr::mutate(
    ID   = as.integer(ID),
    SEX  = as.integer(SEX),
    RACE = as.integer(RACE)
  )
# Error in `dplyr::transmute()`:
# ℹ In argument: `AMT = DOSE_MG`.
# Caused by error:
# ! object 'DOSE_MG' not found
#      ▆
#   1. ├─... %>% ...
#   2. ├─dplyr::mutate(...)
#   3. ├─dplyr::select(., -SCENARIO, -DATASET)
#   4. ├─dplyr::filter(., DATASET == 1)
#   5. ├─global to_nm_dataset(sim_obs_scn16)
#   6. │ └─... %>% ...
#   7. ├─dplyr::transmute(...)
#   8. ├─dplyr:::transmute.data.frame(...)
#   9. │ └─dplyr:::mutate_cols(.data, dots, by)
#  10. │   ├─base::withCallingHandlers(...)
#  11. │   └─dplyr:::mutate_col(dots[[i]], data, mask, new_columns)
#  12. │     └─mask$eval_all_mutate(quo)
#  13. │       └─dplyr (local) eval()
#  14. └─base::.handleSimpleError(...)
#  15.   └─dplyr (local) h(simpleError(msg, call))
#  16.     └─rlang::abort(message, class = error_class, parent = parent, call = error_call)
DOSE_MG <- 100
TVCL <- 0.6;  TVQ <- 1.8;  TVVc <- 20;  TVVp <- 80
k10_typ <- TVCL / TVVc
k12_typ <- TVQ  / TVVc
k21_typ <- TVQ  / TVVp
sum_k    <- k10_typ + k12_typ + k21_typ
beta_typ <- 0.5 * (sum_k - sqrt(sum_k^2 - 4 * k10_typ * k21_typ))
t_half_typ <- log(2) / beta_typ
cat(sprintf("Typical terminal half-life: %.2f h\n", t_half_typ)) #Typical terminal half-life: 141.29 h

hl_mult      <- c(0, 0.05, 0.1, 0.5, 1, 3)
sample_times <- hl_mult * t_half_typ
# Typical terminal half-life: 141.29 h
ds16_01 <- to_nm_dataset(sim_obs_scn16) %>%
  dplyr::filter(DATASET == 1) %>%
  dplyr::select(-SCENARIO, -DATASET) %>%
  dplyr::mutate(
    ID   = as.integer(ID),
    SEX  = as.integer(SEX),
    RACE = as.integer(RACE)
  )
saveRDS(ds16_01, file.path(stage1_dir16, "nm_scn16_ds01.rds"))
ds16_01 <- readRDS(file.path(stage1_dir16, "nm_scn16_ds01.rds"))
library(tidyverse)
ds16_01_60 <-  ds16_01 %>% filter(ID <=60)
saveRDS(ds16_01_60, file.path(stage1_dir16, "nm_scn16_ds01_60.rds"))
true_2cmt_scn16_refexp <- function() {
  ini({
    lTVCL      <- log(0.6)
    lTVQ       <- log(1.8)
    lTVVc      <- log(20)
    lTVVp      <- log(80)
    lTVKA      <- fix(log(0.7))    # KA unidentifiable
    TH_BW_CL   <- 0.75
    TH_CRCL_CL <- 0.5
    TH_BW_VC   <- 1.0
    TH_SEX_VC  <- 0.5

    eta.cl + eta.vc ~ c(0.1, 0.02, 0.1)
    prop.err <- 0.1
  })
  model({
    cl_typ <- exp(lTVCL) * (BW / 70)^TH_BW_CL * (CrCL / 95)^TH_CRCL_CL
    vc_typ <- exp(lTVVc) * (BW / 70)^TH_BW_VC * (1 + TH_SEX_VC * SEX)
    cl     <- cl_typ * exp(eta.cl)
    vc     <- vc_typ * exp(eta.vc)
    q      <- exp(lTVQ)
    vp     <- exp(lTVVp)
    ka     <- exp(lTVKA)

    k10 <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)      = -ka * depot
    d/dt(central)    =  ka * depot - k10 * central - k12 * central + k21 * peripheral
    d/dt(peripheral) =  k12 * central - k21 * peripheral

    cp = central / vc
    cp ~ prop(prop.err)
  })
}
true_2cmt_scn16_lin <- function() {
  ini({
    lTVCL      <- log(0.6)
    lTVQ       <- log(1.8)
    lTVVc      <- log(20)
    lTVVp      <- log(80)
    lTVKA      <- fix(log(0.7))
    TH_BW_CL   <- 0.75
    TH_CRCL_CL <- 0.5
    TH_BW_VC   <- 1.0
    TH_SEX_VC  <- 0.5

    eta.cl + eta.vc ~ c(0.1, 0.02, 0.1)
    prop.err <- 0.1
  })
  model({
    ## BW and CrCL via log() => power on natural scale (matches refexp).
    ## SEX preserves the linear-on-natural-scale parameterisation
    ## (1 + theta*SEX) used in the simulation; log() lifts it to log scale.
    lTVCL_typ <- lTVCL + TH_BW_CL * log(BW / 70) + TH_CRCL_CL * log(CrCL / 95)
    lTVVc_typ <- lTVVc + TH_BW_VC * log(BW / 70) + log(1 + TH_SEX_VC * SEX)
    cl        <- exp(lTVCL_typ + eta.cl)
    vc        <- exp(lTVVc_typ + eta.vc)
    q         <- exp(lTVQ)
    vp        <- exp(lTVVp)
    ka        <- exp(lTVKA)

    k10 <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)      = -ka * depot
    d/dt(central)    =  ka * depot - k10 * central - k12 * central + k21 * peripheral
    d/dt(peripheral) =  k12 * central - k21 * peripheral

    cp = central / vc
    cp ~ prop(prop.err)
  })
}
true_2cmt_scn16_lin <- function() {
  ini({
    lTVCL      <- log(0.6)
    lTVQ       <- log(1.8)
    lTVVc      <- log(20)
    lTVVp      <- log(80)
    lTVKA      <- fix(log(0.7))
    TH_BW_CL   <- 0.75
    TH_CRCL_CL <- 0.5
    TH_BW_VC   <- 1.0
    TH_SEX_VC  <- 0.5

    eta.cl + eta.vc ~ c(0.1, 0.02, 0.1)
    prop.err <- 0.1
  })
  model({
    ## BW and CrCL via log() => power on natural scale (matches refexp).
    ## SEX preserves the linear-on-natural-scale parameterisation
    ## (1 + theta*SEX) used in the simulation; log() lifts it to log scale.
    lTVCL_typ <- lTVCL + TH_BW_CL * log(BW / 70) + TH_CRCL_CL * log(CrCL / 95)
    lTVVc_typ <- lTVVc + TH_BW_VC * log(BW / 70) + log(TH_SEX_VC * SEX) #align with nlmixr2csm
    cl        <- exp(lTVCL_typ + eta.cl)
    vc        <- exp(lTVVc_typ + eta.vc)
    q         <- exp(lTVQ)
    vp        <- exp(lTVVp)
    ka        <- exp(lTVKA)

    k10 <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)      = -ka * depot
    d/dt(central)    =  ka * depot - k10 * central - k12 * central + k21 * peripheral
    d/dt(peripheral) =  k12 * central - k21 * peripheral

    cp = central / vc
    cp ~ prop(prop.err)
  })
}
t_fit_true16_refexp <- system.time(
  fit_true16_refexp <- nlmixr2(true_2cmt_scn16_refexp, ds16_01,
                               est = "focei", control = scm_focei_n)
)
# ℹ parameter labels from comments will be replaced by 'label()'
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities                                                        
# → calculate ∂(f)/∂(η)                                                            
# → calculate ∂(R²)/∂(η)                                                           
# → finding duplicate expressions in inner model...                                
# → optimizing duplicate expressions in inner model...                             
# → finding duplicate expressions in EBE model...                                  
# → optimizing duplicate expressions in EBE model...                               
# → compiling inner model...                                                       
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...                                
# → compiling EBE model...                                                         
# ✔ done
# → compiling events FD model...
# ✔ done
# 
# 
# Timing stopped at: 178.8 4.17 66.57
scm_focei_n <- nlmixr2est::foceiControl(
  sigdig             = 4,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = FALSE,    # SCM doesn't need IPRED / CWRES tables
  covMethod          = "r,s",       # SCM doesn't need cov matrix for LRT
  maxOuterIterations = 2000,
  maxInnerIterations = 2000
)
t_fit_true16_refexp <- system.time(
  fit_true16_refexp <- nlmixr2(true_2cmt_scn16_refexp, ds16_01,
                               est = "focei", control = scm_focei_n)
)
# ℹ parameter labels from comments will be replaced by 'label()'
# calculating covariance matrix
# [====|====|====|====|====|====|====|====|====|====] 100%; 0:01:45 done
saveRDS(fit_true16_refexp,
        file.path(stage1_dir16, "fit_true_scn16_ds01_refexp.rds"))
t_fit_true16_lin <- system.time(
  fit_true16_lin <- nlmixr2(true_2cmt_scn16_lin, ds16_01,
                            est = "focei", control = scm_focei_n)
)
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# [====|====|====|====|====|====|====|====|====|====] 100%; 0:00:00 → calculate sensitivities
# → calculate ∂(f)/∂(η)                                                            
# → calculate ∂(R²)/∂(η)                                                           
# → finding duplicate expressions in inner model...                                
# → optimizing duplicate expressions in inner model...                             
# → finding duplicate expressions in EBE model...                                  
# → optimizing duplicate expressions in EBE model...                               
# → compiling inner model...                                                       
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...                                
# → compiling EBE model...                                                         
# ✔ done
# → compiling events FD model...
# ✔ done
# calculating covariance matrix
# S matrix calculation failed; Switch to R-matrix covariance.:00:10 
# [====|====|====|====|====|====|====|====|====|====] 100%; 0:00:13 done
saveRDS(fit_true16_lin,
        file.path(stage1_dir16, "fit_true_scn16_ds01_lin.rds"))
diag_true16_refexp <- diagnose_fit(fit_true16_refexp)
diag_true16_lin    <- diagnose_fit(fit_true16_lin)
est_true16_refexp <- extract_params_long(fit_true16_refexp)
est_true16_lin    <- extract_params_long(fit_true16_lin)
err_true16_refexp <- rel_err_one(est_true16_refexp, true_params, scenario_id = 16)
err_true16_lin    <- rel_err_one(est_true16_lin,    true_params, scenario_id = 16)
err_true16_compare <- dplyr::bind_rows(
  err_true16_refexp %>% dplyr::mutate(parameterisation = "refexp"),
  err_true16_lin    %>% dplyr::mutate(parameterisation = "lin")
) %>%
  dplyr::select(parameterisation, parameter, true_value,
                estimate, abs_err, rel_err, rel_err_pct)
err_true16_wide <- err_true16_compare %>%
  tidyr::pivot_wider(
    id_cols     = c(parameter, true_value),
    names_from  = parameterisation,
    values_from = c(estimate, rel_err_pct),
    names_glue  = "{.value}_{parameterisation}"
  ) %>%
  dplyr::mutate(
    abs_diff_estimate    = abs(estimate_refexp - estimate_lin),
    abs_diff_rel_err_pct = abs(rel_err_pct_refexp - rel_err_pct_lin)
  )
saveRDS(err_true16_wide,
        file.path(stage1_dir16, "rel_err_true_scn16_ds01_wide.rds"))
.build_part1_row <- function(label, diag_x, t_x) {
  tibble::tibble(
    step             = "true_model_fit",
    parameterisation = label,
    scenario         = 9,
    dataset          = 1,
    converged        = diag_x$converged,
    objf             = diag_x$objf,
    cov_step_ok      = diag_x$cov_ok,
    cond_num         = diag_x$cond_num,
    cond_num_sqrt    = diag_x$cond_num_sqrt,  # NONMEM convention (< 1000 healthy)
    runtime_sec      = unname(t_x["elapsed"])
  )
}
part1_summary16 <- dplyr::bind_rows(
  .build_part1_row("refexp", diag_true16_refexp, t_fit_true16_refexp),
  .build_part1_row("lin",    diag_true16_lin,    t_fit_true16_lin)
) %>%
  dplyr::mutate(scenario = 16)
true_2cmt_scn16_lin <- function() {
  ini({
    lTVCL      <- log(0.6)
    lTVQ       <- log(1.8)
    lTVVc      <- log(20)
    lTVVp      <- log(80)
    lTVKA      <- fix(log(0.7))
    TH_BW_CL   <- 0.75
    TH_CRCL_CL <- 0.5
    TH_BW_VC   <- 1.0
    TH_SEX_VC  <- 0.5

    eta.cl + eta.vc ~ c(0.1, 0.02, 0.1)
    prop.err <- 0.1
  })
  model({
    ## BW and CrCL via log() => power on natural scale (matches refexp).
    ## SEX preserves the linear-on-natural-scale parameterisation
    ## (1 + theta*SEX) used in the simulation; log() lifts it to log scale.
    lTVCL_typ <- lTVCL + TH_BW_CL * log(BW / 70) + TH_CRCL_CL * log(CrCL / 95)
    lTVVc_typ <- lTVVc + TH_BW_VC * log(BW / 70) + log(1 + TH_SEX_VC * SEX) #align with nlmixr2csm
    cl        <- exp(lTVCL_typ + eta.cl)
    vc        <- exp(lTVVc_typ + eta.vc)
    q         <- exp(lTVQ)
    vp        <- exp(lTVVp)
    ka        <- exp(lTVKA)

    k10 <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)      = -ka * depot
    d/dt(central)    =  ka * depot - k10 * central - k12 * central + k21 * peripheral
    d/dt(peripheral) =  k12 * central - k21 * peripheral

    cp = central / vc
    cp ~ prop(prop.err)
  })
}
t_fit_true16_lin <- system.time(
  fit_true16_lin <- nlmixr2(true_2cmt_scn16_lin, ds16_01,
                            est = "focei", control = scm_focei_n)
)
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# [====|====|====|====|====|====|====|====|====|====] 100%; 0:00:00 → calculate sensitivities
# → calculate ∂(f)/∂(η)                                                            
# → calculate ∂(R²)/∂(η)                                                           
# → finding duplicate expressions in inner model...                                
# → optimizing duplicate expressions in inner model...                             
# → finding duplicate expressions in EBE model...                                  
# → optimizing duplicate expressions in EBE model...                               
# → compiling inner model...                                                       
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...                                
# → compiling EBE model...                                                         
# ✔ done
# → compiling events FD model...
# ✔ done
# calculating covariance matrix
# [====|====|====|====|====|====|====|====|====|====] 100%; 0:01:08 done
saveRDS(fit_true16_lin,
        file.path(stage1_dir16, "fit_true_scn16_ds01_lin.rds"))
diag_true16_lin    <- diagnose_fit(fit_true16_lin)
est_true16_lin    <- extract_params_long(fit_true16_lin)
err_true16_lin    <- rel_err_one(est_true16_lin,    true_params, scenario_id = 16)
err_true16_compare <- dplyr::bind_rows(
  err_true16_refexp %>% dplyr::mutate(parameterisation = "refexp"),
  err_true16_lin    %>% dplyr::mutate(parameterisation = "lin")
) %>%
  dplyr::select(parameterisation, parameter, true_value,
                estimate, abs_err, rel_err, rel_err_pct)
err_true16_wide <- err_true16_compare %>%
  tidyr::pivot_wider(
    id_cols     = c(parameter, true_value),
    names_from  = parameterisation,
    values_from = c(estimate, rel_err_pct),
    names_glue  = "{.value}_{parameterisation}"
  ) %>%
  dplyr::mutate(
    abs_diff_estimate    = abs(estimate_refexp - estimate_lin),
    abs_diff_rel_err_pct = abs(rel_err_pct_refexp - rel_err_pct_lin)
  )
diag_summary16 <- dplyr::bind_rows(
  .build_part1_row("refexp", diag_true16_refexp, t_fit_true16_refexp),
  .build_part1_row("lin",    diag_true16_lin,    t_fit_true16_lin)
) %>%
  dplyr::mutate(scenario = 16)
fit_true16_refexp
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#            OBJF       AIC       BIC Log-likelihood Condition#(Cov)
# FOCEi -16641.59 -13309.41 -13243.47       6666.707        282.7779
#       Condition#(Cor)
# FOCEi        3.098302
# 
# ── Time (sec fit_true16_refexp$time): ──
# 
#             setup optimize covariance preprocess postprocess     other
# elapsed 0.0591244 82.12501   105.0896       0.13        0.05 0.3563173
# 
# ── Population Parameters (fit_true16_refexp$parFixed or $parFixedDf): ──
# 
#                    Parameter    Est.       SE    %RSE Back-transformed(95%CI)
# lTVCL                        -0.5123   0.0169   3.299 0.5991 (0.5796, 0.6193)
# lTVQ                          0.6088 0.007157   1.175    1.838 (1.813, 1.864)
# lTVVc                          3.043  0.02257  0.7417    20.97 (20.06, 21.92)
# lTVVp                          4.381 0.004272 0.09751    79.91 (79.24, 80.58)
# lTVKA      KA unidentifiable -0.3567    FIXED   FIXED                 -0.3567
# TH_BW_CL                      0.7059  0.06005   8.507 0.7059 (0.5882, 0.8236)
# TH_CRCL_CL                    0.5226  0.03922   7.505 0.5226 (0.4457, 0.5994)
# TH_BW_VC                       0.899  0.06614   7.356   0.899 (0.7694, 1.029)
# TH_SEX_VC                     0.4666  0.03516   7.535 0.4666 (0.3977, 0.5355)
# prop.err                     0.09906                                  0.09906
# eta.cl                                                                       
# eta.vc                                                                       
#            BSV(CV%) Shrink(SD)%
# lTVCL                          
# lTVQ                           
# lTVVc                          
# lTVVp                          
# lTVKA                          
# TH_BW_CL                       
# TH_CRCL_CL                     
# TH_BW_VC                       
# TH_SEX_VC                      
# prop.err                       
# eta.cl         30.9      1.50% 
# eta.vc         32.8      11.1% 
#  
#   Covariance Type (fit_true16_refexp$covMethod): s
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.151   
#  
# 
#   Full BSV covariance (fit_true16_refexp$omega) 
#     or correlation (fit_true16_refexp$omegaR; diagonals=SDs)
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found (fit_true16_refexp$runInfo):
#    • tolerances (atol/rtol) were increased (after 4 bad solves) for some difficult ODE solving during the optimization. can control with foceiControl(stickyRecalcN=) consider increasing sigdig/atol/rtol changing initial estimates or changing the structural model 
#    • gradient problems with covariance; see $scaleInfo 
#    • using S matrix to calculate covariance, can check sandwich or R matrix with $covRS and $covR 
#    • S matrix had problems solving for some subject and parameters 
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#    • Hessian reset during optimization; (can control by foceiControl(resetHessianAndEta=.)) 
#   Censoring (fit_true16_refexp$censInformation): No censoring
#   Minimization message (fit_true16_refexp$message):  
#     Normal exit from bobyqa 

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

# ---- 1. Side-by-side fit summary -------------------------------------------
cat("OFV refexp = ", round(fit_true16_refexp$objf, 3), "\n")
cat("OFV lin    = ", round(fit_true16_lin$objf,    3), "\n")
cat("dOFV (refexp - lin) = ",
    round(fit_true16_refexp$objf - fit_true16_lin$objf, 3),
    "  (positive => lin found a better minimum)\n\n")

cat("Condition #(Cov)  refexp: ", round(diag_true16_refexp$cond_num, 2), "\n")
cat("Condition #(Cov)  lin   : ", round(diag_true16_lin$cond_num,    2), "\n\n")

# ---- 2. Compare population parameter estimates -----------------------------
pop_refexp <- fit_true16_refexp$parFixedDf %>%
  tibble::rownames_to_column("param") %>%
  dplyr::select(param, est_refexp = Estimate, se_refexp = SE)
pop_lin <- fit_true16_lin$parFixedDf %>%
  tibble::rownames_to_column("param") %>%
  dplyr::select(param, est_lin = Estimate, se_lin = SE)

pop_compare <- dplyr::full_join(pop_refexp, pop_lin, by = "param") %>%
  dplyr::mutate(
    abs_diff = abs(est_refexp - est_lin),
    rel_diff_pct = 100 * abs_diff / pmax(abs(est_refexp), abs(est_lin), 1e-9)
  )
print(pop_compare)
# OFV refexp =  -16641.59 
# OFV lin    =  -16900.84 
# dOFV (refexp - lin) =  259.252   (positive => lin found a better minimum)
# 
# Condition #(Cov)  refexp:  282.78 
# Condition #(Cov)  lin   :  498.63 
# 
#         param  est_refexp   se_refexp     est_lin     se_lin     abs_diff rel_diff_pct
# 1       lTVCL -0.51226955 0.016899117 -0.50564713 0.04243573 0.0066224224   1.29276126
# 2        lTVQ  0.60882809 0.007156739  0.60910723 0.01832525 0.0002791385   0.04582748
# 3       lTVVc  3.04311974 0.022569654  3.05035569 0.05870790 0.0072359476   0.23721652
# 4       lTVVp  4.38084905 0.004271728  4.37993417 0.01290470 0.0009148739   0.02088348
# 5       lTVKA -0.35667494          NA -0.35667494         NA 0.0000000000   0.00000000
# 6    TH_BW_CL  0.70589293 0.060051560  0.65759512 0.17093132 0.0482978102   6.84208724
# 7  TH_CRCL_CL  0.52255265 0.039219570  0.53936531 0.10885736 0.0168126668   3.11712052
# 8    TH_BW_VC  0.89904766 0.066136823  0.88579543 0.20720128 0.0132522350   1.47403031
# 9   TH_SEX_VC  0.46663184 0.035161187  0.45571210 0.12710200 0.0109197356   2.34011799
# 10   prop.err  0.09906355          NA  0.09892471         NA 0.0001388440   0.14015654
# 11     eta.cl          NA          NA          NA         NA           NA           NA
# 12     eta.vc          NA          NA          NA         NA           NA           NA

# ---- 3. Compare omega (BSV) and predicted cp at a reference subject -------
cat("--- Omega (BSV cov matrix) ---\n")
cat("refexp:\n"); print(fit_true16_refexp$omega)
cat("\nlin:\n");    print(fit_true16_lin$omega)

# ---- 4. Inner-loop ETA distribution: do both fits use the same etas? ------
cat("\n--- ETA summary (per-subject EBE) ---\n")
eta_re <- fit_true16_refexp$eta %>% dplyr::select(dplyr::starts_with("eta"))
eta_li <- fit_true16_lin$eta    %>% dplyr::select(dplyr::starts_with("eta"))

cat("refexp eta.cl: mean=", round(mean(eta_re$eta.cl), 4),
    " sd=", round(sd(eta_re$eta.cl), 4), "\n")
cat("lin    eta.cl: mean=", round(mean(eta_li$eta.cl), 4),
    " sd=", round(sd(eta_li$eta.cl), 4), "\n")

cat("refexp eta.vc: mean=", round(mean(eta_re$eta.vc), 4),
    " sd=", round(sd(eta_re$eta.vc), 4), "\n")
cat("lin    eta.vc: mean=", round(mean(eta_li$eta.vc), 4),
    " sd=", round(sd(eta_li$eta.vc), 4), "\n")

# ---- 5. Look at run-info / messages ---------------------------------------
cat("\n--- $runInfo (refexp) ---\n");  print(fit_true16_refexp$runInfo)
cat("\n--- $runInfo (lin) ---\n");      print(fit_true16_lin$runInfo)
# --- Omega (BSV cov matrix) ---
# refexp:
#            eta.cl     eta.vc
# eta.cl 0.09142213 0.01455074
# eta.vc 0.01455074 0.10203501
# 
# lin:
#            eta.cl     eta.vc
# eta.cl 0.09038468 0.01582280
# eta.vc 0.01582280 0.09830706
# 
# --- ETA summary (per-subject EBE) ---
# refexp eta.cl: mean= 0.0019  sd= 0.2978 
# lin    eta.cl: mean= 0.0012  sd= 0.2975 
# refexp eta.vc: mean= 0.0089  sd= 0.2841 
# lin    eta.vc: mean= 0.0081  sd= 0.2817 
# 
# --- $runInfo (refexp) ---
# [1] "tolerances (atol/rtol) were increased (after 4 bad solves) for some difficult ODE solving during the optimization.\ncan control with foceiControl(stickyRecalcN=)\nconsider increasing sigdig/atol/rtol changing initial estimates or changing the structural model"
# [2] "gradient problems with covariance; see $scaleInfo"                                                                                                                                                                                                                  
# [3] "using S matrix to calculate covariance, can check sandwich or R matrix with $covRS and $covR"                                                                                                                                                                       
# [4] "S matrix had problems solving for some subject and parameters"                                                                                                                                                                                                      
# [5] "last objective function was not at minimum, possible problems in optimization"                                                                                                                                                                                      
# [6] "ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.))"                                                                                                                                                                            
# [7] "initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=))"                                                                                                                                                                                    
# [8] "Hessian reset during optimization; (can control by foceiControl(resetHessianAndEta=.))"                                                                                                                                                                             
# 
# --- $runInfo (lin) ---
# [1] "gradient problems with covariance; see $scaleInfo"                                           
# [2] "using S matrix to calculate covariance, can check sandwich or R matrix with $covRS and $covR"
# [3] "last objective function was not at minimum, possible problems in optimization"               
# [4] "ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.))"     
# [5] "initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=))"             

# ---- 6. Algebraic identity check via direct rxode2 prediction --------------
# Take ONE row of fitted predictions; verify cp(refexp) == cp(lin) at the
# *same* (theta, eta).  If they match to machine precision, the models are
# pointwise identical, and the OFV gap is purely a FOCEi inner-loop artifact.

# Grab refexp's EBE for subject 1 + a few observations
ipred_re <- as.data.frame(fit_true16_refexp)
ipred_li <- as.data.frame(fit_true16_lin)

# Same column structure?
cat("refexp cols: ", paste(head(names(ipred_re), 10), collapse=", "), "\n")
cat("lin    cols: ", paste(head(names(ipred_li), 10), collapse=", "), "\n\n")

# Show a few observed rows side by side -- ipred should agree if and only
# if the EBE etas agree, since the predictive equations are identical.
s1_re <- ipred_re %>% dplyr::filter(ID == ID[1]) %>%
  dplyr::select(ID, TIME, DV, IPRED, PRED) %>% head(6)
s1_li <- ipred_li %>% dplyr::filter(ID == ID[1]) %>%
  dplyr::select(ID, TIME, DV, IPRED, PRED) %>% head(6)

cat("--- Subject 1, first 6 observations ---\n")
cat("refexp:\n");  print(s1_re)
cat("\nlin:\n");    print(s1_li)

# Pointwise diff in PRED (population, eta=0) -- this is the cleanest check
# of model algebraic identity:
both <- dplyr::inner_join(
  ipred_re %>% dplyr::select(ID, TIME, PRED_re = PRED, IPRED_re = IPRED),
  ipred_li %>% dplyr::select(ID, TIME, PRED_li = PRED, IPRED_li = IPRED),
  by = c("ID", "TIME")
)
cat("\n--- |PRED_refexp - PRED_lin| stats (eta=0, should be ~0 if algebra matches) ---\n")
cat("max(|PRED diff|):  ", format(max(abs(both$PRED_re - both$PRED_li)), digits=4), "\n")
cat("mean(|PRED diff|): ", format(mean(abs(both$PRED_re - both$PRED_li)), digits=4), "\n")
cat("\n--- |IPRED_refexp - IPRED_lin| stats (with EBE etas, will differ if etas differ) ---\n")
cat("max(|IPRED diff|): ", format(max(abs(both$IPRED_re - both$IPRED_li)), digits=4), "\n")
cat("mean(|IPRED diff|):", format(mean(abs(both$IPRED_re - both$IPRED_li)), digits=4), "\n")
# Error in `as.data.frame.default()`:
# ! cannot coerce class ‘"nlmixr2FitCore"’ to a data.frame
#     ▆
#  1. ├─base::as.data.frame(fit_true16_refexp)
#  2. └─base::as.data.frame.default(fit_true16_refexp)

# ---- 6a. Direct algebraic identity check on the cl/vc formulas -------------
# Use the fitted refexp params and evaluate cl, vc from BOTH algebraic forms
# at the same (BW, CrCL, SEX, eta).  If the two forms are pointwise identical
# the differences must be exactly 0 (modulo floating-point roundoff ~1e-15).
th <- fit_true16_refexp$parFixedDf$Estimate
names(th) <- rownames(fit_true16_refexp$parFixedDf)

# Sample 10 random (BW, CrCL, SEX, eta) triplets from the dataset
set.seed(1)
samp <- ds16_01 %>%
  dplyr::filter(EVID == 0) %>%
  dplyr::slice_sample(n = 10) %>%
  dplyr::select(ID, BW, CrCL, SEX)
samp$eta.cl <- rnorm(10, 0, 0.3)
samp$eta.vc <- rnorm(10, 0, 0.3)

with(samp, {
  cl_refexp <- exp(th["lTVCL"]) * (BW / 70)^th["TH_BW_CL"] *
               (CrCL / 95)^th["TH_CRCL_CL"] * exp(eta.cl)
  cl_lin    <- exp(th["lTVCL"] + th["TH_BW_CL"] * log(BW / 70) +
                   th["TH_CRCL_CL"] * log(CrCL / 95) + eta.cl)

  vc_refexp <- exp(th["lTVVc"]) * (BW / 70)^th["TH_BW_VC"] *
               (1 + th["TH_SEX_VC"] * SEX) * exp(eta.vc)
  vc_lin    <- exp(th["lTVVc"] + th["TH_BW_VC"] * log(BW / 70) +
                   log(1 + th["TH_SEX_VC"] * SEX) + eta.vc)

  cat("max |cl_refexp - cl_lin| = ",
      format(max(abs(cl_refexp - cl_lin)), digits = 4), "\n")
  cat("max |vc_refexp - vc_lin| = ",
      format(max(abs(vc_refexp - vc_lin)), digits = 4), "\n")
})

# ---- 6b. Per-subject EBE etas: how different are they really? -------------
eta_re <- fit_true16_refexp$eta %>%
  dplyr::select(ID, eta.cl_re = eta.cl, eta.vc_re = eta.vc)
eta_li <- fit_true16_lin$eta %>%
  dplyr::select(ID, eta.cl_li = eta.cl, eta.vc_li = eta.vc)

eta_both <- dplyr::inner_join(eta_re, eta_li, by = "ID") %>%
  dplyr::mutate(
    d_cl = eta.cl_re - eta.cl_li,
    d_vc = eta.vc_re - eta.vc_li
  )
cat("\n--- Per-subject EBE eta differences (refexp - lin) ---\n")
cat("eta.cl  diff: max=", round(max(abs(eta_both$d_cl)), 4),
    "  mean=", round(mean(abs(eta_both$d_cl)), 4),
    "  cor(re, li)=", round(cor(eta_both$eta.cl_re, eta_both$eta.cl_li), 4), "\n")
cat("eta.vc  diff: max=", round(max(abs(eta_both$d_vc)), 4),
    "  mean=", round(mean(abs(eta_both$d_vc)), 4),
    "  cor(re, li)=", round(cor(eta_both$eta.vc_re, eta_both$eta.vc_li), 4), "\n")

cat("\nN subjects:", nrow(eta_both), "\n")
cat("If max|diff| were ~0, the inner loops converged to the same mode.\n")
# max |cl_refexp - cl_lin| =  2.22e-16 
# max |vc_refexp - vc_lin| =  1.421e-14 
# 
# --- Per-subject EBE eta differences (refexp - lin) ---
# eta.cl  diff: max= 0.027   mean= 0.0074   cor(re, li)= 0.9995 
# eta.vc  diff: max= 0.037   mean= 0.0047   cor(re, li)= 0.9998 
# 
# N subjects: 300 
# If max|diff| were ~0, the inner loops converged to the same mode.

# ---- Convergence sanity check for both scn16 fits ---------------------------
# nlmixr2 stores the iteration trace in $iter; the final OFV difference and
# gradient magnitudes tell us whether "not at minimum" is cosmetic or real.

cat("================================================================\n")
cat("Convergence diagnostics: refexp vs lin (scenario 16, ds01)\n")
cat("================================================================\n\n")

for (lbl in c("refexp", "lin")) {
  fit <- get(paste0("fit_true16_", lbl))
  cat("--- ", lbl, " ---\n", sep = "")

  # 1. How many outer iterations did BOBYQA do?
  cat("  Outer iterations:    ", fit$message %||% "<no message>", "\n")
  cat("  Final OFV:           ", round(fit$objf, 4), "\n")

  # 2. Best-seen vs final OFV  (the actual cause of the warning)
  if (!is.null(fit$iter) && nrow(fit$iter) > 0) {
    ofv_trace <- fit$iter$objf
    cat("  Best OFV seen:       ", round(min(ofv_trace, na.rm = TRUE), 4), "\n")
    cat("  Final - best OFV:    ", round(fit$objf - min(ofv_trace, na.rm = TRUE), 4),
        "  (cosmetic if < 0.1, mild concern 0.1-1, real concern > 1)\n")
    cat("  N outer iter:        ", nrow(fit$iter), "\n")
  }

  # 3. Final scaled gradients  (covariance struggle = large grads)
  if (!is.null(fit$scaleInfo)) {
    g <- fit$scaleInfo$gradient
    if (!is.null(g)) {
      cat("  Max |grad|:          ", format(max(abs(g)), digits = 3),
          "  (well-converged < 0.01, OK < 0.1, suspect > 1)\n")
      cat("  Median |grad|:       ", format(median(abs(g)), digits = 3), "\n")
    }
  }
  cat("\n")
}

# 4. The OFV gap between the two fits, expressed two ways
cat("--- Direct comparison ---\n")
cat("  OFV(refexp) - OFV(lin) = ",
    round(fit_true16_refexp$objf - fit_true16_lin$objf, 3),
    " (positive => lin converged deeper)\n")
cat("  Average OFV per subject (N=300):\n")
cat("    refexp: ", round(fit_true16_refexp$objf / 300, 4), "\n")
cat("    lin   : ", round(fit_true16_lin$objf    / 300, 4), "\n")
cat("    gap   : ", round((fit_true16_refexp$objf - fit_true16_lin$objf) / 300, 4),
    "  (per-subject log-likelihood difference)\n")
# ================================================================
# Convergence diagnostics: refexp vs lin (scenario 16, ds01)
# ================================================================
# 
# --- refexp ---
#   Outer iterations:     Normal exit from bobyqa 
#   Final OFV:            -16641.59 
# 
# --- lin ---
#   Outer iterations:     Normal exit from bobyqa 
#   Final OFV:            -16900.84 
# 
# --- Direct comparison ---
#   OFV(refexp) - OFV(lin) =  259.252  (positive => lin converged deeper)
#   Average OFV per subject (N=300):
#     refexp:  -55.472 
#     lin   :  -56.3361 
#     gap   :  0.8642   (per-subject log-likelihood difference)

# nlmixr2 stashes the iteration history under different slots depending on
# control settings; let's look at what's actually there.
for (lbl in c("refexp", "lin")) {
  fit <- get(paste0("fit_true16_", lbl))
  cat("--- ", lbl, " ---\n", sep = "")

  # The outer-loop trace is in $parHistData when print=0; that's our best bet
  if (!is.null(fit$parHistData) && nrow(fit$parHistData) > 0) {
    obj_col <- intersect(c("obj","objf","objective","OBJF"), names(fit$parHistData))
    if (length(obj_col)) {
      ofv_trace <- fit$parHistData[[obj_col[1]]]
      cat("  N outer iter recorded: ", length(ofv_trace), "\n")
      cat("  Min OFV in trace:      ", round(min(ofv_trace, na.rm = TRUE), 4), "\n")
      cat("  Final OFV in trace:    ", round(tail(ofv_trace[!is.na(ofv_trace)], 1), 4), "\n")
      cat("  Reported $objf:        ", round(fit$objf, 4), "\n")
      cat("  Final - best (trace):  ", round(tail(ofv_trace[!is.na(ofv_trace)], 1) -
                                              min(ofv_trace, na.rm = TRUE), 4), "\n")
    }
    cat("  parHistData cols: ", paste(names(fit$parHistData), collapse=", "), "\n")
  } else {
    cat("  No parHistData.\n")
  }
  cat("\n")
}
# --- refexp ---
#   N outer iter recorded:  507 
#   Min OFV in trace:       -16885.6 
#   Final OFV in trace:     -16641.58 
#   Reported $objf:         -16641.59 
#   Final - best (trace):   244.0222 
#   parHistData cols:  iter, type, objf, lTVCL, lTVQ, lTVVc, lTVVp, TH_BW_CL, TH_CRCL_CL, TH_BW_VC, TH_SEX_VC, prop.err, o1, o2, o3 
# 
# --- lin ---
#   N outer iter recorded:  624 
#   Min OFV in trace:       -16900.85 
#   Final OFV in trace:     -16900.84 
#   Reported $objf:         -16900.84 
#   Final - best (trace):   0.0068 
#   parHistData cols:  iter, type, objf, lTVCL, lTVQ, lTVVc, lTVVp, TH_BW_CL, TH_CRCL_CL, TH_BW_VC, TH_SEX_VC, prop.err, o1, o2, o3 
# 
focei_grid <- tibble::tribble(
  ~setting,           ~sigdig, ~atol,  ~rtol,  ~stickyRecalcN,
  "baseline",         4,       1e-8,   1e-6,   4,        # matches scm_focei_n
  "fast",             3,       1e-6,   1e-4,   4,        # speed-only
  "no_sticky",        4,       1e-8,   1e-6,   20,       # block auto-relax
  "tight_ode",        4,       1e-10,  1e-8,   4,        # depth via ODE
  "tight_no_sticky",  4,       1e-10,  1e-8,   20,       # belt + braces
  "fast_no_sticky",   3,       1e-6,   1e-4,   20        # fast but protected
)
.bench_one_focei <- function(model_fn, data, ctrl, model_name, setting_name) {
  t0 <- Sys.time()
  fit <- tryCatch(
    nlmixr2(model_fn, data, est = "focei", control = ctrl),
    error = function(e) list(.failed = TRUE, .err = conditionMessage(e))
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  if (isTRUE(fit$.failed)) {
    return(tibble::tibble(
      model = model_name, setting = setting_name,
      elapsed_min = elapsed / 60,
      final_ofv = NA_real_, best_ofv = NA_real_, ofv_gap = NA_real_,
      n_iter = NA_integer_, bad_solves = NA, hessian_reset = NA,
      not_at_min = NA, converged_well = FALSE, error = fit$.err
    ))
  }

  trace_objf <- if (!is.null(fit$parHistData) && nrow(fit$parHistData) > 0)
    fit$parHistData$objf else NA_real_
  best_ofv <- if (any(!is.na(trace_objf))) min(trace_objf, na.rm = TRUE) else NA_real_
  ri <- fit$runInfo %||% character(0)

  tibble::tibble(
    model = model_name, setting = setting_name,
    elapsed_min   = elapsed / 60,
    final_ofv     = fit$objf,
    best_ofv      = best_ofv,
    ofv_gap       = fit$objf - best_ofv,
    n_iter        = if (!is.null(fit$parHistData)) nrow(fit$parHistData) else NA_integer_,
    bad_solves    = any(grepl("bad solves",     ri, fixed = TRUE)),
    hessian_reset = any(grepl("Hessian reset",  ri, fixed = TRUE)),
    not_at_min    = any(grepl("not at minimum", ri, fixed = TRUE)),
    converged_well = !is.na(fit$objf) && !is.na(best_ofv) &&
                     abs(fit$objf - best_ofv) < 0.5,
    error = NA_character_
  )
}
bench_results <- purrr::pmap_dfr(focei_grid, function(setting, sigdig, atol, rtol, stickyRecalcN) {
  cat(sprintf("\n>>> %-18s sigdig=%d atol=%.0e rtol=%.0e sticky=%d\n",
              setting, sigdig, atol, rtol, stickyRecalcN))
  ctrl <- nlmixr2est::foceiControl(
    sigdig             = sigdig,
    outerOpt           = "bobyqa",
    print              = 0,
    calcTables         = FALSE,
    covMethod          = "",
    atol               = atol,
    rtol               = rtol,
    stickyRecalcN      = stickyRecalcN,
    maxOuterIterations = 2000,
    maxInnerIterations = 2000
  )

  res <- dplyr::bind_rows(
    .bench_one_focei(true_2cmt_scn16_lin,    ds16_01_60, ctrl, "lin",    setting),
    .bench_one_focei(true_2cmt_scn16_refexp, ds16_01_60, ctrl, "refexp", setting)
  )
  # Incremental save -- overwritten each cell, so the file always reflects
  # the most recent completed state.
  saveRDS(res, file.path(stage1_dir16, paste0("focei_tune_", setting, ".rds")))
  res
})
# 
# >>> baseline           sigdig=4 atol=1e-08 rtol=1e-06 sticky=4
# Error in `pmap()`:
# ℹ In index: 1.
# Caused by error:
# ! unused argument: 'atol', 'rtol'
#      ▆
#   1. ├─purrr::pmap_dfr(...)
#   2. │ └─purrr::pmap(.l, .f, ...)
#   3. │   └─purrr:::pmap_("list", .l, .f, ..., .progress = .progress)
#   4. │     ├─purrr:::with_indexed_errors(...)
#   5. │     │ └─base::withCallingHandlers(...)
#   6. │     ├─purrr:::call_with_cleanup(...)
#   7. │     └─global .f(...)
#   8. │       └─nlmixr2est::foceiControl(...)
#   9. │         └─base::stop(...)
#  10. └─base::.handleSimpleError(...)
#  11.   └─purrr (local) h(simpleError(msg, call))
#  12.     └─cli::cli_abort(...)
#  13.       └─rlang::abort(...)
focei_grid <- tibble::tribble(
  ~setting,           ~sigdig, ~atol,  ~rtol,  ~stickyRecalcN,
  "baseline",         4,       1e-8,   1e-6,   4,        # matches scm_focei_n
  "fast",             3,       1e-6,   1e-4,   4,        # speed-only
  "no_sticky",        4,       1e-8,   1e-6,   20,       # block auto-relax
  "tight_ode",        4,       1e-10,  1e-8,   4,        # depth via ODE
  "tight_no_sticky",  4,       1e-10,  1e-8,   20,       # belt + braces
  "fast_no_sticky",   3,       1e-6,   1e-4,   20        # fast but protected
)
.bench_one_focei <- function(model_fn, data, ctrl, model_name, setting_name) {
  t0 <- Sys.time()
  fit <- tryCatch(
    nlmixr2(model_fn, data, est = "focei", control = ctrl),
    error = function(e) list(.failed = TRUE, .err = conditionMessage(e))
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  if (isTRUE(fit$.failed)) {
    return(tibble::tibble(
      model = model_name, setting = setting_name,
      elapsed_min = elapsed / 60,
      final_ofv = NA_real_, best_ofv = NA_real_, ofv_gap = NA_real_,
      n_iter = NA_integer_, bad_solves = NA, hessian_reset = NA,
      not_at_min = NA, converged_well = FALSE, error = fit$.err
    ))
  }

  trace_objf <- if (!is.null(fit$parHistData) && nrow(fit$parHistData) > 0)
    fit$parHistData$objf else NA_real_
  best_ofv <- if (any(!is.na(trace_objf))) min(trace_objf, na.rm = TRUE) else NA_real_
  ri <- fit$runInfo %||% character(0)

  tibble::tibble(
    model = model_name, setting = setting_name,
    elapsed_min   = elapsed / 60,
    final_ofv     = fit$objf,
    best_ofv      = best_ofv,
    ofv_gap       = fit$objf - best_ofv,
    n_iter        = if (!is.null(fit$parHistData)) nrow(fit$parHistData) else NA_integer_,
    bad_solves    = any(grepl("bad solves",     ri, fixed = TRUE)),
    hessian_reset = any(grepl("Hessian reset",  ri, fixed = TRUE)),
    not_at_min    = any(grepl("not at minimum", ri, fixed = TRUE)),
    converged_well = !is.na(fit$objf) && !is.na(best_ofv) &&
                     abs(fit$objf - best_ofv) < 0.5,
    error = NA_character_
  )
}
bench_results <- purrr::pmap_dfr(focei_grid, function(setting, sigdig, atol, rtol, stickyRecalcN) {
  cat(sprintf("\n>>> %-18s sigdig=%d atol=%.0e rtol=%.0e sticky=%d\n",
              setting, sigdig, atol, rtol, stickyRecalcN))
  ctrl <- nlmixr2est::foceiControl(
    sigdig             = sigdig,
    outerOpt           = "bobyqa",
    print              = 0,
    calcTables         = FALSE,
    covMethod          = "",
    stickyRecalcN      = stickyRecalcN,
    maxOuterIterations = 2000,
    maxInnerIterations = 2000,
    rxControl          = rxode2::rxControl(atol = atol, rtol = rtol)
  )

  res <- dplyr::bind_rows(
    .bench_one_focei(true_2cmt_scn16_lin,    ds16_01_60, ctrl, "lin",    setting),
    .bench_one_focei(true_2cmt_scn16_refexp, ds16_01_60, ctrl, "refexp", setting)
  )
  # Incremental save -- overwritten each cell, so the file always reflects
  # the most recent completed state.
  saveRDS(res, file.path(stage1_dir16, paste0("focei_tune_", setting, ".rds")))
  res
})
# 
# >>> baseline           sigdig=4 atol=1e-08 rtol=1e-06 sticky=4
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> fast               sigdig=3 atol=1e-06 rtol=1e-04 sticky=4
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> no_sticky          sigdig=4 atol=1e-08 rtol=1e-06 sticky=20
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> tight_ode          sigdig=4 atol=1e-10 rtol=1e-08 sticky=4
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> tight_no_sticky    sigdig=4 atol=1e-10 rtol=1e-08 sticky=20
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> fast_no_sticky     sigdig=3 atol=1e-06 rtol=1e-04 sticky=20
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# R 4.5.3 offline. Waiting to reconnect.
# R 4.5.3 reconnected.
cat("Speed vs convergence depth -- by setting (sorted within model by speed)\n")
# Speed vs convergence depth -- by setting (sorted within model by speed)
bench_summary <- bench_results %>%
  dplyr::select(model, setting, elapsed_min, n_iter, final_ofv, ofv_gap,
                bad_solves, hessian_reset, not_at_min, converged_well) %>%
  dplyr::arrange(model, elapsed_min)
print(bench_summary, n = Inf)
# # A tibble: 12 × 10
#    model  setting         elapsed_min n_iter final_ofv  ofv_gap bad_solves hessian_reset not_at_min converged_well
#    <chr>  <chr>                 <dbl>  <int>     <dbl>    <dbl> <lgl>      <lgl>         <lgl>      <lgl>         
#  1 lin    fast                  0.578   1575    -3447. 0.000350 FALSE      FALSE         TRUE       TRUE          
#  2 lin    baseline              0.689   1500    -3447. 0.00137  FALSE      FALSE         TRUE       TRUE          
#  3 lin    no_sticky             0.705   1500    -3447. 0.00137  FALSE      FALSE         TRUE       TRUE          
#  4 lin    fast_no_sticky        0.846   1575    -3447. 0.000350 FALSE      FALSE         TRUE       TRUE          
#  5 lin    tight_ode             0.860   1587    -3447. 0.000237 FALSE      FALSE         TRUE       TRUE          
#  6 lin    tight_no_sticky       1.05    1587    -3447. 0.000237 FALSE      FALSE         TRUE       TRUE          
#  7 refexp fast                  0.313    615    -3447. 0.00145  FALSE      FALSE         TRUE       TRUE          
#  8 refexp fast_no_sticky        0.397    615    -3447. 0.00145  FALSE      FALSE         TRUE       TRUE          
#  9 refexp no_sticky             0.466    681    -3447. 0.000824 FALSE      FALSE         TRUE       TRUE          
# 10 refexp baseline              0.476    681    -3447. 0.000824 FALSE      FALSE         TRUE       TRUE          
# 11 refexp tight_ode             0.643    786    -3447. 0.000395 FALSE      FALSE         TRUE       TRUE          
# 12 refexp tight_no_sticky       0.670    786    -3447. 0.000395 FALSE      FALSE         TRUE       TRUE          
focei_recommend <- bench_results %>%
  dplyr::filter(converged_well) %>%
  dplyr::group_by(model) %>%
  dplyr::slice_min(elapsed_min, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(model, setting, elapsed_min, ofv_gap, final_ofv, n_iter)
print(focei_recommend)
# # A tibble: 2 × 6
#   model  setting elapsed_min  ofv_gap final_ofv n_iter
#   <chr>  <chr>         <dbl>    <dbl>     <dbl>  <int>
# 1 lin    fast          0.578 0.000350    -3447.   1575
# 2 refexp fast          0.313 0.00145     -3447.    615
cat("\n--- refexp rescue check ---\n")
# 
# --- refexp rescue check ---
refexp_rescues <- bench_results %>%
  dplyr::filter(model == "refexp") %>%
  dplyr::transmute(setting,
                   bad_solves,
                   ofv_gap = round(ofv_gap, 2),
                   converged_well,
                   elapsed_min = round(elapsed_min, 1))
print(refexp_rescues)
# # A tibble: 6 × 5
#   setting         bad_solves ofv_gap converged_well elapsed_min
#   <chr>           <lgl>        <dbl> <lgl>                <dbl>
# 1 baseline        FALSE            0 TRUE                   0.5
# 2 fast            FALSE            0 TRUE                   0.3
# 3 no_sticky       FALSE            0 TRUE                   0.5
# 4 tight_ode       FALSE            0 TRUE                   0.6
# 5 tight_no_sticky FALSE            0 TRUE                   0.7
# 6 fast_no_sticky  FALSE            0 TRUE                   0.4
bench_results_300 <- purrr::pmap_dfr(focei_grid, function(setting, sigdig, atol, rtol, stickyRecalcN) {
  cat(sprintf("\n>>> %-18s sigdig=%d atol=%.0e rtol=%.0e sticky=%d\n",
              setting, sigdig, atol, rtol, stickyRecalcN))
  ctrl <- nlmixr2est::foceiControl(
    sigdig             = sigdig,
    outerOpt           = "bobyqa",
    print              = 0,
    calcTables         = FALSE,
    covMethod          = "",
    stickyRecalcN      = stickyRecalcN,
    maxOuterIterations = 2000,
    maxInnerIterations = 2000,
    rxControl          = rxode2::rxControl(atol = atol, rtol = rtol)
  )

  res300 <- dplyr::bind_rows(
    .bench_one_focei(true_2cmt_scn16_lin,    ds16_01, ctrl, "lin",    setting),
    .bench_one_focei(true_2cmt_scn16_refexp, ds16_01, ctrl, "refexp", setting)
  )
  # Incremental save -- overwritten each cell, so the file always reflects
  # the most recent completed state.
  saveRDS(res300, file.path(stage1_dir16, paste0("focei_tune_", setting, ".rds")))
  res300
})
# 
# >>> baseline           sigdig=4 atol=1e-08 rtol=1e-06 sticky=4
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> fast               sigdig=3 atol=1e-06 rtol=1e-04 sticky=4
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> no_sticky          sigdig=4 atol=1e-08 rtol=1e-06 sticky=20
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> tight_ode          sigdig=4 atol=1e-10 rtol=1e-08 sticky=4
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> tight_no_sticky    sigdig=4 atol=1e-10 rtol=1e-08 sticky=20
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# >>> fast_no_sticky     sigdig=3 atol=1e-06 rtol=1e-04 sticky=20
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
saveRDS(bench_results_300, file.path(stage1_dir16, "focei_tuning_results_300.rds"))
bench_results_300
# # A tibble: 12 × 12
#    model  setting         elapsed_min final_ofv best_ofv  ofv_gap n_iter bad_solves hessian_reset not_at_min converged_well error
#    <chr>  <chr>                 <dbl>     <dbl>    <dbl>    <dbl>  <int> <lgl>      <lgl>         <lgl>      <lgl>          <chr>
#  1 lin    baseline               1.81   -16901.  -16901. 0.00315     600 FALSE      FALSE         TRUE       TRUE           NA   
#  2 refexp baseline               1.92   -16901.  -16901. 0.00282     540 FALSE      FALSE         TRUE       TRUE           NA   
#  3 lin    fast                   1.09   -16901.  -16901. 0.00202     576 FALSE      FALSE         FALSE      TRUE           NA   
#  4 refexp fast                   1.35   -16901.  -16901. 0.00490     498 FALSE      FALSE         TRUE       TRUE           NA   
#  5 lin    no_sticky              1.54   -16901.  -16901. 0.00315     600 FALSE      FALSE         TRUE       TRUE           NA   
#  6 refexp no_sticky              1.94   -16901.  -16901. 0.00282     540 FALSE      FALSE         TRUE       TRUE           NA   
#  7 lin    tight_ode              2.22   -16901.  -16901. 0.000914    747 FALSE      FALSE         TRUE       TRUE           NA   
#  8 refexp tight_ode              2.39   -16901.  -16901. 0.00156     579 FALSE      FALSE         TRUE       TRUE           NA   
#  9 lin    tight_no_sticky        2.27   -16901.  -16901. 0.000914    747 FALSE      FALSE         TRUE       TRUE           NA   
# 10 refexp tight_no_sticky        2.38   -16901.  -16901. 0.00156     579 FALSE      FALSE         TRUE       TRUE           NA   
# 11 lin    fast_no_sticky         1.07   -16901.  -16901. 0.00202     576 FALSE      FALSE         FALSE      TRUE           NA   
# 12 refexp fast_no_sticky         1.31   -16901.  -16901. 0.00490     498 FALSE      FALSE         TRUE       TRUE           NA   
purrr::map_dfr(1:5, function(seed) {
  set.seed(seed)
  ctrl <- nlmixr2est::foceiControl(
    sigdig = 4, outerOpt = "bobyqa", print = 0,
    calcTables = FALSE, covMethod = "", stickyRecalcN = 4,
    maxOuterIterations = 2000, maxInnerIterations = 2000,
    rxControl = rxode2::rxControl(atol = 1e-8, rtol = 1e-6)
  )
  .bench_one_focei(true_2cmt_scn16_refexp, ds16_01, ctrl, "refexp",
                   paste0("seed_", seed))
})
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# # A tibble: 5 × 12
#   model  setting elapsed_min final_ofv best_ofv ofv_gap n_iter bad_solves hessian_reset not_at_min converged_well error
#   <chr>  <chr>         <dbl>     <dbl>    <dbl>   <dbl>  <int> <lgl>      <lgl>         <lgl>      <lgl>          <chr>
# 1 refexp seed_1         1.88   -16901.  -16901. 0.00282    540 FALSE      FALSE         TRUE       TRUE           NA   
# 2 refexp seed_2         1.82   -16901.  -16901. 0.00282    540 FALSE      FALSE         TRUE       TRUE           NA   
# 3 refexp seed_3         1.77   -16901.  -16901. 0.00282    540 FALSE      FALSE         TRUE       TRUE           NA   
# 4 refexp seed_4         1.84   -16901.  -16901. 0.00282    540 FALSE      FALSE         TRUE       TRUE           NA   
# 5 refexp seed_5         1.78   -16901.  -16901. 0.00282    540 FALSE      FALSE         TRUE       TRUE           NA   

# Build the grid baseline ctrl exactly as the grid does
ctrl_grid <- nlmixr2est::foceiControl(
  sigdig             = 4,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = FALSE,
  covMethod          = "",
  stickyRecalcN      = 4,
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl          = rxode2::rxControl(atol = 1e-8, rtol = 1e-6)
)

# Find which top-level arguments differ between scm_focei_n and ctrl_grid
keys <- union(names(scm_focei_n), names(ctrl_grid))
diffs <- purrr::map_dfr(keys, function(k) {
  a <- scm_focei_n[[k]]; b <- ctrl_grid[[k]]
  same <- tryCatch(identical(a, b), error = function(e) FALSE)
  if (same) return(NULL)
  fmt <- function(x) {
    if (is.null(x)) return("<NULL>")
    if (is.atomic(x) && length(x) <= 4) return(paste(format(x), collapse = ","))
    paste0("<", class(x)[1], " len=", length(x), ">")
  }
  tibble::tibble(key = k, scm_focei_n = fmt(a), grid_ctrl = fmt(b))
})

cat("--- Top-level foceiControl keys that DIFFER ---\n")
print(diffs, n = Inf)

# Drill into rxControl specifically (that's the most likely culprit)
cat("\n--- rxControl side-by-side ---\n")
rx_a <- scm_focei_n$rxControl
rx_b <- ctrl_grid$rxControl
rx_keys <- union(names(rx_a), names(rx_b))
rx_diffs <- purrr::map_dfr(rx_keys, function(k) {
  a <- rx_a[[k]]; b <- rx_b[[k]]
  if (identical(a, b)) return(NULL)
  fmt <- function(x) {
    if (is.null(x)) return("<NULL>")
    if (is.atomic(x) && length(x) <= 4) return(paste(format(x), collapse = ","))
    paste0("<", class(x)[1], " len=", length(x), ">")
  }
  tibble::tibble(key = k, scm_focei_n = fmt(a), grid_ctrl = fmt(b))
})
print(rx_diffs, n = Inf)
# --- Top-level foceiControl keys that DIFFER ---
# # A tibble: 4 × 3
#   key              scm_focei_n              grid_ctrl               
#   <chr>            <chr>                    <chr>                   
# 1 iterPrintControl <iterPrintControl len=5> <iterPrintControl len=5>
# 2 covMethod        1                        0                       
# 3 rxControl        <rxControl len=126>      <rxControl len=126>     
# 4 genRxControl     TRUE                     FALSE                   
# 
# --- rxControl side-by-side ---
# # A tibble: 22 × 3
#    key               scm_focei_n  grid_ctrl
#    <chr>             <chr>        <chr>    
#  1 method            2            2        
#  2 atol              5e-07        1e-08    
#  3 rtol              5e-07        1e-06    
#  4 maxsteps          500000       70000    
#  5 covsInterpolation 1            1        
#  6 returnType        0            0        
#  7 sigmaXform        4            4        
#  8 omegaXform        6            6        
#  9 indLinMatExpType  2            2        
# 10 ssAtol            5e-05        1e-08    
# 11 ssRtol            5e-05        1e-06    
# 12 sumType           1            1        
# 13 prodType          1            1        
# 14 atolSens          1.581139e-06 1e-08    
# 15 rtolSens          1.581139e-06 1e-06    
# 16 ssAtolSens        0.0002108483 1e-08    
# 17 ssRtolSens        0.0002108483 1e-06    
# 18 naTimeHandle      1            1        
# 19 naInterpolation   1            1        
# 20 keepInterpolation 2            2        
# 21 linCmtSensType    100          100      
# 22 cvodeLinSolver    1            1        
po
.bench_one_focei_auto <- function(model_fn, data, sigdig, stickyRecalcN,
                                   model_name, setting_name) {
  # NO rxControl arg => genRxControl=TRUE => auto-loose tolerances from sigdig
  ctrl <- nlmixr2est::foceiControl(
    sigdig = sigdig, outerOpt = "bobyqa", print = 0,
    calcTables = FALSE, covMethod = "",
    stickyRecalcN = stickyRecalcN,
    maxOuterIterations = 2000, maxInnerIterations = 2000
  )
  # ... rest of timing logic same as .bench_one_focei
}

# Then: this is the cell that should reproduce the 244-OFV meltdown
.bench_one_focei_auto(true_2cmt_scn16_refexp, ds16_01,
                      sigdig = 4, stickyRecalcN = 4,
                      "refexp", "true_baseline_auto_rx")
# Error:
# ! object 'po' not found

# Option 1: Test the genuine baseline -- no rxControl arg, let sigdig auto-set
# the (loose) ODE tolerances.  This is what scm_focei_n actually does, and it
# is the hypothesised trigger of the Part-1 refexp 244-OFV meltdown.
.bench_one_focei_auto <- function(model_fn, data, sigdig, stickyRecalcN,
                                  model_name, setting_name) {
  ctrl <- nlmixr2est::foceiControl(
    sigdig             = sigdig,
    outerOpt           = "bobyqa",
    print              = 0,
    calcTables         = FALSE,
    covMethod          = "",
    stickyRecalcN      = stickyRecalcN,
    maxOuterIterations = 2000,
    maxInnerIterations = 2000
    # NO rxControl -> genRxControl=TRUE -> auto-loose tolerances from sigdig
  )

  # Sanity-print what rxControl actually got built
  rx <- ctrl$rxControl
  cat(sprintf("    [auto rxControl] atol=%.1e rtol=%.1e atolSens=%.1e rtolSens=%.1e ssAtol=%.1e\n",
              rx$atol, rx$rtol, rx$atolSens, rx$rtolSens, rx$ssAtol))

  t0 <- Sys.time()
  fit <- tryCatch(
    nlmixr2(model_fn, data, est = "focei", control = ctrl),
    error = function(e) list(.failed = TRUE, .err = conditionMessage(e))
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  if (isTRUE(fit$.failed)) {
    return(tibble::tibble(
      model = model_name, setting = setting_name, elapsed_min = elapsed / 60,
      final_ofv = NA_real_, best_ofv = NA_real_, ofv_gap = NA_real_,
      n_iter = NA_integer_, bad_solves = NA, hessian_reset = NA,
      not_at_min = NA, converged_well = FALSE, error = fit$.err
    ))
  }

  trace_objf <- if (!is.null(fit$parHistData) && nrow(fit$parHistData) > 0)
    fit$parHistData$objf else NA_real_
  best_ofv <- if (any(!is.na(trace_objf))) min(trace_objf, na.rm = TRUE) else NA_real_
  ri <- fit$runInfo %||% character(0)

  tibble::tibble(
    model = model_name, setting = setting_name,
    elapsed_min   = elapsed / 60,
    final_ofv     = fit$objf,
    best_ofv      = best_ofv,
    ofv_gap       = fit$objf - best_ofv,
    n_iter        = if (!is.null(fit$parHistData)) nrow(fit$parHistData) else NA_integer_,
    bad_solves    = any(grepl("bad solves",     ri, fixed = TRUE)),
    hessian_reset = any(grepl("Hessian reset",  ri, fixed = TRUE)),
    not_at_min    = any(grepl("not at minimum", ri, fixed = TRUE)),
    converged_well = !is.na(fit$objf) && !is.na(best_ofv) &&
                     abs(fit$objf - best_ofv) < 0.5,
    error = NA_character_
  )
}

# Run BOTH parameterisations with auto-rxControl, the genuine "baseline"
cat("\n=== refexp + auto rxControl (genuine scm_focei_n baseline) ===\n")
bench_auto_refexp <- .bench_one_focei_auto(
  true_2cmt_scn16_refexp, ds16_01, sigdig = 4, stickyRecalcN = 4,
  "refexp", "true_baseline_auto_rx"
)

cat("\n=== lin + auto rxControl ===\n")
bench_auto_lin <- .bench_one_focei_auto(
  true_2cmt_scn16_lin, ds16_01, sigdig = 4, stickyRecalcN = 4,
  "lin", "true_baseline_auto_rx"
)

bench_auto <- dplyr::bind_rows(bench_auto_refexp, bench_auto_lin)
saveRDS(bench_auto, file.path(stage1_dir16, "focei_tune_auto_baseline.rds"))

cat("\n--- Results: auto-rxControl baseline ---\n")
print(bench_auto %>% dplyr::select(model, setting, elapsed_min, n_iter,
                                    final_ofv, best_ofv, ofv_gap,
                                    bad_solves, hessian_reset, not_at_min,
                                    converged_well))
# 
# === refexp + auto rxControl (genuine scm_focei_n baseline) ===
#     [auto rxControl] atol=5.0e-07 rtol=5.0e-07 atolSens=1.6e-06 rtolSens=1.6e-06 ssAtol=5.0e-05
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# 
# === lin + auto rxControl ===
#     [auto rxControl] atol=5.0e-07 rtol=5.0e-07 atolSens=1.6e-06 rtolSens=1.6e-06 ssAtol=5.0e-05
# done
# 
# --- Results: auto-rxControl baseline ---
# # A tibble: 2 × 11
#   model  setting               elapsed_min n_iter final_ofv best_ofv   ofv_gap bad_solves hessian_reset not_at_min converged_well
#   <chr>  <chr>                       <dbl>  <int>     <dbl>    <dbl>     <dbl> <lgl>      <lgl>         <lgl>      <lgl>         
# 1 refexp true_baseline_auto_rx        1.28    486   -16642.  -16886. 244.      TRUE       TRUE          TRUE       FALSE         
# 2 lin    true_baseline_auto_rx        1.37    624   -16901.  -16901.   0.00693 FALSE      FALSE         TRUE       TRUE          
1300/24
# [1] 54.16667
54/7
# [1] 7.714286

# Sanity-check: do IDs 1-100 vs full N=300 give comparable covariate structure?
# Use dataset 1 of scenario 16 (already loaded as ds16_01).
full_300 <- ds16_01 %>% dplyr::distinct(ID, .keep_all = TRUE)
sub_100  <- full_300 %>% dplyr::filter(ID <= 100)

cat("--- Marginal summaries ---\n")
mk_sum <- function(d) {
  d %>%
    dplyr::summarise(
      N        = dplyr::n(),
      BW_mean  = round(mean(BW),    2), BW_sd  = round(sd(BW),  2),
      CrCL_mean= round(mean(CrCL),  2), CrCL_sd= round(sd(CrCL),2),
      pct_F    = round(100 * mean(SEX == 1),  1),  # 1 = female in your coding
      pct_RACE2= round(100 * mean(RACE == 2), 1)
    )
}
cat("Full N=300:\n");  print(mk_sum(full_300))
cat("Subset N=100:\n"); print(mk_sum(sub_100))

cat("\n--- Continuous covariate correlation (BW vs CrCL) ---\n")
cat(sprintf("Full N=300:   r = %.3f\n",
            cor(full_300$BW, full_300$CrCL)))
cat(sprintf("Subset N=100: r = %.3f\n",
            cor(sub_100$BW,  sub_100$CrCL)))

# Stratification check: are IDs 1-100 a random slice or were they sorted by
# SEX/RACE during simulation?  If the SEX/RACE proportions are wildly off,
# the original sampler stratified and subsetting will bias the dataset.
cat("\n--- Stratification check: SEX by ID block ---\n")
full_300 %>%
  dplyr::mutate(block = cut(ID, breaks = c(0, 100, 200, 300),
                            labels = c("1-100","101-200","201-300"))) %>%
  dplyr::group_by(block) %>%
  dplyr::summarise(N = dplyr::n(),
                   pct_F = round(100 * mean(SEX == 1), 1),
                   pct_RACE_nonref = round(100 * mean(RACE != 1), 1),
                   .groups = "drop")
# --- Marginal summaries ---
# Full N=300:
# # A tibble: 1 × 7
#       N BW_mean BW_sd CrCL_mean CrCL_sd pct_F pct_RACE2
#   <int>   <dbl> <dbl>     <dbl>   <dbl> <dbl>     <dbl>
# 1   300    82.4  19.2      106.    39.1    45         0
# Subset N=100:
# # A tibble: 1 × 7
#       N BW_mean BW_sd CrCL_mean CrCL_sd pct_F pct_RACE2
#   <int>   <dbl> <dbl>     <dbl>   <dbl> <dbl>     <dbl>
# 1   100    88.5  19.7      112.    42.8   100         0
# 
# --- Continuous covariate correlation (BW vs CrCL) ---
# Full N=300:   r = 0.590
# Subset N=100: r = 0.616
# 
# --- Stratification check: SEX by ID block ---
# # A tibble: 3 × 4
#   block       N pct_F pct_RACE_nonref
#   <fct>   <int> <dbl>           <dbl>
# 1 1-100     100   100               6
# 2 101-200   100    35               9
# 3 201-300   100     0               6

# What's the TRUE value of TH_SEX_VC across scenarios?
true_params %>%
  dplyr::filter(grepl("SEX", parameter)) %>%
  dplyr::distinct(parameter, true_value, scenario) %>%
  dplyr::arrange(parameter, scenario) %>%
  head(20)
# # A tibble: 16 × 3
#    parameter true_value scenario
#    <chr>          <dbl>    <int>
#  1 VcSEX            0          1
#  2 VcSEX            0.5        2
#  3 VcSEX            0          3
#  4 VcSEX            0.5        4
#  5 VcSEX            0          5
#  6 VcSEX            0.5        6
#  7 VcSEX            0          7
#  8 VcSEX            0.5        8
#  9 VcSEX            0          9
# 10 VcSEX            0.5       10
# 11 VcSEX            0         11
# 12 VcSEX            0.5       12
# 13 VcSEX            0         13
# 14 VcSEX            0.5       14
# 15 VcSEX            0         15
# 16 VcSEX            0.5       16
colnames(ds16_01)
#  [1] "ID"   "TIME" "EVID" "AMT"  "CMT"  "DV"   "BW"   "BMI"  "CrCL" "SEX"  "RACE"
library(nlmixr2scm)
library(nlmixr2utils)
rxode2::rxClean()
library(nlmixr2)
library(rxode2)
library(tidyverse)
virtual_pop <- read_csv("virtual_population_250x300.csv")
# Rows: 75000 Columns: 7
# ── Column specification ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Delimiter: ","
# dbl (7): BW, BMI, CrCL, SEX, RACE, ID, DATASET
# 
# ℹ Use `spec()` to retrieve the full column specification for this data.
# ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
PsN_scenarios <- data.frame(
  scenario = 1:16,
  I_BW_CL = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
  I_CRCL_CL = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
  I_BW_VC = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
  I_SEX_VC = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
)
sample_dataset_etas <- function(n_subj    = 300,
                                omega_var = 0.1,
                                omega_cov = 0.02,
                                rho_low   = ETA_RHO_LOW,
                                rho_high  = ETA_RHO_HIGH,
                                max_tries = 1000L,
                                seed      = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- matrix(c(omega_var, omega_cov,
                    omega_cov, omega_var), 2, 2)
  for (a in seq_len(max_tries)) {
    e   <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = Sigma)
    rho <- stats::cor(e[, 1], e[, 2])
    if (rho >= rho_low && rho <= rho_high) {
      return(list(eta_cl   = e[, 1],
                  eta_vc   = e[, 2],
                  rho      = rho,
                  attempts = a))
    }
  }
  stop("No valid etas after ", max_tries, " attempts.")
}
build_eta_table <- function(scenarios  = PsN_scenarios$scenario,
                            n_datasets = 250L,
                            n_subj     = 300L) {
  out <- vector("list", length(scenarios) * n_datasets)
  k   <- 0L
  for (scn in scenarios) {
    for (ds in seq_len(n_datasets)) {
      k <- k + 1L
      eta <- sample_dataset_etas(
        n_subj = n_subj,
        seed   = 2026L + scn * 1e4L + ds
      )
      out[[k]] <- tibble::tibble(
        SCENARIO     = scn,
        DATASET      = ds,
        SUBJECT      = seq_len(n_subj),
        eta_cl       = eta$eta_cl,
        eta_vc       = eta$eta_vc,
        rho_realised = eta$rho,
        attempts     = eta$attempts
      )
    }
  }
  dplyr::bind_rows(out)
}
eta_table <- build_eta_table()
# Error in `sample_dataset_etas()`:
# ! object 'ETA_RHO_LOW' not found
#     ▆
#  1. └─global build_eta_table()
#  2.   └─global sample_dataset_etas(...)
ETA_RHO_LOW  <- 0.15
ETA_RHO_HIGH <- 0.25
sample_dataset_etas <- function(n_subj    = 300,
                                omega_var = 0.1,
                                omega_cov = 0.02,
                                rho_low   = ETA_RHO_LOW,
                                rho_high  = ETA_RHO_HIGH,
                                max_tries = 1000L,
                                seed      = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- matrix(c(omega_var, omega_cov,
                    omega_cov, omega_var), 2, 2)
  for (a in seq_len(max_tries)) {
    e   <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = Sigma)
    rho <- stats::cor(e[, 1], e[, 2])
    if (rho >= rho_low && rho <= rho_high) {
      return(list(eta_cl   = e[, 1],
                  eta_vc   = e[, 2],
                  rho      = rho,
                  attempts = a))
    }
  }
  stop("No valid etas after ", max_tries, " attempts.")
}
build_eta_table <- function(scenarios  = PsN_scenarios$scenario,
                            n_datasets = 250L,
                            n_subj     = 300L) {
  out <- vector("list", length(scenarios) * n_datasets)
  k   <- 0L
  for (scn in scenarios) {
    for (ds in seq_len(n_datasets)) {
      k <- k + 1L
      eta <- sample_dataset_etas(
        n_subj = n_subj,
        seed   = 2026L + scn * 1e4L + ds
      )
      out[[k]] <- tibble::tibble(
        SCENARIO     = scn,
        DATASET      = ds,
        SUBJECT      = seq_len(n_subj),
        eta_cl       = eta$eta_cl,
        eta_vc       = eta$eta_vc,
        rho_realised = eta$rho,
        attempts     = eta$attempts
      )
    }
  }
  dplyr::bind_rows(out)
}
eta_table <- build_eta_table()
eta_table_qc <- eta_table %>%
  dplyr::distinct(SCENARIO, DATASET, rho_realised, attempts) %>%
  dplyr::summarise(
    min_rho        = min(rho_realised),
    max_rho        = max(rho_realised),
    mean_rho       = mean(rho_realised),
    mean_attempts  = mean(attempts),
    max_attempts   = max(attempts)
  )
sim_mod_fixed <- rxode2::rxode2({
  CL_cov = TVCL * (BW   / BW_REF  )^(TH_BW_CL   * I_BW_CL  ) *
                  (CrCL / CRCL_REF)^(TH_CRCL_CL * I_CRCL_CL)
  Vc_cov = TVVc * (BW   / BW_REF  )^(TH_BW_VC   * I_BW_VC  ) *
                  exp(TH_SEX_VC * SEX * I_SEX_VC)
  CL = CL_cov * exp(eta_cl)
  Vc = Vc_cov * exp(eta_vc)
  Q  = TVQ
  Vp = TVVp
  KA = TVKA
  k10 = CL / Vc
  k12 = Q  / Vc
  k21 = Q  / Vp
  d/dt(depot)      = -KA * depot
  d/dt(central)    =  KA * depot - k10 * central - k12 * central + k21 * peripheral
  d/dt(peripheral) =  k12 * central - k21 * peripheral
  cp = central / Vc
})
simulate_scenario_v2 <- function(scn,
                                 sim_mod       = sim_mod_fixed,
                                 scenarios     = PsN_scenarios,
                                 icov          = iCov_full,
                                 event_table   = ev,
                                 sampling_grid = sample_times,
                                 t_half        = t_half_typ,
                                 etas          = eta_table,
                                 prop_err      = 0.1,
                                 seed_residual = 4242L + scn) {

  scn_row <- scenarios %>% dplyr::filter(scenario == scn)
  if (nrow(scn_row) != 1L) stop("scenario ", scn, " not found.")

  ## Per-subject covariates + pre-sampled etas
  icov_scn <- icov %>%
    dplyr::left_join(
      etas %>%
        dplyr::filter(SCENARIO == scn) %>%
        dplyr::select(DATASET, SUBJECT, eta_cl, eta_vc),
      by = c("DATASET", "SUBJECT")
    ) %>%
    dplyr::select(id, BW, CrCL, SEX, eta_cl, eta_vc)

  ## Scenario-level thetas
  params_scn <- c(
    TVCL       = 0.6,  TVQ       = 1.8,
    TVVc       = 20,   TVVp      = 80,
    TVKA       = 0.7,
    BW_REF     = 70,   CRCL_REF  = 95,
    TH_BW_CL   = 0.75, TH_CRCL_CL = 0.5,
    TH_BW_VC   = 1.0,  TH_SEX_VC  = log(1.5),   # nlmixr2scm convention: exp(theta*SEX); log(1.5) ~ 0.4055 keeps the 1.5x biology of the original design
    I_BW_CL    = scn_row$I_BW_CL,
    I_CRCL_CL  = scn_row$I_CRCL_CL,
    I_BW_VC    = scn_row$I_BW_VC,
    I_SEX_VC   = scn_row$I_SEX_VC
  )

  sim_raw_scn <- rxode2::rxSolve(
    sim_mod,
    params     = params_scn,
    events     = event_table,
    iCov       = icov_scn,
    returnType = "tibble"
  )

  set.seed(seed_residual)
  sim_raw_scn <- sim_raw_scn %>%
    dplyr::mutate(cp_obs = cp * (1 + stats::rnorm(dplyr::n(), 0, prop_err)))

  sim_raw_scn %>%
    dplyr::filter(time %in% sampling_grid) %>%
    dplyr::mutate(
      SCENARIO = scn,
      cp_ipred = cp,
      HL_MULT  = round(time / t_half, 4)
    ) %>%
    ## Drop iCov pass-through eta columns so the join below doesn't create suffixes
    dplyr::select(-dplyr::any_of(c("eta_cl", "eta_vc"))) %>%
    dplyr::left_join(
      icov %>% dplyr::select(id, DATASET, SUBJECT, BMI, RACE),
      by = "id"
    ) %>%
    dplyr::left_join(
      etas %>% dplyr::filter(SCENARIO == scn) %>%
        dplyr::select(DATASET, SUBJECT, eta_cl, eta_vc),
      by = c("DATASET", "SUBJECT")
    ) %>%
    dplyr::select(SCENARIO, DATASET, SUBJECT, HL_MULT, time,
                  BW, BMI, CrCL, SEX, RACE,
                  eta_cl, eta_vc, CL, Vc,
                  cp_ipred, cp_obs)
}
out_dir_v2 <- "simulated_virtual_dataset_eta_filtered"
saveRDS(sim_pkobs_nlmixr2para, file.path(out_dir_v2, "sim_pkobs_nlmixr2para_16scenarios.rds"))
# Error:
# ! object 'sim_pkobs_nlmixr2para' not found
#     ▆
#  1. └─base::saveRDS(sim_pkobs_nlmixr2para, file.path(out_dir_v2, "sim_pkobs_nlmixr2para_16scenarios.rds"))
sim_pkobs_nlmixr2para <- purrr::map(
  .x = PsN_scenarios$scenario,
  .f = function(scn) {
    message(sprintf("Simulating scenario %02d (rejection-sampled etas) ...", scn))
    out <- simulate_scenario_v2(scn)
    saveRDS(out, file.path(out_dir_v2,
                           sprintf("sim_obs_scenario_%02d.rds", scn)))
    out
  },
  .progress = TRUE
)
# Simulating scenario 01 (rejection-sampled etas) ...
# Error in `purrr::map()`:
# ℹ In index: 1.
# Caused by error in `simulate_scenario_v2()`:
# ! object 'iCov_full' not found
#      ▆
#   1. ├─purrr::map(...)
#   2. │ └─purrr:::map_("list", .x, .f, ..., .progress = .progress)
#   3. │   ├─purrr:::with_indexed_errors(...)
#   4. │   │ └─base::withCallingHandlers(...)
#   5. │   ├─purrr:::call_with_cleanup(...)
#   6. │   └─global .f(.x[[i]], ...)
#   7. │     └─global simulate_scenario_v2(scn)
#   8. │       └─... %>% dplyr::select(id, BW, CrCL, SEX, eta_cl, eta_vc)
#   9. ├─dplyr::select(., id, BW, CrCL, SEX, eta_cl, eta_vc)
#  10. ├─dplyr::left_join(...)
#  11. └─base::.handleSimpleError(...)
#  12.   └─purrr (local) h(simpleError(msg, call))
#  13.     └─cli::cli_abort(...)
#  14.       └─rlang::abort(...)
iCov_full <- virtual_pop %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  dplyr::select(id, DATASET, SUBJECT = ID,
                BW, BMI, CrCL, SEX, RACE)
TVCL <- 0.6;  TVQ <- 1.8;  TVVc <- 20;  TVVp <- 80
k10_typ <- TVCL / TVVc
k12_typ <- TVQ  / TVVc
k21_typ <- TVQ  / TVVp
k21_typ <- TVQ  / TVVp
sum_k    <- k10_typ + k12_typ + k21_typ
beta_typ <- 0.5 * (sum_k - sqrt(sum_k^2 - 4 * k10_typ * k21_typ))
beta_typ <- 0.5 * (sum_k - sqrt(sum_k^2 - 4 * k10_typ * k21_typ))
t_half_typ <- log(2) / beta_typ
cat(sprintf("Typical terminal half-life: %.2f h\n", t_half_typ))
# Typical terminal half-life: 141.29 h
hl_mult      <- c(0, 0.05, 0.1, 0.5, 1, 3)
sample_times <- hl_mult * t_half_typ
DOSE_MG <- 100
sim_pkobs_nlmixr2para <- purrr::map(
  .x = PsN_scenarios$scenario,
  .f = function(scn) {
    message(sprintf("Simulating scenario %02d (rejection-sampled etas) ...", scn))
    out <- simulate_scenario_v2(scn)
    saveRDS(out, file.path(out_dir_v2,
                           sprintf("sim_obs_scenario_%02d.rds", scn)))
    out
  },
  .progress = TRUE
)
# Simulating scenario 01 (rejection-sampled etas) ...
# Error in `purrr::map()`:
# ℹ In index: 1.
# Caused by error in `simulate_scenario_v2()`:
# ! object 'ev' not found
#      ▆
#   1. ├─purrr::map(...)
#   2. │ └─purrr:::map_("list", .x, .f, ..., .progress = .progress)
#   3. │   ├─purrr:::with_indexed_errors(...)
#   4. │   │ └─base::withCallingHandlers(...)
#   5. │   ├─purrr:::call_with_cleanup(...)
#   6. │   └─global .f(.x[[i]], ...)
#   7. │     └─global simulate_scenario_v2(scn)
#   8. │       ├─rxode2::rxSolve(...)
#   9. │       └─rxode2:::rxSolve.default(...)
#  10. │         ├─base::setdiff(...)
#  11. │         ├─base::intersect(tolower(names(events)), tolower(names(.ctl$iCov)))
#  12. │         └─base::tolower(names(events))
#  13. └─base::.handleSimpleError(`<fn>`, "object 'ev' not found", base::quote(simulate_scenario_v2(scn)))
#  14.   └─purrr (local) h(simpleError(msg, call))
#  15.     └─cli::cli_abort(...)
#  16.       └─rlang::abort(...)
ev_one <- rxode2::et(amt = DOSE_MG, cmt = "depot", time = 0) %>%
  rxode2::et(time = sample_times) %>%
  as.data.frame()
ev <- ev_one %>%
  dplyr::slice(rep(dplyr::row_number(), n_total)) %>%
  dplyr::mutate(id = rep(seq_len(n_total), each = nrow(ev_one))) %>%
  dplyr::arrange(id, time)
# Error in `dplyr::slice()`:
# ℹ In argument: `rep(dplyr::row_number(), n_total)`.
# Caused by error:
# ! object 'n_total' not found
#      ▆
#   1. ├─... %>% dplyr::arrange(id, time)
#   2. ├─dplyr::arrange(., id, time)
#   3. ├─dplyr::mutate(., id = rep(seq_len(n_total), each = nrow(ev_one)))
#   4. ├─dplyr::slice(., rep(dplyr::row_number(), n_total))
#   5. ├─dplyr:::slice.data.frame(., rep(dplyr::row_number(), n_total))
#   6. │ └─dplyr:::slice_rows(.data, dots, by)
#   7. │   └─dplyr:::slice_eval(mask, dots, error_call = error_call, user_env = user_env)
#   8. │     ├─base::withCallingHandlers(...)
#   9. │     └─mask$eval_all(quo(impl(!!!dots)))
#  10. ├─dplyr (local) impl(rep(dplyr::row_number(), n_total))
#  11. ├─rlang (local) rep(dplyr::row_number(), n_total)
#  12. └─base::.handleSimpleError(...)
#  13.   └─dplyr (local) h(simpleError(msg, call))
#  14.     └─rlang::abort(bullets, call = error_call, parent = cnd)
n_total <- nrow(virtual_pop)
ev <- ev_one %>%
  dplyr::slice(rep(dplyr::row_number(), n_total)) %>%
  dplyr::mutate(id = rep(seq_len(n_total), each = nrow(ev_one))) %>%
  dplyr::arrange(id, time)
sim_pkobs_nlmixr2para <- purrr::map(
  .x = PsN_scenarios$scenario,
  .f = function(scn) {
    message(sprintf("Simulating scenario %02d (rejection-sampled etas) ...", scn))
    out <- simulate_scenario_v2(scn)
    saveRDS(out, file.path(out_dir_v2,
                           sprintf("sim_obs_scenario_%02d.rds", scn)))
    out
  },
  .progress = TRUE
)
# Simulating scenario 01 (rejection-sampled etas) ...
# Simulating scenario 02 (rejection-sampled etas) ...                             
# Simulating scenario 03 (rejection-sampled etas) ...                             
# Simulating scenario 04 (rejection-sampled etas) ...                             
# Simulating scenario 05 (rejection-sampled etas) ...                             
# Simulating scenario 06 (rejection-sampled etas) ...                             
# Simulating scenario 07 (rejection-sampled etas) ...                             
# Simulating scenario 08 (rejection-sampled etas) ...                             
# Simulating scenario 09 (rejection-sampled etas) ...                             
# Simulating scenario 10 (rejection-sampled etas) ...                             
# Simulating scenario 11 (rejection-sampled etas) ...                             
# Simulating scenario 12 (rejection-sampled etas) ...                             
# Simulating scenario 13 (rejection-sampled etas) ...                             
# Simulating scenario 14 (rejection-sampled etas) ...                             
# Simulating scenario 15 (rejection-sampled etas) ...                             
# Simulating scenario 16 (rejection-sampled etas) ...                             
#                                                                                 
names(sim_pkobs_nlmixr2para) <- sprintf("scenario_%02d", PsN_scenarios$scenario)
sim_pkobs_nlmixr2para <- dplyr::bind_rows(sim_pkobs_nlmixr2para)
saveRDS(sim_pkobs_nlmixr2para, file.path(out_dir_v2, "sim_pkobs_nlmixr2para_16scenarios.rds"))
eta_cor_nlmixr2para <- sim_pkobs_nlmixr2para  %>%
  dplyr::distinct(SCENARIO, DATASET, SUBJECT, eta_cl, eta_vc) %>%
  dplyr::group_by(SCENARIO, DATASET) %>%
  dplyr::summarise(rho_eta = stats::cor(eta_cl, eta_vc), .groups = "drop") %>%
  dplyr::mutate(pass = rho_eta >= ETA_RHO_LOW & rho_eta <= ETA_RHO_HIGH)
eta_cor_summary_nlmixr2para <- eta_cor_nlmixr2para %>%
  dplyr::group_by(SCENARIO) %>%
  dplyr::summarise(
    n_datasets = dplyr::n(),
    n_pass     = sum(pass),
    pct_pass   = mean(pass) * 100,
    rho_mean   = mean(rho_eta),
    rho_sd     = stats::sd(rho_eta),
    .groups    = "drop"
  )
build_sister_cohort <- function(n_subj,
                                n_datasets = 250L,
                                tol        = 0.05,
                                pop_seed   = 4242L + n_subj) {
  out_dir_n <- sprintf("simulated_virtual_dataset_eta_filtered_N%d", n_subj)
  if (!dir.exists(out_dir_n)) dir.create(out_dir_n, recursive = TRUE)

  ## --- (a) Bootstrap covariates: same algorithm as the N=300 block above
  set.seed(pop_seed)
  sex_prop_n <- mean(pop$SEX)
  males_p    <- dplyr::filter(pop, SEX == 1)
  females_p  <- dplyr::filter(pop, SEX == 0)
  n_male_n   <- round(n_subj * sex_prop_n)
  n_fem_n    <- n_subj - n_male_n

  datasets_n <- list()
  attempt_n  <- 0L
  while (length(datasets_n) < n_datasets) {
    attempt_n <- attempt_n + 1L
    samp <- dplyr::bind_rows(
      males_p  [sample(nrow(males_p),   n_male_n, replace = TRUE), ],
      females_p[sample(nrow(females_p), n_fem_n,  replace = TRUE), ]
    )
    diff_cor <- abs(stats::cor(samp) - ref_cor)
    if (all(diff_cor[upper.tri(diff_cor)] <= tol)) {
      datasets_n[[length(datasets_n) + 1L]] <- samp
    }
    if (attempt_n > 1e5L) stop("Tolerance too strict at N = ", n_subj)
  }
  for (i in seq_along(datasets_n)) {
    datasets_n[[i]]$ID      <- 1:n_subj
    datasets_n[[i]]$DATASET <- i
  }
  vp_n <- dplyr::bind_rows(datasets_n)
  saveRDS(vp_n, file.path(out_dir_n,
                          sprintf("virtual_population_250x%d.rds", n_subj)))

  ## --- (b) eta table: re-use build_eta_table() at the requested n_subj
  eta_table_n <- build_eta_table(n_datasets = n_datasets, n_subj = n_subj)
  saveRDS(eta_table_n, file.path(out_dir_n, "eta_table_valid.rds"))

  ## --- (c) iCov_n: same schema simulate_scenario_v2() expects from iCov_full
  ##     (id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)
  iCov_n <- vp_n %>%
    dplyr::mutate(SUBJECT = ID,
                  id      = (DATASET - 1L) * n_subj + ID) %>%
    dplyr::select(id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)

  ## --- (d) Run all 16 scenarios via the existing simulate_scenario_v2()
  sim_list_n <- purrr::map(
    .x = PsN_scenarios$scenario,
    .f = function(scn) {
      message(sprintf("[N=%d] Simulating scenario %02d ...", n_subj, scn))
      out <- simulate_scenario_v2(scn, icov = iCov_n, etas = eta_table_n)
      saveRDS(out, file.path(out_dir_n,
                             sprintf("sim_obs_scenario_%02d.rds", scn)))
      out
    },
    .progress = TRUE
  )
  names(sim_list_n) <- sprintf("scenario_%02d", PsN_scenarios$scenario)
  sim_all_n <- dplyr::bind_rows(sim_list_n)
  saveRDS(sim_all_n, file.path(out_dir_n, "sim_obs_all_scenarios.rds"))

  invisible(list(virtual_pop = vp_n,
                 eta_table   = eta_table_n,
                 sim_obs_all = sim_all_n,
                 out_dir     = out_dir_n))
}
cohort_40 <- build_sister_cohort(n_subj = 40L)
# Error in `build_sister_cohort()`:
# ! object 'pop' not found
#     ▆
#  1. └─global build_sister_cohort(n_subj = 40L)
#  2.   └─base::mean(pop$SEX)
cohort_80 <- build_sister_cohort(n_subj = 80L)
# Error in `build_sister_cohort()`:
# ! object 'pop' not found
#     ▆
#  1. └─global build_sister_cohort(n_subj = 80L)
#  2.   └─base::mean(pop$SEX)
demo   <- read_xpt("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/NHANES_dataset/DEMO_L.XPT")   %>% select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3)
# Error in `read_xpt()`:
# ! could not find function "read_xpt"
#     ▆
#  1. ├─... %>% select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3)
#  2. └─dplyr::select(., SEQN, RIDAGEYR, RIAGENDR, RIDRETH3)
library(haven)
demo   <- read_xpt("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/NHANES_dataset/DEMO_L.XPT")   %>% select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3)
bmx    <- read_xpt("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/NHANES_dataset/BMX_L.XPT")    %>% select(SEQN, BMXWT, BMXHT, BMXBMI)
biopro <- read_xpt("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/NHANES_dataset/BIOPRO_L.XPT") %>% select(SEQN, LBXSCR)
nhanes <- demo %>%
  inner_join(bmx, by = "SEQN") %>%
  inner_join(biopro, by = "SEQN")
pop <- nhanes %>%
  filter(RIDAGEYR > 17,
         RIDRETH3 %in% c(3, 6)) %>%   # 3 = Non-Hispanic White, 6 = Non-Hispanic Asian
  mutate(
    SEX  = ifelse(RIAGENDR == 2, 0, 1),       # 0 = female, 1 = male
    RACE = ifelse(RIDRETH3  == 6, 0, 1),      # 0 = Asian, 1 = White
    BW   = BMXWT,
    BMI  = BMXBMI,
    # Cockcroft-Gault CrCL (mL/min)
    CrCL = ((140 - RIDAGEYR) * BW / (72 * LBXSCR)) * ifelse(SEX == 0, 0.85, 1)
  ) %>%
  dplyr::select(BW, BMI, CrCL, SEX, RACE) %>%
  na.omit() %>%
  filter(CrCL > 0, CrCL < 250)
cohort_40 <- build_sister_cohort(n_subj = 40L)
# Error in `build_sister_cohort()`:
# ! object 'ref_cor' not found
#     ▆
#  1. └─global build_sister_cohort(n_subj = 40L)
cohort_40 <- build_sister_cohort(n_subj = 40L)
# Error in `build_sister_cohort()`:
# ! object 'ref_cor' not found
#     ▆
#  1. └─global build_sister_cohort(n_subj = 40L)
build_sister_cohort <- function(n_subj,
                                n_datasets = 250L,
                                tol        = 0.05,
                                pop_seed   = 4242L + n_subj) {
  out_dir_n <- sprintf("simulated_virtual_dataset_eta_filtered_N%d", n_subj)
  if (!dir.exists(out_dir_n)) dir.create(out_dir_n, recursive = TRUE)

  ## --- (a) Bootstrap covariates: same algorithm as the N=300 block above
  set.seed(pop_seed)
  sex_prop_n <- mean(pop$SEX)
  males_p    <- dplyr::filter(pop, SEX == 1)
  females_p  <- dplyr::filter(pop, SEX == 0)
  n_male_n   <- round(n_subj * sex_prop_n)
  n_fem_n    <- n_subj - n_male_n

  datasets_n <- list()
  attempt_n  <- 0L
  while (length(datasets_n) < n_datasets) {
    attempt_n <- attempt_n + 1L
    samp <- dplyr::bind_rows(
      males_p  [sample(nrow(males_p),   n_male_n, replace = TRUE), ],
      females_p[sample(nrow(females_p), n_fem_n,  replace = TRUE), ]
    )
    diff_cor <- abs(stats::cor(samp) - ref_cor)
    if (all(diff_cor[upper.tri(diff_cor)] <= tol)) {
      datasets_n[[length(datasets_n) + 1L]] <- samp
    }
    if (attempt_n > 1e5L) stop("Tolerance too strict at N = ", n_subj)
  }
  for (i in seq_along(datasets_n)) {
    datasets_n[[i]]$ID      <- 1:n_subj
    datasets_n[[i]]$DATASET <- i
  }
  vp_n <- dplyr::bind_rows(datasets_n)
  saveRDS(vp_n, file.path(out_dir_n,
                          sprintf("virtual_population_250x%d.rds", n_subj)))

  ## --- (b) eta table: re-use build_eta_table() at the requested n_subj
  eta_table_n <- build_eta_table(n_datasets = n_datasets, n_subj = n_subj)
  saveRDS(eta_table_n, file.path(out_dir_n, "eta_table_valid.rds"))

  ## --- (c) iCov_n: same schema simulate_scenario_v2() expects from iCov_full
  ##     (id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)
  iCov_n <- vp_n %>%
    dplyr::mutate(SUBJECT = ID,
                  id      = (DATASET - 1L) * n_subj + ID) %>%
    dplyr::select(id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)

  ## --- (d) Run all 16 scenarios via the existing simulate_scenario_v2()
  sim_list_n <- purrr::map(
    .x = PsN_scenarios$scenario,
    .f = function(scn) {
      message(sprintf("[N=%d] Simulating scenario %02d ...", n_subj, scn))
      out <- simulate_scenario_v2(scn, icov = iCov_n, etas = eta_table_n)
      saveRDS(out, file.path(out_dir_n,
                             sprintf("sim_obs_scenario_%02d.rds", scn)))
      out
    },
    .progress = TRUE
  )
  names(sim_list_n) <- sprintf("scenario_%02d", PsN_scenarios$scenario)
  sim_all_n <- dplyr::bind_rows(sim_list_n)
  saveRDS(sim_all_n, file.path(out_dir_n, "sim_obs_all_scenarios.rds"))

  invisible(list(virtual_pop = vp_n,
                 eta_table   = eta_table_n,
                 sim_obs_all = sim_all_n,
                 out_dir     = out_dir_n))
}
cohort_40 <- build_sister_cohort(n_subj = 40L)
# Error in `build_sister_cohort()`:
# ! object 'ref_cor' not found
#     ▆
#  1. └─global build_sister_cohort(n_subj = 40L)
ref_cor <- cor(pop)
cohort_40 <- build_sister_cohort(n_subj = 40L)
# There were 50 or more warnings (use warnings() to see the first 50)
# Error in `if (all(diff_cor[upper.tri(diff_cor)] <= tol)) ...`:
# ! missing value where TRUE/FALSE needed
#     ▆
#  1. └─global build_sister_cohort(n_subj = 40L)
cohort_80 <- build_sister_cohort(n_subj = 80L)
# There were 30 warnings (use warnings() to see them)
# Error in `if (all(diff_cor[upper.tri(diff_cor)] <= tol)) ...`:
# ! missing value where TRUE/FALSE needed
#     ▆
#  1. └─global build_sister_cohort(n_subj = 80L)
warnings()
# Warning messages:
# 1: In stats::cor(samp) : the standard deviation is zero
# 2: In stats::cor(samp) : the standard deviation is zero
# 3: In stats::cor(samp) : the standard deviation is zero
# 4: In stats::cor(samp) : the standard deviation is zero
# 5: In stats::cor(samp) : the standard deviation is zero
# 6: In stats::cor(samp) : the standard deviation is zero
# 7: In stats::cor(samp) : the standard deviation is zero
# 8: In stats::cor(samp) : the standard deviation is zero
# 9: In stats::cor(samp) : the standard deviation is zero
# 10: In stats::cor(samp) : the standard deviation is zero
# 11: In stats::cor(samp) : the standard deviation is zero
# 12: In stats::cor(samp) : the standard deviation is zero
# 13: In stats::cor(samp) : the standard deviation is zero
# 14: In stats::cor(samp) : the standard deviation is zero
# 15: In stats::cor(samp) : the standard deviation is zero
# 16: In stats::cor(samp) : the standard deviation is zero
# 17: In stats::cor(samp) : the standard deviation is zero
# 18: In stats::cor(samp) : the standard deviation is zero
# 19: In stats::cor(samp) : the standard deviation is zero
# 20: In stats::cor(samp) : the standard deviation is zero
# 21: In stats::cor(samp) : the standard deviation is zero
# 22: In stats::cor(samp) : the standard deviation is zero
# 23: In stats::cor(samp) : the standard deviation is zero
# 24: In stats::cor(samp) : the standard deviation is zero
# 25: In stats::cor(samp) : the standard deviation is zero
# 26: In stats::cor(samp) : the standard deviation is zero
# 27: In stats::cor(samp) : the standard deviation is zero
# 28: In stats::cor(samp) : the standard deviation is zero
# 29: In stats::cor(samp) : the standard deviation is zero
# 30: In stats::cor(samp) : the standard deviation is zero
build_sister_cohort <- function(n_subj,
                                n_datasets = 250L,
                                tol        = 0.05,
                                pop_seed   = 4242L + n_subj) {
  out_dir_n <- sprintf("simulated_virtual_dataset_eta_filtered_N%d", n_subj)
  if (!dir.exists(out_dir_n)) dir.create(out_dir_n, recursive = TRUE)

  ## --- (a) Bootstrap covariates
  set.seed(pop_seed)
  sex_prop_n <- mean(pop$SEX)
  males_p    <- dplyr::filter(pop, SEX == 1)
  females_p  <- dplyr::filter(pop, SEX == 0)
  n_male_n   <- round(n_subj * sex_prop_n)
  n_fem_n    <- n_subj - n_male_n

  datasets_n <- list()
  attempt_n  <- 0L
  n_bad_var  <- 0L
  while (length(datasets_n) < n_datasets) {
    attempt_n <- attempt_n + 1L
    samp <- dplyr::bind_rows(
      males_p  [sample(nrow(males_p),   n_male_n, replace = TRUE), ],
      females_p[sample(nrow(females_p), n_fem_n,  replace = TRUE), ]
    )

    # Skip draws with zero-variance columns (e.g. all subjects same RACE),
    # which would make cor() return NA and crash the comparison.
    col_sd <- vapply(samp, stats::sd, numeric(1))
    if (any(is.na(col_sd)) || any(col_sd == 0)) {
      n_bad_var <- n_bad_var + 1L
      if (attempt_n > 1e5L) stop("Tolerance too strict at N = ", n_subj)
      next
    }

    diff_cor <- abs(stats::cor(samp) - ref_cor)
    upper    <- diff_cor[upper.tri(diff_cor)]
    if (anyNA(upper)) {
      if (attempt_n > 1e5L) stop("Tolerance too strict at N = ", n_subj)
      next
    }

    if (all(upper <= tol)) {
      datasets_n[[length(datasets_n) + 1L]] <- samp
    }
    if (attempt_n > 1e5L) stop("Tolerance too strict at N = ", n_subj)
  }
  message(sprintf("[N=%d] kept %d / tried %d (%d zero-var draws skipped)",
                  n_subj, length(datasets_n), attempt_n, n_bad_var))

  for (i in seq_along(datasets_n)) {
    datasets_n[[i]]$ID      <- 1:n_subj
    datasets_n[[i]]$DATASET <- i
  }
  vp_n <- dplyr::bind_rows(datasets_n)
  saveRDS(vp_n, file.path(out_dir_n,
                          sprintf("virtual_population_250x%d.rds", n_subj)))

  ## --- (b) eta table
  eta_table_n <- build_eta_table(n_datasets = n_datasets, n_subj = n_subj)
  saveRDS(eta_table_n, file.path(out_dir_n, "eta_table_valid.rds"))

  ## --- (c) iCov_n
  iCov_n <- vp_n |>
    dplyr::mutate(SUBJECT = ID,
                  id      = (DATASET - 1L) * n_subj + ID) |>
    dplyr::select(id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)

  ## --- (d) Run 16 scenarios
  sim_list_n <- purrr::map(
    .x = PsN_scenarios$scenario,
    .f = function(scn) {
      message(sprintf("[N=%d] Simulating scenario %02d ...", n_subj, scn))
      out <- simulate_scenario_v2(scn, icov = iCov_n, etas = eta_table_n)
      saveRDS(out, file.path(out_dir_n,
                             sprintf("sim_obs_scenario_%02d.rds", scn)))
      out
    },
    .progress = TRUE
  )
  names(sim_list_n) <- sprintf("scenario_%02d", PsN_scenarios$scenario)
  sim_all_n <- dplyr::bind_rows(sim_list_n)
  saveRDS(sim_all_n, file.path(out_dir_n, "sim_obs_all_scenarios.rds"))

  invisible(list(virtual_pop = vp_n,
                 eta_table   = eta_table_n,
                 sim_obs_all = sim_all_n,
                 out_dir     = out_dir_n))
}
cohort_80 <- build_sister_cohort(n_subj = 80L)
# Error in `build_sister_cohort()`:
# ! Tolerance too strict at N = 80
#     ▆
#  1. └─global build_sister_cohort(n_subj = 80L)
cohort_40 <- build_sister_cohort(n_subj = 40L, max_attempts = 1e6L,tol = 0.05 * sqrt(300/40))
# Error in `build_sister_cohort()`:
# ! unused argument (max_attempts = 1000000)
build_sister_cohort <- function(n_subj,
                                n_datasets = 100L,
                                tol        = 0.05,
                                pop_seed   = 4242L + n_subj) {
  out_dir_n <- sprintf("simulated_virtual_dataset_eta_filtered_N%d", n_subj)
  if (!dir.exists(out_dir_n)) dir.create(out_dir_n, recursive = TRUE)

  ## --- (a) Bootstrap covariates: same algorithm as the N=300 block above
  set.seed(pop_seed)
  sex_prop_n <- mean(pop$SEX)
  males_p    <- dplyr::filter(pop, SEX == 1)
  females_p  <- dplyr::filter(pop, SEX == 0)
  n_male_n   <- round(n_subj * sex_prop_n)
  n_fem_n    <- n_subj - n_male_n

  datasets_n <- list()
  attempt_n  <- 0L
  while (length(datasets_n) < n_datasets) {
    attempt_n <- attempt_n + 1L
    samp <- dplyr::bind_rows(
      males_p  [sample(nrow(males_p),   n_male_n, replace = TRUE), ],
      females_p[sample(nrow(females_p), n_fem_n,  replace = TRUE), ]
    )
    diff_cor <- abs(stats::cor(samp) - ref_cor)
    if (all(diff_cor[upper.tri(diff_cor)] <= tol)) {
      datasets_n[[length(datasets_n) + 1L]] <- samp
    }
    if (attempt_n > 1e5L) stop("Tolerance too strict at N = ", n_subj)
  }
  for (i in seq_along(datasets_n)) {
    datasets_n[[i]]$ID      <- 1:n_subj
    datasets_n[[i]]$DATASET <- i
  }
  vp_n <- dplyr::bind_rows(datasets_n)
  saveRDS(vp_n, file.path(out_dir_n,
                          sprintf("virtual_population_250x%d.rds", n_subj)))

  ## --- (b) eta table: re-use build_eta_table() at the requested n_subj
  eta_table_n <- build_eta_table(n_datasets = n_datasets, n_subj = n_subj)
  saveRDS(eta_table_n, file.path(out_dir_n, "eta_table_valid.rds"))

  ## --- (c) iCov_n: same schema simulate_scenario_v2() expects from iCov_full
  ##     (id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)
  iCov_n <- vp_n %>%
    dplyr::mutate(SUBJECT = ID,
                  id      = (DATASET - 1L) * n_subj + ID) %>%
    dplyr::select(id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)

  ## --- (d) Run all 16 scenarios via the existing simulate_scenario_v2()
  sim_list_n <- purrr::map(
    .x = PsN_scenarios$scenario,
    .f = function(scn) {
      message(sprintf("[N=%d] Simulating scenario %02d ...", n_subj, scn))
      out <- simulate_scenario_v2(scn, icov = iCov_n, etas = eta_table_n)
      saveRDS(out, file.path(out_dir_n,
                             sprintf("sim_obs_scenario_%02d.rds", scn)))
      out
    },
    .progress = TRUE
  )
  names(sim_list_n) <- sprintf("scenario_%02d", PsN_scenarios$scenario)
  sim_all_n <- dplyr::bind_rows(sim_list_n)
  saveRDS(sim_all_n, file.path(out_dir_n, "sim_obs_all_scenarios.rds"))

  invisible(list(virtual_pop = vp_n,
                 eta_table   = eta_table_n,
                 sim_obs_all = sim_all_n,
                 out_dir     = out_dir_n))
}
cohort_40 <- build_sister_cohort(n_subj = 40L, tol = 0.05 * sqrt(300/40))
# Warning messages:
# 1: In stats::cor(samp) : the standard deviation is zero
# 2: In stats::cor(samp) : the standard deviation is zero
# 3: In stats::cor(samp) : the standard deviation is zero
# Error in `if (all(diff_cor[upper.tri(diff_cor)] <= tol)) ...`:
# ! missing value where TRUE/FALSE needed
#     ▆
#  1. └─global build_sister_cohort(n_subj = 40L, tol = 0.05 * sqrt(300/40))
cohort_80 <- build_sister_cohort(n_subj = 80L, tol = 0.05 * sqrt(300/80))
# Warning message:
# In stats::cor(samp) : the standard deviation is zero
# Error in `if (all(diff_cor[upper.tri(diff_cor)] <= tol)) ...`:
# ! missing value where TRUE/FALSE needed
#     ▆
#  1. └─global build_sister_cohort(n_subj = 80L, tol = 0.05 * sqrt(300/80))
build_sister_cohort <- function(n_subj,
                                n_datasets = 250L,
                                tol_scale  = 1.0,        # multiply baseline tol
                                tol_base   = 0.05,
                                pop_seed   = 4242L + n_subj) {

  out_dir_n <- sprintf("simulated_virtual_dataset_eta_filtered_N%d", n_subj)
  if (!dir.exists(out_dir_n)) dir.create(out_dir_n, recursive = TRUE)

  ## ---- (a) Stratified bootstrap by SEX x RACE ----------------------------
  set.seed(pop_seed)

  # Target counts per (SEX, RACE) cell, proportional to pop
  cell_props <- pop |>
    dplyr::count(SEX, RACE) |>
    dplyr::mutate(p = n / sum(n),
                  n_target = round(p * n_subj))

  # Fix rounding drift so cells sum to n_subj
  diff <- n_subj - sum(cell_props$n_target)
  if (diff != 0) {
    cell_props$n_target[which.max(cell_props$p)] <-
      cell_props$n_target[which.max(cell_props$p)] + diff
  }

  pop_split <- split(pop, list(pop$SEX, pop$RACE), drop = TRUE)
  cell_key  <- paste(cell_props$SEX, cell_props$RACE, sep = ".")

  draw_one <- function() {
    purrr::map2_dfr(cell_key, cell_props$n_target, function(k, n_k) {
      d <- pop_split[[k]]
      d[sample(nrow(d), n_k, replace = TRUE), , drop = FALSE]
    })
  }

  # Scale tol with 1/sqrt(N) (SD of sample corr ~ 1/sqrt(N))
  tol <- tol_base * tol_scale * sqrt(300 / n_subj)

  datasets_n <- list();  attempt_n <- 0L;  n_skip <- 0L
  while (length(datasets_n) < n_datasets) {
    attempt_n <- attempt_n + 1L
    samp <- draw_one()

    # Safety net: skip any draw with degenerate columns
    col_sd <- vapply(samp, stats::sd, numeric(1))
    if (any(is.na(col_sd)) || any(col_sd == 0)) { n_skip <- n_skip + 1L; next }

    u <- abs(stats::cor(samp) - ref_cor)[upper.tri(ref_cor)]
    if (!anyNA(u) && all(u <= tol)) {
      datasets_n[[length(datasets_n) + 1L]] <- samp
    }
    if (attempt_n > 1e5L) stop(sprintf(
      "Tolerance too strict at N = %d (kept %d / tried %d, %d skipped)",
      n_subj, length(datasets_n), attempt_n, n_skip))
  }
  message(sprintf("[N=%d] tol=%.3f  kept %d / tried %d (%d degenerate skipped)",
                  n_subj, tol, length(datasets_n), attempt_n, n_skip))

  for (i in seq_along(datasets_n)) {
    datasets_n[[i]]$ID      <- seq_len(n_subj)
    datasets_n[[i]]$DATASET <- i
  }
  vp_n <- dplyr::bind_rows(datasets_n)
  saveRDS(vp_n, file.path(out_dir_n,
                          sprintf("virtual_population_250x%d.rds", n_subj)))

  ## ---- (b) eta table (CL/Vc covariance preserved here) -------------------
  eta_table_n <- build_eta_table(n_datasets = n_datasets, n_subj = n_subj)
  saveRDS(eta_table_n, file.path(out_dir_n, "eta_table_valid.rds"))

  ## ---- (c) iCov_n --------------------------------------------------------
  iCov_n <- vp_n |>
    dplyr::mutate(SUBJECT = ID,
                  id      = (DATASET - 1L) * n_subj + ID) |>
    dplyr::select(id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)

  ## ---- (d) 16 scenarios --------------------------------------------------
  sim_list_n <- purrr::map(PsN_scenarios$scenario, function(scn) {
    message(sprintf("[N=%d] scenario %02d ...", n_subj, scn))
    out <- simulate_scenario_v2(scn, icov = iCov_n, etas = eta_table_n)
    saveRDS(out, file.path(out_dir_n,
                           sprintf("sim_obs_scenario_%02d.rds", scn)))
    out
  }, .progress = TRUE)
  names(sim_list_n) <- sprintf("scenario_%02d", PsN_scenarios$scenario)
  sim_all_n <- dplyr::bind_rows(sim_list_n)
  saveRDS(sim_all_n, file.path(out_dir_n, "sim_obs_all_scenarios.rds"))

  invisible(list(virtual_pop = vp_n,
                 eta_table   = eta_table_n,
                 sim_obs_all = sim_all_n,
                 out_dir     = out_dir_n))
}

cohort_40 <- build_sister_cohort(n_subj = 40L)   # tol ≈ 0.137 internally
cohort_80 <- build_sister_cohort(n_subj = 80L)   # tol ≈ 0.097 internally
# [N=40] tol=0.137  kept 250 / tried 2035 (0 degenerate skipped)
# [N=40] scenario 01 ...
# Error in `purrr::map()`:
# ℹ In index: 1.
# Caused by error:
# ! the 'id' in the iCov must have 1 unique match to the event table
#      ▆
#   1. ├─global build_sister_cohort(n_subj = 40L)
#   2. │ └─purrr::map(...)
#   3. │   └─purrr:::map_("list", .x, .f, ..., .progress = .progress)
#   4. │     ├─purrr:::with_indexed_errors(...)
#   5. │     │ └─base::withCallingHandlers(...)
#   6. │     ├─purrr:::call_with_cleanup(...)
#   7. │     └─.f(.x[[i]], ...)
#   8. │       └─global simulate_scenario_v2(scn, icov = iCov_n, etas = eta_table_n)
#   9. │         ├─rxode2::rxSolve(...)
#  10. │         └─rxode2:::rxSolve.default(...)
#  11. │           └─base::tryCatch(...)
#  12. │             └─base (local) tryCatchList(expr, classes, parentenv, handlers)
#  13. │               └─base (local) tryCatchOne(expr, names, parentenv, handlers[[1L]])
#  14. │                 └─value[[3L]](cond)
#  15. │                   └─base::stop(e)
#  16. └─purrr (local) `<fn>`(`<Rcpp::xc>`)
#  17.   └─cli::cli_abort(...)
#  18.     └─rlang::abort(...)
build_sister_cohort <- function(n_subj,
                                n_datasets   = 250L,
                                tol_base     = 0.05,
                                tol_scale    = 1.0,
                                max_attempts = 1e5L,
                                pop_seed     = 4242L + n_subj) {

  out_dir_n <- sprintf("simulated_virtual_dataset_eta_filtered_N%d", n_subj)
  if (!dir.exists(out_dir_n)) dir.create(out_dir_n, recursive = TRUE)

  ## ---- (a) Stratified bootstrap by SEX x RACE -----------------------------
  ##   Fixing the per-cell counts (instead of only SEX totals) removes the
  ##   zero-variance-in-RACE failure mode that crashes pure bootstrap at
  ##   small N, and leaves only within-cell continuous-covariate noise to be
  ##   controlled by the correlation gate below.
  set.seed(pop_seed)

  cell_props <- pop %>%
    dplyr::count(SEX, RACE, name = "n_cell") %>%
    dplyr::mutate(p        = n_cell / sum(n_cell),
                  n_target = round(p * n_subj))

  # Repair rounding drift so cells sum exactly to n_subj
  drift <- n_subj - sum(cell_props$n_target)
  if (drift != 0L) {
    bump <- which.max(cell_props$p)
    cell_props$n_target[bump] <- cell_props$n_target[bump] + drift
  }

  pop_split <- split(pop, interaction(pop$SEX, pop$RACE, drop = TRUE))
  cell_key  <- paste(cell_props$SEX, cell_props$RACE, sep = ".")

  draw_one <- function() {
    purrr::map2_dfr(
      pop_split[cell_key], cell_props$n_target,
      function(d, k) dplyr::slice_sample(d, n = k, replace = TRUE)
    )
  }

  # SD of a sample correlation ~ 1/sqrt(N); scale the gate to match.
  tol <- tol_base * tol_scale * sqrt(300 / n_subj)

  datasets_n <- vector("list", n_datasets)
  n_kept     <- 0L
  attempt_n  <- 0L
  n_skip     <- 0L

  while (n_kept < n_datasets) {
    attempt_n <- attempt_n + 1L
    samp <- draw_one()

    col_sd <- vapply(samp, stats::sd, numeric(1))
    if (any(is.na(col_sd)) || any(col_sd == 0)) {
      n_skip <- n_skip + 1L
    } else {
      u <- abs(stats::cor(samp) - ref_cor)[upper.tri(ref_cor)]
      if (!anyNA(u) && all(u <= tol)) {
        n_kept <- n_kept + 1L
        datasets_n[[n_kept]] <- samp
      }
    }
    if (attempt_n > max_attempts) {
      stop(sprintf(
        "Tolerance too strict at N = %d (kept %d / tried %d, %d skipped)",
        n_subj, n_kept, attempt_n, n_skip
      ))
    }
  }
  message(sprintf(
    "[N=%d] tol=%.3f  kept %d / tried %d (%d degenerate skipped)",
    n_subj, tol, n_kept, attempt_n, n_skip
  ))

  # Tag ID + DATASET and stack into one long tibble
  vp_n <- datasets_n %>%
    purrr::map(function(d) dplyr::mutate(d, ID = seq_len(n_subj))) %>%
    dplyr::bind_rows(.id = "DATASET") %>%
    dplyr::mutate(DATASET = as.integer(DATASET))

  saveRDS(
    vp_n,
    file.path(out_dir_n,
              sprintf("virtual_population_%dx%d.rds", n_datasets, n_subj))
  )

  ## ---- (b) eta table: re-use build_eta_table() at requested n_subj -------
  ##   Preserves CL/Vc covariance (rho in [ETA_RHO_LOW, ETA_RHO_HIGH]).
  eta_table_n <- build_eta_table(n_datasets = n_datasets, n_subj = n_subj)
  saveRDS(eta_table_n, file.path(out_dir_n, "eta_table_valid.rds"))

  ## ---- (c) iCov_n: schema simulate_scenario_v2() expects -----------------
  iCov_n <- vp_n %>%
    dplyr::mutate(SUBJECT = ID,
                  id      = (DATASET - 1L) * n_subj + ID) %>%
    dplyr::select(id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)

  ## ---- (c2) Cohort-sized event table -------------------------------------
  ##   The default `ev` is ev_one replicated 75000x for N=300; it does not
  ##   match cohort-sized iCov, hence the rxSolve "1 unique match" error.
  ##   Rebuild it from the same ev_one template at the cohort's row count.
  n_total_n <- n_datasets * n_subj
  ev_n <- ev_one %>%
    dplyr::slice(rep(dplyr::row_number(), n_total_n)) %>%
    dplyr::mutate(id = rep(seq_len(n_total_n), each = nrow(ev_one))) %>%
    dplyr::arrange(id, time)

  stopifnot(
    nrow(iCov_n) == n_total_n,
    setequal(iCov_n$id, unique(ev_n$id))
  )

  ## ---- (d) Run all 16 scenarios via the existing simulator ---------------
  sim_list_n <- purrr::map(
    .x = PsN_scenarios$scenario,
    .f = function(scn) {
      message(sprintf("[N=%d] Simulating scenario %02d ...", n_subj, scn))
      out <- simulate_scenario_v2(
        scn,
        icov        = iCov_n,
        etas        = eta_table_n,
        event_table = ev_n
      )
      saveRDS(
        out,
        file.path(out_dir_n,
                  sprintf("sim_obs_scenario_%02d.rds", scn))
      )
      out
    },
    .progress = TRUE
  ) %>%
    purrr::set_names(sprintf("scenario_%02d", PsN_scenarios$scenario))

  sim_all_n <- dplyr::bind_rows(sim_list_n)
  saveRDS(sim_all_n, file.path(out_dir_n, "sim_obs_all_scenarios.rds"))

  invisible(list(
    virtual_pop = vp_n,
    eta_table   = eta_table_n,
    sim_obs_all = sim_all_n,
    out_dir     = out_dir_n
  ))
}
cohort_40 <- build_sister_cohort(n_subj = 40L)
# [N=40] tol=0.137  kept 250 / tried 2035 (0 degenerate skipped)
# [N=40] Simulating scenario 01 ...
# [N=40] Simulating scenario 02 ...   6% |  ETA: 18s                              
# [N=40] Simulating scenario 03 ...  12% |  ETA: 18s                              
# [N=40] Simulating scenario 04 ...  19% |  ETA: 16s                              
# [N=40] Simulating scenario 05 ...  25% |  ETA: 16s                              
# [N=40] Simulating scenario 06 ...  31% |  ETA: 15s                              
# [N=40] Simulating scenario 07 ...  38% |  ETA: 14s                              
# [N=40] Simulating scenario 08 ...  44% |  ETA: 13s                              
# [N=40] Simulating scenario 09 ...  50% |  ETA: 12s                              
# [N=40] Simulating scenario 10 ...  56% |  ETA: 11s                              
# [N=40] Simulating scenario 11 ...  62% |  ETA:  9s                              
# [N=40] Simulating scenario 12 ...  69% |  ETA:  8s                              
# [N=40] Simulating scenario 13 ...  75% |  ETA:  6s                              
# [N=40] Simulating scenario 14 ...  81% |  ETA:  5s                              
# [N=40] Simulating scenario 15 ...  88% |  ETA:  3s                              
# [N=40] Simulating scenario 16 ...  94% |  ETA:  2s                              
#                                                                                 
cohort_80 <- build_sister_cohort(n_subj = 80L)
# [N=80] tol=0.097  kept 250 / tried 1809 (0 degenerate skipped)
# [N=80] Simulating scenario 01 ...
# [N=80] Simulating scenario 02 ...   6% |  ETA:  1m                              
# [N=80] Simulating scenario 03 ...  12% |  ETA: 48s                              
# [N=80] Simulating scenario 04 ...  19% |  ETA: 48s                              
# [N=80] Simulating scenario 05 ...  25% |  ETA: 47s                              
# [N=80] Simulating scenario 06 ...  31% |  ETA: 42s                              
# [N=80] Simulating scenario 07 ...  38% |  ETA: 39s                              
# [N=80] Simulating scenario 08 ...  44% |  ETA: 34s                              
# [N=80] Simulating scenario 09 ...  50% |  ETA: 31s                              
# [N=80] Simulating scenario 10 ...  56% |  ETA: 27s                              
# [N=80] Simulating scenario 11 ...  62% |  ETA: 23s                              
# [N=80] Simulating scenario 12 ...  69% |  ETA: 19s                              
# [N=80] Simulating scenario 13 ...  75% |  ETA: 16s                              
# [N=80] Simulating scenario 14 ...  81% |  ETA: 12s                              
# [N=80] Simulating scenario 15 ...  88% |  ETA:  8s                              
# [N=80] Simulating scenario 16 ...  94% |  ETA:  4s                              
#                                                                                 
# R 4.5.3 offline. Waiting to reconnect.
# R 4.5.3 reconnected.
out_dir <- "simulated_virtual_dataset"
PsN_scenarios <- data.frame(
  scenario = 1:16,
  I_BW_CL = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
  I_CRCL_CL = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
  I_BW_VC = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
  I_SEX_VC = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
)
true_params_long <- function(scenarios = PsN_scenarios,
                             theta = c(TVCL = 0.6, TVQ = 1.8, TVVc = 20,
                                       TVVp = 80, TVKA = 0.7,
                                       TH_BW_CL = 0.75, TH_CRCL_CL = 0.5,
                                       TH_BW_VC = 1.0,
                                       ## nlmixr2scm log-additive convention:
                                       ## Vc = TVVc * exp(TH_SEX_VC * SEX).
                                       ## log(1.5) ~ 0.4055 preserves the 1.5x
                                       ## female/male Vc effect of the original
                                       ## PsN-convention design (theta = 0.5).
                                       TH_SEX_VC = log(1.5)),
                             omega = c(var_cl = 0.1, var_vc = 0.1,
                                       cov_cl_vc = 0.02),
                             prop_err = 0.1) {
  scenarios %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      TVCL       = theta[["TVCL"]],
      TVVc       = theta[["TVVc"]],
      TVQ        = theta[["TVQ"]],
      TVVp       = theta[["TVVp"]],
      TVKA       = theta[["TVKA"]],
      CLBW     = theta[["TH_BW_CL"]]   * I_BW_CL,
      CLcrCL   = theta[["TH_CRCL_CL"]] * I_CRCL_CL,
      VcBW     = theta[["TH_BW_VC"]]   * I_BW_VC,
      VcSEX    = theta[["TH_SEX_VC"]]  * I_SEX_VC,
      var_CL   = omega[["var_cl"]], #Variance
      var_Vc   = omega[["var_vc"]], #variance
      cov_VcCL = omega[["cov_cl_vc"]],#covariance
      ResErr   = prop_err #standard deviation
    ) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(
      cols      = c(TVCL, TVVc, TVQ, TVVp, TVKA,
                    CLBW, CLcrCL, VcBW, VcSEX,
                    var_CL, var_Vc, cov_VcCL, ResErr),
      names_to  = "parameter",
      values_to = "true_value"
    ) %>%
    dplyr::select(scenario, parameter, true_value,
                  I_BW_CL, I_CRCL_CL, I_BW_VC, I_SEX_VC)
}
true_params <- true_params_long()
saveRDS(true_params, file.path(out_dir, "true_params_long.rds"))
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
library(nlmixr2scm)
devtools::load_all("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
# ℹ Loading nlmixr2scm
base_2cmt_oral <- function() {
  ini({
    lTVCL <- log(0.6)
    lTVQ  <- log(1.8)
    lTVVc <- log(20)
    lTVVp <- log(80)
    ## KA fixed: the design's first non-zero sample (~7 h) is well past
    ## absorption (KA = 0.7 /h, absorption t1/2 ~ 1 h), so KA cannot be
    ## informatively estimated from the data.  Estimating it inflates
    ## the cov-matrix condition number by ~20x (raw CN ~35000 -> ~1800)
    ## without changing OFV or the SCM decisions.  Fixing matches the
    ## warfarin / Khandelwal convention.
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

    k10 <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)      = -ka * depot
    d/dt(central)    =  ka * depot - k10 * central - k12 * central + k21 * peripheral
    d/dt(peripheral) =  k12 * central - k21 * peripheral

    cp = central / vc
    cp ~ prop(prop.err)
  })
}
t_fit_base <- system.time(
  fit_base <- nlmixr2(base_2cmt_oral, ds01,
                      est = "focei", control = scm_focei)
)
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities                                                        
# → calculate ∂(f)/∂(η)                                                            
# → calculate ∂(R²)/∂(η)                                                           
# → finding duplicate expressions in inner model...                                
# → optimizing duplicate expressions in inner model...                             
# → finding duplicate expressions in EBE model...                                  
# → optimizing duplicate expressions in EBE model...                               
# → compiling inner model...                                                       
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...                                
# → compiling EBE model...                                                         
# ✔ done
# → compiling events FD model...
# ✔ done
# done
package_scm_result <- function(label, scm_res, runtime_sec,
                               true_long = true_params, scenario_id = 9,
                               refit_control = NULL) {
  ## Pick the final fit.  runSCM() returns scm_res$resFwd / scm_res$resBck as
  ## UNNAMED 3-element lists: [[1]] = nlmixr2FitCore fit, [[2]] = step table,
  ## [[3]] = final-selection data.frame.  The previous code looked for
  ## `$finalFit` which never exists, so final_fit was always NULL.  Backward
  ## takes precedence because it represents the post-pruning model.
  .pickFit <- function(x) {
    if (is.null(x)) return(NULL)
    cand <- if (is.list(x) && length(x) >= 1L) x[[1L]] else x
    if (inherits(cand, "nlmixr2FitCore")) cand else NULL
  }
  final_fit <- .pickFit(scm_res$resBck)
  if (is.null(final_fit)) final_fit <- .pickFit(scm_res$resFwd)

  ## Optional: refit final model with full diagnostics (cov, IPRED tables).
  ## When the refit succeeds it REPLACES final_fit so every downstream
  ## extraction (estimates, rel_err, diag) reflects the diagnostic fit.
  final_fit_refit <- NULL
  if (!is.null(final_fit) && !is.null(refit_control)) {
    final_fit_refit <- refit_final_model(final_fit, control = refit_control)
    if (!is.null(final_fit_refit)) final_fit <- final_fit_refit
  }

  ## Selected pairs from summaryTable.  The decision column is `included`
  ## with values "yes"/"no" (forward) and "retained"/"dropped" (backward);
  ## the old code filtered on inFinal/accepted/kept which never exist, so
  ## `selected` collapsed back to the full step_hist.
  selected <- if (!is.null(scm_res$summaryTable)) {
    st <- as.data.frame(scm_res$summaryTable)
    if ("included" %in% colnames(st)) {
      st[st$included %in% c("yes", "retained"), , drop = FALSE]
    } else {
      st
    }
  } else {
    NULL
  }

  final_est <- if (!is.null(final_fit)) extract_params_long(final_fit) else NULL
  rel_err   <- if (!is.null(final_est)) {
    rel_err_one(final_est, true_long, scenario_id)
  } else NULL
  diag      <- if (!is.null(final_fit)) diagnose_fit(final_fit) else NULL
  parFixed  <- if (!is.null(final_fit_refit)) final_fit_refit$parFixedDf else NULL

  list(
    label           = label,
    selected        = selected,
    step_hist       = scm_res$summaryTable,
    final_est       = final_est,
    rel_err         = rel_err,
    diag            = diag,
    parFixed        = parFixed,
    final_fit_refit = final_fit_refit,
    runtime_sec     = runtime_sec,
    raw             = scm_res
  )
}
out_dir <- "simulated_virtual_dataset"
stage1_dir <- file.path(out_dir_v2, "stage1_smoke_scn09_ds01")
ds01 <- readRDS(file.path(stage1_dir, "nm_scn09_ds01.rds"))
runSCM_traced <- function(label, ...) {
  dots      <- list(...)
  n_pairs   <- if (!is.null(dots$pairsVec))   length(dots$pairsVec) else NA
  n_workers <- if (!is.null(dots$workers))    dots$workers          else 1L
  search    <- if (!is.null(dots$searchType)) dots$searchType       else "scm"

  cat(sprintf(
    "\n=== [%s] %s starting %-8s | %d candidate(s), %d worker(s) ===\n",
    label, format(Sys.time(), "%H:%M:%S"), search, n_pairs, n_workers
  ), file = stderr())

  t0  <- Sys.time()
  res <- withCallingHandlers(
    do.call(nlmixr2scm::runSCM, dots),
    message = function(m) {
      txt <- conditionMessage(m)
      ## Tag only step/candidate boundary messages; let cli rules pass through.
      if (grepl("step\\s+\\d+,\\s*candidate|forward search|backward search",
                txt, ignore.case = TRUE)) {
        elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("[%s | +%5.1f min] ",
                    format(Sys.time(), "%H:%M:%S"), elapsed_min),
            file = stderr())
      }
      ## Don't muffle -- let the original cli message print as usual.
    }
  )
  elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  cat(sprintf(
    "=== [%s] %s DONE | elapsed %.1f min (%.0f s) ===\n\n",
    label, format(Sys.time(), "%H:%M:%S"),
    elapsed_s / 60, elapsed_s
  ), file = stderr())

  attr(res, "elapsed_s") <- elapsed_s
  res
}
scm_focei_maxiteration <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl          = rxode2::rxControl(atol = 1e-8, rtol = 1e-6),
  covMethod  = ""         # SCM doesn't need cov matrix for LRT
)
res_fwd_auto_maxiteration <- runSCM_traced(
  label       = "forward",
  data        = ds01,
  fit         = fit_base,
  varsVec    = c("cl", "vc"),
  covarsVec  = "BW",
  catvarsVec = "SEX",
  shapes     = c("power", "lin"),
  searchType  = "forward",
  control     = scm_focei_maxiteration,    # slim control: no tables, no cov, sigdig=4 bobyqa
  saveModels  = FALSE,
  workers     = 3L,           # 4 cores: leave 1 free for OS / Positron
  print       = 100,           # FOCEi iteration progress every 100 iters
  maxRetries = 1L #no retries for this smoke test
)
# 
# === [forward] 11:36:58 starting forward  | NA candidate(s), 3 worker(s) ===
# ── SCM Summary ─────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -16722.83
# ℹ Base model params : 9
# ℹ Search type : forward
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (cl, vc)
# ℹ Covariates : 3
# (BW_power, BW_lin, SEX_1)
# ℹ Total candidates : 6
# ── Categorical covariates ──────────────────────────────────────────────────────────────────────────────────────────────
# ℹ SEX: reference = '0'
# Indicators (1): SEX_1
# ── Relationships to test ───────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ cl [power]
# 2. BW_lin ~ cl [lin]
# 3. BW_power ~ vc [power]
# 4. BW_lin ~ vc [lin]
# 5. SEX_1 ~ cl [cat]
# 6. SEX_1 ~ vc [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Proceed with SCM? [y/N]: 
# 
# [11:37:05 | +  0.1 min] 
# ── starting forward search... ──────────────────────────────────────────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:43:40 | +  6.7 min] 
# → Forward step 1, candidate 1/6: BW_power ~ cl [power]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
#  
#  
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:43:40 | +  6.7 min] 
# → Forward step 1, candidate 2/6: BW_lin ~ cl [lin]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:43:40 | +  6.7 min] 
# → Forward step 1, candidate 3/6: BW_power ~ vc [power]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
#  
#  
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:43:40 | +  6.7 min] 
# → Forward step 1, candidate 4/6: BW_lin ~ vc [lin]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:44:28 | +  7.5 min] 
# → Forward step 1, candidate 5/6: SEX_1 ~ cl [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:44:28 | +  7.5 min] 
# → Forward step 1, candidate 6/6: SEX_1 ~ vc [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
#  
#  
# 
# ── best model at step 1: ───────────────────────────────────────────────────────────────────────────────────────────────
#        step    covar var shape      objf deltObjf       AIC      BIC numParams  qchisqr pchisqr included searchType
# eta.cl    1 BW_power  cl power -16805.34 82.51189 -13479.16 -13429.7        10 3.841459       0      yes    forward
#               covNames covarEffect bsvReduction
# eta.cl cov_BW_power_cl   0.7743659     25.56932
# ℹ Dropping alternative shape(s) for BW~cl: BW_lin
# 
# ── accepted BW_power~cl ──
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:46:43 | +  9.7 min] 
# → Forward step 2, candidate 1/4: BW_power ~ vc [power]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:50:10 | + 13.2 min] 
# → Forward step 2, candidate 2/4: BW_lin ~ vc [lin]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
#  
#  
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:50:10 | + 13.2 min] 
# → Forward step 2, candidate 3/4: SEX_1 ~ cl [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
#  
#  
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:50:10 | + 13.2 min] 
# → Forward step 2, candidate 4/4: SEX_1 ~ vc [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# [11:50:10 | + 13.2 min] 
# ── OFV did not improve, exiting forward search ... ─────────────────────────────────────────────────────────────────────
# 
# [11:50:10 | + 13.2 min] 
# ── forward search complete ──
# 
# ── SCM Step Summary ────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~cl [power]  -16722.830  -16805.342    82.512     0.0000  Added
# Forward    2     BW_power~vc [power]  -16805.342  -16806.748     1.407     0.2356  Not selected
# ── SCM All Candidates ──────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~cl [power]  -16722.830  -16805.342    82.512     0.0000  Added
# Forward    1     BW_lin~cl [lin]      -16722.830  -16740.620    17.790     0.0000  Not selected
# Forward    1     SEX_1~cl [cat]       -16722.830  -16729.462     6.632     0.0100  Not selected
# Forward    1     SEX_1~vc [cat]       -16722.830  -16723.000     0.171     0.6797  Not selected
# Forward    1     BW_lin~vc [lin]      -16722.830  -16722.865     0.035     0.8508  Not selected
# Forward    1     BW_power~vc [power]  -16722.830  -16722.804    -0.025     1.0000  Not selected
# 
# Forward    2     BW_power~vc [power]  -16805.342  -16806.748     1.407     0.2356  Not selected
# Forward    2     SEX_1~vc [cat]       -16805.342  -16806.006     0.665     0.4149  Not selected
# Forward    2     BW_lin~vc [lin]      -16805.342  -16795.536    -9.806     1.0000  Not selected
# Forward    2     SEX_1~cl [cat]       -16805.342  -16805.250    -0.092     1.0000  Not selected
# ── Final model ─────────────────────────────────────────────────────────────────────────────────────────────────────────
# Retained:
# BW_power~cl [power]
# ℹ Final model OFV: -16805.342
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [forward] 11:50:10 DONE | elapsed 13.2 min (792 s) ===
# 
# Warning messages:
# 1: ! BW_power ~ cl: unrealistic OFV on attempt 1/2: p-value underflow (dOFV = 82.352).
# ℹ Retrying with perturbed init. 
# 2: ! BW_power ~ cl: unrealistic OFV after all 2 attempts: p-value underflow (dOFV = 82.512). Accepting best available
#   result (attempt 2/2, dOFV = 82.512). 
# 3: ! BW_power ~ vc: unrealistic OFV on attempt 1/2: OFV increased vs parent (-16722.804 > -16722.83).
# ℹ Retrying with perturbed init. 
# 4: ! BW_power ~ vc: unrealistic OFV after all 2 attempts: OFV increased vs parent (-16722.777 > -16722.83). Accepting best
#   available result (attempt 1/2, dOFV = -0.025). 
# 5: ! SEX_1 ~ vc: unrealistic OFV on attempt 1/2: OFV increased vs parent (-16722.801 > -16722.83).
# ℹ Retrying with perturbed init. 
# 6: ! BW_lin ~ vc: unrealistic OFV on attempt 1/2: OFV increased vs parent (-16795.536 > -16805.342).
# ℹ Retrying with perturbed init. 
# 7: ! BW_lin ~ vc: unrealistic OFV after all 2 attempts: OFV increased vs parent (-16795.525 > -16805.342). Accepting best
#   available result (attempt 1/2, dOFV = -9.806). 
# 8: ! SEX_1 ~ cl: unrealistic OFV on attempt 1/2: OFV increased vs parent (-16805.188 > -16805.342).
# ℹ Retrying with perturbed init. 
# 9: ! SEX_1 ~ cl: unrealistic OFV after all 2 attempts: OFV increased vs parent (-16805.25 > -16805.342). Accepting best
#   available result (attempt 2/2, dOFV = -0.092). 
print(bench_results_300)
# # A tibble: 12 × 12
#    model  setting      elapsed_min final_ofv best_ofv ofv_gap n_iter bad_solves hessian_reset not_at_min converged_well error
#    <chr>  <chr>              <dbl>     <dbl>    <dbl>   <dbl>  <int> <lgl>      <lgl>         <lgl>      <lgl>          <chr>
#  1 lin    baseline            1.81   -16901.  -16901. 3.15e-3    600 FALSE      FALSE         TRUE       TRUE           NA   
#  2 refexp baseline            1.92   -16901.  -16901. 2.82e-3    540 FALSE      FALSE         TRUE       TRUE           NA   
#  3 lin    fast                1.09   -16901.  -16901. 2.02e-3    576 FALSE      FALSE         FALSE      TRUE           NA   
#  4 refexp fast                1.35   -16901.  -16901. 4.90e-3    498 FALSE      FALSE         TRUE       TRUE           NA   
#  5 lin    no_sticky           1.54   -16901.  -16901. 3.15e-3    600 FALSE      FALSE         TRUE       TRUE           NA   
#  6 refexp no_sticky           1.94   -16901.  -16901. 2.82e-3    540 FALSE      FALSE         TRUE       TRUE           NA   
#  7 lin    tight_ode           2.22   -16901.  -16901. 9.14e-4    747 FALSE      FALSE         TRUE       TRUE           NA   
#  8 refexp tight_ode           2.39   -16901.  -16901. 1.56e-3    579 FALSE      FALSE         TRUE       TRUE           NA   
#  9 lin    tight_no_st…        2.27   -16901.  -16901. 9.14e-4    747 FALSE      FALSE         TRUE       TRUE           NA   
# 10 refexp tight_no_st…        2.38   -16901.  -16901. 1.56e-3    579 FALSE      FALSE         TRUE       TRUE           NA   
# 11 lin    fast_no_sti…        1.07   -16901.  -16901. 2.02e-3    576 FALSE      FALSE         FALSE      TRUE           NA   
# 12 refexp fast_no_sti…        1.31   -16901.  -16901. 4.90e-3    498 FALSE      FALSE         TRUE       TRUE           NA   
res_fwd_auto <- readRDS(file.path(stage1_dir, "res_fwd_auto.rds"))
res_fwd_auto
# $summaryTable
#          step    covar var shape      objf      deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType        covNames
# eta.cl      1 BW_power  cl power -16805.30  8.246967e+01 -13479.12 -13429.66        10 3.841459 0.000000e+00      yes    forward cov_BW_power_cl
# eta.cl1     1   BW_lin  cl   lin -16740.62  1.778961e+01 -13414.44 -13364.98        10 3.841459 2.467274e-05       no    forward   cov_BW_lin_cl
# eta.vc      1 BW_power  vc power -16722.77 -5.664084e-02 -13396.59 -13347.13        10 3.841459 1.000000e+00       no    forward cov_BW_power_vc
# eta.vc1     1   BW_lin  vc   lin -16722.86  3.230686e-02 -13396.68 -13347.22        10 3.841459 8.573558e-01       no    forward   cov_BW_lin_vc
# eta.cl2     1    SEX_1  cl   cat -16476.53 -2.462966e+02 -13150.35 -13100.89        10 3.841459 1.000000e+00       no    forward    cov_SEX_1_cl
# eta.vc2     1    SEX_1  vc   cat -16722.83 -4.848361e-03 -13396.65 -13347.19        10 3.841459 1.000000e+00       no    forward    cov_SEX_1_vc
# eta.vc3     2 BW_power  vc power -16546.98 -2.583218e+02 -13218.80 -13163.84        11 3.841459 1.000000e+00       no    forward cov_BW_power_vc
# eta.vc11    2   BW_lin  vc   lin -16795.53 -9.772491e+00 -13467.35 -13412.39        11 3.841459 1.000000e+00       no    forward   cov_BW_lin_vc
# eta.cl3     2    SEX_1  cl   cat -16546.28 -2.590156e+02 -13218.11 -13163.15        11 3.841459 1.000000e+00       no    forward    cov_SEX_1_cl
# eta.vc21    2    SEX_1  vc   cat -16547.14 -2.581562e+02 -13218.96 -13164.01        11 3.841459 1.000000e+00       no    forward    cov_SEX_1_vc
#            covarEffect  bsvReduction
# eta.cl    0.7711465051  2.478980e+01
# eta.cl1   0.0012500000 -3.381735e-13
# eta.vc   -0.0197720250 -3.524291e-01
# eta.vc1  -0.0005022676  5.278545e-09
# eta.cl2   0.0894402138  2.609573e+00
# eta.vc2   0.0214849318 -7.610753e-01
# eta.vc3   0.1248087956  7.898341e-01
# eta.vc11  0.0012500000  3.616738e-01
# eta.cl3   0.0140360840  1.913167e-02
# eta.vc21  0.0330576561 -1.064878e+00
# 
# $resFwd
# $resFwd[[1]]
# ── nlmixr² FOCEi (outer: bobyqa) ──
# 
#           OBJF       AIC       BIC Log-likelihood
# FOCEi -16805.3 -13479.12 -13429.66        6748.56
# 
# ── Time (sec $time): ──
# 
#             setup optimize preprocess postprocess     other
# elapsed 0.0616791 133.6713       0.11        0.06 0.3270012
# 
# ── Population Parameters ($parFixed or $parFixedDf): ──
# 
#                    Est. Back-transformed BSV(CV%) Shrink(SD)%
# lTVCL           -0.4184           0.6581     32.2      1.45% 
# lTVQ             0.6166            1.853                     
# lTVVc             3.056            21.24     31.6      13.0% 
# lTVVp             4.382               80                     
# lTVKA           -0.3567          -0.3567                     
# prop.err         0.1012           0.1012                     
# cov_BW_power_cl  0.7711           0.7711                     
#  
#   Correlations in between subject variability (BSV) matrix:
#     cor:eta.vc,eta.cl 
#            0.174   
#  
# 
#   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#   Information about run found ($runInfo):
#    • last objective function was not at minimum, possible problems in optimization 
#    • ETAs were reset to zero during optimization; (Can control by foceiControl(resetEtaP=.)) 
#    • initial ETAs were nudged; (can control by foceiControl(etaNudge=., etaNudge2=)) 
#   Censoring ($censInformation): No censoring
#   Minimization message ($message):  
#     Normal exit from bobyqa 
# 
# $resFwd[[2]]
#          step    covar var shape      objf      deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType        covNames
# eta.cl      1 BW_power  cl power -16805.30  8.246967e+01 -13479.12 -13429.66        10 3.841459 0.000000e+00      yes    forward cov_BW_power_cl
# eta.cl1     1   BW_lin  cl   lin -16740.62  1.778961e+01 -13414.44 -13364.98        10 3.841459 2.467274e-05       no    forward   cov_BW_lin_cl
# eta.vc      1 BW_power  vc power -16722.77 -5.664084e-02 -13396.59 -13347.13        10 3.841459 1.000000e+00       no    forward cov_BW_power_vc
# eta.vc1     1   BW_lin  vc   lin -16722.86  3.230686e-02 -13396.68 -13347.22        10 3.841459 8.573558e-01       no    forward   cov_BW_lin_vc
# eta.cl2     1    SEX_1  cl   cat -16476.53 -2.462966e+02 -13150.35 -13100.89        10 3.841459 1.000000e+00       no    forward    cov_SEX_1_cl
# eta.vc2     1    SEX_1  vc   cat -16722.83 -4.848361e-03 -13396.65 -13347.19        10 3.841459 1.000000e+00       no    forward    cov_SEX_1_vc
# eta.vc3     2 BW_power  vc power -16546.98 -2.583218e+02 -13218.80 -13163.84        11 3.841459 1.000000e+00       no    forward cov_BW_power_vc
# eta.vc11    2   BW_lin  vc   lin -16795.53 -9.772491e+00 -13467.35 -13412.39        11 3.841459 1.000000e+00       no    forward   cov_BW_lin_vc
# eta.cl3     2    SEX_1  cl   cat -16546.28 -2.590156e+02 -13218.11 -13163.15        11 3.841459 1.000000e+00       no    forward    cov_SEX_1_cl
# eta.vc21    2    SEX_1  vc   cat -16547.14 -2.581562e+02 -13218.96 -13164.01        11 3.841459 1.000000e+00       no    forward    cov_SEX_1_vc
#            covarEffect  bsvReduction
# eta.cl    0.7711465051  2.478980e+01
# eta.cl1   0.0012500000 -3.381735e-13
# eta.vc   -0.0197720250 -3.524291e-01
# eta.vc1  -0.0005022676  5.278545e-09
# eta.cl2   0.0894402138  2.609573e+00
# eta.vc2   0.0214849318 -7.610753e-01
# eta.vc3   0.1248087956  7.898341e-01
# eta.vc11  0.0012500000  3.616738e-01
# eta.cl3   0.0140360840  1.913167e-02
# eta.vc21  0.0330576561 -1.064878e+00
# 
# $resFwd[[3]]
#   var    covar       type center raw_col level has_missing missing_fill missing_check shape    covExpr      init lower upper
# 1  cl BW_power continuous     80      BW  <NA>       FALSE            0          <NA> power log(BW/80) 0.7711465    -5     5
# 
# 
# $resBck
# NULL
# 
# attr(,"elapsed_s")
# [1] 263.7854
scm_focei_maxiteration0 <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  covMethod  = ""         # SCM doesn't need cov matrix for LRT
)
res_fwd_auto_maxiteration0 <- runSCM_traced(
  label       = "forward",
  data        = ds01,
  fit         = fit_base,
  varsVec    = c("cl", "vc"),
  covarsVec  = "BW",
  catvarsVec = "SEX",
  shapes     = c("power", "lin"),
  searchType  = "forward",
  control     = scm_focei_maxiteration0,    # slim control: no tables, no cov, sigdig=4 bobyqa
  saveModels  = FALSE,
  workers     = 3L,           # 4 cores: leave 1 free for OS / Positron
  print       = 100,           # FOCEi iteration progress every 100 iters
  maxRetries = 0L #no retries for this smoke test
)
# 
# === [forward] 12:03:49 starting forward  | NA candidate(s), 3 worker(s) ===
# ── SCM Summary ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ Estimation method : focei
# ℹ Control : foceiControl (user-supplied)
# ℹ Base model OFV : -16722.83
# ℹ Base model params : 9
# ℹ Search type : forward
# ℹ p-value fwd / bck : 0.05 / 0.01
# ℹ Output folder : (none — saveModels = FALSE)
# ℹ Parameters : 2
# (cl, vc)
# ℹ Covariates : 3
# (BW_power, BW_lin, SEX_1)
# ℹ Total candidates : 6
# ── Categorical covariates ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ℹ SEX: reference = '0'
# Indicators (1): SEX_1
# ── Relationships to test ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 1. BW_power ~ cl [power]
# 2. BW_lin ~ cl [lin]
# 3. BW_power ~ vc [power]
# 4. BW_lin ~ vc [lin]
# 5. SEX_1 ~ cl [cat]
# 6. SEX_1 ~ vc [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Proceed with SCM? [y/N]: 
# 
# [12:03:52 | +  0.1 min] 
# ── starting forward search... ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [12:06:29 | +  2.7 min] 
# → Forward step 1, candidate 1/6: BW_power ~ cl [power]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [12:06:29 | +  2.7 min] 
# → Forward step 1, candidate 2/6: BW_lin ~ cl [lin]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [12:06:29 | +  2.7 min] 
# → Forward step 1, candidate 3/6: BW_power ~ vc [power]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [12:06:29 | +  2.7 min] 
# → Forward step 1, candidate 4/6: BW_lin ~ vc [lin]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [12:06:40 | +  2.8 min] 
# → Forward step 1, candidate 5/6: SEX_1 ~ cl [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [12:06:40 | +  2.8 min] 
# → Forward step 1, candidate 6/6: SEX_1 ~ vc [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# 
# ── best model at step 1: ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#        step    covar var shape     objf deltObjf       AIC       BIC numParams  qchisqr pchisqr included searchType        covNames covarEffect bsvReduction
# eta.cl    1 BW_power  cl power -16805.3 82.46967 -13479.12 -13429.66        10 3.841459       0      yes    forward cov_BW_power_cl   0.7711465      24.7898
# ℹ Dropping alternative shape(s) for BW~cl: BW_lin
# 
# ── accepted BW_power~cl ──
# 
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [15:11:49 | +188.0 min] 
# → Forward step 2, candidate 1/4: BW_power ~ vc [power]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [15:12:25 | +188.6 min] 
# → Forward step 2, candidate 2/4: BW_lin ~ vc [lin]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [15:12:25 | +188.6 min] 
# → Forward step 2, candidate 3/4: SEX_1 ~ cl [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [15:12:25 | +188.6 min] 
# → Forward step 2, candidate 4/4: SEX_1 ~ vc [cat]
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# [15:12:25 | +188.6 min] 
# ── OFV did not improve, exiting forward search ... ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 
# [15:12:25 | +188.6 min] 
# ── forward search complete ──
# 
# ── SCM Step Summary ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~cl [power]  -16722.830  -16805.300    82.470     0.0000  Added
# Forward    2     BW_power~vc [power]  -16805.300  -16547.123  -258.176     1.0000  Not selected
# ── SCM All Candidates ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Direction  Step  Relation                Ref OFV         OFV      dOFV    p-value  Decision 
# ------------------------------------------------------------------------------------------- 
# Forward    1     BW_power~cl [power]  -16722.830  -16805.300    82.470     0.0000  Added
# Forward    1     BW_lin~cl [lin]      -16722.830  -16740.619    17.790     0.0000  Not selected
# Forward    1     BW_lin~vc [lin]      -16722.830  -16722.862     0.032     0.8574  Not selected
# Forward    1     BW_power~vc [power]  -16722.830  -16722.773    -0.057     1.0000  Not selected
# Forward    1     SEX_1~cl [cat]       -16722.830  -16475.670  -247.160     1.0000  Not selected
# Forward    1     SEX_1~vc [cat]       -16722.830  -16722.825    -0.005     1.0000  Not selected
# 
# Forward    2     BW_power~vc [power]  -16805.300  -16547.123  -258.176     1.0000  Not selected
# Forward    2     BW_lin~vc [lin]      -16805.300  -16795.527    -9.772     1.0000  Not selected
# Forward    2     SEX_1~cl [cat]       -16805.300  -16546.250  -259.050     1.0000  Not selected
# Forward    2     SEX_1~vc [cat]       -16805.300  -16546.991  -258.309     1.0000  Not selected
# ── Final model ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# Retained:
# BW_power~cl [power]
# ℹ Final model OFV: -16805.3
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ✔ Log files written to c:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm: scm_log.txt, scm_step_summary.csv, scm_all_candidates.csv
# === [forward] 15:12:27 DONE | elapsed 188.6 min (11318 s) ===
# 
# Warning messages:
# 1: ! BW_power ~ cl: unrealistic OFV after all 1 attempt: p-value underflow (dOFV = 82.47). Accepting best available result (attempt 1/1, dOFV = 82.47). 
# 2: ! BW_power ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-16722.773 > -16722.83). Accepting best available result (attempt 1/1, dOFV = -0.057). 
# 3: ! SEX_1 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-16475.67 > -16722.83). Accepting best available result (attempt 1/1, dOFV = -247.16). 
# 4: ! SEX_1 ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-16722.825 > -16722.83). Accepting best available result (attempt 1/1, dOFV = -0.005). 
# 5: ! BW_power ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-16547.123 > -16805.3). Accepting best available result (attempt 1/1, dOFV = -258.176). 
# 6: ! BW_lin ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-16795.527 > -16805.3). Accepting best available result (attempt 1/1, dOFV = -9.772). 
# 7: ! SEX_1 ~ cl: unrealistic OFV after all 1 attempt: OFV increased vs parent (-16546.25 > -16805.3). Accepting best available result (attempt 1/1, dOFV = -259.05). 
# 8: ! SEX_1 ~ vc: unrealistic OFV after all 1 attempt: OFV increased vs parent (-16546.991 > -16805.3). Accepting best available result (attempt 1/1, dOFV = -258.309). 
true_2cmt_scn16_refexp <- function() {
  ini({
    lTVCL      <- log(0.6)
    lTVQ       <- log(1.8)
    lTVVc      <- log(20)
    lTVVp      <- log(80)
    lTVKA      <- fix(log(0.7))    # KA unidentifiable
    TH_BW_CL   <- 0.75
    TH_CRCL_CL <- 0.5
    TH_BW_VC   <- 1.0
    TH_SEX_VC  <- log(1.5)         # nlmixr2scm convention: exp(theta*SEX)

    eta.cl + eta.vc ~ c(0.1, 0.02, 0.1)
    prop.err <- 0.1
  })
  model({
    cl_typ <- exp(lTVCL) * (BW / 70)^TH_BW_CL * (CrCL / 95)^TH_CRCL_CL
    vc_typ <- exp(lTVVc) * (BW / 70)^TH_BW_VC * exp(TH_SEX_VC * SEX)
    cl     <- cl_typ * exp(eta.cl)
    vc     <- vc_typ * exp(eta.vc)
    q      <- exp(lTVQ)
    vp     <- exp(lTVVp)
    ka     <- exp(lTVKA)

    k10 <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)      = -ka * depot
    d/dt(central)    =  ka * depot - k10 * central - k12 * central + k21 * peripheral
    d/dt(peripheral) =  k12 * central - k21 * peripheral

    cp = central / vc
    cp ~ prop(prop.err)
  })
}
true_2cmt_scn16_lin <- function() {
  ini({
    lTVCL      <- log(0.6)
    lTVQ       <- log(1.8)
    lTVVc      <- log(20)
    lTVVp      <- log(80)
    lTVKA      <- fix(log(0.7))
    TH_BW_CL   <- 0.75
    TH_CRCL_CL <- 0.5
    TH_BW_VC   <- 1.0
    TH_SEX_VC  <- log(1.5)         # nlmixr2scm convention

    eta.cl + eta.vc ~ c(0.1, 0.02, 0.1)
    prop.err <- 0.1
  })
  model({
    ## Pure log-additive form: matches `runSCM(shapes = c("power", "cat"))`
    ## verbatim.  BW and CrCL via log() => power on natural scale.  SEX is
    ## the indicator-coded categorical; coefficient enters linearly on the
    ## log scale (= exp(theta*SEX) on natural scale, matching the simulator).
    lTVCL_typ <- lTVCL + TH_BW_CL * log(BW / 70) + TH_CRCL_CL * log(CrCL / 95)
    lTVVc_typ <- lTVVc + TH_BW_VC * log(BW / 70) + TH_SEX_VC * SEX
    cl        <- exp(lTVCL_typ + eta.cl)
    vc        <- exp(lTVVc_typ + eta.vc)
    q         <- exp(lTVQ)
    vp        <- exp(lTVVp)
    ka        <- exp(lTVKA)

    k10 <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)      = -ka * depot
    d/dt(central)    =  ka * depot - k10 * central - k12 * central + k21 * peripheral
    d/dt(peripheral) =  k12 * central - k21 * peripheral

    cp = central / vc
    cp ~ prop(prop.err)
  })
}
out_dir_v2   <- "simulated_virtual_dataset_eta_filtered"
out_dir_v2   <- "simulated_virtual_dataset_eta_filtered"
stage1_dir16 <- file.path(out_dir_v2, "stage1_smoke_scn16_ds01")
if (!dir.exists(stage1_dir16)) dir.create(stage1_dir16, recursive = TRUE)
sim_obs_scn16 <- readRDS(file.path(out_dir_v2, "sim_obs_scenario_16.rds"))
ds16_01 <- to_nm_dataset(sim_obs_scn16) %>%
  dplyr::filter(DATASET == 1) %>%
  dplyr::select(-SCENARIO, -DATASET) %>%
  dplyr::mutate(
    ID   = as.integer(ID),
    SEX  = as.integer(SEX),
    RACE = as.integer(RACE)
  )
saveRDS(ds16_01, file.path(stage1_dir16, "nm_scn16_ds01.rds"))
out_dir_v2_N80   <- "simulated_virtual_dataset_eta_filtered_N80"
if (!dir.exists(stage1_dir16_N80)) dir.create(stage1_dir16_N80, recursive = TRUE)
# Error:
# ! object 'stage1_dir16_N80' not found
#     ▆
#  1. └─base::dir.exists(stage1_dir16_N80)
stage1_dir16_N80 <- file.path(out_dir_v2, "stage1_smoke_scn16_ds01_N80")
out_dir_v2_N80   <- "simulated_virtual_dataset_eta_filtered_N80"
stage1_dir16_N80 <- file.path(out_dir_v2_N80 , "stage1_smoke_scn16_ds01_N80")
if (!dir.exists(stage1_dir16_N80)) dir.create(stage1_dir16_N80, recursive = TRUE)
sim_obs_scn16_N80 <- readRDS(file.path(out_dir_v2_N80, "sim_obs_scenario_16.rds"))
ds16_01_N80 <- to_nm_dataset(sim_obs_scn16_N80) %>%
  dplyr::filter(DATASET == 1) %>%
  dplyr::select(-SCENARIO, -DATASET) %>%
  dplyr::mutate(
    ID   = as.integer(ID),
    SEX  = as.integer(SEX),
    RACE = as.integer(RACE)
  )
saveRDS(ds16_01_N80, file.path(stage1_dir16_N80, "nm_scn16_ds01_N80.rds"))
out_dir_v2_N40   <- "simulated_virtual_dataset_eta_filtered_N40"
stage1_dir16_N40 <- file.path(out_dir_v2_N40 , "stage1_smoke_scn16_ds01_N40")
if (!dir.exists(stage1_dir16_N40)) dir.create(stage1_dir16_N40, recursive = TRUE)
sim_obs_scn16_N40 <- readRDS(file.path(out_dir_v2_N40, "sim_obs_scenario_16.rds"))
ds16_01_N40 <- to_nm_dataset(sim_obs_scn16_N40) %>%
  dplyr::filter(DATASET == 1) %>%
  dplyr::select(-SCENARIO, -DATASET) %>%
  dplyr::mutate(
    ID   = as.integer(ID),
    SEX  = as.integer(SEX),
    RACE = as.integer(RACE)
  )
saveRDS(ds16_01_N40, file.path(stage1_dir16_N40, "nm_scn16_ds01_N40.rds"))
ds16_01_N40 <- readRDS(file.path(stage1_dir16_N40, "nm_scn16_ds01_N40.rds"))
scm_focei_production <- nlmixr2est::foceiControl(
  sigdig             = 3,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = FALSE,
  covMethod          = "",
  stickyRecalcN      = 20,                        # belt-and-braces
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl          = rxode2::rxControl(atol = 1e-6, rtol = 1e-4)
)
t_fit_true16_refexp80 <- system.time(
  fit_true16_refexp <- nlmixr2(true_2cmt_scn16_refexp, ds16_01_N80 ,
                               est = "focei", control = scm_focei_production)
)
t_fit_true16_lin80 <- system.time(
  fit_true16_lin <- nlmixr2(true_2cmt_scn16_lin, ds16_01_N80 ,
                            est = "focei", control = scm_focei_production)
)

t_fit_true16_refexp40 <- system.time(
  fit_true16_refexp <- nlmixr2(true_2cmt_scn16_refexp, ds16_01_N40 ,
                               est = "focei", control = scm_focei_production)
)

t_fit_true16_lin40 <- system.time(
  fit_true16_lin <- nlmixr2(true_2cmt_scn16_lin, ds16_01_N40 ,
                            est = "focei", control = scm_focei_production)
)
# ℹ parameter labels from comments will be replaced by 'label()'
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities                                                        
# → calculate ∂(f)/∂(η)                                                            
# → calculate ∂(R²)/∂(η)                                                           
# → finding duplicate expressions in inner model...                                
# → optimizing duplicate expressions in inner model...                             
# → finding duplicate expressions in EBE model...                                  
# → optimizing duplicate expressions in EBE model...                               
# → compiling inner model...                                                       
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...                                
# → compiling EBE model...                                                         
# ✔ done
# → compiling events FD model...
# ✔ done
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities                                                        
# → calculate ∂(f)/∂(η)                                                            
# → calculate ∂(R²)/∂(η)                                                           
# → finding duplicate expressions in inner model...                                
# → optimizing duplicate expressions in inner model...                             
# → finding duplicate expressions in EBE model...                                  
# → optimizing duplicate expressions in EBE model...                               
# → compiling inner model...                                                       
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...                                
# → compiling EBE model...                                                         
# ✔ done
# → compiling events FD model...
# ✔ done
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
t_fit_true16_refexp300 <- system.time(
  fit_true16_refexp <- nlmixr2(true_2cmt_scn16_refexp, ds16_01,
                               est = "focei", control = scm_focei_production)
)
t_fit_true16_lin300 <- system.time(
  fit_true16_lin <- nlmixr2(true_2cmt_scn16_lin, ds16_01,
                            est = "focei", control = scm_focei_production)
)
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
t_fit_true16_refexp300 <- system.time(
  fit_true16_refexp300 <- nlmixr2(true_2cmt_scn16_refexp, ds16_01,
                                  est = "focei", control = scm_focei_production)
)
t_fit_true16_lin300 <- system.time(
  fit_true16_lin300 <- nlmixr2(true_2cmt_scn16_lin, ds16_01,
                               est = "focei", control = scm_focei_production)
)

t_fit_true16_refexp80 <- system.time(
  fit_true16_refexp80 <- nlmixr2(true_2cmt_scn16_refexp, ds16_01_N80,
                                 est = "focei", control = scm_focei_production)
)
t_fit_true16_lin80 <- system.time(
  fit_true16_lin80 <- nlmixr2(true_2cmt_scn16_lin, ds16_01_N80,
                              est = "focei", control = scm_focei_production)
)

t_fit_true16_refexp40 <- system.time(
  fit_true16_refexp40 <- nlmixr2(true_2cmt_scn16_refexp, ds16_01_N40,
                                 est = "focei", control = scm_focei_production)
)
t_fit_true16_lin40 <- system.time(
  fit_true16_lin40 <- nlmixr2(true_2cmt_scn16_lin, ds16_01_N40,
                              est = "focei", control = scm_focei_production)
)
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
# ℹ parameter labels from comments will be replaced by 'label()'
# done
t_fit_true16_refexp <- t_fit_true16_refexp300
t_fit_true16_lin    <- t_fit_true16_lin300
fit_registry16 <- tibble::tibble(
  parameterisation = c("refexp", "lin", "refexp", "lin", "refexp", "lin"),
  n_subj           = c(300L, 300L, 80L, 80L, 40L, 40L),
  fit              = list(fit_true16_refexp300, fit_true16_lin300,
                          fit_true16_refexp80,  fit_true16_lin80,
                          fit_true16_refexp40,  fit_true16_lin40),
  runtime          = list(t_fit_true16_refexp300, t_fit_true16_lin300,
                          t_fit_true16_refexp80,  t_fit_true16_lin80,
                          t_fit_true16_refexp40,  t_fit_true16_lin40)
)
diag_summary16 <- fit_registry16 %>%
  dplyr::mutate(
    diag = purrr::map(fit, diagnose_fit),
    row  = purrr::pmap(
      list(parameterisation, diag, runtime),
      function(p, d, r) .build_part1_row(p, d, r)
    )
  ) %>%
  dplyr::select(parameterisation, n_subj, row) %>%
  tidyr::unnest(row) %>%
  dplyr::mutate(scenario = 16) %>%
  dplyr::relocate(scenario, n_subj, parameterisation) %>%
  dplyr::arrange(dplyr::desc(n_subj), parameterisation)
# Error in `tidyr::unnest()`:
# ! Can't duplicate names between the affected columns and the original data.
# ✖ These names are duplicated:
#   ℹ `parameterisation`, from `row`.
# ℹ Use `names_sep` to disambiguate using the column name.
# ℹ Or use `names_repair` to specify a repair strategy.
#      ▆
#   1. ├─... %>% ...
#   2. ├─dplyr::arrange(., dplyr::desc(n_subj), parameterisation)
#   3. ├─dplyr::relocate(., scenario, n_subj, parameterisation)
#   4. ├─dplyr::mutate(., scenario = 16)
#   5. ├─tidyr::unnest(., row)
#   6. └─tidyr:::unnest.data.frame(., row)
#   7.   └─tidyr::unpack(...)
#   8.     └─tidyr:::check_outer_inner_duplicate(cols, names(data), error_call = error_call)
#   9.       └─cli::cli_abort(message, call = error_call)
#  10.         └─rlang::abort(...)
diag_summary16 <- fit_registry16 %>%
  dplyr::mutate(
    diag = purrr::map(fit, diagnose_fit),
    row  = purrr::pmap(
      list(parameterisation, diag, runtime),
      function(p, d, r) .build_part1_row(p, d, r)
    )
  ) %>%
  dplyr::select(n_subj, row) %>%
  tidyr::unnest(row) %>%
  dplyr::mutate(scenario = 16) %>%
  dplyr::relocate(scenario, n_subj, parameterisation) %>%
  dplyr::arrange(dplyr::desc(n_subj), parameterisation)
err_true16_compare <- fit_registry16 %>%
  dplyr::mutate(
    est     = purrr::map(fit, extract_params_long),
    rel_err = purrr::map(
      est, rel_err_one,
      true_long   = true_params,
      scenario_id = 16
    )
  ) %>%
  dplyr::select(parameterisation, n_subj, rel_err) %>%
  tidyr::unnest(rel_err) %>%
  dplyr::select(n_subj, parameterisation, parameter, true_value,
                estimate, abs_err, rel_err, rel_err_pct) %>%
  dplyr::arrange(dplyr::desc(n_subj), parameter, parameterisation)
err_true16_wide <- err_true16_compare %>%
  tidyr::pivot_wider(
    id_cols     = c(n_subj, parameter, true_value),
    names_from  = parameterisation,
    values_from = c(estimate, rel_err_pct),
    names_glue  = "{.value}_{parameterisation}"
  ) %>%
  dplyr::mutate(
    abs_diff_estimate    = abs(estimate_refexp - estimate_lin),
    abs_diff_rel_err_pct = abs(rel_err_pct_refexp - rel_err_pct_lin)
  ) %>%
  dplyr::arrange(dplyr::desc(n_subj), parameter)
extract_params_long <- function(fit, includeCov = TRUE) {
  theta <- fit$theta
  om    <- fit$omega

  ## Fixed-effects: keep the "TV" prefix so names line up with true_params_long()
  ## (TVCL, TVVc, TVQ, TVVp, TVKA).  fit$theta still stores them on the LOG
  ## scale (lTVCL etc.), so we exp() back to natural-scale typical values.
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

  ## Covariate effects: read from fit$theta if present, else NA.
  ##   The theta names emitted by runSCM() preserve the original case of the
  ##   data column (e.g. "cov_BW_power_cl" when the data has column "BW"),
  ##   so the regex must match case-insensitively to catch both "BW" and
  ##   "wt" / "bw" parameterisations.
  cov_names_in_fit <- names(theta)
  cov_map <- list(
    CLBW   = c("TH_BW_CL",   grep("^cov_(bw|wt)_power_cl$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    CLcrCL = c("TH_CRCL_CL", grep("^cov_crcl_power_cl$",      cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcBW   = c("TH_BW_VC",   grep("^cov_(bw|wt)_power_vc$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcSEX  = c("TH_SEX_VC",  grep("^cov_sex(_male)?_cat_vc$", cov_names_in_fit, value = TRUE, ignore.case = TRUE))
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
err_true16_compare <- fit_registry16 %>%
  dplyr::mutate(
    est     = purrr::map(fit, extract_params_long),
    rel_err = purrr::map(
      est, rel_err_one,
      true_long   = true_params,
      scenario_id = 16
    )
  ) %>%
  dplyr::select(parameterisation, n_subj, rel_err) %>%
  tidyr::unnest(rel_err) %>%
  dplyr::select(n_subj, parameterisation, parameter, true_value,
                estimate, abs_err, rel_err, rel_err_pct) %>%
  dplyr::arrange(dplyr::desc(n_subj), parameter, parameterisation)
err_true16_wide <- err_true16_compare %>%
  tidyr::pivot_wider(
    id_cols     = c(n_subj, parameter, true_value),
    names_from  = parameterisation,
    values_from = c(estimate, rel_err_pct),
    names_glue  = "{.value}_{parameterisation}"
  ) %>%
  dplyr::mutate(
    abs_diff_estimate    = abs(estimate_refexp - estimate_lin),
    abs_diff_rel_err_pct = abs(rel_err_pct_refexp - rel_err_pct_lin)
  ) %>%
  dplyr::arrange(dplyr::desc(n_subj), parameter)
saveRDS(err_true16_wide, file.path(stage1_dir16, "part1_err_wide.rds"))
err_true16_lin_plot <- err_true16_compare %>%
  dplyr::filter(parameterisation == "lin") %>%
  dplyr::mutate(
    abs_rel_err_pct = abs(rel_err_pct),
    n_label  = factor(n_subj,
                      levels = c(300L, 80L, 40L),
                      labels = c("N = 300", "N = 80", "N = 40")),
    is_worst = parameter == "cov_VcCL"
  )
p_err_true16_lin <- ggplot2::ggplot(
  err_true16_lin_plot,
  ggplot2::aes(x = n_label, y = abs_rel_err_pct,
               group = parameter, colour = is_worst)
) +
  ggplot2::geom_line(ggplot2::aes(linewidth = is_worst), alpha = 0.9) +
  ggplot2::geom_point(ggplot2::aes(size = is_worst)) +
  ggplot2::geom_text(
    data = dplyr::filter(err_true16_lin_plot, n_label == "N = 40"),
    ggplot2::aes(label = parameter),
    hjust       = 0,
    nudge_x     = 0.10,
    size        = 3,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = c(`FALSE` = "grey70", `TRUE` = "#C0392B"),
    guide  = "none"
  ) +
  ggplot2::scale_linewidth_manual(
    values = c(`FALSE` = 0.4, `TRUE` = 1.4),
    guide  = "none"
  ) +
  ggplot2::scale_size_manual(
    values = c(`FALSE` = 1.4, `TRUE` = 3),
    guide  = "none"
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0.02, 0.10))
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title    = "Estimation error grows as cohort shrinks; cov_VcCL is the worst-recovered parameter",
    subtitle = "Lin parameterisation, scenario 16, dataset 01",
    x        = NULL,
    y        = "|Relative error|",
    caption  = "Each line = one model parameter; cov_VcCL highlighted in red."
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.margin         = ggplot2::margin(5.5, 60, 5.5, 5.5),
    panel.grid.minor    = ggplot2::element_blank()
  )
print(p_err_true16_lin)
p_err_true16_lin <- ggplot2::ggplot(
  err_true16_lin_plot,
  ggplot2::aes(x = n_label, y = abs_rel_err_pct,
               group = parameter, colour = is_worst)
) +
  ggplot2::geom_line(ggplot2::aes(linewidth = is_worst), alpha = 0.9) +
  ggplot2::geom_point(ggplot2::aes(size = is_worst)) +
  ggplot2::geom_text(
    data = dplyr::filter(err_true16_lin_plot, n_label == "N = 40"),
    ggplot2::aes(label = parameter),
    hjust       = 0,
    nudge_x     = 0.10,
    size        = 4,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(
    values = c(`FALSE` = "grey70", `TRUE` = "#C0392B"),
    guide  = "none"
  ) +
  ggplot2::scale_linewidth_manual(
    values = c(`FALSE` = 0.4, `TRUE` = 1.4),
    guide  = "none"
  ) +
  ggplot2::scale_size_manual(
    values = c(`FALSE` = 1.4, `TRUE` = 3),
    guide  = "none"
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0.02, 0.10))
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title    = "Estimation error grows as cohort shrinks; cov_VcCL is the worst-recovered parameter",
    subtitle = "Lin parameterisation, scenario 16, dataset 01",
    x        = NULL,
    y        = "|Relative error|",
    caption  = "Each line = one model parameter; cov_VcCL highlighted in red."
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    plot.margin         = ggplot2::margin(5.5, 60, 5.5, 5.5),
    panel.grid.minor    = ggplot2::element_blank()
  )
print(p_err_true16_lin)
err_true16_lin_order <- err_true16_lin_plot %>%
  dplyr::filter(n_label == "N = 40") %>%
  dplyr::arrange(dplyr::desc(abs_rel_err_pct)) %>%
  dplyr::pull(parameter)

err_true16_lin_bars <- err_true16_lin_plot %>%
  dplyr::mutate(parameter = factor(parameter, levels = err_true16_lin_order))

p_err_true16_lin_bar <- ggplot2::ggplot(
  err_true16_lin_bars,
  ggplot2::aes(x = parameter, y = abs_rel_err_pct, fill = n_label)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                    width    = 0.75) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f%%", abs_rel_err_pct)),
    position = ggplot2::position_dodge(width = 0.8),
    vjust    = -0.4,
    size     = 2.7
  ) +
  ggplot2::scale_fill_manual(
    name   = "Cohort size",
    values = c("N = 300" = "#9ECAE1",
               "N = 80"  = "#4292C6",
               "N = 40"  = "#08519C")
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    title    = "Per-parameter |relative error| grows as cohort shrinks",
    subtitle = "Lin parameterisation, scenario 16, dataset 01 (parameters sorted by worst-N error)",
    x        = NULL,
    y        = "|Relative error|"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    axis.text.x         = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid.major.x  = ggplot2::element_blank(),
    panel.grid.minor    = ggplot2::element_blank(),
    legend.position     = "top"
  )

print(p_err_true16_lin_bar)
ggplot2::ggsave(
  file.path(stage1_dir16, "part1_err_lin_by_N_bar.png"),
  p_err_true16_lin_bar,
  width = 9, height = 5, dpi = 150
)
ggplot2::ggsave(
  file.path(stage1_dir16, "part1_err_lin_by_N.png"),
  p_err_true16_lin,
  width = 10, height = 5, dpi = 150
)
ggplot2::ggsave(
  file.path(stage1_dir16, "part1_err_lin_by_N_bar.png"),
  p_err_true16_lin_bar,
  width = 13, height = 5, dpi = 150
)
err_true16_lin_order <- err_true16_lin_plot %>%
  dplyr::filter(n_label == "N = 40") %>%
  dplyr::arrange(dplyr::desc(abs_rel_err_pct)) %>%
  dplyr::pull(parameter)
err_true16_lin_bars <- err_true16_lin_plot %>%
  dplyr::mutate(parameter = factor(parameter, levels = err_true16_lin_order))
p_err_true16_lin_bar <- ggplot2::ggplot(
  err_true16_lin_bars,
  ggplot2::aes(x = parameter, y = abs_rel_err_pct, fill = n_label)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                    width    = 0.75) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f%%", abs_rel_err_pct)),
    position = ggplot2::position_dodge(width = 0.8),
    vjust    = -0.4,
    size     = 3
  ) +
  ggplot2::scale_fill_manual(
    name   = "Cohort size",
    values = c("N = 300" = "#9ECAE1",
               "N = 80"  = "#4292C6",
               "N = 40"  = "#08519C")
  ) +
  ggplot2::scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0, 0.12))
  ) +
  ggplot2::labs(
    title    = "Per-parameter |relative error| grows as cohort shrinks",
    subtitle = "Lin parameterisation, scenario 16, dataset 01 (parameters sorted by worst-N error)",
    x        = NULL,
    y        = "|Relative error|"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title.position = "plot",
    axis.text.x         = ggplot2::element_text(angle = 35, hjust = 1),
    panel.grid.major.x  = ggplot2::element_blank(),
    panel.grid.minor    = ggplot2::element_blank(),
    legend.position     = "top"
  )
ggplot2::ggsave(
  file.path(stage1_dir16, "part1_err_lin_by_N_bar.png"),
  p_err_true16_lin_bar,
  width = 14, height = 5, dpi = 150
)
rxode2::rxControl()
# $scale
# NULL
# 
# $method
# liblsoda 
#        2 
# 
# $atol
# [1] 1e-08
# 
# $rtol
# [1] 1e-06
# 
# $maxsteps
# [1] 70000
# 
# $hmin
# [1] 0
# 
# $hmax
# [1] NA
# 
# $hini
# [1] 0
# 
# $maxordn
# [1] 12
# 
# $maxords
# [1] 5
# 
# $covsInterpolation
# locf 
#    1 
# 
# $addCov
# [1] TRUE
# 
# $returnType
# rxSolve 
#       0 
# 
# $sigma
# NULL
# 
# $sigmaDf
# NULL
# 
# $nCoresRV
# [1] 1
# 
# $sigmaIsChol
# [1] FALSE
# 
# $sigmaSeparation
# [1] "auto"
# 
# $sigmaXform
# identity 
#        4 
# 
# $nDisplayProgress
# [1] 10000
# 
# $amountUnits
# [1] NA
# 
# $timeUnits
# [1] "hours"
# 
# $addDosing
# [1] FALSE
# 
# $stateTrim
# [1] Inf
# 
# $updateObject
# [1] FALSE
# 
# $omega
# NULL
# 
# $omegaDf
# NULL
# 
# $omegaIsChol
# [1] FALSE
# 
# $omegaSeparation
# [1] "auto"
# 
# $omegaXform
# variance 
#        6 
# 
# $nSub
# [1] 1
# 
# $thetaMat
# NULL
# 
# $thetaDf
# NULL
# 
# $thetaIsChol
# [1] FALSE
# 
# $nStud
# [1] 1
# 
# $dfSub
# [1] 0
# 
# $dfObs
# [1] 0
# 
# $seed
# NULL
# 
# $nsim
# NULL
# 
# $minSS
# [1] 10
# 
# $maxSS
# [1] 10000
# 
# $strictSS
# [1] 1
# 
# $infSSstep
# [1] 12
# 
# $istateReset
# [1] TRUE
# 
# $subsetNonmem
# [1] TRUE
# 
# $hmaxSd
# [1] 0
# 
# $maxAtolRtolFactor
# [1] 0.1
# 
# $from
# NULL
# 
# $to
# NULL
# 
# $by
# NULL
# 
# $length.out
# NULL
# 
# $iCov
# NULL
# 
# $keep
# NULL
# 
# $keepF
# character(0)
# 
# $drop
# NULL
# 
# $warnDrop
# [1] TRUE
# 
# $omegaLower
# [1] -Inf
# 
# $omegaUpper
# [1] Inf
# 
# $sigmaLower
# [1] -Inf
# 
# $sigmaUpper
# [1] Inf
# 
# $thetaLower
# [1] -Inf
# 
# $thetaUpper
# [1] Inf
# 
# $indLinPhiM
# [1] 0
# 
# $indLinPhiTol
# [1] 1e-07
# 
# $indLinMatExpType
# expokit 
#       2 
# 
# $indLinMatExpOrder
# [1] 6
# 
# $idFactor
# [1] TRUE
# 
# $mxhnil
# [1] 0
# 
# $hmxi
# [1] 0
# 
# $warnIdSort
# [1] TRUE
# 
# $ssAtol
# [1] 1e-08
# 
# $ssRtol
# [1] 1e-06
# 
# $safeZero
# [1] 1
# 
# $sumType
# pairwise 
#        1 
# 
# $prodType
# long double 
#           1 
# 
# $resample
# NULL
# 
# $resampleID
# [1] TRUE
# 
# $maxwhile
# [1] 100000
# 
# $cores
# [1] 0
# 
# $atolSens
# [1] 1e-08
# 
# $rtolSens
# [1] 1e-06
# 
# $ssAtolSens
# [1] 1e-08
# 
# $ssRtolSens
# [1] 1e-06
# 
# $simVariability
# [1] NA
# 
# $nLlikAlloc
# NULL
# 
# $useStdPow
# [1] 0
# 
# $naTimeHandle
# ignore 
#      1 
# 
# $addlKeepsCov
# [1] FALSE
# 
# $addlDropSs
# [1] TRUE
# 
# $ssAtDoseTime
# [1] TRUE
# 
# $ss2cancelAllPending
# [1] FALSE
# 
# $naInterpolation
# locf 
#    1 
# 
# $keepInterpolation
# na 
#  2 
# 
# $safeLog
# [1] 1
# 
# $safePow
# [1] 1
# 
# $ssSolved
# [1] TRUE
# 
# $linCmtSensType
# auto 
#  100 
# 
# $linCmtSensH
# [1] 1e-04
# 
# $linCmtGillFtol
# [1] 0
# 
# $linCmtGillK
# [1] 20
# 
# $linCmtGillStep
# [1] 4
# 
# $linCmtGillRtol
# [1] 1.490116e-08
# 
# $linCmtShiErr
# [1] 1.490116e-08
# 
# $linCmtShiMax
# [1] 20
# 
# $linCmtScale
# [1] 0 0 0 0 0 0 0
# 
# $linCmtHcmt
# [1] 1
# 
# $linCmtHmeanI
# geometric 
#         2 
# 
# $linCmtHmeanO
# geometric 
#         2 
# 
# $linCmtSuspect
# [1] 1e-06
# 
# $linCmtForwardMax
# [1] 2
# 
# $indOwnAlloc
# [1] -1
# 
# $maxExtra
# [1] 1000
# 
# $tolFactor
# NULL
# 
# $serializeFile
# NULL
# 
# $dense
# [1] FALSE
# 
# $cvodeLinSolver
# dense 
#     1 
# 
# $stiff2
# [1] 0
# 
# $autoSwitchMaxStiff
# [1] 10
# 
# $autoSwitchMaxNonstiff
# [1] 3
# 
# $autoSwitchStiffFirst
# [1] 0
# 
# $autoSwitchNonstifftol
# [1] 0.9
# 
# $autoSwitchStifftol
# [1] 0.9
# 
# $autoSwitchDtfac
# [1] 2
# 
# $autoSwitchSwitchMax
# [1] 5
# 
# $useLinCmt
# [1] TRUE
# 
# $.zeros
# NULL
# 
# attr(,"class")
# [1] "rxControl"
nlmixr2est::foceiControl()
# $maxOuterIterations
# [1] 5000
# 
# $maxInnerIterations
# [1] 1000
# 
# $n1qn1nsim
# [1] 10001
# 
# $iterPrintControl
# $every
# [1] 1
# 
# $ncol
# [1] 8
# 
# $headerEvery
# [1] 10
# 
# $useColor
# [1] TRUE
# 
# $simple
# [1] FALSE
# 
# attr(,"class")
# [1] "iterPrintControl" "list"            
# 
# $lbfgsLmm
# [1] 7
# 
# $lbfgsPgtol
# [1] 0
# 
# $lbfgsFactr
# [1] 4.5036e+10
# 
# $scaleTo
# [1] 1
# 
# $epsilon
# [1] 1e-05
# 
# $derivEps
# [1] 2.980232e-07 2.980232e-07
# 
# $derivMethod
# [1] 3
# 
# $covDerivMethod
# [1] 1
# 
# $covMethod
# [1] 1
# 
# $centralDerivEps
# [1] 2.980232e-07 2.980232e-07
# 
# $eigen
# [1] 1
# 
# $diagXform
# [1] "sqrt"
# 
# $iovXform
# [1] "sd"
# 
# $sumProd
# [1] FALSE
# 
# $optExpression
# [1] TRUE
# 
# $literalFix
# [1] TRUE
# 
# $literalFixRes
# [1] TRUE
# 
# $outerOpt
# [1] 1
# 
# $ci
# [1] 0.95
# 
# $sigdig
# [1] 4
# 
# $sigdigTable
# [1] 4
# 
# $scaleObjective
# [1] 0
# 
# $boundTol
# [1] 0.005
# 
# $calcTables
# [1] TRUE
# 
# $noAbort
# [1] 1
# 
# $interaction
# [1] 1
# 
# $cholSEtol
# [1] 6.055454e-06
# 
# $hessEps
# [1] 6.055454e-06
# 
# $hessEpsLlik
# [1] 6.055454e-06
# 
# $optimHessType
# [1] 1
# 
# $optimHessCovType
# [1] 1
# 
# $cholAccept
# [1] 0.001
# 
# $resetEtaSize
# [1] 1.439531
# 
# $resetThetaSize
# [1] 1.959964
# 
# $resetThetaFinalSize
# [1] 1.439531
# 
# $diagOmegaBoundUpper
# [1] 5
# 
# $diagOmegaBoundLower
# [1] 100
# 
# $cholSEOpt
# [1] 0
# 
# $cholSECov
# [1] 0
# 
# $fo
# [1] 0
# 
# $covTryHarder
# [1] 0
# 
# $outerOptFun
# NULL
# 
# $rhobeg
# [1] 0.2
# 
# $rhoend
# [1] 1e-05
# 
# $npt
# NULL
# 
# $rel.tol
# [1] 1e-05
# 
# $x.tol
# [1] 1e-05
# 
# $eval.max
# [1] 4000
# 
# $iter.max
# [1] 2000
# 
# $innerOpt
# [1] 1
# 
# $abstol
# [1] 1e-05
# 
# $reltol
# [1] 1e-05
# 
# $derivSwitchTol
# [1] 2e-05
# 
# $resetHessianAndEta
# [1] 0
# 
# $stateTrim
# [1] Inf
# 
# $gillK
# [1] 10
# 
# $gillKcov
# [1] 10
# 
# $gillKcovLlik
# [1] 10
# 
# $gillRtol
# [1] 1.490116e-08
# 
# $gillStep
# [1] 4
# 
# $gillStepCov
# [1] 2
# 
# $gillStepCovLlik
# [1] 4.5
# 
# $scaleType
# [1] 2
# 
# $normType
# [1] 1
# 
# $scaleC
# NULL
# 
# $scaleCmin
# [1] 1e-05
# 
# $scaleCmax
# [1] 1e+05
# 
# $scaleC0
# [1] 1e+05
# 
# $outerOptTxt
# [1] "lbfgsb3c"
# 
# $rmatNorm
# [1] 1
# 
# $rmatNormLlik
# [1] 1
# 
# $smatNorm
# [1] 1
# 
# $smatNormLlik
# [1] 1
# 
# $covGillF
# [1] 1
# 
# $optGillF
# [1] 1
# 
# $gillFtol
# [1] 0
# 
# $gillFtolCov
# [1] 0
# 
# $gillFtolCovLlik
# [1] 0
# 
# $covSmall
# [1] 1e-05
# 
# $adjLik
# [1] TRUE
# 
# $gradTrim
# [1] Inf
# 
# $gradCalcCentralSmall
# [1] 1e-04
# 
# $gradCalcCentralLarge
# [1] 10000
# 
# $etaNudge
# [1] 1.131586
# 
# $etaNudge2
# [1] 1.518182
# 
# $maxOdeRecalc
# [1] 5
# 
# $odeRecalcFactor
# [1] 3.162278
# 
# $nRetries
# [1] 3
# 
# $seed
# [1] 42
# 
# $resetThetaCheckPer
# [1] 0.1
# 
# $etaMat
# NULL
# 
# $repeatGillMax
# [1] 1
# 
# $stickyRecalcN
# [1] 4
# 
# $indTolRelax
# [1] TRUE
# 
# $eventType
# [1] 2
# 
# $gradProgressOfvTime
# [1] 10
# 
# $addProp
# [1] "combined2"
# 
# $badSolveObjfAdj
# [1] 100
# 
# $compress
# [1] FALSE
# 
# $rxControl
# $scale
# NULL
# 
# $method
# liblsoda 
#        2 
# 
# $atol
# [1] 5e-07
# 
# $rtol
# [1] 5e-07
# 
# $maxsteps
# [1] 500000
# 
# $hmin
# [1] 0
# 
# $hmax
# [1] NA
# 
# $hini
# [1] 0
# 
# $maxordn
# [1] 12
# 
# $maxords
# [1] 5
# 
# $covsInterpolation
# locf 
#    1 
# 
# $addCov
# [1] TRUE
# 
# $returnType
# rxSolve 
#       0 
# 
# $sigma
# NULL
# 
# $sigmaDf
# NULL
# 
# $nCoresRV
# [1] 1
# 
# $sigmaIsChol
# [1] FALSE
# 
# $sigmaSeparation
# [1] "auto"
# 
# $sigmaXform
# identity 
#        4 
# 
# $nDisplayProgress
# [1] 10000
# 
# $amountUnits
# [1] NA
# 
# $timeUnits
# [1] "hours"
# 
# $addDosing
# [1] FALSE
# 
# $stateTrim
# [1] Inf
# 
# $updateObject
# [1] FALSE
# 
# $omega
# NULL
# 
# $omegaDf
# NULL
# 
# $omegaIsChol
# [1] FALSE
# 
# $omegaSeparation
# [1] "auto"
# 
# $omegaXform
# variance 
#        6 
# 
# $nSub
# [1] 1
# 
# $thetaMat
# NULL
# 
# $thetaDf
# NULL
# 
# $thetaIsChol
# [1] FALSE
# 
# $nStud
# [1] 1
# 
# $dfSub
# [1] 0
# 
# $dfObs
# [1] 0
# 
# $seed
# NULL
# 
# $nsim
# NULL
# 
# $minSS
# [1] 10
# 
# $maxSS
# [1] 10000
# 
# $strictSS
# [1] 1
# 
# $infSSstep
# [1] 12
# 
# $istateReset
# [1] TRUE
# 
# $subsetNonmem
# [1] TRUE
# 
# $hmaxSd
# [1] 0
# 
# $maxAtolRtolFactor
# [1] 0.1
# 
# $from
# NULL
# 
# $to
# NULL
# 
# $by
# NULL
# 
# $length.out
# NULL
# 
# $iCov
# NULL
# 
# $keep
# NULL
# 
# $keepF
# character(0)
# 
# $drop
# NULL
# 
# $warnDrop
# [1] TRUE
# 
# $omegaLower
# [1] -Inf
# 
# $omegaUpper
# [1] Inf
# 
# $sigmaLower
# [1] -Inf
# 
# $sigmaUpper
# [1] Inf
# 
# $thetaLower
# [1] -Inf
# 
# $thetaUpper
# [1] Inf
# 
# $indLinPhiM
# [1] 0
# 
# $indLinPhiTol
# [1] 1e-07
# 
# $indLinMatExpType
# expokit 
#       2 
# 
# $indLinMatExpOrder
# [1] 6
# 
# $idFactor
# [1] TRUE
# 
# $mxhnil
# [1] 0
# 
# $hmxi
# [1] 0
# 
# $warnIdSort
# [1] TRUE
# 
# $ssAtol
# [1] 5e-05
# 
# $ssRtol
# [1] 5e-05
# 
# $safeZero
# [1] 1
# 
# $sumType
# pairwise 
#        1 
# 
# $prodType
# long double 
#           1 
# 
# $resample
# NULL
# 
# $resampleID
# [1] TRUE
# 
# $maxwhile
# [1] 100000
# 
# $cores
# [1] 0
# 
# $atolSens
# [1] 1.581139e-06
# 
# $rtolSens
# [1] 1.581139e-06
# 
# $ssAtolSens
# [1] 0.0002108483
# 
# $ssRtolSens
# [1] 0.0002108483
# 
# $simVariability
# [1] NA
# 
# $nLlikAlloc
# NULL
# 
# $useStdPow
# [1] 0
# 
# $naTimeHandle
# ignore 
#      1 
# 
# $addlKeepsCov
# [1] FALSE
# 
# $addlDropSs
# [1] TRUE
# 
# $ssAtDoseTime
# [1] TRUE
# 
# $ss2cancelAllPending
# [1] FALSE
# 
# $naInterpolation
# locf 
#    1 
# 
# $keepInterpolation
# na 
#  2 
# 
# $safeLog
# [1] 1
# 
# $safePow
# [1] 1
# 
# $ssSolved
# [1] TRUE
# 
# $linCmtSensType
# auto 
#  100 
# 
# $linCmtSensH
# [1] 1e-04
# 
# $linCmtGillFtol
# [1] 0
# 
# $linCmtGillK
# [1] 20
# 
# $linCmtGillStep
# [1] 4
# 
# $linCmtGillRtol
# [1] 1.490116e-08
# 
# $linCmtShiErr
# [1] 1.490116e-08
# 
# $linCmtShiMax
# [1] 20
# 
# $linCmtScale
# [1] 0 0 0 0 0 0 0
# 
# $linCmtHcmt
# [1] 1
# 
# $linCmtHmeanI
# geometric 
#         2 
# 
# $linCmtHmeanO
# geometric 
#         2 
# 
# $linCmtSuspect
# [1] 1e-06
# 
# $linCmtForwardMax
# [1] 2
# 
# $indOwnAlloc
# [1] -1
# 
# $maxExtra
# [1] 1000
# 
# $tolFactor
# NULL
# 
# $serializeFile
# NULL
# 
# $dense
# [1] FALSE
# 
# $cvodeLinSolver
# dense 
#     1 
# 
# $stiff2
# [1] 0
# 
# $autoSwitchMaxStiff
# [1] 10
# 
# $autoSwitchMaxNonstiff
# [1] 3
# 
# $autoSwitchStiffFirst
# [1] 0
# 
# $autoSwitchNonstifftol
# [1] 0.9
# 
# $autoSwitchStifftol
# [1] 0.9
# 
# $autoSwitchDtfac
# [1] 2
# 
# $autoSwitchSwitchMax
# [1] 5
# 
# $useLinCmt
# [1] TRUE
# 
# $.zeros
# NULL
# 
# attr(,"class")
# [1] "rxControl"
# 
# $genRxControl
# [1] TRUE
# 
# $skipCov
# NULL
# 
# $fallbackFD
# [1] FALSE
# 
# $shi21maxOuter
# [1] 0
# 
# $shi21maxInner
# [1] 20
# 
# $shi21maxInnerCov
# [1] 20
# 
# $shi21maxFD
# [1] 20
# 
# $smatPer
# [1] 0.6
# 
# $sdLowerFact
# [1] 0.001
# 
# $zeroGradFirstReset
# [1] TRUE
# 
# $zeroGradRunReset
# [1] TRUE
# 
# $zeroGradBobyqa
# [1] TRUE
# 
# $mceta
# [1] -1
# 
# $nAGQ
# [1] 0
# 
# $agqHi
# [1] Inf
# 
# $agqLow
# [1] -Inf
# 
# $boundedTransform
# [1] TRUE
# 
# attr(,"class")
# [1] "foceiControl"
scm_focei_n <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = "",         # SCM doesn't need cov matrix for LRT
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
)
# Error in `nlmixr2est::foceiControl()`:
# ! argument is missing, with no default
#     ▆
#  1. └─nlmixr2est::foceiControl(...)
#  2.   └─nlmixr2est::.absorbIterPrintControl(...)
scm_focei_n <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = "",         # SCM doesn't need cov matrix for LRT
  maxOuterIterations = 2000,
  maxInnerIterations = 2000
)
t_fit_true16_refexp300_n <- system.time(
  fit_true16_refexp300_n <- nlmixr2(true_2cmt_scn16_refexp, ds16_01,
                                  est = "focei", control = scm_focei_n)
)
# ℹ parameter labels from comments will be replaced by 'label()'
# done
t_fit_true16_lin300_n <- system.time(
  fit_true16_lin300_n <- nlmixr2(true_2cmt_scn16_lin, ds16_01,
                               est = "focei", control = scm_focei_n)
)
# ℹ parameter labels from comments will be replaced by 'label()'
# done
diag_true_refexp300_n <- diagnose_fit(fit_true16_refexp300_n)
diag_true_lin300_n   <- diagnose_fit(fit_true16_lin300_n)
t_fit_base16_N80 <- system.time(
  fit_base16_N80 <- nlmixr2(base_2cmt_oral, ds16_01_N80,
                        est = "focei", control = scm_focei_n)
)
# done
scm_focei_n <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = "r.s",         # SCM doesn't need cov matrix for LRT
  maxOuterIterations = 2000,
  maxInnerIterations = 2000
)
# Error in `match.arg()`:
# ! 'arg' should be one of “r,s”, “r”, “s”
#     ▆
#  1. └─nlmixr2est::foceiControl(...)
#  2.   └─base::match.arg(covMethod)
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
# [09:37:15 | +  0.1 min] 
# ── starting forward search... ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:42:54 | +  5.8 min] 
# → Forward step 1, candidate 1/16: BW_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:42:54 | +  5.8 min] 
# → Forward step 1, candidate 2/16: BW_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:42:54 | +  5.8 min] 
# → Forward step 1, candidate 3/16: BW_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:42:54 | +  5.8 min] 
# → Forward step 1, candidate 4/16: BW_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:42:54 | +  5.8 min] 
# → Forward step 1, candidate 5/16: CrCL_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 6/16: CrCL_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 7/16: CrCL_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 8/16: CrCL_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 9/16: BMI_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 10/16: BMI_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 11/16: BMI_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 12/16: BMI_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 13/16: SEX_1 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 14/16: SEX_1 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 15/16: RACE_0 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:43:05 | +  6.0 min] 
# → Forward step 1, candidate 16/16: RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# 
# ── best model at step 1: ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#        step    covar var shape      objf deltObjf       AIC      BIC numParams  qchisqr     pchisqr included searchType        covNames covarEffect bsvReduction
# eta.vc    1 BW_power  vc power -4419.326 23.40874 -3519.145 -3481.58        10 3.841459 1.30982e-06      yes    forward cov_BW_power_vc    1.099023     43.70168
# ℹ Dropping alternative shape(s) for BW~vc: BW_lin
# 
# ── accepted BW_power~vc ──
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:45 | + 12.6 min] 
# → Forward step 2, candidate 1/14: BW_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:45 | + 12.6 min] 
# → Forward step 2, candidate 2/14: BW_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:45 | + 12.6 min] 
# → Forward step 2, candidate 3/14: CrCL_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:45 | + 12.6 min] 
# → Forward step 2, candidate 4/14: CrCL_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:45 | + 12.6 min] 
# → Forward step 2, candidate 5/14: CrCL_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:45 | + 12.6 min] 
# → Forward step 2, candidate 6/14: CrCL_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:45 | + 12.6 min] 
# → Forward step 2, candidate 7/14: BMI_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:46 | + 12.6 min] 
# → Forward step 2, candidate 8/14: BMI_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:46 | + 12.6 min] 
# → Forward step 2, candidate 9/14: BMI_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:51 | + 12.7 min] 
# → Forward step 2, candidate 10/14: BMI_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:51 | + 12.7 min] 
# → Forward step 2, candidate 11/14: SEX_1 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:51 | + 12.7 min] 
# → Forward step 2, candidate 12/14: SEX_1 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:51 | + 12.7 min] 
# → Forward step 2, candidate 13/14: RACE_0 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:49:51 | + 12.7 min] 
# → Forward step 2, candidate 14/14: RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ── best model at step 2: ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#         step      covar var shape      objf deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType          covNames covarEffect bsvReduction
# eta.cl2    2 CrCL_power  cl power -4449.376 30.05001 -3547.195 -3505.457        11 3.841459 4.210461e-08      yes    forward cov_CrCL_power_cl    0.614296     33.31519
# ℹ Dropping alternative shape(s) for CrCL~cl: CrCL_lin
# 
# ── accepted CrCL_power~cl ──
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:55:58 | + 18.8 min] 
# → Forward step 3, candidate 1/12: BW_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:55:58 | + 18.8 min] 
# → Forward step 3, candidate 2/12: BW_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:55:58 | + 18.8 min] 
# → Forward step 3, candidate 3/12: CrCL_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:55:58 | + 18.8 min] 
# → Forward step 3, candidate 4/12: CrCL_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:55:58 | + 18.8 min] 
# → Forward step 3, candidate 5/12: BMI_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:55:58 | + 18.8 min] 
# → Forward step 3, candidate 6/12: BMI_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:55:58 | + 18.8 min] 
# → Forward step 3, candidate 7/12: BMI_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:55:58 | + 18.8 min] 
# → Forward step 3, candidate 8/12: BMI_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:57:12 | + 20.1 min] 
# → Forward step 3, candidate 9/12: SEX_1 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:57:12 | + 20.1 min] 
# → Forward step 3, candidate 10/12: SEX_1 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:57:12 | + 20.1 min] 
# → Forward step 3, candidate 11/12: RACE_0 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [09:57:12 | + 20.1 min] 
# → Forward step 3, candidate 12/12: RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ── best model at step 3: ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#         step     covar var shape      objf deltObjf       AIC      BIC numParams  qchisqr      pchisqr included searchType         covNames covarEffect bsvReduction
# eta.cl2    3 BMI_power  cl power -4460.772 11.39668 -3556.591 -3510.68        12 3.841459 0.0007357546      yes    forward cov_BMI_power_cl   0.7450663     6.095665
# ℹ Dropping alternative shape(s) for BMI~cl: BMI_lin
# 
# ── accepted BMI_power~cl ──
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:02:14 | + 25.1 min] 
# → Forward step 4, candidate 1/10: BW_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:02:14 | + 25.1 min] 
# → Forward step 4, candidate 2/10: BW_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:02:14 | + 25.1 min] 
# → Forward step 4, candidate 3/10: CrCL_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:03:54 | + 26.8 min] 
# → Forward step 4, candidate 4/10: CrCL_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:03:54 | + 26.8 min] 
# → Forward step 4, candidate 5/10: BMI_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:03:54 | + 26.8 min] 
# → Forward step 4, candidate 6/10: BMI_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:03:54 | + 26.8 min] 
# → Forward step 4, candidate 7/10: SEX_1 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:03:54 | + 26.8 min] 
# → Forward step 4, candidate 8/10: SEX_1 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:03:54 | + 26.8 min] 
# → Forward step 4, candidate 9/10: RACE_0 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:03:54 | + 26.8 min] 
# → Forward step 4, candidate 10/10: RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ── best model at step 4: ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#         step covar var shape      objf deltObjf       AIC       BIC numParams  qchisqr      pchisqr included searchType     covNames covarEffect bsvReduction
# eta.vc4    4 SEX_1  vc   cat -4472.988 12.21602 -3566.807 -3516.722        13 3.841459 0.0004738083      yes    forward cov_SEX_1_vc   0.2795625     15.22219
# 
# ── accepted SEX_1~vc ──
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:11:31 | + 34.4 min] 
# → Forward step 5, candidate 1/9: BW_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:11:31 | + 34.4 min] 
# → Forward step 5, candidate 2/9: BW_lin ~ cl [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:11:31 | + 34.4 min] 
# → Forward step 5, candidate 3/9: CrCL_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:11:31 | + 34.4 min] 
# → Forward step 5, candidate 4/9: CrCL_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:11:31 | + 34.4 min] 
# → Forward step 5, candidate 5/9: BMI_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:11:31 | + 34.4 min] 
# → Forward step 5, candidate 6/9: BMI_lin ~ vc [lin]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:11:41 | + 34.6 min] 
# → Forward step 5, candidate 7/9: SEX_1 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:11:41 | + 34.6 min] 
# → Forward step 5, candidate 8/9: RACE_0 ~ cl [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:11:41 | + 34.6 min] 
# → Forward step 5, candidate 9/9: RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# [10:11:41 | + 34.6 min] 
# ── OFV did not improve, exiting forward search ... ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 
# [10:11:41 | + 34.6 min] 
# ── forward search complete ──
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
# $final_est
# # A tibble: 13 × 2
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
# 13 VcSEX      NA     
# 
# $rel_err
# # A tibble: 13 × 6
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
#  9 VcSEX          0.405  NA      NA        NA          NA       
# 10 var_CL         0.1     0.110   9.94e- 3  9.94e- 2    9.94e+ 0
# 11 var_Vc         0.1     0.0869 -1.31e- 2 -1.31e- 1   -1.31e+ 1
# 12 cov_VcCL       0.02    0.0277  7.68e- 3  3.84e- 1    3.84e+ 1
# 13 ResErr         0.1     0.106   6.26e- 3  6.26e- 2    6.26e+ 0
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
# 
# $raw
# $raw$summaryTable
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
# $raw$resFwd
# $raw$resFwd[[1]]
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
# $raw$resFwd[[2]]
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
# $raw$resFwd[[3]]
#    var      covar        type   center raw_col level has_missing missing_fill missing_check shape                    covExpr      init lower upper
# 3   vc   BW_power  continuous  76.6000      BW  <NA>       FALSE            0          <NA> power               log(BW/76.6) 0.9934590    -5     5
# 5   cl CrCL_power  continuous 100.1423    CrCL  <NA>       FALSE            0          <NA> power log(CrCL/100.142331480428) 0.4420569    -5     5
# 9   cl  BMI_power  continuous  27.1000     BMI  <NA>       FALSE            0          <NA> power              log(BMI/27.1) 0.6544207    -5     5
# 14  vc      SEX_1 categorical       NA     SEX     1       FALSE            0          <NA>   cat                      SEX_1 0.2795625    -5     5
# 
# 
# $raw$resBck
# NULL
# 
# attr(,"elapsed_s")
# [1] 2074.707
log(1.5)
# [1] 0.4054651
extract_params_long <- function(fit, includeCov = TRUE) {
  theta <- fit$theta
  om    <- fit$omega

  ## Fixed-effects: keep the "TV" prefix so names line up with true_params_long()
  ## (TVCL, TVVc, TVQ, TVVp, TVKA).  fit$theta still stores them on the LOG
  ## scale (lTVCL etc.), so we exp() back to natural-scale typical values.
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

  ## Covariate effects: read from fit$theta if present, else NA.
  ##   The theta names emitted by runSCM() preserve the original case of the
  ##   data column AND insert the SHAPE (power/lin) for continuous covariates
  ##   or the LEVEL VALUE (e.g. "_1" for SEX==1, "_2" for RACE==2) for
  ##   categorical covariates -- there is NO literal "_cat_" token.  The
  ##   regexes below mirror these two naming conventions exactly.
  ##
  ##   Naming conventions (verified against runSCM output):
  ##     continuous : cov_<COVAR>_<SHAPE>_<VAR>     e.g. cov_BW_power_cl
  ##     categorical: cov_<COVAR>_<LEVEL>_<VAR>     e.g. cov_SEX_1_vc
  ##
  ##   Extra keys (CLBMI, VcBMI, VcCrCL, VcRACE) are included so that
  ##   collinearity-driven false positives at small N (e.g. BMI->CL stealing
  ##   the slot of BW->CL) are surfaced in the wide compare table instead of
  ##   being silently dropped.
  cov_names_in_fit <- names(theta)
  cov_map <- list(
    # --- True relations in the scenario-16 simulator ---------------------
    CLBW   = c("TH_BW_CL",    grep("^cov_(bw|wt)_power_cl$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    CLcrCL = c("TH_CRCL_CL",  grep("^cov_crcl_power_cl$",      cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcBW   = c("TH_BW_VC",    grep("^cov_(bw|wt)_power_vc$",   cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    VcSEX  = c("TH_SEX_VC",   grep("^cov_sex_[^_]+_vc$",       cov_names_in_fit, value = TRUE, ignore.case = TRUE)),
    # --- Distractors / collinear false positives -------------------------
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
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities                                                        
# → calculate ∂(f)/∂(η)                                                            
# → calculate ∂(R²)/∂(η)                                                           
# → finding duplicate expressions in inner model...                                
# → optimizing duplicate expressions in inner model...                             
# → finding duplicate expressions in EBE model...                                  
# → optimizing duplicate expressions in EBE model...                               
# → compiling inner model...                                                       
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...                                
# → compiling EBE model...                                                         
# ✔ done
# → compiling events FD model...
# ✔ done
# calculating covariance matrix
# [====|====|====|====|====|====|====|====|====|====] 100%; 0:00:38 done
# 
# [10:40:11 | +  1.9 min] 
# ── starting backward search... ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:46:17 | +  8.0 min] 
# → Backward step 1, candidate 1/7: BW_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:46:17 | +  8.0 min] 
# → Backward step 1, candidate 2/7: CrCL_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:48:00 | +  9.7 min] 
# → Backward step 1, candidate 3/7: BW_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:48:00 | +  9.7 min] 
# → Backward step 1, candidate 4/7: BMI_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:48:00 | +  9.7 min] 
# → Backward step 1, candidate 5/7: CrCL_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:48:00 | +  9.7 min] 
# → Backward step 1, candidate 6/7: SEX_1 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# rxode2 5.1.3 using 4 threads (see ?getRxThreads)
#   no cache: create with `rxCreateCache()`
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:48:00 | +  9.7 min] 
# → Backward step 1, candidate 7/7: RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# 
# ── removing covariate at step 1: ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#         step     covar var shape      objf     deltObjf       AIC       BIC numParams  qchisqr pchisqr included searchType         covNames covarEffect bsvReduction
# eta.vc1    1 BMI_power  vc power -4474.449 -0.007083477 -3564.268 -3505.835        15 6.634897       1  dropped   backward cov_BMI_power_vc   0.1259127    -4.932421
# 
# ── removed BMI_power~vc ──
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:52:43 | + 14.4 min] 
# → Backward step 2, candidate 1/6: BW_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:52:43 | + 14.4 min] 
# → Backward step 2, candidate 2/6: CrCL_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:53:10 | + 14.9 min] 
# → Backward step 2, candidate 3/6: BW_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:53:10 | + 14.9 min] 
# → Backward step 2, candidate 4/6: CrCL_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:53:10 | + 14.9 min] 
# → Backward step 2, candidate 5/6: SEX_1 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:53:10 | + 14.9 min] 
# → Backward step 2, candidate 6/6: RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ── removing covariate at step 2: ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#         step      covar var shape      objf  deltObjf       AIC       BIC numParams  qchisqr   pchisqr included searchType          covNames covarEffect bsvReduction
# eta.vc1    2 CrCL_power  vc power -4474.088 0.3605457 -3565.907 -3511.648        14 6.634897 0.5482033  dropped   backward cov_CrCL_power_vc  0.05910443     5.495205
# 
# ── removed CrCL_power~vc ──
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:56:54 | + 18.6 min] 
# → Backward step 3, candidate 1/5: BW_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:56:54 | + 18.6 min] 
# → Backward step 3, candidate 2/5: CrCL_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:56:54 | + 18.6 min] 
# → Backward step 3, candidate 3/5: BW_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:57:35 | + 19.3 min] 
# → Backward step 3, candidate 4/5: SEX_1 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:57:35 | + 19.3 min] 
# → Backward step 3, candidate 5/5: RACE_0 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ── removing covariate at step 3: ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#         step  covar var shape      objf   deltObjf      AIC       BIC numParams  qchisqr pchisqr included searchType      covNames covarEffect bsvReduction
# eta.vc2    3 RACE_0  vc   cat -4474.661 -0.5722571 -3568.48 -3518.394        13 6.634897       1  dropped   backward cov_RACE_0_vc -0.04371761    -4.370627
# 
# ── removed RACE_0~vc ──
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [10:59:45 | + 21.5 min] 
# → Backward step 4, candidate 1/4: BW_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:00:32 | + 22.2 min] 
# → Backward step 4, candidate 2/4: CrCL_power ~ cl [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:00:32 | + 22.2 min] 
# → Backward step 4, candidate 3/4: BW_power ~ vc [power]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# [11:00:32 | + 22.2 min] 
# → Backward step 4, candidate 4/4: SEX_1 ~ vc [cat]
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#  
#  
# → loading into symengine environment...
# → pruning branches (`if`/`else`) of full model...
# ✔ done
# → calculate jacobian
# → calculate sensitivities
# → calculate ∂(f)/∂(η)
# → calculate ∂(R²)/∂(η)
# → finding duplicate expressions in inner model...
# → optimizing duplicate expressions in inner model...
# → finding duplicate expressions in EBE model...
# → optimizing duplicate expressions in EBE model...
# → compiling inner model...
#  
#  
# ✔ done
# → finding duplicate expressions in FD model...
# → optimizing duplicate expressions in FD model...
# → compiling EBE model...
#  
#  
# ✔ done
# → compiling events FD model...
#  
#  
# ✔ done
# [11:00:32 | + 22.2 min] 
# ── All remaining covariates significant (p <= 0.01), backward search complete. ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 
# [11:00:32 | + 22.2 min] 
# ── backward search complete ──
# 
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
# 
# $raw
# $raw$summaryTable
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
# $raw$resFwd
# NULL
# 
# $raw$resBck
# $raw$resBck[[1]]
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
# $raw$resBck[[2]]
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
# 
# attr(,"elapsed_s")
# [1] 1335.51