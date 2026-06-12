##Install remote & install nlmixr2utils and nlmixr2scm from GitHub
#options(repos = c(CRAN = "https://cran.rstudio.com"))
#options(download.file.method = "wininet")
#install.packages("remotes")
#remotes::install_github("kestrel99/nlmixr2utils")
#remotes::install_github("kestrel99/nlmixr2scm")
##updated all R packages. 


##1. Run test suites first for nlmixr2scm
## NOTE: ~/.Rprofile sets options(rxode2.cache.dir = ...);
##       ~/.Renviron sets TMPDIR (must be in .Renviron, not .Rprofile,
##       because R locks tempdir() at start-up).
##       Open them manually only if you need to inspect / edit:
# usethis::edit_r_profile(scope = "user")
# file.edit("C:/Users/LIUYA8J/.Renviron")

library(testthat)
library(nlmixr2scm)
library(nlmixr2utils)
#testthat::test_file("tests/testthat/test-scm.R") #Pass228, warning from nlmixrest[near-singular cov on toy fixture]
#testthat::test_file("tests/testthat/test-parsing.R")#pass29, after fixing the function calling issue


## Run once per R session, before any nlmixr2() fit
##   With dev rxode2 (PR-1072) + nlmixr2est (PR-664) the cross-DLL OpenMP
##   bug on Windows is fixed -- multi-threaded FOCEi is now safe. We just
##   drop any stale compiled DLLs from prior runs and let the packages pick
##   their own thread counts (rxode2 defaults to detectCores()).
rxode2::rxClean()
library(nlmixr2)
library(rxode2)
library(tidyverse)
library(devtools)
library(haven)

#BenchMark Article Journal of Pharmacokinetics and Pharmacodynamics (2019) 46:273–285
## rxode popPK model scenarios-------
scm_2cmt_oral_rx <- function() {
  ini({
    ## Structural fixed effects from the article
    TVCL <- 0.6      # L/h
    TVQ  <- 1.8      # L/h
    TVVc <- 20       # L
    TVVp <- 80       # L
    TVKA <- 0.7      # 1/h

    ## Reference covariate values
    ## Replace with the article specific values if using their exact NHANES covariate dataset
    BW_REF   <- 70
    CRCL_REF <- 95

    ## Covariate-effect parameters
    ## Set to article true values when reproducing their scenarios exactly
    TH_BW_CL   <- 0.75
    TH_CRCL_CL <- 0.5
    TH_BW_VC   <- 1.0
    TH_SEX_VC  <- 0.5

    ## Scenario switches
    ## 0 means covariate relationship absent
    ## 1 means covariate relationship present
    I_BW_CL   <- 0
    I_CRCL_CL <- 0
    I_BW_VC   <- 0
    I_SEX_VC  <- 0

    ## Between subject variability
    eta_vc + eta_cl ~ c(
      0.1,   # variance eta_vc sqrt(0.1)= 31.6% between-subject variability(CV)
      0.02,   # cov eta_vc eta_cl #0.20 pearson's correlation
      0.1    # variance eta_cl
    )

    ## Proportional residual standard deviation  ~10% proportion error (CV)
    prop_err <- 0.1
  })

  model({
    ## Covariate models
    CL_cov <- TVCL *
      (BW / BW_REF)^(TH_BW_CL * I_BW_CL) *
      (CrCL / CRCL_REF)^(TH_CRCL_CL * I_CRCL_CL)

    Vc_cov <- TVVc *
      (BW / BW_REF)^(TH_BW_VC * I_BW_VC) *
      (1 + TH_SEX_VC * SEX * I_SEX_VC) ## Try linearized way & this way; FOCEI should be similar for continous but interesting for categorical

    ## Individual parameters
    CL <- CL_cov * exp(eta_cl)
    Vc <- Vc_cov * exp(eta_vc)
    Q  <- TVQ
    Vp <- TVVp
    KA <- TVKA

    ## Micro-rate constants
    k10 <- CL / Vc
    k12 <- Q / Vc
    k21 <- Q / Vp

    ## Two-compartment oral absorption model
    d/dt(depot)   = -KA * depot
    d/dt(central) =  KA * depot - k10 * central - k12 * central + k21 * peripheral
    d/dt(peripheral) = k12 * central - k21 * peripheral

    ## Plasma concentration

    cp = central / Vc

    ## Observation model

    cp ~ prop(prop_err)
  })
}



##2021-2023 cycle covariate population for SCM simulation--------
#Covariates:
##   BW   = body weight, kg
##   BMI  = body mass index, kg/m2
##   CrCL = creatinine clearance, mL/min
##   SEX  = 0 female, 1 male
##   RACE = 0 Asian, 1 white
## Required NHANES tables:
##   DEMO_*   : demographics, including age, sex, race
##   BMX_*    : body measures, including weight and BMI
##   BIOPRO_* : biochemistry profile, including serum creatinine
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
  filter(CrCL > 0, CrCL < 250)   # remove implausible values

#"original observed" correlation the article preserves within ±0.05.
ref_cor <- cor(pop)
print(round(ref_cor, 3))

## Bootstrap 250 Datasets of 300 Subjects (Preserving Correlations)--------
set.seed(4242)
n_subj      <- 300
n_datasets  <- 250
tol         <- 0.05
sex_prop    <- mean(pop$SEX)

# Stratify by SEX to preserve proportion
males   <- pop %>% filter(SEX == 1)
females <- pop %>% filter(SEX == 0)
n_male  <- round(n_subj * sex_prop)
n_fem   <- n_subj - n_male

datasets <- list()
attempt  <- 0
while (length(datasets) < n_datasets) {
  attempt <- attempt + 1
  samp <- bind_rows(
    males[sample(nrow(males),     n_male, replace = TRUE), ],
    females[sample(nrow(females), n_fem,  replace = TRUE), ]
  )
  
  # Check correlation tolerance against reference
  diff_cor <- abs(cor(samp) - ref_cor)
  if (all(diff_cor[upper.tri(diff_cor)] <= tol)) {
    datasets[[length(datasets) + 1]] <- samp
  }
  if (attempt > 1e5) stop("Tolerance too strict — relax it.")
}
length(datasets)   #250

for (i in seq_along(datasets)) {
  datasets[[i]]$ID      <- 1:n_subj
  datasets[[i]]$DATASET <- i
}

virtual_pop <- bind_rows(datasets)
saveRDS(virtual_pop, "virtual_population_250x300.rds")
write.csv(virtual_pop, "virtual_population_250x300.csv", row.names = FALSE)
virtual_pop <- read_csv("virtual_population_250x300.csv")

#Sanity Check:
#1.correlation structure preseved
print(round(ref_cor, 3))
#all 250*300 patients
round(cor(virtual_pop[, c("BW","BMI","CrCL","SEX","RACE")]),3)
#by dataset
cor_by_ds <- virtual_pop %>%
  group_by(DATASET) %>%
  summarise(
    r_BW_BMI   = cor(BW, BMI),
    r_BW_CrCL  = cor(BW, CrCL),
    r_BMI_RACE = cor(BMI, RACE)
  )
summary(cor_by_ds)

#2 Descriptive Summary 
summary(virtual_pop)

virtual_pop <- read_csv("virtual_population_250x300.csv")


## Creat 16 scenarios with different combinations of covariate effects----------
## BW on CL, CrCL on CL, BW on V, sex on V
PsN_scenarios <- data.frame(
  scenario = 1:16,
  I_BW_CL = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
  I_CRCL_CL = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
  I_BW_VC = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
  I_SEX_VC = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
)

# simulating individual PK parameters with η_CL, η_Vc ~ MVN with correlation 0.2 and concentration 
## ============================================================================
## Simulate using `scm_2cmt_oral_rx()` directly  with rxSolve
##   - rxode2 samples eta_cl, eta_vc internally from the model's omega
##     (var = 0.1 each, cov = 0.02  =>  correlation = 0.2 as requested
##   - Scenario indicators are switched on/off via `ini()` updater
##   - Per-subject covariates (BW, CrCL, SEX) supplied via `iCov`
##   - Sampling at 0, 0.05, 0.1, 0.5, 1, 3 x typical terminal half-life
## ============================================================================

## Example: scenario 1: No covariate effects--------------
theta_s1 <-  PsN_scenarios %>% filter(scenario == 1) %>% select(-scenario)
params_s1 <- c(
  I_BW_CL   = theta_s1$I_BW_CL,
  I_CRCL_CL = theta_s1$I_CRCL_CL,
  I_BW_VC   = theta_s1$I_BW_VC,
  I_SEX_VC  = theta_s1$I_SEX_VC)

## ---- (a) Build the scenario-specific model from the UI function ----------
mod    <- scm_2cmt_oral_rx()
mod_s1 <- mod %>% rxode2::ini(
  I_BW_CL   = theta_s1$I_BW_CL,
  I_CRCL_CL = theta_s1$I_CRCL_CL,
  I_BW_VC   = theta_s1$I_BW_VC,
  I_SEX_VC  = theta_s1$I_SEX_VC
)

## ---- (b) Typical terminal half-life (beta phase of 2-cmt model) ----------
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

## ---- (c) Per-subject covariate table (iCov) ------------------------------
n_total <- nrow(virtual_pop)
iCov <- virtual_pop %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  dplyr::select(id, DATASET, SUBJECT = ID, BW, CrCL, SEX)

## ---- (d) Event table: 100 mg single oral dose + sampling times -----------
##   Build a one-subject template (1 dose + 6 observation times), then
##   replicate it for every subject so each `id` appears with the full
##   dose + sampling schedule. This matches the `id` column in iCov.
DOSE_MG <- 100
ev_one <- rxode2::et(amt = DOSE_MG, cmt = "depot", time = 0) %>%
  rxode2::et(time = sample_times) %>%
  as.data.frame()

ev <- ev_one %>%
  dplyr::slice(rep(dplyr::row_number(), n_total)) %>%
  dplyr::mutate(id = rep(seq_len(n_total), each = nrow(ev_one))) %>%
  dplyr::arrange(id, time)

## ---- (e) Solve ODE -----------------------------------------------------------
##   rxode2 samples eta_cl, eta_vc per subject from the model's omega
##   (variance = 0.1 each, covariance = 0.02  =>  rho [correlation]= 0.2 = covariance/(square root(varianceX*varianceY))).
set.seed(2026)
sim_raw <- rxode2::rxSolve(
  mod_s1,
  events     = ev,
  iCov       = iCov %>% dplyr::select(id, BW, CrCL, SEX),
  returnType = "tibble"
)

## ---- (f) Tidy output ------------------------------------------------------
##   `cp`  = IPRED (no residual error), determistically produced by `cp = central/Vc`
##   `sim` = observation with proportional residual error already applied by
##           rxode2 because the model has `cp ~ prop(prop_err)`.
##   Etas are not returned as columns; recover them from CL/CL_cov and Vc/Vc_cov.
sim_obs1 <- sim_raw %>%
  dplyr::filter(time %in% sample_times) %>%
  dplyr::mutate(
    eta_cl   = log(CL / CL_cov),
    eta_vc   = log(Vc / Vc_cov),
    cp_ipred = cp,
    cp_obs   = sim,
    HL_MULT  = round(time / t_half_typ, 4)
  ) %>%
  dplyr::left_join(
    iCov %>% dplyr::select(id, DATASET, SUBJECT),
    by = "id"
  ) %>%
  dplyr::select(DATASET, SUBJECT, HL_MULT, time,
                BW, CrCL, SEX,
                cp_ipred, cp_obs)


sim_obs_sum1 <- sim_obs1 %>% 
  dplyr::group_by(HL_MULT) %>%
  dplyr::summarise(
    n         = dplyr::n(),
    median_cp = stats::median(cp_ipred),
    p05       = stats::quantile(cp_ipred, 0.05),
    p95       = stats::quantile(cp_ipred, 0.95),
    .groups   = "drop"
  )


#Use this order:
#Simulate scenario 1 only.
#Fit base model to dataset 1.
#Fit base model to all 250 scenario-1 datasets.
#Run runscm() on scenario 1, dataset 1.
#Run runscm() on 10 scenario-1 datasets.
#Run runscm() on all 250 scenario-1 datasets.
#Then expand to scenarios 2 to 16.






















#### Generalized simulation for all 16 scenarios (function + loop)###----------------------
## ============================================================================
##   - Reuses the per-subject covariate table and event grid built above
##   - Per-scenario seed = 2026 + scenario  (independent eta draws)
##   - Output: one tibble per scenario, plus a combined `sim_obs_all`
##   - Saves each scenario to disk as RDS so the 16-scenario set can be
##     reloaded without re-simulating
## ============================================================================

#Strategy1:  use omega metrix to sample etas and residual error automically by the popPKmodel------
             #only half of datasets meet the eta correlation constraint, but it's straightforward and preserves the full variability in the model.
## 16 scenarios: BW on CL, CrCL on CL, BW on V, sex on V
PsN_scenarios <- data.frame(
  scenario = 1:16,
  I_BW_CL = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
  I_CRCL_CL = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
  I_BW_VC = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
  I_SEX_VC = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
)
##(i) Preperation of covariate & ev for simulation-----------
## ---- (a) Expanded iCov: include BMI and RACE for downstream SCM ---------
##   The model itself only references BW / CrCL / SEX, but BMI and RACE must
##   be carried in the analysis dataset so SCM can test them as covariates.
iCov_full <- virtual_pop %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  dplyr::select(id, DATASET, SUBJECT = ID,
                BW, BMI, CrCL, SEX, RACE)

## ---- (b) Typical terminal half-life (beta phase of 2-cmt model) ----------
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
## ---- (c) Event table: 100 mg single oral dose + sampling times -----------
##   Build a one-subject template (1 dose + 6 observation times), then
##   replicate it for every subject so each `id` appears with the full
##   dose + sampling schedule. This matches the `id` column in iCov.
DOSE_MG <- 100
ev_one <- rxode2::et(amt = DOSE_MG, cmt = "depot", time = 0) %>%
  rxode2::et(time = sample_times) %>%
  as.data.frame()

ev <- ev_one %>%
  dplyr::slice(rep(dplyr::row_number(), n_total)) %>%
  dplyr::mutate(id = rep(seq_len(n_total), each = nrow(ev_one))) %>%
  dplyr::arrange(id, time)


## ---- (ii) Scenario-simulation helper ------------------------------------
simulate_scenario <- function(scn,
                              base_mod      = mod,
                              scenarios     = PsN_scenarios,
                              icov          = iCov_full,
                              event_table   = ev,
                              sampling_grid = sample_times,
                              t_half        = t_half_typ,
                              seed          = 2026L + scn) {
  ## (a) switch covariate-effect indicators for this scenario
  scn_row <- scenarios %>% dplyr::filter(scenario == scn)
  if (nrow(scn_row) != 1L)
    stop("scenario ", scn, " not found in `scenarios`.")

  mod_scn <- base_mod %>% rxode2::ini(
    I_BW_CL   = scn_row$I_BW_CL,
    I_CRCL_CL = scn_row$I_CRCL_CL,
    I_BW_VC   = scn_row$I_BW_VC,
    I_SEX_VC  = scn_row$I_SEX_VC
  )

  ## (b) solve (rxode2 samples eta_cl, eta_vc from omega each call)
  set.seed(seed)
  sim_raw_scn <- rxode2::rxSolve(
    mod_scn,
    events     = event_table,
    iCov       = icov %>% dplyr::select(id, BW, CrCL, SEX),
    returnType = "tibble"
  )

  ## (c) tidy: keep IPRED (`cp`), residual-error obs (`sim`), recover etas
  sim_raw_scn %>%
    dplyr::filter(time %in% sampling_grid) %>%
    dplyr::mutate(
      SCENARIO = scn,
      eta_cl   = log(CL / CL_cov),
      eta_vc   = log(Vc / Vc_cov),
      cp_ipred = cp,
      cp_obs   = sim,
      HL_MULT  = round(time / t_half, 4)
    ) %>%
    dplyr::left_join(
      icov %>% dplyr::select(id, DATASET, SUBJECT, BMI, RACE),
      by = "id"
    ) %>%
    dplyr::select(SCENARIO, DATASET, SUBJECT, HL_MULT, time,
                  BW, BMI, CrCL, SEX, RACE,
                  eta_cl, eta_vc, CL, Vc,
                  cp_ipred, cp_obs)
}

## ---- (iii) Run all 16 scenarios -----------------------------------------
##   Each scenario produces ~450k rows (75000 subjects x 6 obs).
##   Save per scenario to keep memory footprint manageable.
out_dir <- "simulated_virtual_dataset" #with covariate & PK sampling date ready for model fitting.
if (!dir.exists(out_dir)) dir.create(out_dir)

sim_obs_list <- purrr::map(
  .x        = PsN_scenarios$scenario,
  .f        = function(scn) {
    message(sprintf("Simulating scenario %02d ...", scn))
    out <- simulate_scenario(scn)
    saveRDS(out, file.path(out_dir, sprintf("sim_obs_scenario_%02d.rds", scn)))
    out
  },
  .progress = TRUE
)
names(sim_obs_list) <- sprintf("scenario_%02d", PsN_scenarios$scenario)

## ---- (iv) Combine and persist the full 7.2M-record dataset --------------
sim_obs_all <- dplyr::bind_rows(sim_obs_list)
saveRDS(sim_obs_all, file.path(out_dir, "sim_obs_all_scenarios.rds"))
#R’s native single-object storage format.
sim_obs_all <- readRDS("sim_obs_all_scenarios.rds")

#write.csv(sim_obs_all, "sim_obs_all_scenarios.csv", row.names = FALSE) too large, didn't excute.

## ============================================================================
## Eta-correlation QC per (SCENARIO, DATASET)
## ----------------------------------------------------------------------------
##   Article constraint: keep only datasets where cor(eta_CL, eta_Vc) is in
##   [0.15, 0.25] (target 0.20). Below we COMPUTE this per dataset and tag
##   datasets that pass / fail. Datasets failing the constraint should be
##   excluded from SSE and SCM analyses (or replaced by oversampling).
## ============================================================================
ETA_RHO_LOW  <- 0.15
ETA_RHO_HIGH <- 0.25

eta_cor_per_dataset <- sim_obs_all %>%
  dplyr::distinct(SCENARIO, DATASET, SUBJECT, eta_cl, eta_vc) %>%
  dplyr::group_by(SCENARIO, DATASET) %>%
  dplyr::summarise(
    n_subj   = dplyr::n(),
    rho_eta  = stats::cor(eta_cl, eta_vc),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(
    pass_eta_cor = rho_eta >= ETA_RHO_LOW & rho_eta <= ETA_RHO_HIGH
  )

## How many datasets per scenario meet the constraint?
eta_cor_summary <- eta_cor_per_dataset %>%
  dplyr::group_by(SCENARIO) %>%
  dplyr::summarise(
    n_datasets   = dplyr::n(),
    n_pass       = sum(pass_eta_cor),
    pct_pass     = mean(pass_eta_cor) * 100,
    rho_mean     = mean(rho_eta),
    rho_sd       = stats::sd(rho_eta),
    .groups      = "drop"
  )

saveRDS(eta_cor_per_dataset, file.path(out_dir, "eta_cor_per_dataset.rds"))
saveRDS(eta_cor_summary,     file.path(out_dir, "eta_cor_summary.rds"))
## Take-home message: ~58%-68% meet the condition of covCL_Vc0.15-0.25




#Strategy2:  use the pre-sampled etas (with guaranteed correlation) and add residual error manually after solving the ODE with fixed etas (no random sampling)----------
# This is more complex but guarantees all datasets meet the eta correlation constraint.
## ============================================================================
## Rejection-sampled etas: guarantee 250 valid datasets per scenario
## ----------------------------------------------------------------------------
##   The internal `rxode2` eta sampling above leaves ~58-67% of datasets in
##   [0.15, 0.25]. Article requires *all 250* per scenario in that band.
##   Strategy:
##     1. Pre-sample (eta_cl, eta_vc) per (SCENARIO, DATASET) of 300 subjects
##        with rejection sampling on cor(eta_cl, eta_vc).
##     2. Run rxode2 with a model that takes etas as inputs (no random draws).
##     3. Add proportional residual error manually (replaces `cp ~ prop()`).
##   Result: 16 x 250 datasets, each with rho in [0.15, 0.25] by construction.
## ============================================================================

## ---- (i) Helper: rejection-sample one dataset of (eta_cl, eta_vc) --------

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

## ---- (ii) Build a (SCENARIO, DATASET, SUBJECT) eta table -----------------
##   16 scenarios x 250 datasets = 4000 mvrnorm calls (~1 s total).
##   Seed mixes scenario x dataset for reproducibility.
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

## QC: every (SCENARIO, DATASET) rho should be in [0.15, 0.25] by construction
eta_table_qc <- eta_table %>%
  dplyr::distinct(SCENARIO, DATASET, rho_realised, attempts) %>%
  dplyr::summarise(
    min_rho        = min(rho_realised),
    max_rho        = max(rho_realised),
    mean_rho       = mean(rho_realised),
    mean_attempts  = mean(attempts),
    max_attempts   = max(attempts)
  )

## ---- (iii) Simulation model with etas as INPUTS (no internal draws) ------
##   Mirrors `scm_2cmt_oral_rx()` exactly, but eta_cl / eta_vc enter through
##   `iCov` instead of being sampled from omega.
sim_mod_fixed <- rxode2::rxode2({
  CL_cov = TVCL * (BW   / BW_REF  )^(TH_BW_CL   * I_BW_CL  ) *
                  (CrCL / CRCL_REF)^(TH_CRCL_CL * I_CRCL_CL)
  Vc_cov = TVVc * (BW   / BW_REF  )^(TH_BW_VC   * I_BW_VC  ) *
                  (1 + TH_SEX_VC * SEX * I_SEX_VC)
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

## ---- (iv) Scenario simulator using the fixed-eta table ------------------
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
    TH_BW_VC   = 1.0,  TH_SEX_VC  = 0.5,
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

## ---- (v) Run all 16 scenarios ------------------------------------------
out_dir_v2 <- "simulated_virtual_dataset_eta_filtered"
if (!dir.exists(out_dir_v2)) dir.create(out_dir_v2)

sim_obs_list_v2 <- purrr::map(
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
names(sim_obs_list_v2) <- sprintf("scenario_%02d", PsN_scenarios$scenario)

sim_obs_all_v2 <- dplyr::bind_rows(sim_obs_list_v2)
saveRDS(sim_obs_all_v2, file.path(out_dir_v2, "sim_obs_all_scenarios.rds"))
saveRDS(eta_table,      file.path(out_dir_v2, "eta_table_valid.rds"))

## ---- (vi) QC: confirm 100% pass rate -----------------------------------
eta_cor_v2 <- sim_obs_all_v2 %>%
  dplyr::distinct(SCENARIO, DATASET, SUBJECT, eta_cl, eta_vc) %>%
  dplyr::group_by(SCENARIO, DATASET) %>%
  dplyr::summarise(rho_eta = stats::cor(eta_cl, eta_vc), .groups = "drop") %>%
  dplyr::mutate(pass = rho_eta >= ETA_RHO_LOW & rho_eta <= ETA_RHO_HIGH)

eta_cor_summary_v2 <- eta_cor_v2 %>%
  dplyr::group_by(SCENARIO) %>%
  dplyr::summarise(
    n_datasets = dplyr::n(),
    n_pass     = sum(pass),
    pct_pass   = mean(pass) * 100,
    rho_mean   = mean(rho_eta),
    rho_sd     = stats::sd(rho_eta),
    .groups    = "drop"
  )
saveRDS(eta_cor_summary_v2, file.path(out_dir_v2, "eta_cor_summary_v2.rds"))



## ============================================================================
## True parameter table per scenario (for RMRSE)----------
## ----------------------------------------------------------------------------
##   RMRSE = sqrt(mean(((estimate - true) / true)^2)) per parameter, per scenario.
##   Long format makes joining estimates (one row per parameter per fit) trivial.
##
##   Mapping (article name -> simulation symbol):
##     CL      <- TVCL              Vc       <- TVVc
##     Q       <- TVQ               Vp       <- TVVp
##     KA      <- TVKA              ResErr   <- prop_err  (proportional SD)
##     CLBW    <- TH_BW_CL  * I_BW_CL
##     CLcrCL  <- TH_CRCL_CL * I_CRCL_CL
##     VcBW    <- TH_BW_VC   * I_BW_VC
##     VcSEX   <- TH_SEX_VC  * I_SEX_VC
##     var_CL  <- omega(eta_cl, eta_cl) = 0.1
##     var_Vc  <- omega(eta_vc, eta_vc) = 0.1
##     cov_VcCL<- omega(eta_vc, eta_cl) = 0.02
## ============================================================================
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
                                       TH_BW_VC = 1.0, TH_SEX_VC = 0.5),
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


true_params <- readRDS(file.path(out_dir, "true_params_long.rds"))
## Convenience wide table (one row per scenario, one column per parameter)
true_params_wide <- true_params %>%
  tidyr::pivot_wider(id_cols = scenario,
                     names_from = parameter,
                     values_from = true_value)



## Use Scenario 9 [beta_WtCl=0.75]to test for nlmixr2 model refitting/runSCM
## Test every feature of runscm(): forward selection, backwald elimination, user-specified covariate; full covariate building
## work on a single dataset with 300 patients first, benchmark the runtime. 


## ============================================================================
## Convert to NONMEM-format dataset for nlmixr2 fitting / runSCM-------
## ----------------------------------------------------------------------------
##   For each (SCENARIO, DATASET, SUBJECT) we create:
##     - one EVID = 1 dose row at TIME = 0 with AMT = DOSE_MG, CMT = "depot"
##     - 6 EVID = 0 observation rows at sample_times with DV = cp_obs
##   ID is unique within a (SCENARIO, DATASET) fit; for a global ID use
##   paste(SCENARIO, DATASET, SUBJECT, sep = "_").
## ============================================================================

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


## ############################################################################
## STAGE 1 SMOKE TEST -- one dataset from scenario 9------------
## ############################################################################
##
## Goal of this stage
## ------------------
## (1) Confirm the simulation output can be re-fitted with the *true* model.
## (2) Quantify estimation bias / precision against `true_params` using relative squared error
## (3) Exercise every feature of `nlmixr2scm::runSCM()`:
##       a. forward selection only             searchType = "forward"
##       b. backward elimination only          searchType = "backward"
##       c. full SCM (forward then backward)   searchType = "scm"
##       d. user-specified candidate pairs     pairsVec = ...
##       e. full covariate building            (b) with all candidates pre-included
## (4) Benchmark wall-clock runtime for every step so the 250-dataset / 16-
##     scenario sweep can be sized.
##
## Why scenario 9
## --------------
## Scenario 9 = the simplest non-null scenario:
##   I_BW_CL = 1, I_CRCL_CL = 0, I_BW_VC = 0, I_SEX_VC = 0
## Only BW->CL has a true effect (TH_BW_CL = 0.75), so SCM should:
##   - keep CLBW (true positive)
##   - reject CLcrCL, VcBW, VcSEX (true negatives at the article p-values)
##
## Article parameters to recover (true_params)
## -------------------------------------------
##   Structural    : CL, Vc, Q, Vp, KA
##   Covariate     : CLBW (only nonzero in scenario 9)
##                   CLcrCL, VcBW, VcSEX (zero -> not estimated by base model)
##   Random effects: var_CL, var_Vc, cov_VcCL
##   Residual      : ResErr
##
## RMRSE
## -----
##   RMRSE_p = sqrt( mean_d ( ((est_p,d - true_p) / true_p)^2 ) )
## With one dataset this collapses to relative error=|est - true| / |true|; we compute it as
## a smoke test, then expand to 250 datasets in stage 2.
## ############################################################################

## ############################################################################
## STAGE 1 SMOKE TEST -- one dataset from scenario 9
## ############################################################################
##
## Analytical plan
## ---------------
## Part 1: Robustness of model refitting (true model)
##   1.1 Load scenario 9, dataset 1.
##   1.2 Convert to NONMEM-style dataset.
##   1.3 Fit the *true* scenario-9 model (no covariates EXCEPT BW->CL).
##   1.4 Extract parameter estimates, relative error vs `true_params`,
##       convergence status, OFV, runtime.
##
## Part 2: runSCM feature tests (base model = no covariates)
##   2.1 Fit base model (no covariates) -- this is the SCM starting point.
##   2.2 SCM forward selection only
##   2.3 SCM backward elimination only
##   2.4 User-specified single covariate relationship (BW on CL)
##   2.5 Full SCM (forward then backward)
##   For each SCM run, capture:
##       - final selected covariates
##       - SCM step history (summaryTable)
##       - final-model parameter estimates (with covariate coefficients)
##       - runtime
##
## Why scenario 9
## --------------
## Scenario 9 = simplest non-null scenario:
##   I_BW_CL = 1, I_CRCL_CL = 0, I_BW_VC = 0, I_SEX_VC = 0
## Only BW->CL has a true effect (TH_BW_CL = 0.75), so SCM should:
##   - keep CLBW (true positive)
##   - reject CLcrCL, VcBW, VcSEX, plus BMI/RACE on CL and Vc (true negatives)
##
## Article parameters (true_params, scenario 9)
## --------------------------------------------
##   Structural    : CL = 0.6, Vc = 20, Q = 1.8, Vp = 80, KA = 0.7
##   Covariate     : CLBW = 0.75   (only nonzero in scenario 9)
##                   CLcrCL = VcBW = VcSEX = 0
##   Random effects: var_CL = 0.1, var_Vc = 0.1, cov_VcCL = 0.02
##   Residual      : ResErr = 0.1
##
## Relative error per parameter (one dataset)
##   rel_err_p = (estimate_p - true_p) / true_p
## (The full RMRSE definition collapses to this when N_dataset = 1.)
## ############################################################################

## ============================================================================
## Part 1: Robustness of model refitting -- TRUE scenario-9 model
## ============================================================================
## ---- 1.1 / 1.2  Build NM-format dataset for SCENARIO = 9, DATASET = 1 ----
out_dir_v2 <- "simulated_virtual_dataset_eta_filtered"
stage1_dir <- file.path(out_dir_v2, "stage1_smoke_scn09_ds01")
if (!dir.exists(stage1_dir)) dir.create(stage1_dir, recursive = TRUE)

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

sim_obs_scn09 <- readRDS(file.path(out_dir_v2, "sim_obs_scenario_09.rds"))
ds01 <- to_nm_dataset(sim_obs_scn09) %>%
  dplyr::filter(DATASET == 1) %>%
  dplyr::select(-SCENARIO, -DATASET) %>%
  dplyr::mutate(
    ID   = as.integer(ID),
    SEX  = as.integer(SEX),
    RACE = as.integer(RACE)
  )
stopifnot(
  dplyr::n_distinct(ds01$ID) == 300,
  sum(ds01$EVID == 1) == 300,
  sum(ds01$EVID == 0) == 300 * length(sample_times)
)
saveRDS(ds01, file.path(stage1_dir, "nm_scn09_ds01.rds"))

out_dir_v2 <- "simulated_virtual_dataset_eta_filtered"
stage1_dir <- file.path(out_dir_v2, "stage1_smoke_scn09_ds01")
ds01 <- readRDS(file.path(stage1_dir, "nm_scn09_ds01.rds"))

## ---- 1.3  TRUE model: BW on CL covariate baked in -----------------------
##   Same structure as scm_2cmt_oral_rx() with the scenario-9 indicators
##   hard-wired. Used to re-estimate parameters and assess bias.
##   - BW reference value: 70 kg (matches simulation)
##   - Covariate term:  (BW / 70) ^ TH_BW_CL  with TH_BW_CL estimated.
##   - Fixed effects on log scale; covariate coefficient on natural scale.
true_2cmt_scn09_refexp <- function() {
  ini({
    lTVCL    <- log(0.6)
    lTVQ     <- log(1.8)
    lTVVc    <- log(20)
    lTVVp    <- log(80)
    lTVKA    <- log(0.7)
    TH_BW_CL <- 0.75    # power exponent for BW on CL

    eta.cl + eta.vc ~ c(
      0.1,
      0.02,
      0.1
    )

    prop.err <- 0.1
  })
  model({
    cl_typ <- exp(lTVCL) * (BW / 70)^TH_BW_CL  
    cl     <- cl_typ * exp(eta.cl)
    vc     <- exp(lTVVc + eta.vc)
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


## ---- 1.3 (b)  TRUE model -- LINEAR-ON-LOG-SCALE parameterisation --------
##   Algebraically identical to true_2cmt_scn09_refexp:
##     refexp:  cl_typ = exp(lTVCL) * (BW/70)^TH_BW_CL ;  cl = cl_typ*exp(eta.cl)
##     lin   :  lTVCL_typ = lTVCL + TH_BW_CL*log(BW/70) ; cl = exp(lTVCL_typ + eta.cl)
##   Both expand to  exp( lTVCL + TH_BW_CL*log(BW/70) + eta.cl ).
##   The `lin` form is the more idiomatic nlmixr2 parameterisation: FOCEI
##   gradients are computed directly on log-scale, which is usually more
##   numerically stable.  Fitting BOTH provides a sanity check on the
##   simulation-then-estimation loop -- estimates should match to within
##   numerical tolerance.  Disagreement flags a FOCEI-stability issue
##   in one of the parameterisations.
true_2cmt_scn09_lin <- function() {
  ini({
    lTVCL    <- log(0.6)
    lTVQ     <- log(1.8)
    lTVVc    <- log(20)
    lTVVp    <- log(80)
    lTVKA    <- log(0.7)
    TH_BW_CL <- 0.75    # log-scale slope of BW on CL (== power exponent)

    eta.cl + eta.vc ~ c(
      0.1,
      0.02,
      0.1
    )

    prop.err <- 0.1
  })
  model({
    lTVCL_typ <- lTVCL + TH_BW_CL * log(BW / 70)   # linear on log scale
    cl        <- exp(lTVCL_typ + eta.cl)
    vc        <- exp(lTVVc + eta.vc)
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

## ---- 1.4  Fit BOTH parameterisations + extract estimates / relative error
##   Strategy:
##     1. Fit refexp form (centre-then-power) on full ds01
##     2. Fit lin    form (linear-on-log)    on full ds01  -- canonical
##     3. Compare estimates side-by-side as a virtual-cohort sanity check
##
##   `fit_true` / `t_fit_true` are aliased to the LIN fit so the existing
##   downstream code (overall_runtime, etc.) works unchanged.

stage1_focei <- nlmixr2est::foceiControl(
  print      = 0,        # silent: per-iteration log is noisy in Positron
  calcTables = TRUE,     # keep IPRED / CWRES / NPDE tables
  covMethod  = "r,s",     # sandwich SEs (the article default)
  sigdig     = 4,
  foceiControl(outerOpt = "bobyqa")
)

t_fit_true_refexp <- system.time(
  fit_true_refexp <- nlmixr2(true_2cmt_scn09_refexp, ds01,
                             est = "focei", control = stage1_focei)
)
saveRDS(fit_true_refexp,
        file.path(stage1_dir, "fit_true_scn09_ds01_refexp.rds"))

t_fit_true_lin <- system.time(
  fit_true_lin <- nlmixr2(true_2cmt_scn09_lin, ds01,
                          est = "focei", control = stage1_focei)
)
saveRDS(fit_true_lin,
        file.path(stage1_dir, "fit_true_scn09_ds01_lin.rds"))


## Convergence diagnostics
diagnose_fit <- function(fit) {
  list(
    converged = !is.null(fit$objf) && is.finite(fit$objf),
    objf      = if (!is.null(fit$objf)) fit$objf else NA_real_,
    cond_num  = if (!is.null(fit$conditionNumber)) fit$conditionNumber else NA_real_,
    cov_ok    = isTRUE(!is.null(fit$cov) && all(is.finite(diag(fit$cov))))
  )
}
diag_true_refexp <- diagnose_fit(fit_true_refexp) #-16805
diag_true_lin    <- diagnose_fit(fit_true_lin) #-16345- not global mininum

## ---- Benchmark: parameterisation × foceiControl × n_reps -----------------
##   Cross both parameterisations (refexp, lin) with three control settings
##   to see whether the `sigdig = 3` convergence break is unique to `lin`
##   or also bites `refexp`. 6 configs × 3 reps = 18 fits (~30 min total).

benchmark_fit <- function(model, data, control,
                          par_name, setting_name, n_reps = 3L) {
  purrr::map_dfr(seq_len(n_reps), function(rep) {
    rxode2::rxClean()                          # clean state per rep
    t <- system.time(
      fit <- tryCatch(
        nlmixr2(model, data, est = "focei", control = control),
        error = function(e) NULL
      )
    )
    tibble::tibble(
      par_name  = par_name,
      setting   = setting_name,
      label     = paste(par_name, setting_name, sep = "_"),
      rep       = rep,
      elapsed_s = unname(t["elapsed"]),
      objf      = if (is.null(fit)) NA_real_ else fit$objf,
      converged = !is.null(fit) && is.finite(fit$objf)
    )
  })
}

bench_models <- list(
  refexp = true_2cmt_scn09_refexp,
  lin    = true_2cmt_scn09_lin
)

bench_controls <- list(
  sigdig3        = nlmixr2est::foceiControl(sigdig = 3, print = 0,
                                            calcTables = FALSE),
  sigdig4        = nlmixr2est::foceiControl(sigdig = 4, print = 0,
                                            calcTables = FALSE),
  sigdig4_bobyqa = nlmixr2est::foceiControl(sigdig = 4, outerOpt = "bobyqa",
                                            print = 0, calcTables = FALSE)
)

bench_grid <- tidyr::expand_grid(
  par_name     = names(bench_models),
  setting_name = names(bench_controls)
)

t_bench <- system.time({
  bench_results <- purrr::pmap_dfr(
    bench_grid,
    function(par_name, setting_name) {
      benchmark_fit(
        model        = bench_models[[par_name]],
        data         = ds01,
        control      = bench_controls[[setting_name]],
        par_name     = par_name,
        setting_name = setting_name
      )
    }
  )
})
saveRDS(bench_results, file.path(stage1_dir, "bench_results.rds"))
bench_results <- readRDS(file.path(stage1_dir, "bench_results.rds"))
## Per-(parameterisation, setting) summary
bench_summary <- bench_results |>
  dplyr::group_by(par_name, setting) |>
  dplyr::summarise(
    n          = dplyr::n(),
    mean_s     = mean(elapsed_s),
    median_s   = median(elapsed_s),
    sd_s       = stats::sd(elapsed_s),
    objf_med   = median(objf),
    n_converged = sum(converged),
    .groups    = "drop"
  ) |>
  dplyr::arrange(par_name, setting)
saveRDS(bench_summary, file.path(stage1_dir, "bench_summary.rds"))

##Also test another optimizer (lbfgsb3c) and combine with the original
##sigdig3 / sigdig4 / sigdig4_bobyqa runs into one unified bench_summary.
##
##  Reuses the existing `benchmark_fit()` defined above -- no need to clone it.
##  The original `bench_results` is reloaded from disk so this block is
##  idempotent: re-running it appends fresh lbfgsb3c rows without re-running
##  the 6 already-completed configurations.

bench_controls_lbfgs <- list(
  sigdig3_lbfgsb3c = nlmixr2est::foceiControl(sigdig = 3, outerOpt = "lbfgsb3c",
                                              print = 0, calcTables = FALSE),
  sigdig4_lbfgsb3c = nlmixr2est::foceiControl(sigdig = 4, outerOpt = "lbfgsb3c",
                                              print = 0, calcTables = FALSE)
)

bench_grid_lbfgs <- tidyr::expand_grid(
  par_name     = names(bench_models),
  setting_name = names(bench_controls_lbfgs)
)

t_bench_lbfgs <- system.time({
  bench_results_lbfgs <- purrr::pmap_dfr(
    bench_grid_lbfgs,
    function(par_name, setting_name) {
      benchmark_fit(
        model        = bench_models[[par_name]],
        data         = ds01,
        control      = bench_controls_lbfgs[[setting_name]],
        par_name     = par_name,
        setting_name = setting_name
      )
    }
  )
})
saveRDS(bench_results_lbfgs,
        file.path(stage1_dir, "bench_results_lbfgs.rds"))

## ---- Combine with the original 3-setting run ---------------------------
bench_results_all <- dplyr::bind_rows(
  readRDS(file.path(stage1_dir, "bench_results.rds")),
  bench_results_lbfgs
)
saveRDS(bench_results_all,
        file.path(stage1_dir, "bench_results_all.rds"))

## ---- One unified summary across all 5 settings -------------------------
bench_summary_all <- bench_results_all |>
  dplyr::group_by(par_name, setting) |>
  dplyr::summarise(
    n           = dplyr::n(),
    mean_s      = mean(elapsed_s),
    median_s    = stats::median(elapsed_s),
    sd_s        = stats::sd(elapsed_s),
    objf_med    = stats::median(objf),
    n_converged = sum(converged),
    .groups     = "drop"
  ) |>
  dplyr::arrange(par_name, setting)
saveRDS(bench_summary_all,
        file.path(stage1_dir, "bench_summary_all.rds"))

## Wide view: parameterisation on cols, setting on rows
bench_wide <- bench_summary_all |>
  dplyr::select(par_name, setting, median_s, objf_med) |>
  tidyr::pivot_wider(
    id_cols     = setting,
    names_from  = par_name,
    values_from = c(median_s, objf_med),
    names_glue  = "{par_name}_{.value}"
  )


## Map fit$theta + fit$omega -> article parameter names
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

  ## Covariate effects: read from fit$theta if present, else NA
  cov_names_in_fit <- names(theta)
  cov_map <- list(
    CLBW   = c("TH_BW_CL", grep("^cov_(bw|wt)_power_cl$",   cov_names_in_fit, value = TRUE)),
    CLcrCL = c("TH_CRCL_CL", grep("^cov_crcl_power_cl$",    cov_names_in_fit, value = TRUE)),
    VcBW   = c("TH_BW_VC", grep("^cov_(bw|wt)_power_vc$",   cov_names_in_fit, value = TRUE)),
    VcSEX  = c("TH_SEX_VC", grep("^cov_sex(_male)?_cat_vc$", cov_names_in_fit, value = TRUE))
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

est_true_refexp <- extract_params_long(fit_true_refexp)
est_true_lin    <- extract_params_long(fit_true_lin)

## Relative error per parameter (one dataset)
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

err_true_refexp <- rel_err_one(est_true_refexp, true_params, scenario_id = 9)
err_true_lin    <- rel_err_one(est_true_lin,    true_params, scenario_id = 9)

## Long-format side-by-side comparison
err_true_compare <- dplyr::bind_rows(
  err_true_refexp %>% dplyr::mutate(parameterisation = "refexp"),
  err_true_lin    %>% dplyr::mutate(parameterisation = "lin")
) %>%
  dplyr::select(parameterisation, parameter, true_value,
                estimate, abs_err, rel_err, rel_err_pct)

## Wide view: estimate + rel_err_pct per parameter, side-by-side
err_true_wide <- err_true_compare %>%
  tidyr::pivot_wider(
    id_cols     = c(parameter, true_value),
    names_from  = parameterisation,
    values_from = c(estimate, rel_err_pct),
    names_glue  = "{.value}_{parameterisation}"
  ) %>%
  dplyr::mutate(
    abs_diff_estimate     = abs(estimate_refexp - estimate_lin),
    abs_diff_rel_err_pct  = abs(rel_err_pct_refexp - rel_err_pct_lin)
  )

saveRDS(err_true_refexp,
        file.path(stage1_dir, "rel_err_true_scn09_ds01_refexp.rds"))
saveRDS(err_true_lin,
        file.path(stage1_dir, "rel_err_true_scn09_ds01_lin.rds"))
saveRDS(err_true_compare,
        file.path(stage1_dir, "rel_err_true_scn09_ds01_compare.rds"))
saveRDS(err_true_wide,
        file.path(stage1_dir, "rel_err_true_scn09_ds01_wide.rds"))


## ---- Per-model summary tibble (one row per parameterisation) -----------
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
    runtime_sec      = unname(t_x["elapsed"])
  )
}
part1_summary <- dplyr::bind_rows(
  .build_part1_row("refexp", diag_true_refexp, t_fit_true_refexp),
  .build_part1_row("lin",    diag_true_lin,    t_fit_true_lin)
)
saveRDS(part1_summary, file.path(stage1_dir, "part1_summary.rds"))

## ============================================================================
## Part 2: runSCM feature tests -- BASE (no-covariate) model as starting point
## ============================================================================
scm_focei <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = ""         # SCM doesn't need cov matrix for LRT
)
## ---- 2.1  Base model (no covariates) ------------------------------------
base_2cmt_oral <- function() {
  ini({
    lTVCL <- log(0.6)
    lTVQ  <- log(1.8)
    lTVVc <- log(20)
    lTVVp <- log(80)
    lTVKA <- log(0.7)

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
saveRDS(fit_base, file.path(stage1_dir, "fit_base.rds"))


## ---- Helper: package one runSCM result for downstream comparison --------
##   Returns:
##     - selected:    accepted (var, covar, shape) pairs
##     - step_hist:   full SCM step history (forward + backward summary)
##     - final_est:   extract_params_long() of the final model
##     - rel_err:     relative error vs true scenario-9 parameters
##     - runtime_sec: wall-clock seconds for the runSCM call
##     - diag:        convergence diagnostics for the final fit
package_scm_result <- function(label, scm_res, runtime_sec,
                               true_long = true_params, scenario_id = 9) {
  ## Pick the final fit (backward if available, else forward, else base)
  final_fit <-
    if (!is.null(scm_res$resBck) && !is.null(scm_res$resBck$finalFit)) {
      scm_res$resBck$finalFit
    } else if (!is.null(scm_res$resFwd) && !is.null(scm_res$resFwd$finalFit)) {
      scm_res$resFwd$finalFit
    } else {
      NULL
    }

  ## Selected pairs from summaryTable (rows where the relation is in final model)
  selected <- if (!is.null(scm_res$summaryTable)) {
    st <- as.data.frame(scm_res$summaryTable)
    keep_col <- intersect(c("inFinal", "accepted", "kept"), colnames(st))
    if (length(keep_col) >= 1L) st[as.logical(st[[keep_col[1]]]), , drop = FALSE]
    else st
  } else {
    NULL
  }

  final_est <- if (!is.null(final_fit)) extract_params_long(final_fit) else NULL
  rel_err   <- if (!is.null(final_est)) {
    rel_err_one(final_est, true_long, scenario_id)
  } else NULL
  diag      <- if (!is.null(final_fit)) diagnose_fit(final_fit) else NULL

  list(
    label        = label,
    selected     = selected,
    step_hist    = scm_res$summaryTable,
    final_est    = final_est,
    rel_err      = rel_err,
    diag         = diag,
    runtime_sec  = runtime_sec,
    raw          = scm_res
  )
}



## ---- Candidate covariate-parameter pairs --------------------------------
candidate_pairs_full <- list(
  list(var = "cl", covar = "BW",   shapes = "power"),
  list(var = "cl", covar = "CrCL", shapes = "power"),
  list(var = "cl", covar = "BMI",  shapes = "power"),
  list(var = "cl", covar = "RACE", shapes = "cat"),
  list(var = "cl", covar = "SEX",  shapes = "cat"),
  list(var = "vc", covar = "BW",   shapes = "power"),
  list(var = "vc", covar = "CrCL", shapes = "power"),
  list(var = "vc", covar = "BMI",  shapes = "power"),
  list(var = "vc", covar = "RACE", shapes = "cat"),
  list(var = "vc", covar = "SEX",  shapes = "cat")
)
candidate_pairs_test <- list(
  list(var = "cl", covar = "BW",   shapes = "power"),
  list(var = "cl", covar = "CrCL", shapes = "power"),
 list(var = "vc", covar = "BW",   shapes = "power"),
  list(var = "vc", covar = "CrCL", shapes = "power")
)

## ---- Wrapper: runSCM with wall-clock progress tracking -----------------
##   Adds three things on top of nlmixr2scm::runSCM():
##     1. Pre-flight banner: HH:MM start, N candidates, N workers, search type.
##     2. Live progress -- every internal cli message about a step / candidate
##        is prefixed with [HH:MM:SS | +M.m min] so you can see how long each
##        candidate fit took even when print = 0 hides FOCEi iteration logs.
##     3. Post-flight banner: HH:MM end, total elapsed in minutes + seconds.
##   All other arguments pass through to runSCM() unchanged.
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

## ---- 2.2  Forward selection only ----------------------------------------
res_fwd <- runSCM_traced(
  label       = "forward",
  fit         = fit_base,
  pairsVec    = candidate_pairs_test,
  searchType  = "forward",
  control     = scm_focei,    # slim control: no tables, no cov, sigdig=4 bobyqa
  saveModels  = FALSE,
  workers     = 3L,           # 4 cores: leave 1 free for OS / Positron
  print       = 100,           # FOCEi iteration progress every 100 iters
  maxRetries = 0L #no retries for this smoke test
)
t_fwd     <- attr(res_fwd, "elapsed_s")
test_fwd  <- package_scm_result("forward_only", res_fwd, t_fwd)



## ---- 2.3  Backward elimination only -------------------------------------
##   Start with all 10 candidate relations included, then prune.
res_bck <- runSCM_traced(
  label              = "backward",
  fit                = fit_base,
  pairsVec           = candidate_pairs_test,
  catvarsVec         = c("SEX", "RACE"),
  searchType         = "backward",
  includedRelations  = candidate_pairs_test,
  control            = scm_focei,
  saveModels         = FALSE,
  workers            = 3L,
  print              = 100,
  maxRetries = 0L
)
t_bck     <- attr(res_bck, "elapsed_s")
test_bck  <- package_scm_result("backward_only", res_bck, t_bck)


## ---- 2.4  User-specified single covariate relationship (BW on CL) -------
res_user <- runSCM_traced(
  label       = "user",
  fit         = fit_base,
  pairsVec    = list(list(var = "cl", covar = "BW", shapes = "power")),
  searchType  = "scm",
  control     = scm_focei,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100
  maxRetries = 0L
)
t_user     <- attr(res_user, "elapsed_s")
test_user  <- package_scm_result("user_specified_BWonCL", res_user, t_user)


## ---- 2.5  Full SCM (forward then backward) ------------------------------
res_full <- runSCM_traced(
  label       = "full_scm",
  fit         = fit_base,
  pairsVec    = candidate_pairs_full,
  catvarsVec  = c("SEX", "RACE"),
  searchType  = "scm",
  control     = scm_focei,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries = 0L
)
t_full     <- attr(res_full, "elapsed_s")
test_full  <- package_scm_result("full_scm", res_full, t_full)


## ---- Aggregate Part 2 results -------------------------------------------
scm_tests <- list(
  forward_only          = test_fwd,
  backward_only         = test_bck,
  user_specified_BWonCL = test_user,
  full_scm              = test_full
)
saveRDS(scm_tests, file.path(stage1_dir, "scm_tests.rds"))

## Per-test summary tibble
part2_summary <- purrr::map_dfr(scm_tests, function(x) {
  tibble::tibble(
    test           = x$label,
    n_selected     = if (is.null(x$selected)) NA_integer_ else nrow(x$selected),
    converged      = if (is.null(x$diag)) NA else x$diag$converged,
    objf           = if (is.null(x$diag)) NA_real_ else x$diag$objf,
    cov_step_ok    = if (is.null(x$diag)) NA else x$diag$cov_ok,
    cond_num       = if (is.null(x$diag)) NA_real_ else x$diag$cond_num,
    runtime_sec    = x$runtime_sec
  )
})
saveRDS(part2_summary, file.path(stage1_dir, "part2_summary.rds"))


## ============================================================================
## Stage-1 reporting
## ============================================================================
overall_runtime <- dplyr::bind_rows(
  tibble::tibble(step = "fit_true_model_refexp",
                 runtime_sec = unname(t_fit_true_refexp["elapsed"])),
  tibble::tibble(step = "fit_true_model_lin",
                 runtime_sec = unname(t_fit_true_lin["elapsed"])),
  tibble::tibble(step = "fit_base_model",
                 runtime_sec = unname(t_fit_base["elapsed"])),
  part2_summary %>% dplyr::transmute(step = paste0("scm_", test), runtime_sec)
) %>%
  dplyr::mutate(
    minutes        = runtime_sec / 60,
    proj_250_min   = minutes * 250,             # one scenario, all datasets
    proj_4000_hr   = minutes * 250 * 16 / 60    # all 16 scenarios x 250 datasets
  )
saveRDS(overall_runtime, file.path(stage1_dir, "overall_runtime.rds"))

cat("\n========== STAGE 1 SMOKE TEST -- summary ==========\n\n")
cat("Part 1 (true-model robustness, refexp + lin parameterisations):\n")
print(part1_summary)
cat("\nPart 1 relative error per parameter (long, refexp vs lin):\n")
print(err_true_compare)
cat("\nPart 1 estimate side-by-side (wide):\n")
print(err_true_wide)
cat("\nPart 1 extended -- all 6 configs (long, baseline + lbfgsb3c):\n")
print(err_true_compare_all, n = Inf)
cat("\nPart 1 extended -- per-config OFV / convergence:\n")
print(err_true_summary_all)
cat("\nPart 1 extended -- estimate side-by-side (wide, all configs):\n")
print(err_true_wide_all, n = Inf, width = Inf)
cat("\nPart 2 (runSCM feature tests):\n"); print(part2_summary)
cat("\nOverall runtime + projections:\n"); print(overall_runtime)