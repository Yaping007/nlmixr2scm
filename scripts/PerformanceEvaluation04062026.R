library(nlmixr2)
library(rxode2)
library(tidyverse)
library(devtools)
library(nhanesA)
library(haven)

#options(download.file.method = "wininet")
#options(repos = c(CRAN = "https://cloud.r-project.org"))
#install.packages("nhanesA")

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
## Creat 16 scenarios with different combinations of covariate effects----------
## BW on CL, CrCL on CL, BW on V, sex on V
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
## True parameter table per scenario (for RMRSE)
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

## Convenience wide table (one row per scenario, one column per parameter)
true_params_wide <- true_params %>%
  tidyr::pivot_wider(id_cols = scenario,
                     names_from = parameter,
                     values_from = true_value)


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

## Subset of valid datasets for downstream SSE/SCM analyses
valid_datasets <- eta_cor_per_dataset %>%
  dplyr::filter(pass_eta_cor) %>%
  dplyr::select(SCENARIO, DATASET)

sim_obs_valid <- sim_obs_all %>%
  dplyr::semi_join(valid_datasets, by = c("SCENARIO", "DATASET"))

saveRDS(sim_obs_valid, file.path(out_dir, "sim_obs_valid.rds"))
## Take-home message: ~58%-68% meet the condition of covCL_Vc0.15-0.25





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
## Convert to NONMEM-format dataset for nlmixr2 fitting / runSCM
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

## Example: NM-format dataset for scenario 1 only
nm_scenario_01 <- to_nm_dataset(sim_obs_list$scenario_01)
saveRDS(nm_scenario_01, file.path(out_dir, "nm_scenario_01.rds")) 
