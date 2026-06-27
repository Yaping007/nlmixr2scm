##Install remote & install nlmixr2utils and nlmixr2scm from GitHub
#options(repos = c(CRAN = "https://cran.rstudio.com"))
options(download.file.method = "wininet")
#install.packages("remotes")
#remotes::install_github("kestrel99/nlmixr2utils")
remotes::install_github("kestrel99/nlmixr2scm")
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
devtools::load_all("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
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


##1.Virtual Patients_2021-2023 cycle NHANES covariate population for SCM simulation--------
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


#1.2Virtual patients (Strategy2) use the pre-sampled etas (with guaranteed correlation) and add residual error manually after solving the ODE with fixed etas (no random sampling)----------
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
PsN_scenarios <- data.frame(
  scenario = 1:16,
  I_BW_CL = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1),
  I_CRCL_CL = c(0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1),
  I_BW_VC = c(0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1),
  I_SEX_VC = c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
)

#(i) Preperation of covariate & ev for simulation-----------
## ---- (a) Expanded iCov: include BMI and RACE for downstream SCM ---------
##   The model itself only references BW / CrCL / SEX, but BMI and RACE must
##   be carried in the analysis dataset so SCM can test them as covariates.
n_total <- nrow(virtual_pop)
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


## ---- (i) Helper: rejection-sample one dataset of (eta_cl, eta_vc) --------
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
##
##   SEX convention: nlmixr2scm log-additive form `exp(theta * SEX)`, which
##   matches what `runSCM(shapes = "cat")` recovers.  PsN-style
##   `(1 + theta * SEX)` was used in an earlier draft and produces a
##   systematic ~19% bias in VcSEX recovery because nlmixr2scm fits the
##   log-additive form.  With `TH_SEX_VC = log(1.5)` the biological effect
##   is unchanged from the original design (female Vc / male Vc = 1.5x when
##   I_SEX_VC = 1).
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

## ---- (v) Run all 16 scenarios ------------------------------------------
out_dir_v2 <- "simulated_virtual_dataset_eta_filtered"
if (!dir.exists(out_dir_v2)) dir.create(out_dir_v2)

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
names(sim_pkobs_nlmixr2para) <- sprintf("scenario_%02d", PsN_scenarios$scenario)

sim_pkobs_nlmixr2para <- dplyr::bind_rows(sim_pkobs_nlmixr2para)
saveRDS(sim_pkobs_nlmixr2para, file.path(out_dir_v2, "sim_pkobs_nlmixr2para_16scenarios.rds"))

##This is the previous simulated datasets with 16 scenerios and (1+beta*sex) convention, which didnot align with nlmixr2 parameterization.
saveRDS(sim_obs_all_v2, file.path(out_dir_v2, "sim_obs_all_scenarios.rds"))

## ---- (vi) QC: confirm 100% pass rate -----------------------------------
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
  ) ##100% pass rate


## ============================================================================
# Sister cohorts: N=40 and N=80, same pipeline as N=300 ----------------------
## ----------------------------------------------------------------------------
##   Re-runs the exact bootstrap + eta-table + scenario-simulator pipeline at
##   smaller N.  Calls the same helpers as N=300 (sample_dataset_etas,
##   build_eta_table, sim_mod_fixed, simulate_scenario_v2), only the sample
##   size differs.  Output schema (columns + row layout) is identical to
##   sim_obs_all_scenarios.rds at N=300.
##   Covariate correlation (BW, BMI, CrCL, SEX, RACE): preserved by the same
##                                                     rejection bootstrap
##                                                     against ref_cor at
##                                                     tol = 0.05.
##   eta_cl / eta_vc covariance: preserved by sample_dataset_etas() (rho in
##                              [ETA_RHO_LOW, ETA_RHO_HIGH] by construction).

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

## --- Build the two sister cohorts -----------------------------------------
##   Wall-clock: ~3-5 min for N=40, ~6-10 min for N=80.
##   Tolerance is auto-scaled by sqrt(300 / n_subj) internally; bump
##   `tol_scale` if a cohort can't reach n_datasets within `max_attempts`.
cohort_40 <- build_sister_cohort(n_subj = 40L)
cohort_80 <- build_sister_cohort(n_subj = 80L)



## ============================================================================
# True parameter table per scenario (for RMRSE)----------
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

true_params <- readRDS(file.path(out_dir, "true_params_long.rds"))

## Use Scenario 9 [beta_WtCl=0.75]to test for nlmixr2 model refitting/runSCM
## Test every feature of runscm(): forward selection, backwald elimination, user-specified covariate; full covariate building
## work on a single dataset with 300 patients first, benchmark the runtime. 


## ============================================================================
# Convert to NONMEM-format dataset for nlmixr2 fitting / runSCM-------
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
#STAGE 1 SMOKE TEST -- one dataset from scenario 9--------
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
##       - convergence status (CN, successful minimization)
##
## Why scenario 9
## --------------
## Scenario 9 = simplest non-null scenario:
##   I_BW_CL = 1, I_CRCL_CL = 0, I_BW_VC = 0, I_SEX_VC = 0
## Only BW->CL has a true effect (TH_BW_CL = 0.75), so SCM should:
##   - keep CLBW (true positive)
##   - reject CLcrCL, VcBW, VcSEX, plus BMI/RACE on CL and Vc (true negatives)
## Relative error per parameter (one dataset)
##   rel_err_p = (estimate_p - true_p) / true_p
## (The full RMRSE definition collapses to this when N_dataset = 1.)
## ############################################################################

## ============================================================================
#Part 1: Robustness of model refitting -- TRUE scenario-9 model--------
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

## --- Build sister cohorts (N=40, N=80) ------------------------------------
##   Defined further up alongside the N=300 simulator; called here because
##   to_nm_dataset() and DOSE_MG must be in scope first.
##   Wall-clock: ~3-5 min for N=40, ~6-10 min for N=80.
cohort_40 <- build_sister_cohort(n_subj = 40L)
cohort_80 <- build_sister_cohort(n_subj = 80L)

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

## ---- 1.3 (a)  TRUE model: BW on CL covariate baked in (exponential as in referenced paper) -----------------------
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
    ## KA is structurally near-unidentifiable with Khandelwal's sampling
    ## design (first non-zero sample at ~7 h, absorption t1/2 ~ 1 h, so
    ## the absorption phase is invisible).  Estimating it inflates the
    ## cov-matrix condition number by ~20x without changing OFV.  Fix at
    ## its true value to match common warfarin-modelling practice.
    lTVKA    <- fix(log(0.7))
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
    ## KA fixed at true value -- see comment in true_2cmt_scn09_refexp().
    lTVKA    <- fix(log(0.7))
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
##   converged     : TRUE iff optimizer reported success (fit$convergence == 0)
##                   AND objf is finite.  fit$convergence is the canonical
##                   optimizer status flag (0 = success); the earlier
##                   is.finite(objf)-only test silently treated non-converged
##                   fits as converged whenever OFV happened to be finite.
##   cond_num      : raw condition number = lambda_max / lambda_min of the
##                   parameter cov matrix.  This is what nlmixr2 stores in
##                   fit$conditionNumber*.  Populated only when covMethod
##                   != "" was used at fit time.
##   cond_num_sqrt : sqrt(cond_num) = sqrt(lambda_max / lambda_min).  This
##                   is the NONMEM / Beal convention reported by $COV with
##                   PRINT=E.  Most PMX papers (incl. Khandelwal 2019)
##                   quote thresholds against this form -- e.g. < 1000 for
##                   a well-conditioned model.  Roughly: sqrt(35492) ~ 188.
##   cov_ok        : TRUE iff fit$cov is a finite-diagonal matrix.  Will be
##                   FALSE when covMethod = "" was used (cov not computed).
##   message       : optimizer exit message, e.g. "Normal exit from bobyqa".
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


## ---- Per-model summary tibble (diagnostic parameters) -----------
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
part1_summary <- dplyr::bind_rows(
  .build_part1_row("refexp", diag_true_refexp, t_fit_true_refexp),
  .build_part1_row("lin",    diag_true_lin,    t_fit_true_lin)
)
saveRDS(part1_summary, file.path(stage1_dir, "part1_summary.rds"))

## ============================================================================
#Part 2: runSCM feature tests -- BASE (no-covariate) model as starting point---------
## ============================================================================
scm_focei <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = ""         # SCM doesn't need cov matrix for LRT
)

## ---- Control for the post-SCM diagnostic refit ---------------------------

## ---- Helper: refit a final SCM-selected model with full diagnostics ------
##   Reuses the model UI baked into the SCM-final fit (which already includes
##   all retained covariate relations) and the same training data, but with
##   a richer foceiControl.  Output is a fully-instrumented nlmixr2 fit with
##     fit$cov, fit$conditionNumber, fit$parFixedDf (SE/%RSE/CI), etc.
##   Errors are caught and returned as NULL with a warning so a bad refit
##   doesn't blow up downstream packaging.
refit_final_model <- function(final_fit, control = final_focei) {
  if (is.null(final_fit)) return(NULL)
  tryCatch({
    nlmixr2(final_fit$ui, nlme::getData(final_fit),
            est = final_fit$est, control = control)
  }, error = function(e) {
    warning("refit_final_model() failed: ", conditionMessage(e),
            call. = FALSE)
    NULL
  })
}
## ---- 2.1  Base model (no covariates) ------------------------------------
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
saveRDS(fit_base, file.path(stage1_dir, "fit_base.rds"))

fit_base <- readRDS(file.path(stage1_dir, "fit_base.rds"))

t_fit_base_cov <- system.time(
  fit_base_cov <- nlmixr2(base_2cmt_oral, ds01,
                      est = "focei", control = final_focei)
)
saveRDS(fit_base_cov, file.path(stage1_dir, "fit_base_cov.rds"))
fit_base_cov <- readRDS(file.path(stage1_dir, "fit_base_cov.rds"))


## ---- Helper: package one runSCM result for downstream comparison --------
##   Returns:
##     - selected:        accepted (var, covar, shape) pairs
##     - step_hist:       full SCM step history (forward + backward summary)
##     - final_est:       extract_params_long() of the final model
##     - rel_err:         relative error vs true scenario parameters
##     - runtime_sec:     wall-clock seconds for the runSCM call
##     - diag:            convergence diagnostics for the final fit
##     - parFixed:        nlmixr2 parFixedDf (Estimate, SE, %RSE, CI) when a
##                         refit was performed; NULL otherwise
##     - refit_done:      TRUE iff refit_control was supplied AND the refit
##                         succeeded; final_fit is then the refit object
##
##   refit_control:
##     When NULL (default) the SCM-final fit (covMethod = "", no tables) is
##     used as-is -- fast, but SE / condition number are unavailable.
##     When non-NULL (typically `final_focei`) the SCM-final model is refit
##     ONCE with the supplied control, producing full diagnostics.  The
##     refitted fit REPLACES final_fit (we don't keep the pre-refit copy --
##     it would just be ~800 KB of duplicate object per .rds).
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
  ## We previously kept the refit under a separate `final_fit_refit` slot,
  ## but it was always identical to `final_fit` post-swap -- pure duplication.
  ## Now we just keep a boolean flag indicating whether the swap happened.
  refit_done <- FALSE
  if (!is.null(final_fit) && !is.null(refit_control)) {
    refit <- refit_final_model(final_fit, control = refit_control)
    if (!is.null(refit)) {
      final_fit  <- refit
      refit_done <- TRUE
    }
  }

  ## --- Final-model covariate relations -------------------------------------
  ## The previous `selected` filter looked at `summaryTable$included` (any
  ## "yes"/"retained" row in the trace) which double-counts: backward listed
  ## each retained relation once per backward step, and full-scm listed every
  ## forward-accepted PLUS every backward-retained relation.  The only
  ## reliable source-of-truth for "what's in the final model" is final_fit's
  ## own theta vector, since runSCM names covariate coefficients
  ##     cov_<COVAR>_<SHAPE-or-LEVEL>_<VAR>
  ## (e.g. cov_BW_power_cl, cov_SEX_1_vc).
  .parse_cov_theta <- function(nms) {
    hits <- grep("^cov_", nms, value = TRUE)
    if (length(hits) == 0L) {
      return(tibble::tibble(theta_name = character(),
                            covar      = character(),
                            shape      = character(),
                            var        = character()))
    }
    parts <- strsplit(sub("^cov_", "", hits), "_", fixed = TRUE)
    tibble::tibble(
      theta_name = hits,
      covar      = vapply(parts, `[`, character(1), 1L),
      shape      = vapply(parts, function(p) paste(p[-c(1L, length(p))],
                                                   collapse = "_"),
                          character(1)),
      var        = vapply(parts, function(p) p[length(p)], character(1))
    )
  }

  selected <- if (!is.null(final_fit)) {
    .parse_cov_theta(names(final_fit$theta)) %>%
      dplyr::mutate(estimate = unname(final_fit$theta[theta_name])) %>%
      dplyr::select(var, covar, shape, theta_name, estimate)
  } else {
    NULL
  }

  final_est <- if (!is.null(final_fit)) extract_params_long(final_fit) else NULL
  rel_err   <- if (!is.null(final_est)) {
    rel_err_one(final_est, true_long, scenario_id)
  } else NULL
  diag      <- if (!is.null(final_fit)) diagnose_fit(final_fit) else NULL
  ## parFixedDf is only populated when covMethod was non-empty; that's exactly
  ## the refit_done case, since the screening control sets covMethod = "".
  parFixed  <- if (refit_done) final_fit$parFixedDf else NULL

  ## --- Packaged result -----------------------------------------------------
  ## `step_hist` is the ONLY view of the search trace we keep.  Previously we
  ## also stored `raw = scm_res`, which contained `$summaryTable` (== step_hist)
  ## and `$resFwd[[2]]` / `$resBck[[2]]` (also == step_hist), so the same data
  ## was being persisted three times per saved .rds.  Drop `raw` entirely and
  ## promote `final_fit` to a top-level slot so downstream code can call e.g.
  ## `extract_params_long(test$final_fit)` without spelunking into `$raw`.
  list(
    label       = label,
    selected    = selected,
    step_hist   = scm_res$summaryTable,
    final_fit   = final_fit,
    final_est   = final_est,
    rel_err     = rel_err,
    diag        = diag,
    parFixed    = parFixed,
    refit_done  = refit_done,
    runtime_sec = runtime_sec
  )
}


## ---- Candidate contiuous covariate-parameter pairs --------------------------------
out_dir <- "simulated_virtual_dataset" #with covariate & PK sampling date ready for model fitting.
true_params <- readRDS(file.path(out_dir, "true_params_long.rds"))
out_dir_v2 <- "simulated_virtual_dataset_eta_filtered"
stage1_dir <- file.path(out_dir_v2, "stage1_smoke_scn09_ds01")
fit_base <- readRDS(file.path(stage1_dir, "fit_base.rds"))
fit_base_cov <- readRDS(file.path(stage1_dir, "fit_base_cov.rds"))

scm_focei <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
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

## ---- 2.2  Forward selection only + explicit pair ----------------------------------------
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
saveRDS(res_fwd, file.path(stage1_dir, "res_fwd.rds"))    # idempotent recovery
test_fwd  <- package_scm_result("forward_only", res_fwd, t_fwd)
saveRDS(test_fwd, file.path(stage1_dir, "test_fwd.rds"))
test_fwd <- readRDS(file.path(stage1_dir, "test_fwd.rds"))

## ---- 2.2.1  Forward selection only + auto-generated pair ----------------------------------------
res_fwd_auto <- runSCM_traced(
  label       = "forward",
  data        = ds01,
  fit         = fit_base,
  varsVec    = c("cl", "vc"),
  covarsVec  = "BW",
  catvarsVec = "SEX",
  shapes     = c("power", "lin"),
  searchType  = "forward",
  control     = scm_focei,    # slim control: no tables, no cov, sigdig=4 bobyqa
  saveModels  = FALSE,
  workers     = 3L,           # 4 cores: leave 1 free for OS / Positron
  print       = 100,           # FOCEi iteration progress every 100 iters
  maxRetries = 0L #no retries for this smoke test
)
t_fwd_auto     <- attr(res_fwd_auto, "elapsed_s")
saveRDS(res_fwd_auto, file.path(stage1_dir, "res_fwd_auto.rds"))    # idempotent recovery
res_fwd_auto <- readRDS(file.path(stage1_dir, "res_fwd_auto.rds"))

test_fwd_auto  <- package_scm_result("forward_only", res_fwd_auto, t_fwd_auto)
saveRDS(test_fwd_auto, file.path(stage1_dir, "test_fwd_auto.rds"))
test_fwd_auto <- readRDS(file.path(stage1_dir, "test_fwd_auto.rds"))

res_fwd_auto_cov <- runSCM_traced(
  label       = "forward",
  data        = ds01,
  fit         =  fit_base_cov,
  varsVec    = c("cl", "vc"),
  covarsVec  = "BW",
  catvarsVec = "SEX",
  shapes     = c("power", "lin"),
  searchType  = "forward",
  control     = final_focei,    # slim control: no tables, no cov, sigdig=4 bobyqa
  saveModels  = FALSE,
  workers     = 3L,           # 4 cores: leave 1 free for OS / Positron
  print       = 100,           # FOCEi iteration progress every 100 iters
  maxRetries = 0L #no retries for this smoke test
)

t_fwd_auto_cov     <- attr(res_fwd_auto_cov, "elapsed_s")
test_fwd_auto_cov  <- package_scm_result("forward_only_cov", res_fwd_auto_cov, t_fwd_auto_cov)
saveRDS(test_fwd_auto_cov, file.path(stage1_dir, "test_fwd_auto_cov.rds"))
test_fwd_auto_cov <- readRDS(file.path(stage1_dir, "test_fwd_auto_cov.rds"))


scm_focei_maxiteration <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl      = rxode2::rxControl(atol = 1e-8, rtol = 1e-6),
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
t_fwd_auto_maxiteration     <- attr(res_fwd_auto_maxiteration, "elapsed_s")
saveRDS(res_fwd_auto_maxiteration, file.path(stage1_dir, "res_fwd_auto_maxiteration.rds"))    # idempotent recovery
test_fwd_auto_maxiteration  <- package_scm_result("forward_only", res_fwd_auto_maxiteration, t_fwd_auto_maxiteration)
saveRDS(test_fwd_auto_maxiteration, file.path(stage1_dir, "test_fwd_auto_maxiteration.rds"))
test_fwd_auto_maxiteration <- readRDS(file.path(stage1_dir, "test_fwd_auto_maxiteration.rds"))

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
t_fwd_auto_maxiteration     <- attr(res_fwd_auto_maxiteration, "elapsed_s")
saveRDS(res_fwd_auto_maxiteration, file.path(stage1_dir, "res_fwd_auto_maxiteration.rds"))    # idempotent recovery
test_fwd_auto_maxiteration  <- package_scm_result("forward_only", res_fwd_auto_maxiteration, t_fwd_auto_maxiteration)
saveRDS(test_fwd_auto_maxiteration, file.path(stage1_dir, "test_fwd_auto_maxiteration.rds"))
test_fwd_auto_maxiteration <- readRDS(file.path(stage1_dir, "test_fwd_auto_maxiteration.rds"))





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
saveRDS(res_bck, file.path(stage1_dir, "res_bck.rds"))
test_bck  <- package_scm_result("backward_only", res_bck, t_bck)
saveRDS(test_bck, file.path(stage1_dir, "test_bck.rds"))


## ---- 2.4  User-specified single covariate relationship (BW on CL) -------
res_user <- runSCM_traced(
  label       = "user",
  fit         = fit_base,
  pairsVec    = list(list(var = "cl", covar = "BW", shapes = "power")),
  searchType  = "scm",
  control     = scm_focei,
  saveModels  = FALSE,
  workers     = 1L,
  print       = 100,
  maxRetries = 0L
)
t_user     <- attr(res_user, "elapsed_s")
saveRDS(res_user, file.path(stage1_dir, "res_user.rds"))
test_user  <- package_scm_result("user_specified_BWonCL", res_user, t_user)
saveRDS(test_user, file.path(stage1_dir, "test_user.rds"))


## ---- 2.5  Full SCM (forward then backward) ------------------------------
res_full <- runSCM_traced(
  label       = "full_scm",
  fit         = fit_base,
  pairsVec    = candidate_pairs_test,
  catvarsVec  = c("SEX", "RACE"),
  searchType  = "scm",
  control     = scm_focei,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries = 0L
)
t_full     <- attr(res_full, "elapsed_s")
saveRDS(res_full, file.path(stage1_dir, "res_full.rds"))
test_full  <- package_scm_result("full_scm", res_full, t_full)
saveRDS(test_full, file.path(stage1_dir, "test_full.rds"))

res_full_cov <- runSCM_traced(
  label       = "full_scm",
  fit         = fit_base_cov,
  pairsVec    = candidate_pairs_test,
  catvarsVec  = c("SEX", "RACE"),
  searchType  = "scm",
  control     = final_focei,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries = 0L
)

t_full_cov     <- attr(res_full_cov , "elapsed_s") #runtime 17min
saveRDS(res_full_cov , file.path(stage1_dir, "res_full_cov.rds"))
test_full_cov  <- package_scm_result("full_scm_cov", res_full_cov, t_full_cov)
saveRDS(test_full_cov, file.path(stage1_dir, "test_full_cov.rds"))


## ============================================================================
## ---- 2.6  Categorical-covariate smoke test (BW continuous + SEX cat) ----
## ============================================================================
##   Sections 2.2-2.5 covered only continuous covariates.  This block re-runs
##   the four search modes on a 2-candidate set:
##     - cl ~ BW   (power, TRUE positive  -- TH_BW_CL = 0.75 in scenario 9)
##     - cl ~ SEX  (cat,   TRUE negative  -- no SEX-on-CL effect in scenario 9)
##
##   Expectation:
##     forward      : keeps BW, rejects SEX
##     backward     : starting from both, drops SEX, retains BW
##     user-only    : single-relation fit on SEX, formally tested and rejected
##     full SCM     : forward picks BW; backward leaves it in -> final = {BW}
##
##   `catvarsVec = "SEX"` is REQUIRED -- it triggers .makeSCMData() which
##   builds the SEX_1 indicator column the package's `cat` shape consumes.
## ============================================================================

out_dir_v2 <- "simulated_virtual_dataset_eta_filtered"
stage1_dir <- file.path(out_dir_v2, "stage1_smoke_scn09_ds01")
fit_base <- readRDS(file.path(stage1_dir, "fit_base.rds"))
scm_focei_n <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = "",         # SCM doesn't need cov matrix for LRT
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
)


candidate_pairs_shape_cattest <- list(
  ## All four continuous shapes for BW on CL.  In scenario 9 the truth is
  ## power (TH_BW_CL = 0.75), so SCM should pick BW_power over BW_lin /
  ## BW_log / BW_identity by lower OFV.  All four shapes will dwarf base
  ## OFV because BW carries strong info; the interesting comparison is
  ## *between* the four shapes (which one minimises OFV).
  list(var = "cl", covar = "BW",  shapes = "power"),
  list(var = "cl", covar = "BW",  shapes = "log"),
  ## Categorical: SEX is 0/1.  Use "cat" -- "power" would give 0^theta
  ## for SEX = 0 (undefined / 0).  catvarsVec = "SEX" passed to runSCM
  ## triggers .makeSCMData() to build the SEX_1 indicator column.
  list(var = "cl", covar = "SEX", shapes = "cat")
)
## ---- 2.6a  Forward selection only ---------------------------------------
res_fwd_cat <- runSCM_traced(
  label       = "forward_cat",
  fit         = fit_base,
  pairsVec    = candidate_pairs_shape_cattest,
  catvarsVec  = "SEX",
  searchType  = "forward",
  control     = scm_focei,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
t_fwd_cat    <- attr(res_fwd_cat, "elapsed_s")
saveRDS(res_fwd_cat, file.path(stage1_dir, "res_fwd_cat.rds"))
test_fwd_cat <- package_scm_result("forward_cat", res_fwd_cat, t_fwd_cat)
saveRDS(test_fwd_cat, file.path(stage1_dir, "test_fwd_cat.rds"))


## ---- 2.6b  Backward elimination only ------------------------------------
res_bck_cat <- runSCM_traced(
  label              = "backward_cat",
  fit                = fit_base,
  pairsVec           = candidate_pairs_shape_cattest,
  catvarsVec         = "SEX",
  searchType         = "backward",
  includedRelations  = candidate_pairs_shape_cattest,
  control            = scm_focei,
  saveModels         = FALSE,
  workers            = 3L,
  print              = 100,
  maxRetries         = 0L
)
t_bck_cat    <- attr(res_bck_cat, "elapsed_s")
saveRDS(res_bck_cat, file.path(stage1_dir, "res_bck_cat.rds"))
test_bck_cat <- package_scm_result("backward_cat", res_bck_cat, t_bck_cat)
saveRDS(test_bck_cat, file.path(stage1_dir, "test_bck_cat.rds"))
test_bck_cat <- readRDS(file.path(stage1_dir, "test_bck_cat.rds"))

## ---- 2.6c  User-specified single relation (SEX on CL only) --------------
##   Tests that the categorical pipeline correctly fits and rejects a true
##   no-effect categorical relation, even when it's the *only* candidate.
res_user_cat <- runSCM_traced(
  label       = "user_cat",
  fit         = fit_base,
  pairsVec    = list(list(var = "cl", covar = "SEX", shapes = "cat")),
  catvarsVec  = "SEX",
  searchType  = "scm",
  control     = scm_focei,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
t_user_cat    <- attr(res_user_cat, "elapsed_s")
saveRDS(res_user_cat, file.path(stage1_dir, "res_user_cat.rds"))
test_user_cat <- package_scm_result("user_specified_SEXonCL", res_user_cat, t_user_cat)
saveRDS(test_user_cat, file.path(stage1_dir, "test_user_cat.rds"))


## ---- 2.6d  Full SCM (forward then backward) -----------------------------
res_full_cat <- runSCM_traced(
  label       = "full_scm_cat",
  fit         = fit_base,
  pairsVec    = candidate_pairs_shape_cattest,
  catvarsVec  = "SEX",
  searchType  = "scm",
  control     = scm_focei,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries  = 0L
)
t_full_cat    <- attr(res_full_cat, "elapsed_s")
saveRDS(res_full_cat, file.path(stage1_dir, "res_full_cat.rds"))
test_full_cat <- package_scm_result("full_scm_cat", res_full_cat, t_full_cat)
saveRDS(test_full_cat, file.path(stage1_dir, "test_full_cat.rds"))


## ############################################################################
# STAGE 1 -- Scenario 16: all four covariate effects active -------------------
## ############################################################################
##   Truth (Khandelwal 2019, scenario 16):
##     BW   on CL  (power)  TH_BW_CL   = 0.75
##     CrCL on CL  (power)  TH_CRCL_CL = 0.50
##     BW   on Vc  (power)  TH_BW_VC   = 1.00
##     SEX  on Vc  (cat) TH_SEX_VC  = 0.50
##   All four are TRUE positives.  The SCM smoke test must keep these and
##   reject the false-positive distractors (CrCL on Vc, SEX on CL).
##
##   Parallel structure to the scenario-9 block above:
##     Part 1 -- true-model robustness (refexp vs lin)
##     Part 2 -- runSCM feature tests (forward / backward / user / full)
##   All fits use scm_focei_n: bobyqa, sigdig 4, maxOuter/Inner = 2000,
##   covMethod = "" (LRT only needs OFV).  Re-fit with final_focei outside
##   this block if SE / %RSE / parFixedDf are needed for any final model.
## ############################################################################
## ============================================================================
# Part0: foceiControl tuning grid (scenario 16, ds01) ----------------------
## ============================================================================
##   Goal: find the fastest foceiControl setting that ALSO converges deeply.
##   In Part 1, refexp walked 244 OFV units away from its own best minimum
##   (bad_solves triggered -> ODE tolerances auto-relaxed -> optimizer
##   wandered in noise).  Lin converged cleanly (gap = 0.007).  The grid
##   probes the four knobs that govern this trade-off:
##
##     sigdig         -- outer-loop convergence threshold (lower => slower
##                       but deeper)
##     atol / rtol    -- ODE solver tolerances (tighter => slower but more
##                       stable gradients)
##     stickyRecalcN  -- # of bad solves before nlmixr2 auto-loosens
##                       tolerances (low => quick relax; high => protect
##                       precision at cost of robustness on hard subjects)
##
##   We treat convergence as a HARD CONSTRAINT and minimise time under it:
##
##       converged_well := |final_ofv - best_ofv_in_trace| < 0.5
##
##   Reported alongside: bad_solves (whether sticky-relax fired), not_at_min
##   (cosmetic final-pass discrepancy), and n_iter (outer-loop budget used).

focei_grid <- tibble::tribble(
  ~setting,           ~sigdig, ~atol,  ~rtol,  ~stickyRecalcN,
  "baseline",         4,       1e-8,   1e-6,   4,        # matches scm_focei_n
  "fast",             3,       1e-6,   1e-4,   4,        # speed-only
  "no_sticky",        4,       1e-8,   1e-6,   20,       # block auto-relax
  "tight_ode",        4,       1e-10,  1e-8,   4,        # depth via ODE
  "tight_no_sticky",  4,       1e-10,  1e-8,   20,       # belt + braces
  "fast_no_sticky",   3,       1e-6,   1e-4,   20        # fast but protected
)

## ---- Per-fit diagnostic helper --------------------------------------------
##   Captures the wall-clock cost AND the four convergence signals needed to
##   judge whether a setting is acceptable.  All counted in seconds; convert
##   to minutes at display time.
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

## ---- Run the grid: 6 settings x 2 parameterisations = 12 fits -------------
##   Wall-clock ~ 30-40 min on a modern laptop.  Each cell saved
##   incrementally so a crash mid-grid doesn't lose prior results.
##
##   NB: atol / rtol are ODE-solver knobs and must be routed through
##   rxode2::rxControl(); foceiControl() does NOT accept them directly in
##   the installed nlmixr2est version.  stickyRecalcN stays on foceiControl
##   because it governs the FOCEi-level "bad solves" auto-relax behaviour.
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

saveRDS(bench_results, file.path(stage1_dir16, "focei_tuning_results.rds"))

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
saveRDS(bench_results_300, file.path(stage1_dir16, "focei_tuning_results_300.rds"))


## ---- Prep ---------------------------------------------------------------

## NM-format dataset for SCENARIO = 16, DATASET = 1 (n=300)----
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
ds16_01 <- readRDS(file.path(stage1_dir16, "nm_scn16_ds01.rds"))

## NM-format dataset for SCENARIO = 16, DATASET = 1 (n=80)-------
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
ds16_01_N80 <- readRDS(file.path(stage1_dir16_N80, "nm_scn16_ds01_N80.rds"))


## NM-format dataset for SCENARIO = 16, DATASET = 1 (n=40)---------
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


## ============================================================================
# Part 1: True-model robustness (refexp vs lin) ------------------------------
## ============================================================================
##   refexp: parameterised on natural scale, exp() at fit boundaries
##   lin   : parameterised log-additively (canonical nlmixr2 form)
##   Both encode the same simulation algebra.  Disagreement between fits flags
##   a FOCEI-stability issue in one of the parameterisations.


scm_focei_n <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = "r,s",         # SCM doesn't need cov matrix for LRT
  maxOuterIterations = 2000,
  maxInnerIterations = 2000
)


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
    TH_SEX_VC  <- log(1.5) #0.405        # nlmixr2scm convention

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


## ---- Fit both parameterisations at three cohort sizes -----------------
##   Each fit gets its own variable name (no clobbering across N).  All six
##   results then flow through one registry tibble for diagnostics + estimates.

t_fit_true16_refexp300_n <- system.time(
  fit_true16_refexp300_n <- nlmixr2(true_2cmt_scn16_refexp, ds16_01,
                                  est = "focei", control = scm_focei_n)
)
t_fit_true16_lin300_n <- system.time(
  fit_true16_lin300_n <- nlmixr2(true_2cmt_scn16_lin, ds16_01,
                               est = "focei", control = scm_focei_n)
)
diag_true_refexp300_n <- diagnose_fit(fit_true16_refexp300_n) #-16805
diag_true_lin300_n   <- diagnose_fit(fit_true16_lin300_n) #-16345- not global mininum-become the same again. 

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

## Aliases for the downstream runtime-projection block (uses N=300 as ref).
t_fit_true16_refexp <- t_fit_true16_refexp300
t_fit_true16_lin    <- t_fit_true16_lin300

## ---- Single registry drives diagnostics + estimates -------------------
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

## ---- Per-fit diagnostics (one row per (N, parameterisation)) ----------
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
saveRDS(diag_summary16, file.path(stage1_dir16, "part1_summary.rds"))

# Alias for the reporting block at the bottom of the script.
part1_summary16 <- diag_summary16


## ---- Estimates -> long; rel-err vs scenario-16 truth ------------------
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
saveRDS(err_true16_compare, file.path(stage1_dir16, "part1_err_long.rds"))


## ---- Wide compare: refexp vs lin at each N (parameter x N grid) -------
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
##


## ---- Visualisation: |rel-err| vs cohort size, lin parameterisation -----
##   Two messages in one panel:
##     1) Every parameter's |rel-err| grows as N shrinks (left -> right).
##     2) cov_VcCL sits well above the rest at every N.
##   Design: one faint grey line per parameter (background reference) +
##   one bold red line for cov_VcCL (the headline finding).  Parameter
##   labels at the right edge (N = 40) identify the gray spaghetti without
##   a legend.
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
ggplot2::ggsave(
  file.path(stage1_dir16, "part1_err_lin_by_N.png"),
  p_err_true16_lin,
  width = 10, height = 5, dpi = 150
)


## ---- Side-by-side bars: |rel-err| per parameter, grouped by N ----------
##   The line chart above answers "which parameter degrades fastest?"; this
##   bar version answers "for parameter X, how big is the error at each N?".
##   Parameters are ordered from worst to best at N=40 so the eye lands on
##   cov_VcCL at the left; bars within a parameter are coloured by cohort
##   size with the largest cohort lightest (intuitive "more data = less
##   error" reading).
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

print(p_err_true16_lin_bar)
ggplot2::ggsave(
  file.path(stage1_dir16, "part1_err_lin_by_N_bar.png"),
  p_err_true16_lin_bar,
  width = 14, height = 5, dpi = 150
)


# Part 2: runSCM feature tests -- BASE (no-covariate) model -------------------
## ============================================================================
##   The base model `base_2cmt_oral` is reused unchanged (covariate-free);
##   only the data differs from the scenario-9 fit so we refit on ds16_01.

base_2cmt_oral <- function() {
  ini({
    lTVCL <- log(0.6)
    lTVQ  <- log(1.8)
    lTVVc <- log(20)
    lTVVp <- log(80)
    ## KA fixed: the design's first non-zero sample (~7 h) is well past
    ## absorption (KA = 0.7 /h, absorption t1/2 ~ 1 h), so KA cannot be
    ## informatively estimated from the data. 
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
t_fit_base16_N80 <- system.time(
  fit_base16_N80 <- nlmixr2(base_2cmt_oral, ds16_01_N80,
                            est = "focei", control = scm_focei_n)
)
saveRDS(fit_base16_N80, file.path(stage1_dir16_N80, "fit_base_N80.rds"))
fit_base16_N80 <- readRDS(file.path(stage1_dir16_N80, "fit_base_N80.rds"))


## ============================================================================
##   Part 2 (N=80): runSCM feature tests under a hard covariate scenario
## ----------------------------------------------------------------------------
##   Why N=80?  Article scenario 16 has 4 true relations (BW->CL, CrCL->CL,
##   BW->Vc, SEX->Vc).  We deliberately stress-test runSCM here with:
##     * smaller N (less data -> looser LRT, more borderline OFV swings)
##     * BMI in the covariate pool (strongly correlated with BW), to make
##       sure the algorithm doesn't lock in a redundant BMI relation when
##       a BW one already wins
##     * RACE in the covariate pool (3-level categorical, including a rare
##       level), to probe categorical-on-Vc handling
##   varsVec / covarsVec / catvarsVec drive the FULL Cartesian product when
##   runSCM is asked to auto-generate the candidate pool.  The user-curated
##   `candidate_pairs_test_corContinusous_cat` is used wherever the search
##   needs an explicit included/test pair list (backward + user).
## ============================================================================

candidate_pairs_test_corContinusous_cat <- list(
  list(var = "cl", covar = "BW",   shapes = "power"),
  list(var = "cl", covar = "CrCL", shapes = "power"),
  list(var = "vc", covar = "BW",   shapes = "power"),
  list(var = "vc", covar = "BMI",  shapes = "power"),
  list(var = "vc", covar = "CrCL", shapes = "power"),
  list(var = "vc", covar = "SEX",  shapes = "cat"),
  list(var = "vc", covar = "RACE", shapes = "cat")
)

# Knobs reused across all four runs.  Defined once so they stay in sync.
scm16_vars       <- c("cl", "vc")
scm16_covars     <- c("BW", "CrCL", "BMI")
scm16_catvars    <- c("SEX", "RACE")
scm16_shapes     <- c("power", "lin")

scm_focei_n <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = "r,s",         # SCM doesn't need cov matrix for LRT
  maxOuterIterations = 2000,
  maxInnerIterations = 2000
)

## ---- 2.2.1  Forward selection (auto-generated Cartesian product) -------
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
t16_N80_fwd_auto    <- attr(res16_N80_fwd_auto, "elapsed_s") #34.6min 
saveRDS(res16_N80_fwd_auto, file.path(stage1_dir16_N80, "res_fwd_auto.rds"))
test16_N80_fwd_auto <- package_scm_result(
  "scn16_N80_forward_auto", res16_N80_fwd_auto, t16_N80_fwd_auto,
  scenario_id = 16
)

saveRDS(test16_N80_fwd_auto, file.path(stage1_dir16_N80, "test_fwd_auto.rds"))


## ---- 2.3  Backward elimination only ------------------------------------
##   Start with the full curated pool included, then prune.  Article truth
##   keeps 4 (BW~cl, CrCL~cl, BW~vc, SEX~vc); the extra BMI~vc / CrCL~vc /
##   RACE~vc relations should all be dropped.  This probes:
##     * collinearity handling (BMI vs BW on Vc)
##     * weak true negative pruning (CrCL on Vc, RACE on Vc)
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
t16_N80_bck    <- attr(res16_N80_bck, "elapsed_s") # 22.3min
saveRDS(res16_N80_bck, file.path(stage1_dir16_N80, "res_bck.rds"))
test16_N80_bck <- package_scm_result(
  "scn16_N80_backward_only", res16_N80_bck, t16_N80_bck,
  scenario_id = 16
)
saveRDS(test16_N80_bck, file.path(stage1_dir16_N80, "test_bck.rds"))


## ---- 2.4  User-specified true relations  --------------
##   Largest-magnitude true positive (TH_BW_VC = 1.0); smoke-tests the
##   user-specified search path on the N=80 cohort.
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

t16_N80_user    <- attr(res16_N80_user, "elapsed_s") # 9.6 mins
saveRDS(res16_N80_user, file.path(stage1_dir16_N80, "res_user.rds"))
test16_N80_user <- package_scm_result(
  "scn16_N80_user_4tr", res16_N80_user, t16_N80_user,
  scenario_id = 16
)
saveRDS(test16_N80_user, file.path(stage1_dir16_N80, "test_user.rds"))

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

t16_N80_user_wr    <- attr(res16_N80_user_wr, "elapsed_s") # 9.6 mins
saveRDS(res16_N80_user_wr, file.path(stage1_dir16_N80, "res_user_wr.rds"))
test16_N80_user_wr <- package_scm_result(
  "scn16_N80_user_2wr", res16_N80_user_wr, t16_N80_user_wr,
  scenario_id = 16
)
saveRDS(test16_N80_user_wr, file.path(stage1_dir16_N80, "test_user.rds"))

## ---- 2.5  Full SCM (forward then backward, full pool incl. BMI/RACE) ---
##   End state should retain the 4 true relations and drop BMI/CrCL/RACE on
##   Vc.  At N = 80 some borderline relations may flicker; runSCM should
##   still converge without errors and produce a plausible final model.

scm_focei_production <- nlmixr2est::foceiControl(
  sigdig             = 3,
  outerOpt           = "bobyqa",
  print              = 0,
  calcTables         = FALSE,
  covMethod          = "r,s",
  stickyRecalcN      = 20,                      
  maxOuterIterations = 2000,
  maxInnerIterations = 2000,
  rxControl          = rxode2::rxControl(atol = 1e-8, rtol = 1e-6) #Even with this setting, dofv still can be negative
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
t16_N80_full    <- attr(res16_N80_full, "elapsed_s") #31.8min
saveRDS(res16_N80_full, file.path(stage1_dir16_N80, "res_full.rds"))
test16_N80_full <- package_scm_result(
  "scn16_N80_full_scm", res16_N80_full, t16_N80_full,
  scenario_id = 16
)
saveRDS(test16_N80_full, file.path(stage1_dir16_N80, "test_full.rds"))

## ---- Aggregate Part 2 --------------------------------------------------
scm_tests16_N80 <- list(
  forward_auto  = test16_N80_fwd_auto,
  backward_only = test16_N80_bck,
  user_4true   = test16_N80_user,
  user_2true   =test16_N80_user_wr,
  full_scm      = test16_N80_full
)

part2_summary16_N80 <- purrr::map_dfr(scm_tests16_N80, function(x) {
  tibble::tibble(
    test          = x$label,
    n_selected    = if (is.null(x$selected)) NA_integer_ else nrow(x$selected),
    converged     = if (is.null(x$diag)) NA else x$diag$converged,
    objf          = if (is.null(x$diag)) NA_real_ else x$diag$objf,
    cov_step_ok   = if (is.null(x$diag)) NA else x$diag$cov_ok,
    cond_num      = if (is.null(x$diag)) NA_real_ else x$diag$cond_num,
    cond_num_sqrt = if (is.null(x$diag)) NA_real_ else x$diag$cond_num_sqrt,
    runtime_min   = (x$runtime_sec)/60
  )
})
saveRDS(part2_summary16_N80, file.path(stage1_dir16_N80, "part2_summary.rds"))




##Try fit-refit approach to save time -----------------
## ============================================================================
##   Part 2 (FAST path): two-tier control + linCmt rewrite of the model
## ----------------------------------------------------------------------------
##   Stacks three speed levers on top of the original Part 2 block above.
##   The originals are untouched so wall-time, OFV and selections are
##   directly comparable via `part2_speedup16_N80` at the bottom.
##
##     1. covMethod="" + calcTables=FALSE during SCM screening (LRT only
##        needs OFV).  The surviving model is refit ONCE with covMethod="r,s"
##        + calcTables=TRUE via package_scm_result(refit_control=...), which
##        in turn calls refit_final_model() defined above (~L1292).
##
##     2. Loose tolerances on the screening control (atol=1e-6, rtol=1e-4).
##        Combined with stickyRecalcN=20, this keeps the per-subject solver
##        precision UNIFORM (no silent per-subject relax cascade), which
##        also addresses the "dofv can still be negative" symptom from
##        scm_focei_production.  Tight tolerances (atol=1e-8, rtol=1e-6)
##        are restored for the single final refit.
##
##     3. linCmt() replaces the 3-state ODE in all three models
##        (base_2cmt_oral, true_2cmt_scn16_refexp, true_2cmt_scn16_lin).
##        rxode2 auto-detects 2-cmt oral from {ka, cl, vc, q, vp}; the
##        analytical Jacobian eliminates ALL solver-noise diagnostics.
##        Empirical speedup on 1-3 cmt oral PK is ~5-15x per fit.
## ============================================================================

# ---- (1) Two-tier foceiControls --------------------------------------------
##   _screen: applied to every SCM candidate.  No cov, no tables, loose tols.
##   _final:  applied to ONE refit per SCM mode, via refit_control= in
##            package_scm_result().  Tight tols + full diagnostics.

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


# ---- (2) linCmt() rewrites --------------------------------------------------
##   `ini()` blocks are byte-identical to the ODE versions; only the model
##   body changes.  rxode2 auto-selects the analytical 2-cmt oral macro
##   because {ka, cl, vc, q, vp} are all bound and there is no `d/dt(...)`.
##   Dose CMT mapping (cmt=1 -> depot, cmt=2 -> central) matches the data
##   produced by `to_nm_dataset()`.

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

true_2cmt_scn16_refexp_linCmt <- function() {
  ini({
    lTVCL      <- log(0.6)
    lTVQ       <- log(1.8)
    lTVVc      <- log(20)
    lTVVp      <- log(80)
    lTVKA      <- fix(log(0.7))
    TH_BW_CL   <- 0.75
    TH_CRCL_CL <- 0.5
    TH_BW_VC   <- 1.0
    TH_SEX_VC  <- log(1.5)

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
    cp     <- linCmt()
    cp ~ prop(prop.err)
  })
}

true_2cmt_scn16_lin_linCmt <- function() {
  ini({
    lTVCL      <- log(0.6)
    lTVQ       <- log(1.8)
    lTVVc      <- log(20)
    lTVVp      <- log(80)
    lTVKA      <- fix(log(0.7))
    TH_BW_CL   <- 0.75
    TH_CRCL_CL <- 0.5
    TH_BW_VC   <- 1.0
    TH_SEX_VC  <- log(1.5)

    eta.cl + eta.vc ~ c(0.1, 0.02, 0.1)
    prop.err <- 0.1
  })
  model({
    lTVCL_typ <- lTVCL + TH_BW_CL * log(BW / 70) + TH_CRCL_CL * log(CrCL / 95)
    lTVVc_typ <- lTVVc + TH_BW_VC * log(BW / 70) + TH_SEX_VC * SEX
    cl        <- exp(lTVCL_typ + eta.cl)
    vc        <- exp(lTVVc_typ + eta.vc)
    q         <- exp(lTVQ)
    vp        <- exp(lTVVp)
    ka        <- exp(lTVKA)
    cp        <- linCmt()
    cp ~ prop(prop.err)
  })
}


# ---- (3) Refit the BASE model on N=80 using the linCmt body ------------------
##   This becomes the parent fit for every fast SCM run below.

t_fit_base16_N80_linCmt <- system.time(
  fit_base16_N80_linCmt <- nlmixr2(base_2cmt_oral_linCmt, ds16_01_N80,
                                    est = "focei", control = scm_focei_screen)
)
saveRDS(fit_base16_N80_linCmt,
        file.path(stage1_dir16_N80, "fit_base_N80_linCmt.rds"))

fit_base16_N80_linCmt <- readRDS(file.path(stage1_dir16_N80, "fit_base_N80_linCmt.rds"))
# ---- (4) Four SCM modes: screen with fast control, refit-once with final ----
##   Each runSCM_traced() call mirrors the corresponding ODE call above; only
##   `fit` and `control` differ.  The final cov + tables come from passing
##   refit_control = scm_focei_final into package_scm_result().

## 4a. Forward selection (auto-generated Cartesian product)
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
t16_N80_fwd_auto_fast <- attr(res16_N80_fwd_auto_fast, "elapsed_s") #19.2min
saveRDS(res16_N80_fwd_auto_fast,
        file.path(stage1_dir16_N80, "res_fwd_auto_fast.rds"))
test16_N80_fwd_auto_fast <- package_scm_result(
  "scn16_N80_forward_auto_fast",
  res16_N80_fwd_auto_fast, t16_N80_fwd_auto_fast,
  scenario_id   = 16,
  refit_control = scm_focei_final
)
saveRDS(test16_N80_fwd_auto_fast,
        file.path(stage1_dir16_N80, "test_fwd_auto_fast.rds"))

## 4b. Backward elimination only (curated pool incl. BMI / RACE distractors)
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
t16_N80_bck_fast <- attr(res16_N80_bck_fast, "elapsed_s")
saveRDS(res16_N80_bck_fast,
        file.path(stage1_dir16_N80, "res_bck_fast.rds"))
test16_N80_bck_fast <- package_scm_result(
  "scn16_N80_backward_only_fast",
  res16_N80_bck_fast, t16_N80_bck_fast,
  scenario_id   = 16,
  refit_control = scm_focei_final
)
saveRDS(test16_N80_bck_fast,
        file.path(stage1_dir16_N80, "test_bck_fast.rds"))

## 4c. User-specified true relations (4 true pairs)
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
t16_N80_user_fast <- attr(res16_N80_user_fast, "elapsed_s")
saveRDS(res16_N80_user_fast,
        file.path(stage1_dir16_N80, "res_user_fast.rds"))
test16_N80_user_fast <- package_scm_result(
  "scn16_N80_user_4tr_fast",
  res16_N80_user_fast, t16_N80_user_fast,
  scenario_id   = 16,
  refit_control = scm_focei_final
)
saveRDS(test16_N80_user_fast,
        file.path(stage1_dir16_N80, "test_user_fast.rds"))

## 4d. Full SCM (forward then backward, full Cartesian pool)
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
t16_N80_full_fast <- attr(res16_N80_full_fast, "elapsed_s")
saveRDS(res16_N80_full_fast,
        file.path(stage1_dir16_N80, "res_full_fast.rds"))
test16_N80_full_fast <- package_scm_result(
  "scn16_N80_full_scm_fast",
  res16_N80_full_fast, t16_N80_full_fast,
  scenario_id   = 16,
  refit_control = scm_focei_final
)
saveRDS(test16_N80_full_fast,
        file.path(stage1_dir16_N80, "test_full_fast.rds"))


# ---- (5) Fast-path summary (parallel to part2_summary16_N80) ----------------
scm_tests16_N80_fast <- list(
  forward_auto  = test16_N80_fwd_auto_fast,
  backward_only = test16_N80_bck_fast,
  user_4true    = test16_N80_user_fast,
  full_scm      = test16_N80_full_fast
)

part2_summary16_N80_fast <- purrr::imap_dfr(scm_tests16_N80_fast, function(x, mode) {
  tibble::tibble(
    mode          = mode,
    test          = x$label,
    n_selected    = if (is.null(x$selected)) NA_integer_ else nrow(x$selected),
    converged     = if (is.null(x$diag)) NA else x$diag$converged,
    objf          = if (is.null(x$diag)) NA_real_ else x$diag$objf,
    cov_step_ok   = if (is.null(x$diag)) NA else x$diag$cov_ok,
    cond_num      = if (is.null(x$diag)) NA_real_ else x$diag$cond_num,
    cond_num_sqrt = if (is.null(x$diag)) NA_real_ else x$diag$cond_num_sqrt,
    runtime_min   = x$runtime_sec / 60
  )
})
saveRDS(part2_summary16_N80_fast,
        file.path(stage1_dir16_N80, "part2_summary_fast.rds"))


# ---- (6) Side-by-side speedup comparison ------------------------------------
##   Joins original (ODE + cov-on-every-step) against fast (linCmt + 2-tier).
##   part2_summary16_N80 has 5 rows; we mirror only the 4 modes that the
##   fast block runs.  `speedup` > 1 means the fast path is faster.

.orig_mode_map <- tibble::tribble(
  ~mode,            ~test,
  "forward_auto",   "scn16_N80_forward_auto",
  "backward_only",  "scn16_N80_backward_only",
  "user_4true",     "scn16_N80_user_4tr",
  "full_scm",       "scn16_N80_full_scm"
)

part2_speedup16_N80 <- dplyr::full_join(
  part2_summary16_N80 %>%
    dplyr::inner_join(.orig_mode_map, by = "test") %>%
    dplyr::select(mode,
                  n_sel_orig    = n_selected,
                  objf_orig     = objf,
                  cond_num_orig = cond_num,
                  runtime_orig  = runtime_min),
  part2_summary16_N80_fast %>%
    dplyr::select(mode,
                  n_sel_fast    = n_selected,
                  objf_fast     = objf,
                  cond_num_fast = cond_num,
                  runtime_fast  = runtime_min),
  by = "mode"
) %>%
  dplyr::mutate(
    speedup        = runtime_orig / runtime_fast,
    dObjf          = objf_fast - objf_orig,
    same_n_sel     = n_sel_orig == n_sel_fast
  ) %>%
  dplyr::arrange(dplyr::desc(speedup))

saveRDS(part2_speedup16_N80,
        file.path(stage1_dir16_N80, "part2_speedup.rds"))




