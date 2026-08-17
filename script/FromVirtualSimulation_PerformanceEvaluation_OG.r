##Install remote & install nlmixr2utils and nlmixr2scm from GitHub
remotes::install_github("kestrel99/nlmixr2utils")
remotes::install_github("kestrel99/nlmixr2scm")

options(download.file.method = "libcurl", url.method = "libcurl")
install.packages("nlmixr2",
  repos = c("https://nlmixr2.r-universe.dev",
            "https://cloud.r-project.org"))

##1. Run test suites first for nlmixr2scm
## NOTE: ~/.Rprofile sets options(rxode2.cache.dir = ...);
##       ~/.Renviron sets TMPDIR (must be in .Renviron, not .Rprofile,
##       because R locks tempdir() at start-up).
##       Open them manually only if you need to inspect / edit:
# usethis::edit_r_profile(scope = "user")
# file.edit("C:/Users/LIUYA8J/.Renviron")

library(testthat)
library(remotes)
library(nlmixr2utils)
library(nlmixr2scm)
#devtools::load_all("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
#testthat::test_file("tests/testthat/test-scm.R") #Pass228, warning from nlmixrest[near-singular cov on toy fixture]
#testthat::test_file("tests/testthat/test-parsing.R") #pass29, after fixing the function calling issue

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
TVCL <- 0.6;  TVQ <- 1.8;  TVVc <- 20;  TVVp <- 80 ; Ka <- 0.7
k10_typ <- TVCL / TVVc
k12_typ <- TVQ  / TVVc
k21_typ <- TVQ  / TVVp
sum_k    <- k10_typ + k12_typ + k21_typ
beta_typ <- 0.5 * (sum_k - sqrt(sum_k^2 - 4 * k10_typ * k21_typ))
t_half_typ <- log(2) / beta_typ
t_half_ab <- log(2) / Ka

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

out_dir <- "simulated_virtual_dataset"
true_params <- readRDS(file.path(out_dir, "true_params_long.rds"))

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

## ============================================================================
## Convergence diagnostics
##   converged    : TRUE iff optimizer reported success (fit$convergence == 0)
##                  AND objf is finite.  fit$convergence is the canonical
##                  optimizer status flag (0 = success); the earlier
##                  is.finite(objf)-only test silently treated non-converged
##                  fits as converged whenever OFV happened to be finite.
##   cond_num_cor : raw condition number of the *correlation* matrix of fixed
##                  effects (lambda_max / lambda_min of cor).  This is the
##                  pharma-standard "Condition#(Cor)" reported by nlmixr2
##                  since 2.1.4.  Populated only when covMethod != "" was
##                  used at fit time -- or after a post-hoc getVarCov() call
##                  that triggers .setCov().  Common threshold: < 1000.
##   cov_ok       : TRUE iff fit$cov is a finite-diagonal matrix.  Will be
##                  FALSE when covMethod = "" was used (cov not computed).
##   message      : optimizer exit message, e.g. "Normal exit from bobyqa".#
## ============================================================================
diagnose_fit <- function(fit) {
  if (is.null(fit)) {
    return(list(converged = NA, objf = NA_real_,
                cond_num_cor = NA_real_,
                cov_ok = NA, message = NA_character_))
  }
  conv_code <- if (!is.null(fit$convergence)) fit$convergence else NA_integer_
  cn_cor    <- fit$conditionNumberCor
  cn_cor_val <- if (is.null(cn_cor)) NA_real_ else as.numeric(cn_cor)
  list(
    converged    = isTRUE(conv_code == 0L) &&
                     !is.null(fit$objf) && is.finite(fit$objf),
    objf         = if (!is.null(fit$objf)) fit$objf else NA_real_,
    cond_num_cor = cn_cor_val,
    cov_ok       = isTRUE(!is.null(fit$cov) && all(is.finite(diag(fit$cov)))),
    message      = if (!is.null(fit$message)) as.character(fit$message) else NA_character_
  )
}


## ============================================================================
## Map fit$theta + fit$omega -> article parameter names
## Extract fixed-effects, random-effects, residual error, and covariate effects into a long tibble
## ============================================================================
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
  ## The regexes below mirror these two naming conventions exactly.
  ##   Naming conventions (verified against runSCM output):
  ##     continuous : cov_<COVAR>_<SHAPE>_<VAR>     e.g. cov_BW_power_cl
  ##     categorical: cov_<COVAR>_<LEVEL>_<VAR>     e.g. cov_SEX_1_vc
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

## ============================================================================
## Calculate Relative error per parameter (one dataset)
## ============================================================================
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

## ============================================================================
## ---- Helper: package one runSCM result for downstream comparison --------
##   Returns:
##     - selected:        accepted (var, covar, shape) pairs
##     - step_hist:       full SCM step history (forward + backward summary)
##     - final_est:       extract_params_long() of the final model
##     - rel_err:         relative error vs true scenario parameters
##     - runtime_sec:     wall-clock seconds for the runSCM call
##     - diag:            convergence diagnostics for the final fit
##     - parFixed:        nlmixr2 parFixedDf (Estimate, SE, %RSE, CI) when
##                         the refit succeeded; else NULL
##     - cov_done:        TRUE iff the final_ctrl refit succeeded and the
##                         resulting fit has a populated $cov + parFixedDf.
##                         FALSE when (a) final_ctrl = NULL, or (b) the refit
##                         hit an error.  When FALSE: cond_num_cor stays NA
##                         and PowerCN becomes unavailable for that dataset.
##
##   Covariance computation (the refit path -- restored after the failed
##   getVarCov() experiment).  Why a full refit instead of just cov:
##     * getVarCov() / .setCov() freezes theta + ETAs at the SCM-final point
##       (maxOuterIterations=0L, maxInnerIterations=0L) and only recomputes
##       the Hessian.  But the screening control stops on sigdig=3 with
##       atol=1e-6 / rtol=1e-4, so theta is NOT at a tight-tol stationary
##       point.  The numerical Hessian at a non-stationary theta is ill-
##       conditioned -> R-matrix inversion fails -> fallback to S-only cov
##       -> SEs are 5-50x too small (S inflated by the systematic non-zero
##       gradient).  No cov-time tolerance tightening can fix this -- only
##       moving theta to a tight-tol stationary point does.
##     * The refit (scm_focei_final: sigdig=4, atol=1e-8, rtol=1e-6,
##       covMethod="r,s") re-converges theta + ETAs at tight tols and lands
##       on a stationary point where R is well-conditioned.  Theta shifts
##       are tiny (<3% on covariate effects in pilot ds01) but the cov
##       quality is qualitatively different.
##     * Cost: ~30-60s/ds, negligible against the ~15 min SCM screening.
##
##   final_ctrl:
##     * NULL (default): skip the refit entirely.  selected / step_hist /
##       final_est / rel_err / diag still populate from the screening-tol
##       fit; parFixed = NULL and cov_done = FALSE.
##     * scm_focei_final (or similar): trigger the refit.  Recommended
##       whenever uncertainty / cond_num_cor / PowerCN are needed.
## ============================================================================
package_scm_result <- function(label, scm_res, runtime_sec,
                               true_long = true_params, scenario_id = 9,
                               final_ctrl = NULL) {
  ## Pick the final fit.
  .pickFit <- function(x) {
    if (is.null(x)) return(NULL)
    cand <- if (is.list(x) && length(x) >= 1L) x[[1L]] else x
    if (inherits(cand, "nlmixr2FitCore")) cand else NULL
  }
  final_fit <- .pickFit(scm_res$resBck)
  if (is.null(final_fit)) final_fit <- .pickFit(scm_res$resFwd)

  ## Optional refit with the diagnostic control.  On success the refit
  ## REPLACES final_fit so every downstream extraction (estimates, rel_err,
  ## diag, parFixed) reflects the tight-tol stationary point.
  cov_done <- FALSE
  if (!is.null(final_fit) && !is.null(final_ctrl)) {
    refit <- tryCatch(
      nlmixr2(final_fit$ui, nlme::getData(final_fit),
              est = final_fit$est, control = final_ctrl),
      error = function(e) {
        warning("package_scm_result(): refit failed: ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(refit)) {
      final_fit <- refit
      cov_done  <- !is.null(refit$cov) && all(is.finite(diag(refit$cov)))
    }
  }

  ## --- Final-model covariate relations -
  ##The only reliable source-of-truth for "what's in the final model" is final_fit's
  ## own theta vector, since runSCM names covariate coefficients
  ##     cov_<COVAR>_<SHAPE-or-LEVEL>_<VAR> (e.g. cov_BW_power_cl, cov_SEX_1_vc).
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
  ## parFixedDf is only populated when the refit succeeded.
  parFixed  <- if (cov_done) final_fit$parFixedDf else NULL

  ## --- Packaged result -----------------------------------------------------
  ## `step_hist` is the ONLY view of the search trace we keep.
  list(
    label       = label,
    selected    = selected,
    step_hist   = scm_res$summaryTable,
    final_fit   = final_fit,
    final_est   = final_est,
    rel_err     = rel_err,
    diag        = diag,
    parFixed    = parFixed,
    cov_done    = cov_done,
    runtime_sec = runtime_sec
  )
}



## ============================================================================
## ---- Wrapper: runSCM with elapsed-time stash ---------------------------
##   Minimal pass-through to nlmixr2scm::runSCM(): times the call and stashes
##   the elapsed seconds on the returned object as attr(res, "elapsed_s").
##   `label` is accepted for backward compatibility with existing call sites
##   (Part 2 one-shots, the per-ds driver) but is otherwise unused.  
## ============================================================================
runSCM_traced <- function(label, ...) {
  t0  <- Sys.time()
  res <- nlmixr2scm::runSCM(...)
  attr(res, "elapsed_s") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  res
}




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

## ---- DATA Preparation ---------------------------------------------------------------
## NM-format dataset for SCENARIO = 16, DATASET = 1 (n=300)----
out_dir_v2   <- "simulated_virtual_dataset_eta_filtered"
stage1_dir16 <- file.path(out_dir_v2, "stage1_smoke_scn16_ds01")
if (!dir.exists(stage1_dir16)) dir.create(stage1_dir16, recursive = TRUE)

sim_obs_scn16 <- readRDS(file.path(out_dir_v2, "sim_obs_scenario_16.rds"))
ds16_01 <- to_nm_dataset(sim_obs_scn16) %>%
  dplyr::filter(DATASET == 1) %>% #getDATASET1
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
    cond_num_cor  = if (is.null(x$diag)) NA_real_ else x$diag$cond_num_cor,
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
##        needs OFV).  The surviving model gets a post-hoc covariance via
##        nlme::getVarCov() called inside package_scm_result() -- this
##        invokes nlmixr2est's .setCov() method (no re-estimation,
##        ~5-20s/ds), populating $cov, $parFixedDf, $conditionNumberCor.
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
##   _final:  tight tols + covMethod="r,s" + calcTables=TRUE.  Fed to
##            package_scm_result(final_ctrl = scm_focei_final) for the
##            single post-SCM refit that gives trustworthy SE / cond_num_cor
##            (see package_scm_result docstring for why a refit is needed).

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
# ---- (4) Four SCM modes: screen with fast control, cov via getVarCov() ----
##   Each runSCM_traced() call mirrors the corresponding ODE call above; only
##   `fit` and `control` differ.  Post-hoc cov + parFixedDf are attached by
##   package_scm_result() via nlme::getVarCov() on the SCM-final fit.

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
  scenario_id   = 16
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
  scenario_id   = 16
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
  scenario_id   = 16
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
  scenario_id   = 16
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
    cond_num_cor  = if (is.null(x$diag)) NA_real_ else x$diag$cond_num_cor,
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
                  n_sel_orig        = n_selected,
                  objf_orig         = objf,
                  cond_num_cor_orig = cond_num_cor,
                  runtime_orig      = runtime_min),
  part2_summary16_N80_fast %>%
    dplyr::select(mode,
                  n_sel_fast        = n_selected,
                  objf_fast         = objf,
                  cond_num_cor_fast = cond_num_cor,
                  runtime_fast      = runtime_min),
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


## ============================================================================
# Part 3: Operating-characteristics pilot (4 datasets, scenario 16, N=80)-----------

##   Wraps the fast-route full-SCM pipeline used above (single dataset) into
##   a per-dataset driver, loops 4 replicates of scenario 16, then computes
##   six OC artefacts:
##     - oc_pilot_power       : Power, PowerCN, PowerMinSuc
##     - oc_pilot_relpower    : relative power at k = 1..Ntrue
##     - oc_pilot_rmrse_uncond: unconditional RMRSE %
##     - oc_pilot_rmrse_cond  : conditional RMRSE % (exact-match runs only)
##     - oc_pilot_timing      : per-dataset + aggregate wall-clock
##     - oc_pilot_per_ds      : audit log (one row per dataset)
##
##   Helpers below are intentionally scenario-agnostic so the same code drops
##   onto the 100 x 16 scale-up.  See the SCALE-UP STUB at the end.
## ============================================================================


# ---- (P3.1) OC helpers (reusable across scenarios) --------------------------

##  Canonical-shape normalizer.  runSCM emits per-LEVEL theta names for
##  categorical covariates (e.g. cov_SEX_1_vc, cov_RACE_2_vc); the parser
##  in package_scm_result() leaves the level digit in `shape`.  For
##  semantic matching against truth (where the relation is one entity)
##  we collapse any all-digit shape to "cat".  This also folds a 3-level
##  cat that selected multiple levels into ONE relation row.
.canon_shape_tbl <- function(rel) {
  rel %>%
    dplyr::mutate(shape = ifelse(grepl("^[0-9]+$", shape), "cat", shape)) %>%
    dplyr::distinct(var, covar, shape)
}

##  Build the truth tibble of "what runSCM should have selected" for a given
##  scenario.  Reads true_long, keeps rows with non-zero true_value AND a
##  parameter name that maps to a (var, covar, shape) triple.  Defaults
##  cover the article scenarios (scn 1-16); pass `shape_map` if you add
##  custom covariate effects later.
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

##  Compare a single run's `selected` tibble against truth.  Returns a
##  list with: n_true_hit, n_false_pos, exact_match (set equality on
##  (var, covar, shape) after canonicalization).
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

##  Compute Power, PowerCN (cond_num_cor < cn_cor_cut), PowerMinSuc
##  (converged == TRUE), plus a relative-power tibble (k = 1..n_true).
##  cond_num_cor is the correlation-matrix condition number (lambda_max /
##  lambda_min) from fit$conditionNumberCor.  Convention: < 1000 is "well
##  conditioned" (FDA / pharma pop-PK threshold).  We do NOT take a sqrt
##  any more -- the old cond_num_sqrt was sqrt(cond_num_cov), and Cor and
##  Cov condition numbers are on different scales (Cor is more sensitive
##  to high parameter correlations, Cov is dominated by scale differences).
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

##  Compute unconditional + conditional RMRSE / MARE for a parameter set.
##    pop_params: vector of population-parameter names (always estimated).
##                Unconditional denom = runs with finite estimate.
##    cov_params: vector of true covariate-effect parameter names.
##                Unconditional denom = runs where that cov was SELECTED
##                (and therefore estimated, => non-NA estimate).
##    Conditional denom for BOTH groups = runs with exact_match == TRUE.
##
##  Metrics returned per parameter:
##    RMRSE_pct: 100 * sqrt(mean(rel_err^2))   -- mean accuracy + variance
##                                                penalty, sensitive to outliers
##    MARE_pct:  100 * median(|rel_err|)       -- robust central tendency,
##                                                outlier-resistant
##    n_used:    # of finite contributions
##
##  Rationale: with N = 4 (pilot) or even N = 250 (scale-up) draws, a single
##  numerically-unstable estimate (e.g. cov_VcCL near the noise floor) can
##  dominate RMRSE.  MARE gives a complementary view that's stable to that
##  single bad draw.  Report both side-by-side.
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


# ---- (P3.2) Per-dataset driver ----------------------------------------------
##  Filters the long sim file to one DATASET, fits base with linCmt, runs full
##  SCM via the fast screening control, then attaches a post-hoc covariance to
##  the SCM-final model via nlme::getVarCov() (inside package_scm_result --
##  no re-estimation, ~5-20s/ds vs ~60s/ds for the old refit path).  Wrapped
##  in tryCatch so a single bad dataset returns a stub list (no $test) instead
##  of halting the loop.  `confirm = FALSE` is forced -- no interactive prompts.
##
##  Resilience design:
##    * `true_long` is threaded explicitly so an unloaded `true_params` in
##      the calling env can't cause silent NULL returns from package_scm_result.
##    * SCM result `res_i` is saved to disk IMMEDIATELY after runSCM_traced,
##      before any packaging.  Even if package_scm_result errors later,
##      the (expensive) SCM artefact stays recoverable.
##    * Error handler echoes the exception via warning() AND writes a
##      per-dataset error log to <save_dir>/<ds_tag>_ERROR.txt.
##
##  Resume support (two-tier cache):
##    * Tier 1 -- packaged result.  If `test_full_fast_<ds_tag>.rds` exists
##      and neither `force_rerun` nor `force_repackage` is TRUE, the driver
##      returns the cached `test_*` object untouched (~ms cost).
##    * Tier 2 -- SCM screening result.  If `res_full_fast_<ds_tag>.rds`
##      exists (and tier 1 was skipped or bypassed via `force_repackage`),
##      the driver skips base-fit + SCM screening and runs only
##      package_scm_result() (~30-60s for the tight-tol refit).  This is
##      the supported recovery path when package_scm_result() crashed on
##      a previous run, or when `scm_focei_final` settings have changed
##      and only the cov/SE step needs to be redone.
##    * Tier 3 -- full pipeline.  Both caches missing, or `force_rerun =
##      TRUE`.  Runs base fit + SCM screening + packaging end-to-end.
##
##    Flags:
##      - `force_rerun = TRUE`     bypasses BOTH caches -> tier 3.
##      - `force_repackage = TRUE` ignores stale `test_*` but still reuses
##        cached `res_*` -> tier 2 (saves the ~15 min SCM screening).
##        When `res_*` is also missing, falls through to tier 3.
##
##    On tier-3 success the driver unlinks any stale
##    `<ds_tag>_ERROR.txt` from a prior failed run so the log doesn't
##    accumulate misleading sentinels.
##
##  Covariance:
##    * The driver passes `final_ctrl = scm_focei_final` to
##      package_scm_result(), which re-fits the SCM-final model ONCE with
##      tight tols + covMethod="r,s".  Cost: ~30-60s/ds (negligible vs the
##      ~15 min SCM screening).  Yields trustworthy parFixed / cond_num_cor.
##    * Why a full refit rather than getVarCov() on the screening-tol fit:
##      the screening optimum is not a tight-tol stationary point, so the
##      numerical Hessian there is ill-conditioned and the cov calculation
##      falls back to a deflated S-only matrix.  See package_scm_result()
##      docstring for the full explanation.
##
##  saveModels = FALSE inside runSCM:
##    * The package's per-accepted-step audit trail (scm_<ds_tag>/...rds)
##      is disabled to keep the working-set disk small and avoid OneDrive
##      sync overhead on each step.  The only persisted artefact per ds is
##      the packaged `test_full_fast_<ds_tag>.rds`.
##
##  Parallelism:
##    * `workers` controls INNER parallelism: candidates within one SCM step.
##      For the pilot (sequential outer loop) workers = 3L is the default.
##      For the scale-up under outer parallelism (future_pmap), set workers
##      = 1L to avoid nested-future thrash.

run_one_dataset_scn16_N80 <- function(ds_id, save_dir,
                                       sim_long     = sim_obs_scn16_N80,
                                       base_fn      = base_2cmt_oral_linCmt,
                                       screen_ctrl  = scm_focei_screen,
                                       vars_vec     = scm16_vars,
                                       covars_vec   = scm16_covars,
                                       catvars_vec  = scm16_catvars,
                                       shapes_vec   = scm16_shapes,
                                       true_long    = true_params,
                                       workers          = 3L,
                                       keep_res         = TRUE,
                                       force_rerun      = FALSE,
                                       force_repackage  = FALSE,
                                       scenario_id      = 16L) {
  ds_tag    <- sprintf("ds%02d", ds_id)
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  test_path <- file.path(save_dir, sprintf("test_full_fast_%s.rds", ds_tag))
  res_path  <- file.path(save_dir, sprintf("res_full_fast_%s.rds",  ds_tag))
  err_path  <- file.path(save_dir, sprintf("%s_ERROR.txt",          ds_tag))

  ## -- Tier 1: packaged result already on disk -> return immediately.
  ##    force_rerun trumps everything; force_repackage skips tier 1 so the
  ##    driver falls through to tier 2 (repackage from cached res_*).
  if (!force_rerun && !force_repackage && file.exists(test_path)) {
    message(sprintf(">>> [%s] tier-1 cached, skipping (%s)",
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

  ## -- Tier 2: SCM screening cached but packaging missing or invalidated.
  ##    Re-run only package_scm_result() on the cached res_*.  This is the
  ##    supported recovery path for a crashed packaging step OR a deliberate
  ##    re-package after scm_focei_final has been tightened/changed.
  ##    Wrapped in its own tryCatch so a refit crash here doesn't lose the
  ##    expensive res_* artefact -- the ERROR.txt is written and the user
  ##    can retry after fixing helpers or final_ctrl.
  if (!force_rerun && file.exists(res_path)) {
    message(sprintf(">>> [%s] tier-2 re-packaging from cached %s",
                    ds_tag, basename(res_path)))
    return(tryCatch({
      res_i     <- readRDS(res_path)
      t_scm_sec <- suppressWarnings(as.numeric(attr(res_i, "elapsed_s")))
      if (!is.finite(t_scm_sec)) {
        cli::cli_warn("[{ds_tag}] cached res_* missing `elapsed_s` attr; t_scm_sec set to NA.")
        t_scm_sec <- NA_real_
      }
      t_pkg <- system.time(
        test_i <- package_scm_result(
          label       = sprintf("scn%02d_%s_N80_full_fast", scenario_id, ds_tag),
          scm_res     = res_i,
          runtime_sec = t_scm_sec,
          true_long   = true_long,
          scenario_id = scenario_id,
          final_ctrl  = scm_focei_final
        )
      )
      saveRDS(test_i, test_path)
      if (file.exists(err_path)) unlink(err_path)  # stale sentinel from prior fail
      message(sprintf("<<< [%s] tier-2 done; refit %.1fs, cov=%s",
                      ds_tag, as.numeric(t_pkg["elapsed"]),
                      if (isTRUE(test_i$cov_done)) "ok" else "FAIL"))
      list(
        ds_id       = ds_id,
        ds_tag      = ds_tag,
        test        = test_i,
        t_base_sec  = NA_real_,
        t_scm_sec   = t_scm_sec,
        t_total_sec = t_scm_sec,
        resumed     = TRUE
      )
    }, error = function(e) {
      msg <- conditionMessage(e)
      writeLines(c(format(Sys.time()), "tier-2 repackage failed:", msg), err_path)
      warning(sprintf("[%s] tier-2 FAILED: %s (see %s)", ds_tag, msg, err_path),
              call. = FALSE, immediate. = TRUE)
      list(ds_id = ds_id, ds_tag = ds_tag, error = msg, res = NULL,
           resumed = TRUE)
    }))
  }

  ## -- Tier 3: full pipeline (base fit + SCM screening + package).
  message(sprintf("\n>>> [%s] tier-3 starting full pipeline at %s",
                  ds_tag, format(Sys.time(), "%H:%M:%S")))

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
      saveModels  = FALSE,       # skip per-step audit trail (OneDrive sync tax)
      workers     = workers,
      print       = 0,           # silence per-iteration FOCEi noise
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

    ## package_scm_result() refits the SCM-final model once with scm_focei_final
    ## (tight tols + covMethod="r,s") to obtain trustworthy SE / cond_num_cor.
    ## Returns cov_done = FALSE if the refit hits an error; downstream OC code
    ## handles that by filtering on cond_num_cor < cn_cor_cut.
    test_i <- package_scm_result(
      label         = sprintf("scn%02d_%s_N80_full_fast", scenario_id, ds_tag),
      scm_res       = res_i,
      runtime_sec   = t_scm_sec,
      true_long     = true_long,         # explicit -- no lazy global lookup
      scenario_id   = scenario_id,
      final_ctrl    = scm_focei_final    # tight-tol refit for cov + SE
    )

    saveRDS(test_i, test_path)
    if (file.exists(err_path)) unlink(err_path)  # clear stale sentinel on recovery

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


# ---- (P3.3) Run the 4-dataset pilot -----------------------------------------
true_params <- readRDS(file.path(out_dir, "true_params_long.rds"))
out_dir <- "simulated_virtual_dataset"
out_dir_v2_N80   <- "simulated_virtual_dataset_eta_filtered_N80"

stage1_dir16_N80 <- file.path(out_dir_v2_N80 , "stage1_smoke_scn16_ds01_N80")
sim_obs_scn16_N80 <- readRDS(file.path(out_dir_v2_N80, "sim_obs_scenario_16.rds")) #250 files N=80

stage1_pilot_scn16_N80 <- file.path(out_dir_v2_N80,
                                     "stage1_pilot_scn16_N80")

## Pilot output dir for the corrected pipeline (refit restored, saveModels=FALSE).
## The _m01 copy contains results from the failed getVarCov() experiment and
## is kept for reference / diff comparisons.
stage1_pilot_scn16_N80 <- file.path(out_dir_v2_N80,
                                     "stage1_pilot_scn16_N80_m02")

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

DOSE_MG <- 100
scm16_vars       <- c("cl", "vc")
scm16_covars     <- c("BW", "CrCL", "BMI")
scm16_catvars    <- c("SEX", "RACE")
scm16_shapes     <- c("power", "lin")


t_pilot_total <- system.time(
  scm_pilot_scn16_N80 <- purrr::map(
    1:4,
    run_one_dataset_scn16_N80,
    save_dir = stage1_pilot_scn16_N80
  ) %>%
    rlang::set_names(sprintf("ds%02d", 1:4))
)
res_full_fast_ds01_T02 <- readRDS("simulated_virtual_dataset_eta_filtered_N80/stage1_pilot_scn16_N80_m02/res_full_fast_ds01.rds")
## Loop over cached SCM screening results.
pilot_dir <- "simulated_virtual_dataset_eta_filtered_N80/stage1_pilot_scn16_N80_m02"
res_files <- list.files(pilot_dir,
                        pattern = "^res_full_fast_ds\\d+\\.rds$",
                        full.names = TRUE)

scm_pilot_scn16_N80 <- purrr::map(res_files, function(rp) {
  ds_tag    <- sub("^res_full_fast_(ds\\d+)\\.rds$", "\\1", basename(rp))
  ds_id     <- as.integer(sub("^ds", "", ds_tag))
  test_path <- file.path(pilot_dir, sprintf("test_full_fast_%s.rds", ds_tag))

  message(sprintf(">>> [%s] loading SCM cache + refitting", ds_tag))
  res_i  <- readRDS(rp)
  t_scm  <- as.numeric(attr(res_i, "elapsed_s"))

  t_pkg <- system.time(
    test_i <- package_scm_result(
      label       = sprintf("scn16_%s_N80_full_fast", ds_tag),
      scm_res     = res_i,
      runtime_sec = t_scm,
      true_long   = true_params,
      scenario_id = 16,
      final_ctrl  = scm_focei_final     # tight-tol refit for cov + SE
    )
  )

  saveRDS(test_i, test_path)
  message(sprintf("<<< [%s] refit %.1fs, cov_done=%s, cn_cor=%.1f",
                  ds_tag, t_pkg["elapsed"],
                  isTRUE(test_i$cov_done),
                  if (is.null(test_i$diag$cond_num_cor)) NA_real_
                  else test_i$diag$cond_num_cor))

  ## Mimic the driver's return shape so downstream OC code works unchanged.
  list(
    ds_id       = ds_id,
    ds_tag      = ds_tag,
    test        = test_i,
    t_base_sec  = NA_real_,                       # base fit not re-timed here
    t_scm_sec   = t_scm,
    t_total_sec = t_scm,
    resumed     = FALSE
  )
}) %>% rlang::set_names(sprintf("ds%02d", seq_along(res_files)))



## Loop summary: fresh runs vs cache hits vs failures.
.summarize_pilot <- function(pilot_list) {
  n_total   <- length(pilot_list)
  n_resumed <- sum(vapply(pilot_list, function(x) isTRUE(x$resumed), logical(1)))
  n_failed  <- sum(vapply(pilot_list,
                          function(x) !is.null(x$error) && is.null(x$test),
                          logical(1)))
  n_fresh   <- n_total - n_resumed - n_failed
  message(sprintf(
    "\n=== pilot loop done: %d fresh + %d resumed + %d failed / %d total ===",
    n_fresh, n_resumed, n_failed, n_total
  ))
  invisible(list(fresh = n_fresh, resumed = n_resumed,
                 failed = n_failed, total = n_total))
}
.summarize_pilot(scm_pilot_scn16_N80)


# ---- (P3.3.5) One-shot migration of cached pilot files ---------------------
##   Earlier pilot runs produced `test_full_fast_<ds_tag>.rds` whose `$diag`
##   slot has the OLD schema (cond_num, cond_num_sqrt) and an obsolete
##   `$refit_done` flag.  After the cond_num_cor + getVarCov refactor those
##   names no longer match what compute_power_block() / oc_pilot_per_ds /
##   compute_rmrse_block expect.  Rather than redo the ~17-min/ds SCM, this
##   helper re-runs diagnose_fit() on each cached fit (filling cond_num_cor)
##   and replaces refit_done with cov_done.  In-place; idempotent (safe to
##   call on already-migrated files).
##
##   Note: cached `$final_fit` came from the old refit path which used
##   covMethod="r,s", so $cov / $conditionNumberCor are already populated.
##   We only need to recompute the per-test $diag and reset the new flag.
.rediag_pilot_cache <- function(save_dir,
                                pattern = "^test_full_fast_ds\\d+\\.rds$") {
  if (!dir.exists(save_dir)) {
    message(sprintf("  no pilot cache dir: %s", save_dir))
    return(invisible(0L))
  }
  files <- list.files(save_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    message(sprintf("  no cached pilot files matching %s in %s",
                    pattern, save_dir))
    return(invisible(0L))
  }
  n_changed <- 0L
  for (f in files) {
    t <- readRDS(f)
    needs_diag  <- is.null(t$diag) || is.null(t$diag$cond_num_cor)
    needs_flag  <- is.null(t$cov_done)
    needs_pfix  <- !is.null(t$final_fit) &&
                    is.null(t$parFixed)  &&
                    !is.null(t$final_fit$parFixedDf)
    if (!(needs_diag || needs_flag || needs_pfix)) next  # already current

    if (!is.null(t$final_fit)) {
      ## Cov already attached from old refit -- diagnose_fit reads it.
      t$diag     <- diagnose_fit(t$final_fit)
      t$cov_done <- !is.null(t$final_fit$cov) &&
                     all(is.finite(diag(t$final_fit$cov)))
      t$parFixed <- if (isTRUE(t$cov_done)) t$final_fit$parFixedDf else NULL
    } else {
      t$diag     <- diagnose_fit(NULL)
      t$cov_done <- FALSE
      t$parFixed <- NULL
    }
    ## Drop the obsolete field name so the schema is clean.
    t$refit_done <- NULL
    saveRDS(t, f)
    n_changed <- n_changed + 1L
  }
  message(sprintf("  re-diagnosed %d / %d cached pilot files in %s",
                  n_changed, length(files), basename(save_dir)))
  invisible(n_changed)
}
.rediag_pilot_cache(stage1_pilot_scn16_N80)

## After migration the in-memory `scm_pilot_scn16_N80` may still hold the
## OLD-schema list from the resume path (run_one_dataset_scn16_N80() does
## `readRDS(test_path)` at the top of every cached run).  Re-load the
## migrated copies so the downstream oc_pilot_per_ds / compute_power_block
## chain sees cond_num_cor + cov_done immediately.
for (ds_id in seq_along(scm_pilot_scn16_N80)) {
  ds_tag <- sprintf("ds%02d", ds_id)
  tp <- file.path(stage1_pilot_scn16_N80,
                  sprintf("test_full_fast_%s.rds", ds_tag))
  if (file.exists(tp)) scm_pilot_scn16_N80[[ds_id]]$test <- readRDS(tp)
}


# ---- (P3.4) Per-dataset audit tibble ----------------------------------------
##   `p$test` is NULL when the driver hit its error handler (it still returns
##   a list carrying $error and any salvaged $res); skip those rows.
true_params <- readRDS(file.path(out_dir, "true_params_long.rds"))
true_rel_scn16 <- extract_true_relations(true_params, scenario_id = 16)

oc_pilot_per_ds <- purrr::map_dfr(scm_pilot_scn16_N80, function(p) {
  if (is.null(p) || is.null(p$test)) return(NULL)
  tt  <- p$test
  mat <- match_selected_to_truth(tt$selected, true_rel_scn16)
  tibble::tibble(
    dataset_id    = p$ds_id,
    ds_tag        = p$ds_tag,
    t_base_min    = p$t_base_sec / 60,
    t_scm_min     = p$t_scm_sec  / 60,
    t_total_min   = p$t_total_sec / 60,
    converged     = isTRUE(tt$diag$converged),
    cov_ok        = isTRUE(tt$diag$cov_ok),
    cond_num_cor  = tt$diag$cond_num_cor,
    n_selected    = if (is.null(tt$selected)) NA_integer_ else nrow(tt$selected),
    n_true_hit    = mat$n_true_hit,
    n_false_pos   = mat$n_false_pos,
    exact_match   = mat$exact_match,
    selected      = list(tt$selected),
    final_est     = list(tt$final_est)
  )
})


# ---- (P3.5) Six OC artefacts ------------------------------------------------
n_true_scn16 <- nrow(true_rel_scn16)   # = 4 for scenario 16


## Pop params for RMRSE (TVKA dropped: fixed in model)
pop_params_rmrse <- c("TVCL", "TVVc", "TVQ", "TVVp",
                      "var_CL", "var_Vc", "cov_VcCL", "ResErr")
cov_params_rmrse <- c("CLBW", "CLcrCL", "VcBW", "VcSEX")

power_out <- compute_power_block(oc_pilot_per_ds,
                                 n_true      = n_true_scn16,
                                 cn_cor_cut  = 1000)
oc_pilot_power    <- power_out$power
oc_pilot_relpower <- power_out$rel_power

rmrse_out <- compute_rmrse_block(
  per_ds      = oc_pilot_per_ds,
  true_long   = true_params,
  scenario_id = 16,
  pop_params  = pop_params_rmrse,
  cov_params  = cov_params_rmrse
)
oc_pilot_rmrse_uncond <- rmrse_out$rmrse_uncond
oc_pilot_rmrse_cond   <- rmrse_out$rmrse_cond

## Timing summary + scale-up projection
oc_pilot_timing <- tibble::tibble(
  metric = c("median_min_per_ds", "max_min_per_ds", "total_pilot_min",
             "projected_hours_100x16_seq",
             "projected_hours_100x16_workers4"),
  value  = c(
    stats::median(oc_pilot_per_ds$t_total_min, na.rm = TRUE),
    max(oc_pilot_per_ds$t_total_min,           na.rm = TRUE),
    sum(oc_pilot_per_ds$t_total_min,           na.rm = TRUE),
    stats::median(oc_pilot_per_ds$t_total_min, na.rm = TRUE) * 100 * 16 / 60,
    stats::median(oc_pilot_per_ds$t_total_min, na.rm = TRUE) * 100 * 16 / 60 / 4
  )
)


# ---- (P3.6) Persist OC artefacts --------------------------------------------
dir.create(stage1_pilot_scn16_N80, showWarnings = FALSE, recursive = TRUE)
saveRDS(oc_pilot_per_ds,
        file.path(stage1_pilot_scn16_N80, "oc_pilot_per_ds.rds"))
saveRDS(oc_pilot_power,
        file.path(stage1_pilot_scn16_N80, "oc_pilot_power.rds"))
saveRDS(oc_pilot_relpower,
        file.path(stage1_pilot_scn16_N80, "oc_pilot_relpower.rds"))
saveRDS(oc_pilot_rmrse_uncond,
        file.path(stage1_pilot_scn16_N80, "oc_pilot_rmrse_uncond.rds"))
saveRDS(oc_pilot_rmrse_cond,
        file.path(stage1_pilot_scn16_N80, "oc_pilot_rmrse_cond.rds"))
saveRDS(oc_pilot_timing,
        file.path(stage1_pilot_scn16_N80, "oc_pilot_timing.rds"))


# ---- SCALE-UP SCAFFOLD (NOT EXECUTED) --------------------------------------
##  Future-ready template for the 16-scenario x 100-dataset study.
##  Layered to keep behavior identical to the pilot when reused, but flips
##  three knobs that matter at scale:
##
##    1. Resume:        force_rerun = FALSE (default).  Reruns of any ds
##                       that already has its test_*.rds are instant.
##    2. Refit policy:  "full" keeps cov / SEs (current pilot setting).
##                       Switch to "skip" later to save ~30-60s per ds
##                       if PowerCN / parFixed are not needed.
##    3. Parallelism:   OUTER via future_pmap across (scenario x dataset);
##                       INNER via runSCM workers = 1L to prevent nested
##                       future plan thrashing on Windows.
##
##  Disk hygiene:
##    * `keep_res = FALSE` drops the candidate trail (~50-200 MB per run);
##       only the slim ~1-5 MB test_*.rds is kept per replicate.
##    * At 1600 runs that's the difference between ~3 GB and ~500 GB of
##       persistence.
##
##  Reproducibility:
##    * `furrr_options(seed = TRUE)` gives reproducible RNG per worker.
##    * `packages = ...` must list every package the worker R session needs
##       (auto-detection sometimes misses S4 methods).

if (FALSE) {  # guard: do not execute via source()

  library(future)
  library(furrr)

  ## Pick worker count: leave 1 physical core free for OS / IDE.
  n_workers_outer <- max(parallel::detectCores(logical = FALSE) - 1L, 1L)
  message(sprintf("outer workers = %d (inner = 1)", n_workers_outer))

  ## Plan grid: every (scenario, dataset) combination.
  plan_grid <- tidyr::expand_grid(scenario_id = 1:16, dataset_id = 1:100)

  ## Activate outer plan; revert when done.
  oplan <- future::plan(future::multisession, workers = n_workers_outer)
  on.exit(future::plan(oplan), add = TRUE)

  ## Per-job invocation.  All driver args are explicit so workers don't have
  ## to chase globals; furrr will serialize what's needed.
  scaleup_results <- furrr::future_pmap(
    plan_grid,
    function(scenario_id, dataset_id) {
      sim_long <- readRDS(file.path(
        "simulated_virtual_dataset_eta_filtered_N80",
        sprintf("sim_obs_scenario_%02d.rds", scenario_id)
      ))
      save_dir <- file.path(
        "simulated_virtual_dataset_eta_filtered_N80",
        sprintf("stage1_full_scn%02d_N80", scenario_id)
      )
      run_one_dataset_scn16_N80(
        ds_id        = dataset_id,
        save_dir     = save_dir,
        sim_long     = sim_long,
        workers      = 1L,         # disable inner parallelism under outer
        keep_res     = FALSE,      # drop heavy candidate trails
        force_rerun  = FALSE,      # resume previously-completed datasets
        scenario_id  = scenario_id
      )
    },
    .options = furrr::furrr_options(
      seed     = TRUE,
      globals  = TRUE,
      packages = c("nlmixr2", "nlmixr2est", "nlmixr2scm",
                   "rxode2", "dplyr", "purrr", "tibble", "tidyr")
    )
  )

  ## Optional: aggregate across scenarios using compute_power_block and
  ## compute_rmrse_block helpers, looping over scenario_id.
}
