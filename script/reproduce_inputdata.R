##   OUTPUT (canonical folder names: sim_obs_N40 / sim_obs_N80)
##     <OUT_ROOT>/sim_obs_N40/sim_obs_scenario_%02d.rds   (16 files)
##     <OUT_ROOT>/sim_obs_N80/sim_obs_scenario_%02d.rds   (16 files)
##   plus, in each folder:
##     virtual_population_<N_DATASETS>x<N>.rds
##     eta_table_valid.rds
##     sim_obs_all_scenarios.rdstasets so the
##   large `.rds` files never need to be transferred / pushed to GitHub.
##
##   Extracts steps 1-7 of PerformanceEvaluation04062026.R:
##     1. NHANES -> `pop`   (covariate population)
##     2. ref_cor <- cor(pop)
##     3. PsN_scenarios     (16-scenario design)
##     4. covariate + event-table prep (t_half_typ, ev_one)
##     5. eta helpers (sample_dataset_etas, build_eta_table)
##     6. simulation model + simulator (sim_mod_fixed, simulate_scenario_v2)
##     7. build_sister_cohort() + the N40 / N80 calls
##
##   OUTPUT (per cohort, 16 files each + aggregates):
##     <OUT_ROOT>/sim_obs_N40/sim_obs_scenario_%02d.rds
##     <OUT_ROOT>/sim_obs_N80/sim_obs_scenario_%02d.rds
##   plus virtual_population_*.rds, eta_table_valid.rds, sim_obs_all_scenarios.rds
##
##
## ----------------------------------------------------------------------------
## HOW TO RUN
## ----------------------------------------------------------------------------
##   Prerequisites:
##     * R packages: dplyr, haven, purrr, tibble, tidyr, rxode2, MASS
##     * The 3 NHANES XPT files in NHANES_DIR:
##         DEMO_L.XPT, BMX_L.XPT, BIOPRO_L.XPT
##
##   Configurable env vars (all optional -- sensible defaults below):
##     NHANES_DIR   folder holding the 3 XPT files
##                  (default: local OneDrive NHANES_dataset path)
##     OUT_ROOT     where sim_obs_N40/ and sim_obs_N80/ are written
##                  (default: "." = repo root)
##     N_DATASETS   datasets per scenario (default: 250)
##
##   -- Command line (from the repo root) ---------------------------------------
##     Rscript script/reproduce_inputdata.R
##
##     # override output location + dataset count:
##     OUT_ROOT=Inputdataset N_DATASETS=250 \
##       Rscript script/reproduce_inputdata.R
##
##     # point at a different NHANES folder:
##     NHANES_DIR=/path/to/NHANES_dataset \
##       Rscript script/reproduce_inputdata.R
##
##   -- Windows PowerShell ------------------------------------------------------
##     $env:OUT_ROOT   = "Inputdataset"
##     $env:N_DATASETS = "250"
##     Rscript script/reproduce_inputdata.R
##
##   -- Interactively in Positron / RStudio -------------------------------------
##     Sys.setenv(NHANES_DIR = "path/to/NHANES_dataset")
##     Sys.setenv(OUT_ROOT   = "Inputdataset")   # optional
##     source("script/reproduce_inputdata.R")
##
##   Runtime note: N_DATASETS=250 x 16 scenarios x 2 cohorts is compute-heavy;
##   set N_DATASETS lower (e.g. 5) for a quick smoke test first.
##   USAGE (from repo root):
##     Rscript script/reproduce_inputdata.R
##   or interactively: source("script/reproduce_inputdata.R")
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(rxode2)
  library(MASS)   # mvrnorm
})

## ---- Config -----------------------------------------------------------------
NHANES_DIR <- Sys.getenv(
  "NHANES_DIR",
  unset = "C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/NHANES_dataset"
)
OUT_ROOT   <- Sys.getenv("OUT_ROOT", unset = ".")
N_DATASETS <- as.integer(Sys.getenv("N_DATASETS", unset = "250"))

xpt_paths <- file.path(NHANES_DIR, c("DEMO_L.XPT", "BMX_L.XPT", "BIOPRO_L.XPT"))
if (!all(file.exists(xpt_paths))) {
  stop("Missing NHANES XPT file(s):\n  ",
       paste(xpt_paths[!file.exists(xpt_paths)], collapse = "\n  "),
       "\nSet NHANES_DIR (env var or top of this script).")
}

## ============================================================================
## STEP 1-2: NHANES -> pop -> ref_cor
## ============================================================================
demo   <- read_xpt(xpt_paths[1]) %>% dplyr::select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3)
bmx    <- read_xpt(xpt_paths[2]) %>% dplyr::select(SEQN, BMXWT, BMXHT, BMXBMI)
biopro <- read_xpt(xpt_paths[3]) %>% dplyr::select(SEQN, LBXSCR)

nhanes <- demo %>%
  dplyr::inner_join(bmx, by = "SEQN") %>%
  dplyr::inner_join(biopro, by = "SEQN")

pop <- nhanes %>%
  dplyr::filter(RIDAGEYR > 17, RIDRETH3 %in% c(3, 6)) %>%
  dplyr::mutate(
    SEX  = ifelse(RIAGENDR == 2, 0, 1),
    RACE = ifelse(RIDRETH3  == 6, 0, 1),
    BW   = BMXWT,
    BMI  = BMXBMI,
    CrCL = ((140 - RIDAGEYR) * BW / (72 * LBXSCR)) * ifelse(SEX == 0, 0.85, 1)
  ) %>%
  dplyr::select(BW, BMI, CrCL, SEX, RACE) %>%
  na.omit() %>%
  dplyr::filter(CrCL > 0, CrCL < 250)

ref_cor <- cor(pop)
message("ref_cor:"); print(round(ref_cor, 3))

## ============================================================================
## STEP 3: scenario design
## ============================================================================
PsN_scenarios <- data.frame(
  scenario = 1:16,
  I_BW_CL   = c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1),
  I_CRCL_CL = c(0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1),
  I_BW_VC   = c(0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1),
  I_SEX_VC  = c(0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1)
)

## ============================================================================
## STEP 4: covariate constants, half-life, sampling grid, event template
## ============================================================================
TVCL <- 0.6; TVQ <- 1.8; TVVc <- 20; TVVp <- 80; Ka <- 0.7
k10_typ <- TVCL / TVVc; k12_typ <- TVQ / TVVc; k21_typ <- TVQ / TVVp
sum_k    <- k10_typ + k12_typ + k21_typ
beta_typ <- 0.5 * (sum_k - sqrt(sum_k^2 - 4 * k10_typ * k21_typ))
t_half_typ <- log(2) / beta_typ

hl_mult      <- c(0, 0.05, 0.1, 0.5, 1, 3)
sample_times <- hl_mult * t_half_typ

DOSE_MG <- 100
ev_one <- rxode2::et(amt = DOSE_MG, cmt = "depot", time = 0) %>%
  rxode2::et(time = sample_times) %>%
  as.data.frame()

## ============================================================================
## STEP 5: eta helpers
## ============================================================================
ETA_RHO_LOW  <- 0.15
ETA_RHO_HIGH <- 0.25

sample_dataset_etas <- function(n_subj = 300, omega_var = 0.1, omega_cov = 0.02,
                                rho_low = ETA_RHO_LOW, rho_high = ETA_RHO_HIGH,
                                max_tries = 1000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- matrix(c(omega_var, omega_cov, omega_cov, omega_var), 2, 2)
  for (a in seq_len(max_tries)) {
    e   <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = Sigma)
    rho <- stats::cor(e[, 1], e[, 2])
    if (rho >= rho_low && rho <= rho_high) {
      return(list(eta_cl = e[, 1], eta_vc = e[, 2], rho = rho, attempts = a))
    }
  }
  stop("No valid etas after ", max_tries, " attempts.")
}

build_eta_table <- function(scenarios = PsN_scenarios$scenario,
                            n_datasets = 250L, n_subj = 300L) {
  out <- vector("list", length(scenarios) * n_datasets); k <- 0L
  for (scn in scenarios) for (ds in seq_len(n_datasets)) {
    k <- k + 1L
    eta <- sample_dataset_etas(n_subj = n_subj, seed = 2026L + scn * 1e4L + ds)
    out[[k]] <- tibble::tibble(
      SCENARIO = scn, DATASET = ds, SUBJECT = seq_len(n_subj),
      eta_cl = eta$eta_cl, eta_vc = eta$eta_vc,
      rho_realised = eta$rho, attempts = eta$attempts
    )
  }
  dplyr::bind_rows(out)
}

## ============================================================================
## STEP 6: simulation model + scenario simulator (etas as inputs)
## ============================================================================
sim_mod_fixed <- rxode2::rxode2({
  CL_cov = TVCL * (BW   / BW_REF  )^(TH_BW_CL   * I_BW_CL  ) *
                  (CrCL / CRCL_REF)^(TH_CRCL_CL * I_CRCL_CL)
  Vc_cov = TVVc * (BW   / BW_REF  )^(TH_BW_VC   * I_BW_VC  ) *
                  exp(TH_SEX_VC * SEX * I_SEX_VC)
  CL = CL_cov * exp(eta_cl)
  Vc = Vc_cov * exp(eta_vc)
  Q  = TVQ; Vp = TVVp; KA = TVKA
  k10 = CL / Vc; k12 = Q / Vc; k21 = Q / Vp
  d/dt(depot)      = -KA * depot
  d/dt(central)    =  KA * depot - k10 * central - k12 * central + k21 * peripheral
  d/dt(peripheral) =  k12 * central - k21 * peripheral
  cp = central / Vc
})

simulate_scenario_v2 <- function(scn, sim_mod = sim_mod_fixed,
                                 scenarios = PsN_scenarios, icov, etas,
                                 event_table, sampling_grid = sample_times,
                                 t_half = t_half_typ, prop_err = 0.1,
                                 seed_residual = 4242L + scn) {
  scn_row <- scenarios %>% dplyr::filter(scenario == scn)
  if (nrow(scn_row) != 1L) stop("scenario ", scn, " not found.")

  icov_scn <- icov %>%
    dplyr::left_join(
      etas %>% dplyr::filter(SCENARIO == scn) %>%
        dplyr::select(DATASET, SUBJECT, eta_cl, eta_vc),
      by = c("DATASET", "SUBJECT")
    ) %>%
    dplyr::select(id, BW, CrCL, SEX, eta_cl, eta_vc)

  params_scn <- c(
    TVCL = 0.6, TVQ = 1.8, TVVc = 20, TVVp = 80, TVKA = 0.7,
    BW_REF = 70, CRCL_REF = 95,
    TH_BW_CL = 0.75, TH_CRCL_CL = 0.5, TH_BW_VC = 1.0, TH_SEX_VC = log(1.5),
    I_BW_CL = scn_row$I_BW_CL, I_CRCL_CL = scn_row$I_CRCL_CL,
    I_BW_VC = scn_row$I_BW_VC, I_SEX_VC = scn_row$I_SEX_VC
  )

  sim_raw_scn <- rxode2::rxSolve(sim_mod, params = params_scn,
                                 events = event_table, iCov = icov_scn,
                                 returnType = "tibble")

  set.seed(seed_residual)
  sim_raw_scn <- sim_raw_scn %>%
    dplyr::mutate(cp_obs = cp * (1 + stats::rnorm(dplyr::n(), 0, prop_err)))

  sim_raw_scn %>%
    dplyr::filter(time %in% sampling_grid) %>%
    dplyr::mutate(SCENARIO = scn, cp_ipred = cp,
                  HL_MULT = round(time / t_half, 4)) %>%
    dplyr::select(-dplyr::any_of(c("eta_cl", "eta_vc"))) %>%
    dplyr::left_join(
      icov %>% dplyr::select(id, DATASET, SUBJECT, BMI, RACE), by = "id") %>%
    dplyr::left_join(
      etas %>% dplyr::filter(SCENARIO == scn) %>%
        dplyr::select(DATASET, SUBJECT, eta_cl, eta_vc),
      by = c("DATASET", "SUBJECT")) %>%
    dplyr::select(SCENARIO, DATASET, SUBJECT, HL_MULT, time,
                  BW, BMI, CrCL, SEX, RACE, eta_cl, eta_vc, CL, Vc,
                  cp_ipred, cp_obs)
}

## ============================================================================
## STEP 7: sister-cohort builder (bootstrap + eta table + simulate 16 scenarios)
## ============================================================================
build_sister_cohort <- function(n_subj, n_datasets = N_DATASETS,
                                tol_base = 0.05, tol_scale = 1.0,
                                max_attempts = 1e5L, pop_seed = 4242L + n_subj,
                                out_root = OUT_ROOT) {

  ## Folder name matches the Remote Explorer layout: sim_obs_N40 / sim_obs_N80
  out_dir_n <- file.path(out_root, sprintf("sim_obs_N%d", n_subj))
  if (!dir.exists(out_dir_n)) dir.create(out_dir_n, recursive = TRUE)

  ## (a) Stratified bootstrap by SEX x RACE, correlation-gated
  set.seed(pop_seed)
  cell_props <- pop %>%
    dplyr::count(SEX, RACE, name = "n_cell") %>%
    dplyr::mutate(p = n_cell / sum(n_cell), n_target = round(p * n_subj))
  drift <- n_subj - sum(cell_props$n_target)
  if (drift != 0L) {
    bump <- which.max(cell_props$p)
    cell_props$n_target[bump] <- cell_props$n_target[bump] + drift
  }
  pop_split <- split(pop, interaction(pop$SEX, pop$RACE, drop = TRUE))
  cell_key  <- paste(cell_props$SEX, cell_props$RACE, sep = ".")

  draw_one <- function() {
    purrr::map2_dfr(pop_split[cell_key], cell_props$n_target,
                    function(d, k) dplyr::slice_sample(d, n = k, replace = TRUE))
  }

  tol <- tol_base * tol_scale * sqrt(300 / n_subj)
  datasets_n <- vector("list", n_datasets)
  n_kept <- 0L; attempt_n <- 0L; n_skip <- 0L
  while (n_kept < n_datasets) {
    attempt_n <- attempt_n + 1L
    samp <- draw_one()
    col_sd <- vapply(samp, stats::sd, numeric(1))
    if (any(is.na(col_sd)) || any(col_sd == 0)) {
      n_skip <- n_skip + 1L
    } else {
      u <- abs(stats::cor(samp) - ref_cor)[upper.tri(ref_cor)]
      if (!anyNA(u) && all(u <= tol)) {
        n_kept <- n_kept + 1L; datasets_n[[n_kept]] <- samp
      }
    }
    if (attempt_n > max_attempts)
      stop(sprintf("Tolerance too strict at N=%d (kept %d / tried %d)",
                   n_subj, n_kept, attempt_n))
  }
  message(sprintf("[N=%d] tol=%.3f kept %d / tried %d (%d skipped)",
                  n_subj, tol, n_kept, attempt_n, n_skip))

  vp_n <- datasets_n %>%
    purrr::map(function(d) dplyr::mutate(d, ID = seq_len(n_subj))) %>%
    dplyr::bind_rows(.id = "DATASET") %>%
    dplyr::mutate(DATASET = as.integer(DATASET))
  saveRDS(vp_n, file.path(out_dir_n,
          sprintf("virtual_population_%dx%d.rds", n_datasets, n_subj)))

  ## (b) eta table
  eta_table_n <- build_eta_table(n_datasets = n_datasets, n_subj = n_subj)
  saveRDS(eta_table_n, file.path(out_dir_n, "eta_table_valid.rds"))

  ## (c) iCov + cohort-sized event table
  iCov_n <- vp_n %>%
    dplyr::mutate(SUBJECT = ID, id = (DATASET - 1L) * n_subj + ID) %>%
    dplyr::select(id, DATASET, SUBJECT, BW, BMI, CrCL, SEX, RACE)
  n_total_n <- n_datasets * n_subj
  ev_n <- ev_one %>%
    dplyr::slice(rep(dplyr::row_number(), n_total_n)) %>%
    dplyr::mutate(id = rep(seq_len(n_total_n), each = nrow(ev_one))) %>%
    dplyr::arrange(id, time)
  stopifnot(nrow(iCov_n) == n_total_n,
            setequal(iCov_n$id, unique(ev_n$id)))

  ## (d) simulate all 16 scenarios
  sim_list_n <- purrr::map(
    PsN_scenarios$scenario,
    function(scn) {
      message(sprintf("[N=%d] scenario %02d ...", n_subj, scn))
      out <- simulate_scenario_v2(scn, icov = iCov_n, etas = eta_table_n,
                                  event_table = ev_n)
      saveRDS(out, file.path(out_dir_n,
              sprintf("sim_obs_scenario_%02d.rds", scn)))
      out
    }, .progress = TRUE
  ) %>% purrr::set_names(sprintf("scenario_%02d", PsN_scenarios$scenario))

  sim_all_n <- dplyr::bind_rows(sim_list_n)
  saveRDS(sim_all_n, file.path(out_dir_n, "sim_obs_all_scenarios.rds"))

  invisible(list(virtual_pop = vp_n, eta_table = eta_table_n,
                 sim_obs_all = sim_all_n, out_dir = out_dir_n))
}

## ============================================================================
## Build both sister cohorts
## ============================================================================
cohort_40 <- build_sister_cohort(n_subj = 40L)
cohort_80 <- build_sister_cohort(n_subj = 80L)

message("\nDone. Wrote:")
message("  ", cohort_40$out_dir)
message("  ", cohort_80$out_dir)
