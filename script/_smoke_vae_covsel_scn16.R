# ============================================================================
# _smoke_vae_covsel_scn16.R   (TEST B -- VAE one-run covariate selection)
# ----------------------------------------------------------------------------
# GOAL
#   Test the VAE paper's headline claim ("one run is all you need"): VAE does
#   SIMULTANEOUS parameter estimation AND covariate selection in a SINGLE
#   nlmixr2(..., est = "vae") fit via its BICc-ELBO penalty.
#   *** runSCM / runSCM_traced are NOT used here. ***
#
# DESIGN (approach A -- how VAE actually works)
#   VAE runs its OWN covariate mechanism against the raw covariate COLUMNS in
#   the data.  It does NOT select from a hand-built candidate-theta superset.
#   So we hand VAE a PLAIN BASE MODEL (no covariate terms) plus a dataset that
#   carries every candidate covariate column, and let VAE auto-promote the
#   winners.  The promoted `beta_<PARAM>_<COV>` terms ARE the selected model.
#
#   TRUE scn16 relationships (what VAE should recover):
#     cl ~ BW   (power), cl ~ CrCL (power), vc ~ BW (power),
#     vc ~ SEX  (cat)   <- "cat", not "lin", per runSCM convention
#
# CRCL CASE-ALIAS (blocker fix)
#   VAE uppercases CrCL -> CRCL when building beta_lTVCL_CRCL; the solver then
#   needs a column named CRCL.  We add an uppercase CRCL alias (= CrCL).
#
# OUTPUT (schema 2.1, with a $covsel block)
#   output/smoke_vae_scm/scn16_<structure>/covsel/res_ds001.{rds,fit.rds,meta.json}
#
# RUN
#   Rscript script/_smoke_vae_covsel_scn16.R              # linCmt (default)
#   Rscript script/_smoke_vae_covsel_scn16.R --structure ode
# ============================================================================

suppressPackageStartupMessages({
  library(nlmixr2); library(nlmixr2est); library(rxode2)
  library(dplyr);   library(tibble);     library(nlmixr2scm)
})

.script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(fa)) return(normalizePath(dirname(fa[1]), mustWork = FALSE))
  if (dir.exists("script")) return("script")
  getwd()
})()

source(file.path(.script_dir, "refit_helpers.R"))
source(file.path(.script_dir, "true_model_factory.R"))
source(file.path(.script_dir, "scm_bench_helpers.R"))
source(file.path(.script_dir, "output_schema.R"))
source(file.path(.script_dir, "estimator_factory.R"))

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x

# ---- SCM covariate axes (runSCM convention: cont -> power/lin; cat -> "cat")
scm_bench_vars   <- c("cl", "vc")
scm_bench_covars <- c("BW", "CrCL", "BMI")
scm_bench_cats   <- c("SEX", "RACE")

parse_args <- function(argv) {
  opts <- list(N = 80L, scenario = 16L, dataset = 1L,
               structure = "linCmt",
               out_root = "output/smoke_vae_scm",
               input_root = "Inputdataset",
               true_params_path = "Inputdataset/true_params_long.rds")
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[i]; val <- function() { i <<- i + 1L; argv[i] }
    switch(a,
      "--N"           = { opts$N          <- as.integer(val()) },
      "--scenario"    = { opts$scenario   <- as.integer(val()) },
      "--dataset"     = { opts$dataset    <- as.integer(val()) },
      "--structure"   = { opts$structure  <- val() },
      "--out_root"    = { opts$out_root   <- val() },
      "--input_root"  = { opts$input_root <- val() },
      "--true_params" = { opts$true_params_path <- val() },
      stop(sprintf("Unknown arg: %s", a))
    )
    i <- i + 1L
  }
  opts
}

opts <- if (!interactive()) parse_args(commandArgs(trailingOnly = TRUE)) else parse_args(character(0))

sim_path <- file.path(opts$input_root, sprintf("sim_obs_N%d", opts$N),
                      sprintf("sim_obs_scenario_%02d.rds", opts$scenario))
stopifnot(file.exists(sim_path), file.exists(opts$true_params_path))

save_dir <- file.path(opts$out_root, sprintf("scn%02d_%s", opts$scenario, opts$structure), "covsel")
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
out_rds  <- file.path(save_dir, sprintf("res_ds%03d.rds", opts$dataset))

sim_all     <- readRDS(sim_path)
true_params <- readRDS(opts$true_params_path)

ds_i <- to_nm_dataset(sim_all) |>
  dplyr::filter(DATASET == opts$dataset) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE)) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::mutate(CRCL = CrCL)   # uppercase alias VAE requires

base_mod <- switch(opts$structure,
  "linCmt" = base_2cmt_oral_linCmt,
  "ode"    = base_2cmt_oral_ode,
  stop("Unknown structure: ", opts$structure)
)

vae_ctrl <- nlmixr2est::vaeControl(
  covariateSelection = TRUE,
  covMethod          = "r,s",
  calcTables         = TRUE,
  print              = 0
)

seed <- seed_for_dataset(opts$dataset); seed_all(seed, "vae")

message(sprintf(">>> [vae|covsel] N%d scn%02d ds%03d struct=%s start %s",
                opts$N, opts$scenario, opts$dataset, opts$structure,
                format(Sys.time(), "%H:%M:%S")))

t0  <- Sys.time()
fit <- tryCatch(
  nlmixr2(base_mod, ds_i, est = "vae", control = vae_ctrl),
  error = function(e) { message("VAE covsel ERROR: ", conditionMessage(e)); NULL }
)
runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

rec <- assemble_common(
  fit         = fit,
  true_mod    = base_mod,
  scenario_id = opts$scenario,
  true_params = true_params,
  runtime_sec = runtime_sec,
  status      = if (is.null(fit)) "error" else "ok",
  estimator   = "vae"
)

# ---- $covsel block: parse VAE's promoted beta_<PARAM>_<COV> terms ----------
.param_to_var <- c(lTVCL = "cl", lTVVc = "vc")
.norm_covar   <- function(x) ifelse(toupper(x) == "CRCL", "CrCL", x)

.parse_beta_theta <- function(nms) {
  hits <- grep("^beta_", nms, value = TRUE)
  if (length(hits) == 0L)
    return(tibble::tibble(theta_name = character(), var = character(), covar = character()))
  m <- regmatches(hits, regexec("^beta_(lTV[[:alnum:]]+)_(.+)$", hits))
  param <- vapply(m, function(x) if (length(x) == 3L) x[2] else NA_character_, character(1))
  covar <- vapply(m, function(x) if (length(x) == 3L) x[3] else NA_character_, character(1))
  tibble::tibble(theta_name = hits,
                 var   = unname(.param_to_var[param]),
                 covar = .norm_covar(covar))
}

selected <- if (!is.null(fit)) {
  th <- tryCatch(fit$theta, error = function(e) NULL)
  bt <- .parse_beta_theta(names(th))
  if (nrow(bt)) {
    bt |>
      dplyr::mutate(estimate = unname(th[theta_name])) |>
      dplyr::filter(!is.na(var)) |>
      dplyr::select(var, covar, theta_name, estimate) |>
      dplyr::arrange(dplyr::desc(abs(estimate)))
  } else NULL
} else NULL

true_set <- tibble::tibble(
  var   = c("cl",    "cl",    "vc",    "vc"),
  covar = c("BW",    "CrCL",  "BW",    "SEX"),
  shape = c("power", "power", "power", "cat")
)

rec$covsel <- list(selected = selected, true_set = true_set)

# ---- VAE-ONLY patch: backfill covariate estimates into rec$rel_err ---------
# assemble_common() / extract_params_long() only know the runSCM/FOCEi covariate
# theta names (TH_*, cov_*), so VAE's beta_<PARAM>_<COV> terms never join and
# appear as NA in rec$rel_err.  We DO NOT modify assemble_common (that stays the
# FOCEi path for the HPCE cohort*scenario runs).  Map VAE's selected betas onto
# the true_params parameter labels here and refill just those rows.
#
# Comparability: for power (continuous) covariates beta is the exponent on
# log(cov/ref) and is centring-invariant, so it is directly comparable to the
# true value regardless of VAE's data-driven reference.  For SEX the beta is the
# log-effect and is likewise directly comparable.  (Structural intercepts lTVCL/
# lTVVc still absorb VAE's reference shift -- expected, not a bug.)
.vae_true_param <- function(var, covar) {
  key <- paste(var, toupper(covar), sep = "|")
  unname(c(
    "cl|BW"   = "CLBW",
    "cl|CRCL" = "CLcrCL",
    "vc|BW"   = "VcBW",
    "vc|SEX"  = "VcSEX"
  )[key])
}

if (!is.null(selected) && !is.null(rec$rel_err)) {
  sel_named <- selected |>
    dplyr::mutate(parameter = .vae_true_param(var, covar)) |>
    dplyr::filter(!is.na(parameter)) |>
    dplyr::select(parameter, vae_estimate = estimate)

  rec$rel_err <- rec$rel_err |>
    dplyr::left_join(sel_named, by = "parameter") |>
    dplyr::mutate(
      estimate    = dplyr::coalesce(estimate, vae_estimate),
      abs_err     = ifelse(is.na(estimate) | is.na(true_value),
                           NA_real_, estimate - true_value),
      rel_err     = ifelse(is.na(estimate) | is.na(true_value) | true_value == 0,
                           NA_real_, (estimate - true_value) / true_value),
      rel_err_pct = rel_err * 100
    ) |>
    dplyr::select(-vae_estimate)
}

rec$sample_N    <- opts$N
rec$scenario_id <- opts$scenario
rec$dataset_id  <- opts$dataset
rec$estimator   <- "vae"
rec$outer_opt   <- NA_character_

write_fit_sidecar(rec, out_rds, fit = fit)

message(sprintf("<<< [vae|covsel] done in %.1fs  objf=%s  converged=%s",
                runtime_sec, format(rec$objf %||% NA), format(rec$converged %||% NA)))

if (!is.null(selected)) {
  message("\nPromoted covariate terms (beta_*, |est| desc):")
  print(selected, n = nrow(selected))

  cmp <- dplyr::full_join(
    dplyr::mutate(dplyr::select(true_set, var, covar), in_true = TRUE),
    dplyr::mutate(dplyr::select(selected, var, covar), in_vae  = TRUE),
    by = c("var", "covar")
  ) |>
    dplyr::distinct() |>
    dplyr::mutate(
      in_true = tidyr::replace_na(in_true, FALSE),
      in_vae  = tidyr::replace_na(in_vae,  FALSE),
      verdict = dplyr::case_when(
        in_true &  in_vae ~ "TP (correct)",
        in_true & !in_vae ~ "FN (missed)",
       !in_true &  in_vae ~ "FP (spurious)"
      )
    )
  message("\nSelection vs. truth (covariate x parameter):")
  print(cmp, n = nrow(cmp))

  n_tp <- sum(cmp$verdict == "TP (correct)")
  n_fn <- sum(cmp$verdict == "FN (missed)")
  n_fp <- sum(cmp$verdict == "FP (spurious)")
  message(sprintf("TP=%d  FN=%d  FP=%d  -> %s",
                  n_tp, n_fn, n_fp,
                  if (n_fn == 0L && n_fp == 0L) "PASS (exact recovery)" else "CHECK"))
} else {
  message("\nNo beta_* covariate terms found on the fit (selection empty or fit failed).")
}

message(sprintf("\nSaved: %s", out_rds))
