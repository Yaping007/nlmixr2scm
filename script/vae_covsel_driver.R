# ============================================================================
# vae_covsel_driver.R   (HPCE VAE one-run covariate-selection driver)
# ----------------------------------------------------------------------------
# GENERALIZED from script/_smoke_vae_covsel_scn16.R for the HPCE pilot grid:
#   3 cohorts (N40/N80/N300) x 16 scenarios x 2 structures (linCmt/ode).
#
# WHAT IT DOES
#   Runs VAE's SINGLE-fit simultaneous parameter estimation + covariate
#   selection (nlmixr2(..., est = "vae"), covariateSelection = TRUE) on one
#   (N, scenario, structure, dataset) cell, and writes a schema-2.1 record.
#   *** runSCM / runSCM_traced are NOT used -- selection is internal to VAE. ***
#
# KEY DIFFERENCE FROM THE SMOKE FILE
#   The smoke file hard-codes the scn16 truth (all four covariates). Here the
#   expected relationship set (`true_set`) is DERIVED PER SCENARIO from the
#   PsN_scenarios indicator columns (I_BW_CL, I_CRCL_CL, I_BW_VC, I_SEX_VC),
#   so the TP/FN/FP verdict and the rel_err backfill are correct for all 16
#   scenarios (scenario 1 = null model, scenario 16 = all four active).
#
# PREFLIGHT
#   VAE needs the `torch` R package (+ a working libtorch backend). On an
#   offline compute node these may be missing, so we check up front and, on
#   failure, write a <name>_ERROR.txt sidecar and exit non-zero (loud fail,
#   matching the existing bench ERROR-file convention).
#
# REUSED INFRASTRUCTURE (unchanged; assemble_common stays the FOCEi path)
#   refit_helpers.R      : to_nm_dataset, PsN_scenarios
#   true_model_factory.R : (loaded for parity with smoke file)
#   scm_bench_helpers.R  : base_2cmt_oral_linCmt / base_2cmt_oral_ode
#   output_schema.R      : assemble_common (schema 2.1), write_fit_sidecar
#   estimator_factory.R  : seed_for_dataset, seed_all
#
# OUTPUT (schema 2.1, with a $covsel block)
#   <out_root>/N<N>/scn<SS>_<structure>/covsel/res_ds<DDD>.{rds,fit.rds,meta.json}
#
# RUN
#   Rscript script/vae_covsel_driver.R --N 80 --scenario 16 \
#       --structure linCmt --dataset 1 --out_root output/vae_covsel_pilot
# ============================================================================

suppressPackageStartupMessages({
  library(nlmixr2); library(nlmixr2est); library(rxode2)
  library(dplyr);   library(tibble);     library(nlmixr2scm)
})

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x

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

# ---- CLI -------------------------------------------------------------------
parse_args <- function(argv) {
  opts <- list(N = 80L, scenario = 16L, dataset = 1L,
               structure = "linCmt",
               out_root = "output/vae_covsel_pilot",
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

# ---- Output paths (needed early so preflight can write an ERROR sidecar) ---
save_dir <- file.path(opts$out_root, sprintf("N%d", opts$N),
                      sprintf("scn%02d_%s", opts$scenario, opts$structure), "covsel")
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
out_rds  <- file.path(save_dir, sprintf("res_ds%03d.rds", opts$dataset))
err_txt  <- file.path(save_dir, sprintf("res_ds%03d_ERROR.txt", opts$dataset))

.write_error <- function(msg) {
  writeLines(c(
    sprintf("VAE covsel FAILED  %s", format(Sys.time())),
    sprintf("cell: N%d scn%02d %s ds%03d",
            opts$N, opts$scenario, opts$structure, opts$dataset),
    "", msg
  ), err_txt)
  message("ERROR sidecar written: ", err_txt)
}

# ---- Torch preflight -------------------------------------------------------
# VAE (LSTM encoder) requires the `torch` R package AND a usable libtorch
# backend (torch::install_torch() downloads it -- may be absent offline).
.torch_ok <- function() {
  if (!requireNamespace("torch", quietly = TRUE)) {
    return(list(ok = FALSE, why = "R package 'torch' is not installed on this node."))
  }
  installed <- tryCatch(isTRUE(torch::torch_is_installed()),
                        error = function(e) FALSE)
  if (!installed) {
    return(list(ok = FALSE, why = paste(
      "'torch' is installed but the libtorch backend is not available",
      "(torch::torch_is_installed() == FALSE). Run torch::install_torch()",
      "once on an ONLINE login node to download libtorch.")))
  }
  list(ok = TRUE, why = "")
}

pf <- .torch_ok()
if (!pf$ok) {
  .write_error(paste0("TORCH PREFLIGHT FAILED\n", pf$why))
  stop(pf$why, call. = FALSE)
}

# ---- Inputs ----------------------------------------------------------------
sim_path <- file.path(opts$input_root, sprintf("sim_obs_N%d", opts$N),
                      sprintf("sim_obs_scenario_%02d.rds", opts$scenario))
if (!file.exists(sim_path))              { .write_error(sprintf("missing sim file: %s", sim_path)); stop("missing sim file", call. = FALSE) }
if (!file.exists(opts$true_params_path)) { .write_error(sprintf("missing true_params: %s", opts$true_params_path)); stop("missing true_params", call. = FALSE) }

sim_all     <- readRDS(sim_path)
true_params <- readRDS(opts$true_params_path)

ds_i <- to_nm_dataset(sim_all) |>
  dplyr::filter(DATASET == opts$dataset) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE)) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::mutate(CRCL = CrCL)   # uppercase alias VAE requires

if (nrow(ds_i) == 0L) {
  .write_error(sprintf("no rows for dataset %d in %s (max DATASET = %s)",
                       opts$dataset, sim_path,
                       tryCatch(max(sim_all$DATASET), error = function(e) NA)))
  stop("empty dataset slice", call. = FALSE)
}

# ---- Per-scenario TRUE covariate set (from PsN_scenarios indicators) -------
# I_BW_CL -> cl~BW (power); I_CRCL_CL -> cl~CrCL (power);
# I_BW_VC -> vc~BW (power); I_SEX_VC -> vc~SEX (cat).
.true_set_for <- function(scenario_id) {
  scn <- PsN_scenarios[PsN_scenarios$scenario == scenario_id, , drop = FALSE]
  if (nrow(scn) != 1L)
    stop(sprintf("scenario %s not found in PsN_scenarios", scenario_id))
  rows <- list()
  if (isTRUE(scn$I_BW_CL   == 1)) rows[[length(rows)+1L]] <- tibble::tibble(var="cl", covar="BW",   shape="power")
  if (isTRUE(scn$I_CRCL_CL == 1)) rows[[length(rows)+1L]] <- tibble::tibble(var="cl", covar="CrCL", shape="power")
  if (isTRUE(scn$I_BW_VC   == 1)) rows[[length(rows)+1L]] <- tibble::tibble(var="vc", covar="BW",   shape="power")
  if (isTRUE(scn$I_SEX_VC  == 1)) rows[[length(rows)+1L]] <- tibble::tibble(var="vc", covar="SEX",  shape="cat")
  if (length(rows) == 0L)
    return(tibble::tibble(var=character(), covar=character(), shape=character()))
  dplyr::bind_rows(rows)
}
true_set <- .true_set_for(opts$scenario)

# ---- Base model (NO covariate terms; VAE auto-selects) ---------------------
base_mod <- switch(opts$structure,
  "linCmt" = base_2cmt_oral_linCmt,
  "ode"    = base_2cmt_oral_ode,
  stop("Unknown structure: ", opts$structure)
)

# ---- VAE control: covariate selection ON (the "one run") -------------------
vae_ctrl <- nlmixr2est::vaeControl(
  covariateSelection = TRUE,
  covMethod          = "r,s",
  calcTables         = TRUE,
  print              = 0
)

# ---- Fit -------------------------------------------------------------------
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

if (is.null(fit)) {
  .write_error("nlmixr2(est = 'vae') returned NULL (fit error; see stderr).")
}

# ---- Schema-2.1 record -----------------------------------------------------
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

rec$covsel <- list(selected = selected, true_set = true_set)

# ---- VAE-ONLY patch: backfill covariate estimates into rec$rel_err ---------
# VAE's beta_<PARAM>_<COV> names never match the FOCEi labels in true_params
# (CLBW/CLcrCL/VcBW/VcSEX), so those rows are NA out of assemble_common. We map
# them here (comparability: power betas are centring-invariant exponents; SEX
# beta is the log-effect). assemble_common itself is untouched (FOCEi path).
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

  if (nrow(sel_named)) {
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
}

# ---- Identity keys ---------------------------------------------------------
rec$sample_N    <- opts$N
rec$scenario_id <- opts$scenario
rec$dataset_id  <- opts$dataset
rec$estimator   <- "vae"
rec$structure   <- opts$structure
rec$outer_opt   <- NA_character_

write_fit_sidecar(rec, out_rds, fit = fit)

# On a successful (non-error) write, remove any stale ERROR sidecar left by a
# previous failed attempt at this same cell, so the directory reflects reality.
if (!is.null(fit) && file.exists(err_txt)) unlink(err_txt)

# ---- Console report --------------------------------------------------------
message(sprintf("<<< [vae|covsel] done in %.1fs  objf=%s  converged=%s",
                runtime_sec, format(rec$objf %||% NA), format(rec$converged %||% NA)))

sel_pairs <- if (!is.null(selected)) {
  dplyr::select(selected, var, covar)
} else {
  tibble::tibble(var = character(), covar = character())
}

cmp <- dplyr::full_join(
  dplyr::mutate(dplyr::select(true_set, var, covar), in_true = TRUE),
  dplyr::mutate(sel_pairs, in_vae = TRUE),
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

if (nrow(cmp)) {
  message("\nSelection vs. truth (covariate x parameter):")
  print(cmp, n = nrow(cmp))
  n_tp <- sum(cmp$verdict == "TP (correct)")
  n_fn <- sum(cmp$verdict == "FN (missed)")
  n_fp <- sum(cmp$verdict == "FP (spurious)")
  message(sprintf("TP=%d  FN=%d  FP=%d  -> %s",
                  n_tp, n_fn, n_fp,
                  if (n_fn == 0L && n_fp == 0L) "PASS (exact recovery)" else "CHECK"))
} else {
  message("\nScenario ", opts$scenario, " is the null model (no true covariates) ",
          "and VAE promoted none -> PASS.")
}

message(sprintf("\nSaved: %s", out_rds))
