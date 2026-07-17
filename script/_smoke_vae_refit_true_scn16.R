# ============================================================================
# _smoke_vae_refit_true_scn16.R   (TASK 1 -- VAE fixed-structure true refit)
# ----------------------------------------------------------------------------
# GOAL
#   Test VAE's POPULATION-ONLY / fixed-structure mode: supply the TRUE covariate
#   model and estimate its population parameters + covariate coefficients WITHOUT
#   any selection.  This is the "refit the true model" task (Task 1), the direct
#   VAE analogue of the FOCEi refit bench -- NOT the one-run covariate selection
#   (that is _smoke_vae_covsel_scn16.R / vae_covsel_driver.R).
#
#   nlmixr2(true_mod, ds, est = "vae",
#           control = vaeControl(covariateSelection = FALSE))
#
# WHY A SEPARATE (FLATTENED) MODEL
#   make_true_model() emits covariate terms on intermediate `_typ` lines
#   (lTVCL_typ <- lTVCL + TH_BW_CL*log(BW/70); cl <- exp(lTVCL_typ + eta.cl)).
#   VAE errors ("cannot find parameter 'NA'") when covariates sit on a `_typ`
#   line, so here we build a FLATTENED twin with every covariate inlined
#   directly inside exp(...).  The FOCEi factory is left untouched.
#
#   The flattened model keeps the SAME theta names as the factory
#   (TH_BW_CL / TH_CRCL_CL / TH_BW_VC / TH_SEX_VC), so extract_params_long()
#   maps them onto the true labels (CLBW / CLcrCL / VcBW / VcSEX) and
#   assemble_common() populates rel_err with NO beta backfill needed.
#
# TRUTH (scn16): all four active
#   cl ~ BW (power, 0.75), cl ~ CrCL (power, 0.50),
#   vc ~ BW (power, 1.00), vc ~ SEX (log-shift, log(1.5) = 0.405)
#
# CRCL CASE-ALIAS: VAE uppercases CrCL -> CRCL, so add a CRCL alias column.
#
# OUTPUT (schema 2.1)
#   output/smoke_vae_refit/scn<SS>_<structure>/refit/res_ds<DDD>.{rds,fit.rds,meta.json}
#
# RUN
#   Rscript script/_smoke_vae_refit_true_scn16.R                 # linCmt (default)
#   Rscript script/_smoke_vae_refit_true_scn16.R --structure ode
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
               out_root = "output/smoke_vae_refit",
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

# ---- Inputs ----------------------------------------------------------------
sim_path <- file.path(opts$input_root, sprintf("sim_obs_N%d", opts$N),
                      sprintf("sim_obs_scenario_%02d.rds", opts$scenario))
stopifnot(file.exists(sim_path), file.exists(opts$true_params_path))

save_dir <- file.path(opts$out_root, sprintf("scn%02d_%s", opts$scenario, opts$structure), "refit")
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

# ---- FOCEi true model (for attrs + model_type in assemble_common) ----------
true_mod <- make_true_model(opts$scenario, boundary = "none",
                            structure = opts$structure)

# ---- Flattened VAE twin: covariates inlined inside exp() -------------------
# Same theta names + true structure as make_true_model(), but no `_typ` lines.
# Covariate terms are conditional on the PsN_scenarios indicators so the twin
# tracks make_true_model() for ALL scenarios (scn1 = null ... scn16 = all four).
.build_vae_true_model <- function(scenario_id, structure) {
  scn <- PsN_scenarios[PsN_scenarios$scenario == scenario_id, , drop = FALSE]
  if (nrow(scn) != 1L) stop("scenario ", scenario_id, " not in PsN_scenarios")

  ini_lines <- c(
    "lTVCL <- log(0.6)",
    "lTVQ  <- log(1.8)",
    "lTVVc <- log(20)",
    "lTVVp <- log(80)",
    "lTVKA <- fix(log(0.7))"
  )
  cl_terms <- character(0); vc_terms <- character(0)
  if (isTRUE(scn$I_BW_CL   == 1)) { ini_lines <- c(ini_lines, "TH_BW_CL   <- 0.5");        cl_terms <- c(cl_terms, "TH_BW_CL * log(BW / 70)") }
  if (isTRUE(scn$I_CRCL_CL == 1)) { ini_lines <- c(ini_lines, "TH_CRCL_CL <- 0.5");        cl_terms <- c(cl_terms, "TH_CRCL_CL * log(CRCL / 95)") }
  if (isTRUE(scn$I_BW_VC   == 1)) { ini_lines <- c(ini_lines, "TH_BW_VC   <- 0.5");        vc_terms <- c(vc_terms, "TH_BW_VC * log(BW / 70)") }
  if (isTRUE(scn$I_SEX_VC  == 1)) { ini_lines <- c(ini_lines, "TH_SEX_VC  <- 0.405");      vc_terms <- c(vc_terms, "TH_SEX_VC * SEX") }
  ini_lines <- c(ini_lines,
                 "eta.cl + eta.vc ~ c(0.1, 0.02, 0.1)",
                 "prop.err <- 0.1")

  cl_rhs <- paste(c("lTVCL", cl_terms), collapse = " + ")
  vc_rhs <- paste(c("lTVVc", vc_terms), collapse = " + ")
  param_lines <- c(
    sprintf("cl <- exp(%s + eta.cl)", cl_rhs),
    sprintf("vc <- exp(%s + eta.vc)", vc_rhs),
    "q  <- exp(lTVQ)",
    "vp <- exp(lTVVp)",
    "ka <- exp(lTVKA)"
  )
  struct_lines <- switch(structure,
    "linCmt" = c("cp <- linCmt()", "cp ~ prop(prop.err)"),
    "ode" = c(
      "d/dt(depot)   <- -ka * depot",
      "d/dt(central) <-  ka * depot - (cl / vc) * central - (q / vc) * central + (q / vp) * periph",
      "d/dt(periph)  <-  (q / vc) * central - (q / vp) * periph",
      "cp <- central / vc",
      "cp ~ prop(prop.err)"
    ),
    stop("Unknown structure: ", structure)
  )
  fn_text <- paste0(
    "function() {\n",
    "  ini({\n    ", paste(ini_lines, collapse = "\n    "), "\n  })\n",
    "  model({\n    ", paste(c(param_lines, struct_lines), collapse = "\n    "), "\n  })\n",
    "}"
  )
  eval(parse(text = fn_text), envir = globalenv())
}

vae_mod <- .build_vae_true_model(opts$scenario, opts$structure)

# ---- VAE control: fixed structure (selection OFF) --------------------------
vae_ctrl <- nlmixr2est::vaeControl(
  covariateSelection = FALSE,   # <-- Task 1: estimate the GIVEN structure only
  covMethod          = "r,s",
  calcTables         = TRUE,
  print              = 0
)

# ---- Fit -------------------------------------------------------------------
seed <- seed_for_dataset(opts$dataset); seed_all(seed, "vae")

message(sprintf(">>> [vae|refit-true] N%d scn%02d ds%03d struct=%s start %s",
                opts$N, opts$scenario, opts$dataset, opts$structure,
                format(Sys.time(), "%H:%M:%S")))

t0  <- Sys.time()
fit <- tryCatch(
  nlmixr2(vae_mod, ds_i, est = "vae", control = vae_ctrl),
  error = function(e) { message("VAE refit ERROR: ", conditionMessage(e)); NULL }
)
runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# ---- Schema-2.1 record -----------------------------------------------------
# true_mod (the FOCEi factory UI) supplies attributes / model_type; the fit is
# the VAE twin's. Theta names match, so rel_err populates directly.
rec <- assemble_common(
  fit         = fit,
  true_mod    = true_mod,
  scenario_id = opts$scenario,
  true_params = true_params,
  runtime_sec = runtime_sec,
  status      = if (is.null(fit)) "error" else "ok",
  estimator   = "vae"
)

rec$sample_N    <- opts$N
rec$scenario_id <- opts$scenario
rec$dataset_id  <- opts$dataset
rec$estimator   <- "vae"
rec$structure   <- opts$structure
rec$outer_opt   <- NA_character_
rec$vae_mode    <- "refit_true_fixed"   # tag: fixed-structure, not selection

write_fit_sidecar(rec, out_rds, fit = fit)

# ---- Console report --------------------------------------------------------
message(sprintf("<<< [vae|refit-true] done in %.1fs  objf=%s  converged=%s",
                runtime_sec, format(rec$objf %||% NA), format(rec$converged %||% NA)))

if (!is.null(rec$rel_err)) {
  message("\nrel_err (estimate vs truth):")
  print(rec$rel_err, n = nrow(rec$rel_err))
  cov_rows <- rec$rel_err |>
    dplyr::filter(parameter %in% c("CLBW", "CLcrCL", "VcBW", "VcSEX"),
                  !is.na(rel_err_pct))
  if (nrow(cov_rows)) {
    message("\nCovariate coefficients (rel_err %):")
    print(dplyr::select(cov_rows, parameter, true_value, estimate, rel_err_pct))
  }
} else {
  message("\nNo rel_err (fit failed).")
}

message(sprintf("\nSaved: %s", out_rds))
