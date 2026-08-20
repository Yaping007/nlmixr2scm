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

# Location-robust helper sourcing: the five helpers live in the top-level
# `script/` dir, but this driver may itself sit in `script/` OR a subfolder
# (e.g. `script/hpce_vae_covsel/`).  Search `.script_dir`, its parent, and
# `script/` so the driver works from either location without edits.
.source_helper <- function(f) {
  cand <- unique(c(file.path(.script_dir, f),
                   file.path(dirname(.script_dir), f),
                   file.path("script", f)))
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("helper not found: ", f,
                         " (looked in: ", paste(cand, collapse = ", "), ")")
  source(hit[1])
}
invisible(lapply(c("refit_helpers.R", "true_model_factory.R",
                   "scm_bench_helpers.R", "output_schema.R",
                   "estimator_factory.R"), .source_helper))

# ---- CLI -------------------------------------------------------------------
parse_args <- function(argv) {
  opts <- list(N = 80L, scenario = 16L, dataset = 1L,
               structure = "linCmt",
               out_root = "output/vae_covsel_pilot",
               input_root = "Inputdataset",
               true_params_path = "Inputdataset/true_params_long.rds",
               # Continuous-covariate shape menu (space-separated). Default
               # "power lin" = competing shapes; --shapes power = power-only.
               shapes = "power lin",
               # PR #921 correlation-aware covariate selection. "on" (default)
               # enables covSelectColinear; "off" reproduces the pre-PR engine.
               colinear = "on",
               # Covariate-colinearity clustering cut (mechanism A). Raised to
               # 0.95 so BW-BMI (r ~ 0.88) never clusters -> isolates the
               # Vc-CL parameter-correlation refinement (mechanism B).
               colinear_cut = 0.95)
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
      "--shapes"      = { opts$shapes      <- val() },
      "--colinear"    = { opts$colinear    <- val() },
      "--colinear_cut" = { opts$colinear_cut <- as.numeric(val()) },
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

# ---- Engine note -----------------------------------------------------------
# nlmixr2est's est = "vae" is a NATIVE C++/Armadillo implementation (LSTM
# encoder with hand-derived analytic backward; decoder = rxode2 solve). It does
# NOT use the `torch` R package or a libtorch backend -- so there is no torch
# preflight. The only optional extra is the `L0Learn` package, used solely for
# large covariate searches (covSelectMethod = "l0learn"/"auto"); with the four
# candidate covariates here the exact branch-and-bound runs and L0Learn is not
# required.

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
# Search space PINNED to the runSCM / PsN candidate set (nlmixr2est 7.0.2).
# This is SCENARIO-INDEPENDENT: every scenario 1..16 searches the SAME menu of
# candidates (cl,vc x BW,CrCL,BMI in {power,lin} + SEX,RACE cat = 16 candidates);
# only the TRUE active subset differs per scenario, and that is captured wholly
# by `true_set` (derived above from PsN_scenarios).  Naming the covariates with
# fixCov=TRUE restricts the search to exactly these five; SEX/RACE = TRUE marks
# them eligible (auto "cat").  covCenter pins PsN's physiological anchors
# (BW=70, CrCL=95); BMI keeps the data median (covCenterType="median"), matching
# PsN (BMI at median) and SCM (median centering).  The power exponent is
# centring-invariant, so coefficients stay comparable to truth.
# Continuous-covariate shape menu from --shapes (space-separated). Default
# c("power","lin") = competing shapes; c("power") = the power-only covariate
# space. SEX/RACE stay categorical (TRUE -> auto "cat").
.cont_shapes <- strsplit(trimws(opts$shapes %||% "power lin"), "\\s+")[[1]]
# PR #921 A/B toggle: "on" -> covSelectColinear = TRUE (correlation-aware
# selection); "off" -> FALSE (pre-PR per-dimension engine).
.colinear_on <- identical(tolower(trimws(opts$colinear %||% "on")), "on")
.colinear_cut <- suppressWarnings(as.numeric(opts$colinear_cut %||% 0.95))
if (!is.finite(.colinear_cut)) .colinear_cut <- 0.95
vae_ctrl <- nlmixr2est::vaeControl(
  covariateSelection = TRUE,
  covSelectColinear  = .colinear_on,
  covSelectColinearCut = .colinear_cut,
  shapes = list(
    BW   = .cont_shapes,
    CrCL = .cont_shapes,
    BMI  = .cont_shapes,
    SEX  = TRUE,       # categorical -> auto "cat"
    RACE = TRUE,       # categorical -> auto "cat"
    fixCov = TRUE      # search ONLY these five covariates (match SCM/PsN)
  ),
  covCenter          = c(BW = 70, CrCL = 95),  # PsN physiological anchors
  covCenterType      = "median",               # BMI -> median (SCM/PsN default)
  catCutoff          = 0.05,
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

# ---- $covsel block: parse VAE's promoted beta coefficients -----------------
# nlmixr2est 7.0.2 renamed the promoted covariate coefficients from the 7.0.1
# underscore form (beta_lTVCL_CRCL) to a DOT-separated form:
#   beta.<param>.<cov>.<shape>   e.g. beta.lTVCL.CRCL.power / beta.lTVVc.BW.lin
#   beta.<param>.<cov>[.<level>] for categoricals   e.g. beta.lTVVc.SEX
# The parser below matches BOTH separators (^beta[._]) and captures the shape.
.param_to_var <- c(lTVCL = "cl", lTVVc = "vc", cl = "cl", vc = "vc")
.norm_covar   <- function(x) ifelse(toupper(x) == "CRCL", "CrCL", x)
.shape_tokens <- c("power", "lin", "log", "identity", "center",
                   "hockey", "hockeyLow", "hockeyHi")

.parse_beta_theta <- function(nms) {
  empty <- tibble::tibble(theta_name = character(), var = character(),
                          covar = character(), shape = character())
  hits <- grep("^beta[._]", nms, value = TRUE)
  if (length(hits) == 0L) return(empty)
  rows <- lapply(hits, function(nm) {
    toks <- strsplit(sub("^beta[._]", "", nm), "[._]")[[1]]
    if (length(toks) < 2L) return(NULL)
    param <- toks[1]
    last  <- toks[length(toks)]
    if (last %in% .shape_tokens) {
      shape <- if (grepl("^hockey", last)) "hockey" else last
      covar <- paste(toks[-c(1L, length(toks))], collapse = ".")
    } else {
      shape <- "cat"; covar <- paste(toks[-1L], collapse = ".")
    }
    tibble::tibble(theta_name = nm,
                   var   = unname(.param_to_var[param]),
                   covar = .norm_covar(covar), shape = shape)
  })
  out <- dplyr::bind_rows(rows)
  if (is.null(out) || nrow(out) == 0L) empty else out
}

# Version-robust source of the promoted beta_* coefficients.  The updated
# nlmixr2est (>= 7.0.0) returns NULL for `fit$theta`, but the estimated fixed
# effects -- INCLUDING the beta_<PARAM>_<COV> covariate terms -- are still in
# the parFixed table (rownames = parameter, column "Estimate").  We therefore
# read a NAMED estimate vector from, in order of preference:
#   1. rec$parFixed  (already assembled, carries the beta_ rows)   [preferred]
#   2. fit$parFixedDf / fit$parFixed  (same table off the fit)
#   3. fit$theta     (legacy path; NULL on the new est, hence the fallback)
# Returns a named numeric vector (names = parameter) or NULL.
.named_theta_estimates <- function(rec, fit) {
  # 1/2: any parFixed-style data.frame with an Estimate column + rownames
  for (pf in list(rec$parFixed,
                  tryCatch(fit$parFixedDf, error = function(e) NULL),
                  tryCatch(fit$parFixed,   error = function(e) NULL))) {
    if (is.data.frame(pf) && nrow(pf) &&
        "Estimate" %in% colnames(pf) && !is.null(rownames(pf))) {
      v <- suppressWarnings(as.numeric(pf[["Estimate"]]))
      names(v) <- rownames(pf)
      return(v)
    }
  }
  # 3: legacy named theta vector
  th <- tryCatch(fit$theta, error = function(e) NULL)
  if (!is.null(th) && length(th)) return(th)
  NULL
}

th_est   <- .named_theta_estimates(rec, fit)
selected <- if (!is.null(fit) && !is.null(th_est)) {
  bt <- .parse_beta_theta(names(th_est))
  if (nrow(bt)) {
    bt |>
      dplyr::mutate(estimate = unname(th_est[theta_name])) |>
      dplyr::filter(!is.na(var)) |>
      dplyr::select(var, covar, shape, theta_name, estimate) |>
      dplyr::arrange(dplyr::desc(abs(estimate)))
  } else NULL
} else NULL

# Loud diagnostic: a converged VAE that promoted NOTHING is plausible only for
# the null scenario (scn 1).  If beta terms are absent for a scenario that HAS
# true covariates, the extraction path is broken -- warn so it is not silently
# scored as all-FN (the exact failure mode of the fit$theta -> NULL regression).
if (!is.null(fit)) {
  n_true_here <- if (!is.null(true_set)) nrow(true_set) else 0L
  n_sel_here  <- if (!is.null(selected)) nrow(selected) else 0L
  if (n_sel_here == 0L && n_true_here > 0L) {
    warning(sprintf(
      paste0("[vae|covsel] extracted 0 promoted beta.* terms but scenario %s ",
             "has %d true covariate(s). parFixed rownames = {%s}. ",
             "Selection scoring will be all-FN -- check the nlmixr2est parameter API."),
      opts$scenario, n_true_here,
      paste(utils::head(rownames(rec$parFixed), 20L), collapse = ", ")),
      call. = FALSE)
  }
}

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
# PR #921 A/B arm tag: "on" (covSelectColinear=TRUE) vs "off" (pre-PR engine).
rec$colinear    <- if (isTRUE(.colinear_on)) "on" else "off"
rec$colinear_cut <- .colinear_cut

write_fit_sidecar(rec, out_rds, fit = fit)

# On a successful (non-error) write, remove any stale ERROR sidecar left by a
# previous failed attempt at this same cell, so the directory reflects reality.
if (!is.null(fit) && file.exists(err_txt)) unlink(err_txt)

# ---- Console report --------------------------------------------------------
message(sprintf("<<< [vae|covsel] done in %.1fs  objf=%s  converged=%s",
                runtime_sec, format(rec$objf %||% NA), format(rec$converged %||% NA)))

sel_pairs <- if (!is.null(selected)) {
  dplyr::select(selected, var, covar, shape)
} else {
  tibble::tibble(var = character(), covar = character(), shape = character())
}

# Shape-aware scoring (like runSCM / PsN): a relationship is TP only when the
# (var, covar) pair AND its functional shape agree with the truth.  A right pair
# on the wrong shape scores FP (wrong-shape term) + FN (unmet true term).
.norm_shape <- function(s) {
  s <- tolower(as.character(s)); s[grepl("^hockey", s)] <- "hockey"
  s[s %in% c("", "na")] <- NA_character_; s
}
.shape_match <- function(a, b) is.na(a) | is.na(b) | (a == b)

tru <- dplyr::transmute(true_set,
                        var = as.character(var), covar = as.character(covar),
                        shape = .norm_shape(shape), .ti = dplyr::row_number())
sel <- dplyr::transmute(sel_pairs,
                        var = as.character(var), covar = as.character(covar),
                        shape = .norm_shape(shape), .si = dplyr::row_number())
j <- dplyr::inner_join(sel, tru, by = c("var", "covar"),
                       suffix = c(".sel", ".tru"),
                       relationship = "many-to-many") |>
  dplyr::filter(.shape_match(shape.sel, shape.tru))
hit_sel <- unique(j$.si); hit_tru <- unique(j$.ti)
cmp <- dplyr::bind_rows(
  tru |> dplyr::filter(.ti %in% hit_tru) |>
    dplyr::transmute(var, covar, shape, verdict = "TP (correct)"),
  tru |> dplyr::filter(!.ti %in% hit_tru) |>
    dplyr::transmute(var, covar, shape, verdict = "FN (missed)"),
  sel |> dplyr::filter(!.si %in% hit_sel) |>
    dplyr::transmute(var, covar, shape, verdict = "FP (spurious)")
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
