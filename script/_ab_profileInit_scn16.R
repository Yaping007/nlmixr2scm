# ============================================================================
# _ab_profileInit_scn16.R  --  does the 1-D profile warm-start rescue SCM?
# ----------------------------------------------------------------------------
# A/B test of runSCM(profileInit=...) on the MOST COMPLEX scenario (scn16,
# all four true covariates), across THREE estimator cells -- the gradient
# optimizers the warm-start was designed to rescue, plus the analytic IRLS
# variant as a reference.
#
#   scn16 truth (Khandelwal 2019): BW->CL (0.75), CrCL->CL (0.5),
#                                  BW->Vc (1.0), SEX->Vc (log 1.5 ~ 0.405)
#   cells     : irlsfoceif + lbfgsb3c   (analytic gradient, IRLS mu-ref)
#               focei      + nlminb     (FD gradient, can stall at 0)
#               focei      + lbfgsb3c   (FD gradient, can stall at 0)
#   arm A     : profileInit = FALSE     (current production behaviour)
#   arm B     : profileInit = TRUE      (1-D FOCEi profile warm-start)
#
# Everything is PRODUCTION wiring reused verbatim from the estimator factory +
# helpers -- no new estimation logic:
#   make_est_control (screen/final tiers) / nlmixr_est_name  (estimator_factory)
#   make_true_model                                          (true_model_factory)
#   base_2cmt_oral_linCmt / runSCM_traced / to_nm_dataset    (scm_bench_helpers)
#   assemble_common / write_fit_sidecar / SCHEMA_VERSION     (output_schema)
#   package_scm_schema21                                     (scm_bench_helpers)
#
# Output (schema 2.1) per cell x arm:
#   output/ab_profileInit/scn16/<cell>/<arm>/res_ds001.rds       (record)
#   output/ab_profileInit/scn16/<cell>/<arm>/res_ds001.fit.rds   (refit winner)
#   output/ab_profileInit/scn16/<cell>/<arm>/res_ds001.meta.json (manifest)
#
# Usage (project console, dev stack via .Rprofile):
#   source("script/_ab_profileInit_scn16.R")
# ============================================================================
suppressPackageStartupMessages({ library(nlmixr2); library(rxode2); library(dplyr) })

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

sd <- "script"
source(file.path(sd, "refit_helpers.R"),          chdir = FALSE)
source(file.path(sd, "true_model_factory.R"),      chdir = FALSE)
source(file.path(sd, "bench_refit_estimators.R"),  chdir = FALSE)  # to_nm_dataset, DOSE_MG
source(file.path(sd, "output_schema.R"),           chdir = FALSE)  # assemble_common, sidecar
source(file.path(sd, "estimator_factory.R"),       chdir = FALSE)  # make_est_control, dispatch
source(file.path(sd, "scm_bench_helpers.R"),       chdir = FALSE)  # base model + packaging

# ---- parallelism (env-driven on HPCE; sensible local default) -------------
#   AB_WORKERS  = number of future multisession workers (SCM candidate fan-out)
#   AB_RXTHREADS= rxode2 ODE-solve threads per worker
# Keep AB_WORKERS x AB_RXTHREADS <= cores requested (#BSUB -n) to avoid oversub.
AB_WORKERS   <- suppressWarnings(as.integer(Sys.getenv("AB_WORKERS",   unset = "3")))
AB_RXTHREADS <- suppressWarnings(as.integer(Sys.getenv("AB_RXTHREADS", unset = "2")))
if (is.na(AB_WORKERS)   || AB_WORKERS   < 1L) AB_WORKERS   <- 3L
if (is.na(AB_RXTHREADS) || AB_RXTHREADS < 1L) AB_RXTHREADS <- 2L

# ---- fixed scenario coordinates -------------------------------------------
N <- 80L; scn <- 16L; ds <- 1L
master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scn))
cat("master:", master_rds, " exists:", file.exists(master_rds), "\n")

# ---- run-scope filters (NULL = full grid) ---------------------------------
# Override these before source() to focus a single cell/arm/structure, e.g.:
#   RUN_CELLS <- "irlsfoceif_lbfgsb3c"; RUN_ARMS <- "profileOn"
#   RUN_STRUCTURE <- "ode"   # "linCmt" (default) or "ode" (analytic gradient)
#
# On HPCE these are supplied via environment variables (AB_CELL / AB_ARM /
# AB_STRUCTURE), one (cell x arm x structure) per LSF array task, so the
# session-variable form need not be edited. Precedence: session var > env > default.
.env_or <- function(nm, default = NULL) {
  v <- Sys.getenv(nm, unset = "")
  if (nzchar(v)) v else default
}
if (!exists("RUN_CELLS"))     RUN_CELLS     <- .env_or("AB_CELL",      NULL)
if (!exists("RUN_ARMS"))      RUN_ARMS      <- .env_or("AB_ARM",       NULL)
if (!exists("RUN_STRUCTURE")) RUN_STRUCTURE <- .env_or("AB_STRUCTURE", "linCmt")
stopifnot(RUN_STRUCTURE %in% c("linCmt", "ode"))

# structure-tagged output dir so linCmt and ode runs never overwrite
out_root <- file.path("output", "ab_profileInit",
                      sprintf("scn%d_%s", scn, RUN_STRUCTURE))

# ---- true parameter table (reuse if already in session) -------------------
if (!exists("true_params")) {
  if (exists("true_params_long") && is.function(true_params_long)) {
    true_params <- true_params_long()
  } else {
    stop("true_params not found and true_params_long() unavailable; ",
         "source true_model_factory/PerformanceEvaluation first.", call. = FALSE)
  }
}

# ---- true model (carries structure + bounds attrs for assemble_common) ----
true_mod <- make_true_model(scn, boundary = "narrow", structure = RUN_STRUCTURE)

# base (covariate-free) model matching the chosen structure
base_model <- if (identical(RUN_STRUCTURE, "ode")) base_2cmt_oral_ode else base_2cmt_oral_linCmt

# ---- production scn16 candidate space (verbatim) --------------------------
scm16_vars    <- c("cl", "vc")
scm16_covars  <- c("BW", "CrCL", "BMI")
scm16_catvars <- c("SEX", "RACE")
scm16_shapes  <- c("power", "lin")

# ---- estimator cells (search estimator; refit uses cell's own tier=final) -
cells <- tibble::tribble(
  ~cell,                  ~estimator,    ~outer_opt,
  "irlsfoceif_lbfgsb3c",  "irlsfoceif",  "lbfgsb3c",
  "focei_nlminb",         "focei",       "nlminb",
  "focei_lbfgsb3c",       "focei",       "lbfgsb3c"
)

# ---- load + NM-shape dataset (production path) ----------------------------
sim_slice <- load_scenario_dataset(master_rds, scenario_id = scn, dataset_id = ds)
ds16_01 <- to_nm_dataset(sim_slice) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))
cat("dataset rows:", nrow(ds16_01), " subjects:", dplyr::n_distinct(ds16_01$ID), "\n")

# ---- one (cell, arm) forward-SCM run --------------------------------------
#   base fit (covariate-free, screen control) + forward SCM + tight-tol refit,
#   packaged into a schema-2.1 record and written to the arm's output dir.
run_cell_arm <- function(cell, estimator, outer_opt, profileInit) {
  arm    <- if (isTRUE(profileInit)) "profileOn" else "profileOff"
  lbl    <- sprintf("scn%d_%s_%s", scn, cell, arm)
  cat(sprintf("\n=== CELL/ARM: %s ===\n", lbl))

  screen_ctrl <- make_est_control(estimator, outer_opt, tier = "screen")$ctrl
  final_ctrl  <- make_est_control(estimator, outer_opt, tier = "final")$ctrl
  est_disp    <- nlmixr_est_name(estimator)

  # base (covariate-free) fit with the cell's screen control
  t_base <- Sys.time()
  fit_base <- nlmixr2(base_model, ds16_01,
                      est = est_disp, control = screen_ctrl)
  base_sec <- as.numeric(difftime(Sys.time(), t_base, units = "secs"))
  cat("  base OFV:", round(as.numeric(fit_base$objf), 3),
      sprintf("(%.0fs)\n", base_sec))

  # forward SCM (screen control) with/without profileInit warm-start
  res <- runSCM_traced(
    label       = lbl,
    data        = ds16_01,
    fit         = fit_base,
    varsVec     = scm16_vars,
    covarsVec   = scm16_covars,
    catvarsVec  = scm16_catvars,
    shapes      = scm16_shapes,
    searchType  = "forward",
    control     = screen_ctrl,
    profileInit = profileInit,
    saveModels  = FALSE,
    workers     = AB_WORKERS,     # profileInit runs INSIDE the parallel candidate map (.plap) -- safe with workers>1
    print       = 100,
    maxRetries  = 0L
  )
  scm_sec <- attr(res, "elapsed_s")

  out_rds <- file.path(out_root, cell, arm, "res_ds001.rds")
  package_scm_schema21(
    scm_res         = res,
    true_mod        = true_mod,
    scenario_id     = scn,
    true_params     = true_params,
    estimator       = estimator,
    outer_opt       = outer_opt,
    refit_ctrl      = final_ctrl,
    refit_estimator = estimator,      # refit with the cell's OWN final settings
    runtimes        = list(base_sec = base_sec, scm_sec = scm_sec),
    identity        = list(sample_N = N, dataset_id = ds, boundary = "narrow"),
    out_rds         = out_rds
  )
}

# ---- run all cells x both arms (honouring RUN_CELLS / RUN_ARMS filters) ----
run_cells <- if (is.null(RUN_CELLS)) cells$cell else intersect(cells$cell, RUN_CELLS)
run_arms  <- if (is.null(RUN_ARMS))  c("profileOff", "profileOn") else RUN_ARMS
cat(sprintf("\n[scope] cells: %s | arms: %s\n",
            paste(run_cells, collapse = ", "), paste(run_arms, collapse = ", ")))

old_rxt <- rxode2::getRxThreads()
rxode2::setRxThreads(AB_RXTHREADS)  # AB_WORKERS x AB_RXTHREADS -- keep <= requested cores
on.exit(rxode2::setRxThreads(old_rxt), add = TRUE)

results <- list()
for (i in seq_len(nrow(cells))) {
  if (!(cells$cell[i] %in% run_cells)) next
  for (arm in run_arms) {
    pinit <- identical(arm, "profileOn")
    key   <- sprintf("%s.%s", cells$cell[i], arm)
    results[[key]] <- run_cell_arm(cells$cell[i], cells$estimator[i],
                                   cells$outer_opt[i], pinit)
  }
}

# ---- side-by-side comparison ----------------------------------------------
n_true_recovered <- function(rec) {
  s <- rec$scm$selected
  if (is.null(s) || nrow(s) == 0L) return(0L)
  key       <- paste(tolower(s$var), toupper(s$covar), sep = "~")
  true_keys <- c("cl~BW", "cl~CRCL", "vc~BW", "vc~SEX")
  sum(toupper(key) %in% toupper(true_keys))
}

sel_tbl <- function(rec, key) {
  s <- rec$scm$selected
  if (is.null(s) || nrow(s) == 0L)
    return(tibble::tibble(cell_arm = key, theta_name = NA_character_,
                          var = NA, covar = NA, shape = NA, estimate = NA_real_))
  dplyr::mutate(s, cell_arm = key) |>
    dplyr::select(cell_arm, theta_name, var, covar, shape, estimate)
}

cat("\n================ SELECTED COVARIATES (final model) ================\n")
selected_cmp <- dplyr::bind_rows(
  Map(sel_tbl, results, names(results))
)
print(as.data.frame(selected_cmp), row.names = FALSE)

# truth for scn16 (nonzero relations only)
truth16 <- true_params |>
  dplyr::filter(scenario == scn,
                parameter %in% c("CLBW", "CLcrCL", "VcBW", "VcSEX"),
                true_value != 0) |>
  dplyr::select(parameter, true_value)
cat("\n---- truth (scn16 nonzero covariate effects) ----\n")
print(as.data.frame(truth16), row.names = FALSE)

cat("\n================ SUMMARY ================\n")
summary_tbl <- tibble::tibble(
  cell_arm      = names(results),
  schema        = vapply(results, function(r) r$schema_version %||% NA_character_, character(1)),
  n_selected    = vapply(results, function(r) nrow(r$scm$selected %||% data.frame()), integer(1)),
  n_true_of_4   = vapply(results, n_true_recovered, integer(1)),
  final_OFV     = vapply(results, function(r) round(r$objf %||% NA_real_, 2), numeric(1)),
  cond_num      = vapply(results, function(r) round(r$cond_num_cor %||% NA_real_, 1), numeric(1)),
  cov_done      = vapply(results, function(r) isTRUE(r$scm$cov_done), logical(1)),
  base_s        = vapply(results, function(r) round(r$runtime$base_sec %||% NA_real_, 0), numeric(1)),
  scm_s         = vapply(results, function(r) round(r$runtime$scm_sec %||% NA_real_, 0), numeric(1)),
  refit_s       = vapply(results, function(r) round(r$runtime$refit_sec %||% NA_real_, 0), numeric(1))
)
print(as.data.frame(summary_tbl), row.names = FALSE)

# ---- sidecar sanity assertions --------------------------------------------
cat("\n================ SIDECAR CHECK ================\n")
for (i in seq_len(nrow(cells))) {
  if (!(cells$cell[i] %in% run_cells)) next
  for (arm in run_arms) {
    stem <- file.path(out_root, cells$cell[i], arm, "res_ds001")
    ok <- c(rds  = file.exists(paste0(stem, ".rds")),
            fit  = file.exists(paste0(stem, ".fit.rds")),
            meta = file.exists(paste0(stem, ".meta.json")))
    cat(sprintf("  %-24s %-11s  rds:%s fit:%s meta:%s\n",
                cells$cell[i], arm, ok["rds"], ok["fit"], ok["meta"]))
  }
}
stopifnot(all(summary_tbl$schema == SCHEMA_VERSION))

cat("\nInterpretation:\n")
cat("  * n_true_of_4 is the SCM 'power' proxy on this single dataset.\n")
cat("  * If profileOn recovers more true relations (or gets closer estimates)\n")
cat("    than profileOff for a gradient cell, the 1-D profile warm-start\n")
cat("    rescues the optimizer from stalling at the zero-effect point.\n")

invisible(list(results = results,
               selected = selected_cmp, summary = summary_tbl))
