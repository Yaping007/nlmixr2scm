# ============================================================================
# _smoke_drivers_lincmt.R  --  linCmt() twin of the Phase-3 driver smoke
# ----------------------------------------------------------------------------
# Recreates the ORIGINAL linCmt()-based true-model fit (the run behind the
# 7-row "Diagnostic comparison" figure) but on the CURRENT reporting stack:
#   * structure = "linCmt"           (closed-form, NOT ODE)
#   * schema 2.1 record + 3 sidecars (via assemble_common / write_fit_sidecar)
#   * estimator-aware diagnostics    (diagnose_fit_table3)
#   * direct foceif/ifoceif dispatch, optExpression=TRUE
#
# The estimator / optimizer settings are IDENTICAL to the ODE smoke -- the ONLY
# difference is the model structure passed to make_true_model(). This is the
# apples-to-apples linCmt vs ODE comparison.
#
# Combo (mirrors the figure EXACTLY -- 6 FOCEi cells + SAEM):
#   focei      x bobyqa
#   focei      x nlminb
#   focei      x lbfgsb3c
#   foceif     x nlminb
#   foceif     x lbfgsb3c
#   irlsfoceif x lbfgsb3c
#   saem       x NA        (post-SAEM focei+bobyqa cov refit, honest label)
#
# The raw nlmixr2 fit for every cell is KEPT: written to <cell>/res_dsDDD.fit.rds
# by write_fit_sidecar(). Reload any cell with readRDS(".../res_dsDDD.fit.rds").
#
# Usage (project console, dev stack via .Rprofile -> nlmixr2est 6.2.0):
#   source("script/_smoke_drivers_lincmt.R")
# ============================================================================
suppressPackageStartupMessages({ library(nlmixr2); library(rxode2) })

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

sd <- "script"
source(file.path(sd, "refit_helpers.R"),          chdir = FALSE)
source(file.path(sd, "true_model_factory.R"),     chdir = FALSE)
source(file.path(sd, "estimator_factory.R"),      chdir = FALSE)
source(file.path(sd, "output_schema.R"),          chdir = FALSE)
source(file.path(sd, "bench_refit_estimators.R"), chdir = FALSE)  # to_nm_dataset()

N <- 80L; scn <- 16L; ds <- 1L; bnd <- "narrow"
STRUCTURE <- "linCmt"

cells <- tibble::tribble(
  ~estimator,     ~outer_opt,
  "focei",        "bobyqa",
  "focei",        "nlminb",
  "focei",        "lbfgsb3c",
  "foceif",       "nlminb",
  "foceif",       "lbfgsb3c",
  "irlsfoceif",   "lbfgsb3c",
  "saem",         NA_character_
)

master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scn))
out_root <- file.path("output", "smoke_lincmt", sprintf("scn%02d", scn))
cat("master:", master_rds, "exists:", file.exists(master_rds), "\n")
cat("structure:", STRUCTURE, " out_root:", out_root, "\n\n")

# ---- load + NM-shape the dataset once (shared by every cell) ---------------
sim_slice <- load_scenario_dataset(master_rds, scenario_id = scn, dataset_id = ds)
ds_nm <- to_nm_dataset(sim_slice) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))
cat("dataset rows:", nrow(ds_nm), " subjects:", dplyr::n_distinct(ds_nm$ID), "\n\n")

# ---- linCmt() true model (same boundary as ODE, structure = "linCmt") ------
true_mod <- make_true_model(scn, boundary = bnd, structure = STRUCTURE)

# ----------------------------------------------------------------------------
# fit_lincmt_one(): the linCmt twin of bench_refit_one(). Same estimator /
# optimizer settings (make_est_control), same SAEM->focei+bobyqa cov refit,
# same schema-2.1 assembly (assemble_common) + sidecars (write_fit_sidecar).
# The raw nlmixr2 fit is persisted to <cell>/res_dsDDD.fit.rds by
# write_fit_sidecar(); this returns just the slim schema-2.1 record.
# ----------------------------------------------------------------------------
fit_lincmt_one <- function(estimator, outer_opt, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_rds <- file.path(out_dir, sprintf("res_ds%03d.rds", ds))

  cfg          <- make_est_control(estimator, outer_opt = outer_opt, tier = "final")
  est_dispatch <- nlmixr_est_name(estimator)

  t0  <- Sys.time()
  fit <- nlmixr2(true_mod, ds_nm, est = est_dispatch, control = cfg$ctrl)
  t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  # Post-SAEM covariance refit (focei + bobyqa; honest "focei" label).
  refit_estimator <- estimator
  t_refit <- NA_real_
  if (isTRUE(cfg$needs_post_refit)) {
    cat("  ... post-SAEM focei cov refit (bobyqa)\n")
    t1 <- Sys.time()
    fit_refit <- tryCatch(
      nlmixr2(true_mod, ds_nm, est = nlmixr_est_name("focei"),
              control = make_saem_refit_control()),
      error = function(e) { message("post-refit failed: ", conditionMessage(e)); NULL })
    t_refit <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
    if (!is.null(fit_refit)) { fit <- fit_refit; refit_estimator <- "focei" }
  }

  common <- assemble_common(
    fit = fit, true_mod = true_mod, scenario_id = scn,
    true_params = true_params, runtime_sec = t_fit,
    status = "ok", estimator = refit_estimator)

  record <- c(list(
    sample_N          = N,
    scenario_id       = as.integer(scn),
    dataset_id        = as.integer(ds),
    estimator         = estimator,
    outer_opt         = if (is.na(outer_opt)) NA_character_ else as.character(outer_opt),
    refit_estimator   = refit_estimator,
    boundary          = bnd,
    refit_runtime_sec = t_refit,
    timestamp         = Sys.time()
  ), common)

  # KEEP the raw nlmixr2 fit alongside the slim schema record (on disk).
  write_fit_sidecar(record, out_rds, fit = fit)
  record
}

# ---- run every cell --------------------------------------------------------
rows <- vector("list", nrow(cells))
for (i in seq_len(nrow(cells))) {
  est <- cells$estimator[i]; opt <- cells$outer_opt[i]
  cell <- if (is.na(opt)) paste0(est, "_NA") else paste(est, opt, sep = "_")
  cat(sprintf("=== [%d/%d] %s (linCmt) ===\n", i, nrow(cells), cell))

  out_dir <- file.path(out_root, cell)
  t0 <- Sys.time()
  res <- tryCatch(fit_lincmt_one(est, opt, out_dir),
                  error = function(e) { message("  cell failed: ", conditionMessage(e)); NULL })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  rds     <- file.path(out_dir, sprintf("res_ds%03d.rds", ds))
  fit_rds <- sub("\\.rds$", ".fit.rds", rds)
  meta    <- sub("\\.rds$", ".meta.json", rds)

  d <- if (!is.null(res)) res$diag_t3 else NULL
  rows[[i]] <- tibble::tibble(
    cell           = cell,
    runtime_s      = round(wall, 0),
    schema_ok      = !is.null(res) && identical(res$schema_version, "2.1"),
    model_type     = res$model_type %||% NA_character_,
    final_est_gone = !is.null(res) && is.null(res$final_est),
    has_rds        = file.exists(rds),
    has_fit        = file.exists(fit_rds),
    has_meta       = file.exists(meta),
    refit_est      = res$refit_estimator %||% NA_character_,
    conv_code      = res$diag$convergence_code %||% NA_integer_,
    ofv            = round(res$objf %||% NA_real_, 2),
    cond_num       = round(res$cond_num_cor %||% NA_real_, 1),
    cn_source      = res$diag$cond_num_cor_source %||% NA_character_,
    msg            = res$diag$message %||% NA_character_,
    est_family     = d$est_family %||% NA_character_,
    cov_source     = d$cov_source %||% NA_character_,
    status         = res$status %||% "error"
  )
}
smoke <- dplyr::bind_rows(rows)

# ---- assertions ------------------------------------------------------------
cat("\n================ SCHEMA-2.1 ASSERTIONS (linCmt) ================\n")
asrt <- smoke |>
  dplyr::transmute(
    cell, model_type, schema_ok, final_est_gone,
    sidecars = has_rds & has_fit & has_meta,
    ok = status == "ok")
print(as.data.frame(asrt), row.names = FALSE)
all_pass <- all(smoke$schema_ok & smoke$final_est_gone &
                smoke$has_rds & smoke$has_fit & smoke$has_meta &
                smoke$status == "ok")
cat("\nALL CELLS PASS:", all_pass, "\n")

# ---- diagnostic comparison (recreates the figure exactly) ------------------
cat(sprintf("\n===== Diagnostic comparison -- scn%02d/ds%02d/N%d (linCmt) =====\n",
            scn, ds, N))
diag_tbl <- smoke |>
  dplyr::transmute(
    Cell      = cell,
    runtime   = paste0(runtime_s, "s"),
    conv_code,
    OFV       = ofv,
    cond_num,
    `CN source` = cn_source,
    msg)
print(as.data.frame(diag_tbl), row.names = FALSE)

invisible(smoke)
