# ============================================================================
# _smoke_drivers_ode.R  --  Phase 3 end-to-end driver smoke (ODE)
# ----------------------------------------------------------------------------
# Reproduces the 7-cell diagnostic comparison from the old linCmt run, but on
# the NEW ODE infrastructure (structure="ode", schema 2.1, estimator-aware
# diagnostics, direct foceif/ifoceif dispatch, optExpression=TRUE). Runs the
# REAL driver bench_refit_one() for each cell at scn16/ds1/N80, then asserts
# every written record is schema-2.1 with its three sidecars, and prints a
# diagnostic-comparison table in the same shape as the linCmt figure.
#
# Combo (FOCEi family only -- mirrors the linCmt table without SAEM):
#   focei      x bobyqa
#   focei      x nlminb
#   focei      x lbfgsb3c
#   foceif     x nlminb
#   foceif     x lbfgsb3c
#   irlsfoceif x lbfgsb3c
#
# Usage (in the project R console, dev stack via .Rprofile):
#   source("script/_smoke_drivers_ode.R")
# ============================================================================
suppressPackageStartupMessages({ library(nlmixr2); library(rxode2) })

sd <- "script"
source(file.path(sd, "refit_helpers.R"),          chdir = FALSE)
source(file.path(sd, "true_model_factory.R"),     chdir = FALSE)
source(file.path(sd, "estimator_factory.R"),      chdir = FALSE)
source(file.path(sd, "output_schema.R"),          chdir = FALSE)
source(file.path(sd, "bench_refit_estimators.R"), chdir = FALSE)

N <- 80L; scn <- 16L; ds <- 1L; bnd <- "narrow"

cells <- tibble::tribble(
  ~estimator,     ~outer_opt,
  "focei",        "bobyqa",
  "focei",        "nlminb",
  "focei",        "lbfgsb3c",
  "foceif",       "nlminb",
  "foceif",       "lbfgsb3c",
  "irlsfoceif",   "lbfgsb3c"
)

master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scn))
out_root <- file.path("output", "smoke_ode", sprintf("scn%02d", scn))
cat("master:", master_rds, "exists:", file.exists(master_rds), "\n")
cat("out_root:", out_root, "\n\n")

# ---- run every cell through the real driver --------------------------------
rows <- vector("list", nrow(cells))
for (i in seq_len(nrow(cells))) {
  est <- cells$estimator[i]; opt <- cells$outer_opt[i]
  cell <- if (is.na(opt)) paste0(est, "_NA") else paste(est, opt, sep = "_")
  cat(sprintf("=== [%d/%d] %s ===\n", i, nrow(cells), cell))

  out_dir <- file.path(out_root, cell)
  t0 <- Sys.time()
  res <- tryCatch(
    bench_refit_one(N = N, scenario_id = scn, dataset_id = ds,
                    estimator = est, outer_opt = opt, boundary = bnd,
                    master_rds = master_rds, out_dir = out_dir,
                    overwrite = TRUE),
    error = function(e) { message("  cell failed: ", conditionMessage(e)); NULL })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  rds  <- file.path(out_dir, sprintf("res_ds%03d.rds", ds))
  fit_rds <- sub("\\.rds$", ".fit.rds", rds)
  meta <- sub("\\.rds$", ".meta.json", rds)

  d  <- if (!is.null(res)) res$diag_t3 else NULL
  rows[[i]] <- tibble::tibble(
    cell            = cell,
    runtime_s       = round(wall, 0),
    schema_ok       = !is.null(res) && identical(res$schema_version, "2.1"),
    model_type      = res$model_type %||% NA_character_,
    final_est_gone  = !is.null(res) && is.null(res$final_est),
    has_rds         = file.exists(rds),
    has_fit         = file.exists(fit_rds),
    has_meta        = file.exists(meta),
    refit_est       = res$refit_estimator %||% NA_character_,
    conv_code       = res$diag$convergence_code %||% NA_integer_,
    ofv             = round(res$objf %||% NA_real_, 2),
    cond_num        = round(res$cond_num_cor %||% NA_real_, 1),
    cn_source       = res$diag$cond_num_cor_source %||% NA_character_,
    est_family      = d$est_family %||% NA_character_,
    cov_source      = d$cov_source %||% NA_character_,
    rnd_err         = d$rnd_err,
    zero_grad       = d$zero_grad,
    se_source       = attr(res$parFixed, "se_source") %||% "nlmixr2",
    status          = res$status %||% "error"
  )
}
smoke <- dplyr::bind_rows(rows)

# ---- assertions ------------------------------------------------------------
cat("\n================ SCHEMA-2.1 ASSERTIONS ================\n")
asrt <- smoke |>
  dplyr::transmute(
    cell,
    schema_ok, final_est_gone,
    sidecars = has_rds & has_fit & has_meta,
    ok = status == "ok"
  )
print(as.data.frame(asrt))
all_pass <- all(smoke$schema_ok & smoke$final_est_gone &
                smoke$has_rds & smoke$has_fit & smoke$has_meta &
                smoke$status == "ok")
cat("\nALL CELLS PASS:", all_pass, "\n")

# ---- diagnostic comparison (linCmt-figure shape, now ODE) ------------------
cat("\n========= Diagnostic comparison -- scn16/ds01/N80 (ODE) =========\n")
diag_tbl <- smoke |>
  dplyr::transmute(
    Cell = cell,
    runtime = paste0(runtime_s, "s"),
    conv_code, OFV = ofv, cond_num, CN_source = cn_source,
    est_family, cov_source, rnd_err, zero_grad
  )
print(as.data.frame(diag_tbl), row.names = FALSE)

invisible(smoke)
