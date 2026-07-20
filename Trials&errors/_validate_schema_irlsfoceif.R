# ---- validate schema 2.0 on a real irlsfoceif fit --------------------------
# Standalone: fits the ODE true model to scn16/ds1/N80 with irlsfoceif,
# assembles the schema-2.0 record, writes the sidecar, and prints the record
# shape so we can eyeball the output scheme before wiring it into the array.
suppressPackageStartupMessages({
  library(nlmixr2); library(rxode2)
})

sd <- "script"
source(file.path(sd, "refit_helpers.R"),      chdir = FALSE)
source(file.path(sd, "true_model_factory.R"), chdir = FALSE)
source(file.path(sd, "estimator_factory.R"),  chdir = FALSE)
source(file.path(sd, "output_schema.R"),      chdir = FALSE)

N <- 80L; scn <- 16L; ds <- 1L; bnd <- "narrow"; est <- "irlsfoceif"; opt <- "lbfgsb3c"

master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scn))
cat("master:", master_rds, "exists:", file.exists(master_rds), "\n")

sim_slice <- load_scenario_dataset(master_rds, scn, ds)
ds_nm <- to_nm_dataset(sim_slice) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))

true_mod <- make_true_model(scn, boundary = bnd, structure = "ode")
cat("model_type:", fit_model_type(true_mod), "\n")

cfg <- make_est_control(est, outer_opt = opt, tier = "final")

est_dispatch <- nlmixr_est_name(est)
cat("grid label:", est, " -> nlmixr2 est:", est_dispatch, "\n")

t0  <- Sys.time()
fit <- nlmixr2(true_mod, ds_nm, est = est_dispatch, control = cfg$ctrl)
t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("fit done in %.1f s; objf=%.3f\n", t_fit, tryCatch(fit$objf, error=function(e) NA)))

common <- assemble_common(fit, true_mod, scn, true_params,
                          runtime_sec = t_fit, status = "ok")
result <- c(list(sample_N = N, scenario_id = scn, dataset_id = ds,
                 estimator = est, outer_opt = opt, refit_estimator = est,
                 boundary = bnd, refit_runtime_sec = NA_real_,
                 timestamp = Sys.time()),
            common)

out_dir <- file.path("output", "schema_validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_rds <- file.path(out_dir, "irlsfoceif_scn16_ds01_N80.rds")
sc <- write_fit_sidecar(result, out_rds, fit = fit)

cat("\n================ SCHEMA-2.1 RECORD ================\n")
cat("top-level keys:\n"); cat(names(result), sep = ", "); cat("\n\n")
cat("schema_version:", result$schema_version, " model_type:", result$model_type, "\n")
cat("objf:", round(result$objf, 3), " converged:", result$converged,
    " cond_num_cor:", round(result$cond_num_cor, 2), "\n")
cat("status:", result$status, " diag$cov_ok:", result$diag$cov_ok, "\n")
cat("diag_t3 flags: min_suc=", result$diag_t3$min_suc,
    " est_bnd=", result$diag_t3$est_bnd,
    " cov_step=", result$diag_t3$cov_step,
    " rnd_err=", result$diag_t3$rnd_err,
    " zero_grad=", result$diag_t3$zero_grad,
    " phys_bnd=", result$diag_t3$phys_bnd, "\n\n")
cat("--- parFixed (SE source:", attr(result$parFixed, "se_source") %||% "nlmixr2", ") ---\n")
print(result$parFixed[, c("Estimate", "SE", "%RSE", "Back-transformed")])
cat("\n--- cov matrix present:", !is.null(result$cov),
    " dim:", paste(dim(result$cov), collapse = "x"), "---\n")
cat("\n--- rel_err (estimate vs true) ---\n")
print(as.data.frame(result$rel_err)[, c("parameter","true_value","estimate","rel_err_pct")])
cat("\nsidecars written:\n  rds :", sc$rds, "\n  fit :", sc$fit,
    "\n  meta:", sc$meta, "\n")
cat("\n--- .meta.json contents ---\n")
cat(readLines(sc$meta), sep = "\n"); cat("\n")
