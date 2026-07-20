# ---- try est="ifoceif" DIRECTLY --------------------------------------------
# Now that the console loads nlmixr2est 6.2.0 (scratch lib), the ifoceif alias
# dispatches (getValidNlmixrControl(., "ifoceif") -> class ifoceiControl). This
# script calls est="ifoceif" directly and compares the result to our proven
# est="focei" + muModel="irls" recipe, to confirm they are numerically the same
# algorithm. optExpression=FALSE is still REQUIRED to avoid the parallel-CSE
# daemon deadlock on the ODE analytic-gradient build.
suppressPackageStartupMessages({ library(nlmixr2); library(rxode2) })

sd <- "script"
source(file.path(sd, "refit_helpers.R"),      chdir = FALSE)
source(file.path(sd, "true_model_factory.R"), chdir = FALSE)
source(file.path(sd, "estimator_factory.R"),  chdir = FALSE)
source(file.path(sd, "output_schema.R"),      chdir = FALSE)

N <- 80L; scn <- 16L; ds <- 1L; bnd <- "narrow"; opt <- "lbfgsb3c"

master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scn))
cat("master:", master_rds, "exists:", file.exists(master_rds), "\n")

sim_slice <- load_scenario_dataset(master_rds, scn, ds)
ds_nm <- to_nm_dataset(sim_slice) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))

true_mod <- make_true_model(scn, boundary = bnd, structure = "ode")

# Build the SAME control we use for irlsfoceif, but hand it to est="ifoceif".
# make_est_control("irlsfoceif") already sets interaction=TRUE, muModel="irls",
# optExpression=FALSE, fallbackFD=TRUE. The ifoceif alias presets the first two
# itself; passing them again is harmless and keeps the deadlock guard.
cfg <- make_est_control("irlsfoceif", outer_opt = opt, tier = "final")

cat("\n>>> calling est = \"ifoceif\" directly ...\n")
t0  <- Sys.time()
fit_direct <- nlmixr2(true_mod, ds_nm, est = "ifoceif", control = cfg$ctrl)
t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("ifoceif direct done in %.1f s; objf=%.4f\n",
            t_fit, tryCatch(fit_direct$objf, error = function(e) NA)))

# ---- compare to the focei+irls recipe (from the earlier saved fit) ----------
prev_fit_rds <- "output/schema_validation/irlsfoceif_scn16_ds01_N80.fit.rds"
cat("\n================ ifoceif  vs  focei+irls ================\n")
cat(sprintf("%-14s %14s %14s %14s\n", "source", "est_dispatch", "objf", "cond_num_cor"))
cat(sprintf("%-14s %14s %14.4f %14.3f\n", "direct", "ifoceif",
            fit_direct$objf,
            tryCatch(kappa(cov2cor(fit_direct$cov), exact = TRUE),
                     error = function(e) NA)))
if (file.exists(prev_fit_rds)) {
  fp <- readRDS(prev_fit_rds)
  cat(sprintf("%-14s %14s %14.4f %14.3f\n", "recipe", "focei+irls",
              fp$objf,
              tryCatch(kappa(cov2cor(fp$cov), exact = TRUE),
                       error = function(e) NA)))
  cat(sprintf("\ndelta objf (direct - recipe): %.3e\n", fit_direct$objf - fp$objf))
  th_d <- fit_direct$theta; th_r <- fp$theta
  common <- intersect(names(th_d), names(th_r))
  cat("max |theta diff| :", max(abs(th_d[common] - th_r[common])), "\n")
} else {
  cat("(no saved focei+irls fit to compare;",
      "run _validate_schema_irlsfoceif.R first)\n")
}

cat("\n>>> fit header:\n"); print(fit_direct)
