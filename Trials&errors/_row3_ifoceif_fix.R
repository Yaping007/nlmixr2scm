# ---- run ONLY the ifoceif optExpression=FALSE fix (row 3, standalone) -------
# Use this after manually stopping the hung default run and RESTARTING R
# (restart clears any orphaned CSE daemon processes from the hung build).
#
# Run:  source("script/_row3_ifoceif_fix.R")
suppressPackageStartupMessages({ library(nlmixr2); library(rxode2) })

## Belt-and-suspenders: forbid any parallel symbolic/expression daemons at the
## global level, so even if optExpression were on nothing can deadlock.
options(rxode2.cores = 1L)
try(rxode2::setRxThreads(1L), silent = TRUE)
Sys.setenv(RXODE2_PRUNE_CORES = "1")

sd <- "script"
source(file.path(sd, "refit_helpers.R"),      chdir = FALSE)
source(file.path(sd, "true_model_factory.R"), chdir = FALSE)
source(file.path(sd, "estimator_factory.R"),  chdir = FALSE)

N <- 80L; scn <- 16L; ds <- 1L; bnd <- "narrow"

master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scn))
sim_slice <- load_scenario_dataset(master_rds, scn, ds)
ds_nm <- to_nm_dataset(sim_slice) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))
true_mod <- make_true_model(scn, boundary = bnd, structure = "ode")

ctrl_fix <- foceiControl(
  interaction   = TRUE,
  muModel       = "irls",
  sigdig        = 5,
  outerOpt      = "lbfgsb3c",
  print         = 0,
  calcTables    = FALSE,
  covMethod     = "",
  optExpression = FALSE,     # <-- the fix: no CSE, no parallel daemons
  fallbackFD    = TRUE,      # guard: degrade to FD if sens eqns fail
  rxControl     = rxControl(atol = 1e-8, rtol = 1e-6, cores = 1L)
)

cat("===== ifoceif optExpression=FALSE (single-core build) =====\n")
cat("  (est='focei' + interaction=TRUE + muModel='irls'  ==  ifoceif)\n")
t0  <- Sys.time()
fit <- nlmixr2(true_mod, ds_nm, est = "focei", control = ctrl_fix)
dt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("\nCOMPLETED in %.1f s   objf=%.3f\n", dt,
            tryCatch(fit$objf, error = function(e) NA_real_)))
print(fit)
