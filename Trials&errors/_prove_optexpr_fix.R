# ---- prove the optExpression=FALSE fix for the ifoceif ODE build hang -------
# Times three configs on scn16/ds1/N80, each wrapped in a hard time cap so a
# deadlocked build aborts instead of hanging for hours. The point is to show:
#   (1) plain focei (no analytic gradient)      -> baseline, always completes
#   (2) ifoceif default (optExpression=TRUE)    -> the hang (capped, expected FAIL)
#   (3) ifoceif optExpression=FALSE             -> the fix (should COMPLETE)
#
# Run:  source("script/_prove_optexpr_fix.R")
suppressPackageStartupMessages({ library(nlmixr2); library(rxode2) })

sd <- "script"
source(file.path(sd, "refit_helpers.R"),      chdir = FALSE)
source(file.path(sd, "true_model_factory.R"), chdir = FALSE)
source(file.path(sd, "estimator_factory.R"),  chdir = FALSE)

N <- 80L; scn <- 16L; ds <- 1L; bnd <- "narrow"
CAP_SEC <- 600  # hard per-fit wall-clock cap (10 min); a hung build aborts here

master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scn))
sim_slice <- load_scenario_dataset(master_rds, scn, ds)
ds_nm <- to_nm_dataset(sim_slice) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))
true_mod <- make_true_model(scn, boundary = bnd, structure = "ode")

# --- helper: run one fit under a hard time cap, capture time + status/objf ---
timed_fit <- function(label, est, ctrl) {
  cat(sprintf("\n===== %-28s =====\n", label))
  t0 <- Sys.time()
  res <- tryCatch({
    setTimeLimit(elapsed = CAP_SEC, transient = TRUE)
    on.exit(setTimeLimit(), add = TRUE)  # clear the cap afterwards
    fit <- nlmixr2(true_mod, ds_nm, est = est, control = ctrl)
    list(ok = TRUE, objf = tryCatch(fit$objf, error = function(e) NA_real_))
  }, error = function(e) list(ok = FALSE, objf = NA_real_, msg = conditionMessage(e)))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (isTRUE(res$ok)) {
    cat(sprintf("  COMPLETED in %.1f s   objf=%.3f\n", dt, res$objf))
  } else {
    cat(sprintf("  ABORTED after %.1f s  (%s)\n", dt,
                if (dt >= CAP_SEC - 5) "HIT TIME CAP -> hang confirmed" else res$msg))
  }
  data.frame(label = label, seconds = round(dt, 1),
             completed = res$ok, objf = round(res$objf, 3))
}

# (1) baseline: plain focei, numeric gradient, derivative-free outer
ctrl_focei <- foceiControl(
  sigdig = 4, outerOpt = "bobyqa", print = 0, calcTables = FALSE,
  covMethod = "", rxControl = rxControl(atol = 1e-8, rtol = 1e-6)
)

# (2) ifoceif DEFAULT: analytic gradient + parallel CSE (the hang). Capped.
ctrl_ifoceif_default <- foceiControl(
  interaction = TRUE, muModel = "irls", sigdig = 5, outerOpt = "lbfgsb3c",
  print = 0, calcTables = FALSE, covMethod = "",
  optExpression = TRUE,                       # <-- the parallel CSE that deadlocks
  rxControl = rxControl(atol = 1e-8, rtol = 1e-6)
)

# (3) ifoceif FIX: same, but optExpression=FALSE -> no daemons, no deadlock
ctrl_ifoceif_fix <- foceiControl(
  interaction = TRUE, muModel = "irls", sigdig = 5, outerOpt = "lbfgsb3c",
  print = 0, calcTables = FALSE, covMethod = "",
  optExpression = FALSE,                      # <-- the fix
  fallbackFD    = TRUE,                       # guard: degrade if sens eqns fail
  rxControl = rxControl(atol = 1e-8, rtol = 1e-6)
)

summary_tbl <- rbind(
  timed_fit("focei_bobyqa (baseline)",       "focei",   ctrl_focei),
  timed_fit("ifoceif default (optExpr=T)",   "ifoceif", ctrl_ifoceif_default),
  timed_fit("ifoceif FIX (optExpr=F)",       "ifoceif", ctrl_ifoceif_fix)
)

cat("\n================ SUMMARY ================\n")
print(summary_tbl, row.names = FALSE)
cat(sprintf("\n(time cap per fit: %d s)\n", CAP_SEC))
