# ============================================================================
# _test_saem_standalone.R  --  isolated SAEM tuning bench (ODE)
# ----------------------------------------------------------------------------
# Standalone probe (NOT part of the driver smoke) to find a SAEM setup that
# gives us the three things we care about -- condition number, convergence
# status, OFV -- as FAST as possible, WITHOUT corrupting the fit.
#
# What each variant tests (scn16 / ds1 / N80, ODE, single-core, TIGHT tol):
#
#   fast     : f-SAEM         nBurn=400 nEm=600 fast=TRUE
#   fast_lik : f-SAEM + logLik=TRUE (standalone SAEM-GQ OFV, no FOCEi refit)
#
# Baseline (standard SAEM 500/1000) is DROPPED here -- it took >30 min and we
# already treat it as the production reference. This probe only measures the
# f-SAEM variants at 400/600 (a safer count than the package default 200/300).
#
# Settings held CONSTANT across variants (from our trials & errors):
#   * atol=1e-8, rtol=1e-6   -- tight tol is NON-NEGOTIABLE. A loose-tol
#                               experiment moved OFV 44 units and flipped
#                               covariate decisions. Never relax.
#   * cores = 1L             -- HPCE thread accounting + reproducibility.
#   * print = 1L             -- OBSERVABLE. print=0 caused the "is it frozen?"
#                               panic (500+1000 silent iterations).
#   * covMethod = "linFim"   -- cheap linearized FIM -> condition number.
#   * covFull  = TRUE        -- full theta+Omega covariance.
#
# The point: MEASURE the f-SAEM speedup and confirm the OFV / condition number
# do not drift vs the 500/1000 baseline. We reject iteration cuts for STANDARD
# SAEM, but f-SAEM (Karimi/Lavielle/Moulines 2020) reaches the same target with
# a smarter MH kernel -- so shorter runs may be legitimate HERE. Verify, don't
# assume.
#
# Usage (project console, dev stack via .Rprofile -> nlmixr2est 6.2.0):
#   source("script/_test_saem_standalone.R")
# ============================================================================
suppressPackageStartupMessages({ library(nlmixr2); library(rxode2) })

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

sd <- "script"
source(file.path(sd, "true_model_factory.R"),     chdir = FALSE)
source(file.path(sd, "refit_helpers.R"),          chdir = FALSE)
source(file.path(sd, "bench_refit_estimators.R"), chdir = FALSE)  # to_nm_dataset()

N <- 80L; scn <- 16L; ds <- 1L; bnd <- "narrow"

master_rds <- file.path("Inputdataset", sprintf("sim_obs_N%d", N),
                        sprintf("sim_obs_scenario_%02d.rds", scn))
cat("master:", master_rds, "exists:", file.exists(master_rds), "\n")

# ---- build true model + dataset -------------------------------------------
true_mod  <- make_true_model(scn, boundary = bnd, structure = "ode")
sim_slice <- load_scenario_dataset(master_rds, scenario_id = scn, dataset_id = ds)
# NM-shape exactly as the driver does: to_nm_dataset() then drop the TIME==0
# placeholder observation row and the SCENARIO/DATASET bookkeeping columns.
ds_nm <- to_nm_dataset(sim_slice) |>
  dplyr::filter(!(EVID == 0L & TIME == 0)) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX), RACE = as.integer(RACE))
cat("dataset rows:", nrow(ds_nm), " subjects:", dplyr::n_distinct(ds_nm$ID), "\n\n")

# ---- condition number from a fit (native slot, else linFim fallback) ------
.cond_num <- function(fit) {
  cn <- tryCatch(as.numeric(fit$conditionNumberCor), error = function(e) NA_real_)
  src <- "native"
  if (!isTRUE(is.finite(cn))) {
    cn <- tryCatch({
      cm  <- fit$cov
      evc <- eigen(cov2cor(cm), only.values = TRUE)$values
      max(evc) / min(evc)
    }, error = function(e) NA_real_)
    src <- "fallback"
  }
  list(cond = cn, source = src)
}

# ---- SAFE objf accessor ----------------------------------------------------
# For a logLik=FALSE SAEM fit there is NO stored objective. Touching $objf
# triggers a LAZY on-the-fly computation: "Calculating -2LL by Gaussian
# quadrature" -> on this bounded-covariate ODE model the GQ FAILS -> it falls
# back to setOfv(obj,"focei") -> internal nlmixr2(fit,...) dispatch hits a
# `nlmixr2FitCoreSilent` object with no method -> UseMethod() error that
# aborts the whole script. `%||%` does NOT catch that (it only guards NULL/NA,
# and the error is thrown, not returned). Wrap it so a lazy-objf failure
# degrades to NA instead of killing the probe.
.safe_objf <- function(fit) {
  tryCatch(as.numeric(fit$objf), error = function(e) {
    message("  (objf unavailable: lazy GQ/FOCEi objective failed -- ",
            "expected for logLik=FALSE SAEM on this model)")
    NA_real_
  })
}

# ---- run one SAEM variant, time it, tabulate ------------------------------
run_variant <- function(label, ctrl) {
  cat(sprintf("=== %s ===\n", label))
  t0  <- Sys.time()
  fit <- tryCatch(nlmixr2(true_mod, ds_nm, est = "saem", control = ctrl),
                  error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL })
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(fit)) {
    return(tibble::tibble(variant = label, wall_s = round(wall, 1),
                          objf = NA_real_, cond_num = NA_real_,
                          cn_source = NA_character_, cov_method = NA_character_,
                          converged = NA))
  }
  # cond# + cov come from the SAEM cov step -- these are the real convergence
  # signal, independent of the (fragile) lazy objf. objf is best-effort only.
  cn      <- .cond_num(fit)
  has_cov <- tryCatch(!is.null(fit$cov), error = function(e) FALSE)
  tibble::tibble(
    variant    = label,
    wall_s     = round(wall, 1),
    objf       = round(.safe_objf(fit), 2),
    cond_num   = round(cn$cond, 2),
    cn_source  = cn$source,
    cov_method = tryCatch(paste(fit$env$covMethod, collapse = ","),
                          error = function(e) NA_character_),
    converged  = isTRUE(is.finite(cn$cond)) && has_cov
  )
}

rxc <- rxode2::rxControl(atol = 1e-8, rtol = 1e-6, cores = 1L)

# f-SAEM (fast=TRUE) crashed with `Mat::operator(): index out of bounds` at the
# first SA iteration when bounded-transform covariate thetas (rxBoundedTr.*)
# were present. Isolate the cause with a 3-way matrix:
#   std          : standard SAEM, bounded transform ON  (does plain SAEM work?)
#   fast_nobnd   : f-SAEM,        bounded transform OFF  (does dropping bnd fix fast?)
#   std_nobnd    : standard SAEM, bounded transform OFF  (control for the bnd change)
ctrl_std <- nlmixr2est::saemControl(
  nBurn = 200L, nEm = 300L, nmc = 3L, fast = FALSE, print = 50L,
  covMethod = "linFim", covFull = TRUE, logLik = FALSE,
  calcTables = TRUE, rxControl = rxc
)
ctrl_fast_nobnd <- nlmixr2est::saemControl(
  nBurn = 200L, nEm = 300L, nmc = 3L, fast = TRUE, print = 50L,
  boundedTransform = FALSE,
  covMethod = "linFim", covFull = TRUE, logLik = FALSE,
  calcTables = TRUE, rxControl = rxc
)
ctrl_std_nobnd <- nlmixr2est::saemControl(
  nBurn = 200L, nEm = 300L, nmc = 3L, fast = FALSE, print = 50L,
  boundedTransform = FALSE,
  covMethod = "linFim", covFull = TRUE, logLik = FALSE,
  calcTables = TRUE, rxControl = rxc
)

variants <- list(
  std        = ctrl_std,
  fast_nobnd = ctrl_fast_nobnd,
  std_nobnd  = ctrl_std_nobnd
)

rows <- lapply(names(variants), function(nm) run_variant(nm, variants[[nm]]))
res  <- dplyr::bind_rows(rows)

cat("\n================ SAEM VARIANT COMPARISON ================\n")
print(as.data.frame(res), row.names = FALSE)

# ---- verdict --------------------------------------------------------------
cat("\n--- interpretation (which cell survived?) ---\n")
ok <- function(v) {
  r <- res[res$variant == v, ]
  if (nrow(r) && isTRUE(r$converged))
    sprintf("OK  (%.0fs, cond#=%.1f, objf=%s)",
            r$wall_s, r$cond_num,
            if (is.finite(r$objf)) sprintf("%.2f", r$objf) else "NA (needs refit)")
  else "CRASHED / no cov"
}
cat("std        (standard SAEM, bnd ON) :", ok("std"),        "\n")
cat("fast_nobnd (f-SAEM,        bnd OFF):", ok("fast_nobnd"), "\n")
cat("std_nobnd  (standard SAEM, bnd OFF):", ok("std_nobnd"),  "\n")
cat("\nDiagnosis:\n",
    "* if std OK        -> plain SAEM is fine; the crash is the f-SAEM kernel.\n",
    "* if fast_nobnd OK -> f-SAEM works ONLY without bounded transform (bug is\n",
    "                      fast + rxBoundedTr covariate thetas).\n",
    "* if std OK but both nobnd crash -> boundedTransform=FALSE is the problem.\n",
    "For production: prefer whichever converges with cond# + OFV. If only 'std'\n",
    "survives, keep standard SAEM (fast=FALSE) -- the f-SAEM speedup is unusable\n",
    "on this bounded-covariate model until the upstream Armadillo bug is fixed.\n",
    sep = "")

invisible(res)
