library(testthat)
library(nlmixr2scm)
devtools::load_all("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
library(nlmixr2utils)
library(nlmixr2)
library(rxode2)
library(tidyverse)

out_dir_v2 <- "simulated_virtual_dataset_eta_filtered"
stage1_dir <- file.path(out_dir_v2, "stage1_smoke_scn09_ds01")
ds01 <- readRDS(file.path(stage1_dir, "nm_scn09_ds01.rds"))

scm_focei <- nlmixr2est::foceiControl(
  sigdig     = 4,
  outerOpt   = "bobyqa",
  print      = 0,
  calcTables = FALSE,     # SCM doesn't need IPRED/CWRES tables
  covMethod  = ""         # SCM doesn't need cov matrix for LRT
)


base_2cmt_oral <- function() {
  ini({
    lTVCL <- log(0.6)
    lTVQ  <- log(1.8)
    lTVVc <- log(20)
    lTVVp <- log(80)
    ## KA fixed: the design's first non-zero sample (~7 h) is well past
    ## absorption (KA = 0.7 /h, absorption t1/2 ~ 1 h), so KA cannot be
    ## informatively estimated from the data.  Estimating it inflates
    ## the cov-matrix condition number by ~20x (raw CN ~35000 -> ~1800)
    ## without changing OFV or the SCM decisions.  Fixing matches the
    ## warfarin / Khandelwal convention.
    lTVKA <- fix(log(0.7))

    eta.cl + eta.vc ~ c(
      0.1,
      0.02,
      0.1
    )

    prop.err <- 0.1
  })
  model({
    cl <- exp(lTVCL + eta.cl)
    vc <- exp(lTVVc + eta.vc)
    q  <- exp(lTVQ)
    vp <- exp(lTVVp)
    ka <- exp(lTVKA)

    k10 <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)      = -ka * depot
    d/dt(central)    =  ka * depot - k10 * central - k12 * central + k21 * peripheral
    d/dt(peripheral) =  k12 * central - k21 * peripheral

    cp = central / vc
    cp ~ prop(prop.err)
  })
}
t_fit_base <- system.time(
  fit_base <- nlmixr2(base_2cmt_oral, ds01,
                      est = "focei", control = scm_focei)
)

candidate_pairs_test <- list(
  list(var = "cl", covar = "BW",   shapes = "power"),
  list(var = "cl", covar = "CrCL", shapes = "power"),
 list(var = "vc", covar = "BW",   shapes = "power"),
  list(var = "vc", covar = "CrCL", shapes = "power")
)

## ---- Wrapper: runSCM with wall-clock progress tracking -----------------
runSCM_traced <- function(label, ...) {
  dots      <- list(...)
  n_pairs   <- if (!is.null(dots$pairsVec))   length(dots$pairsVec) else NA
  n_workers <- if (!is.null(dots$workers))    dots$workers          else 1L
  search    <- if (!is.null(dots$searchType)) dots$searchType       else "scm"

  cat(sprintf(
    "\n=== [%s] %s starting %-8s | %d candidate(s), %d worker(s) ===\n",
    label, format(Sys.time(), "%H:%M:%S"), search, n_pairs, n_workers
  ), file = stderr())

  t0  <- Sys.time()
  res <- withCallingHandlers(
    do.call(nlmixr2scm::runSCM, dots),
    message = function(m) {
      txt <- conditionMessage(m)
      ## Tag only step/candidate boundary messages; let cli rules pass through.
      if (grepl("step\\s+\\d+,\\s*candidate|forward search|backward search",
                txt, ignore.case = TRUE)) {
        elapsed_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
        cat(sprintf("[%s | +%5.1f min] ",
                    format(Sys.time(), "%H:%M:%S"), elapsed_min),
            file = stderr())
      }
      ## Don't muffle -- let the original cli message print as usual.
    }
  )
  elapsed_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  cat(sprintf(
    "=== [%s] %s DONE | elapsed %.1f min (%.0f s) ===\n\n",
    label, format(Sys.time(), "%H:%M:%S"),
    elapsed_s / 60, elapsed_s
  ), file = stderr())

  attr(res, "elapsed_s") <- elapsed_s
  res
}

res_full <- runSCM_traced(
  label       = "full_scm",
  fit         = fit_base,
  pairsVec    = candidate_pairs_test,
  catvarsVec  = c("SEX", "RACE"),
  searchType  = "scm",
  control     = scm_focei,
  saveModels  = FALSE,
  workers     = 3L,
  print       = 100,
  maxRetries = 0L
)
