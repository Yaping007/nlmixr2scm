# ============================================================================
# _test_fix1_irlsfoceif_scn16.R
# ----------------------------------------------------------------------------
# Focused A/B: does irlsfoceif / lbfgsb3c WORK on ODE candidate models when it
# is given a REAL warm start (profileInit-style bobyqa 1-D frozen profile, the
# same thing runSCM does), with Fix 1 (sensitivity-matched atolSens/rtolSens)
# turned ON vs OFF?
#
# This is the missing regime from the earlier cold-start trace: there every
# candidate started at theta=0 and the gradient optimizers stalled. runSCM does
# NOT start cold -- profileInit hands each candidate a bobyqa-profiled non-zero
# init. In THAT regime the only thing left to differentiate irlsfoceif is the
# analytic-outer-gradient -> FD fallback (fd_switch), which is exactly what Fix
# 1 targets.
#
# Design (scenario 16 / N80 / dataset 1, ODE base):
#   For each candidate covariate:
#     step A  warm start  : freeze all structural thetas at the base fit, free
#                           ONLY the new covariate theta, fit focei/bobyqa 1-D
#                           -> profiled init (== runSCM .profileCovInit).
#     step B  irlsfoceif  : fit the FULL candidate model, lbfgsb3c outer,
#                           started at the profiled init, TWICE:
#                             fix1 = FALSE  (loose sens tols, default)
#                             fix1 = TRUE   (sens-matched tols)
#   Compare: dOFV, theta_hat, MOVED/FROZEN, conv, fd_switch, runtime.
#
# PASS for "Fix 1 rescues irlsfoceif" =
#   fix1=TRUE gives fd_switch=FALSE AND a clean, MOVED, sensible-dOFV fit,
#   while fix1=FALSE shows fd_switch=TRUE / stall / Hessian poisoning.
#
# Usage:  source("script/_test_fix1_irlsfoceif_scn16.R")  (repo root)
# ============================================================================

suppressPackageStartupMessages({
  library(nlmixr2)
  library(nlmixr2est)
  library(rxode2)
  library(dplyr)
  library(tibble)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

OPTS <- list(N = 80L, scenario = 16L, dataset = 1L,
             input_root = "Inputdataset",
             structure = "ode",
             cov_bound = c(-5, 5),
             out_rds = "output/diag/fix1_irlsfoceif_scn16_ds01.rds")

.script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(fa)) return(normalizePath(dirname(fa[1]), mustWork = FALSE))
  "script"
})()
source(file.path(.script_dir, "scm_bench_helpers.R"))
source(file.path(.script_dir, "estimator_factory.R"))
if (!exists("DOSE_MG")) DOSE_MG <- 100

# ---- data ------------------------------------------------------------------
sim_path <- file.path(OPTS$input_root, sprintf("sim_obs_N%d", OPTS$N),
                      sprintf("sim_obs_scenario_%02d.rds", OPTS$scenario))
stopifnot(file.exists(sim_path))
ds_i <- to_nm_dataset(readRDS(sim_path)) |>
  dplyr::filter(DATASET == OPTS$dataset) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID), SEX = as.integer(SEX),
                RACE = as.integer(RACE))

# ---- candidate catalogue ---------------------------------------------------
candidates <- tibble::tribble(
  ~key,            ~target, ~theta_name,        ~term,             ~truth,
  "cl~BW.power",   "cl",    "cov_BW_power_cl",   "log(BW / 70)",    TRUE,
  "cl~CrCL.power", "cl",    "cov_CrCL_power_cl", "log(CrCL / 95)",  TRUE,
  "vc~BW.power",   "vc",    "cov_BW_power_vc",   "log(BW / 70)",    TRUE,
  "vc~SEX.lin",    "vc",    "cov_SEX_lin_vc",    "SEX",             TRUE,
  "cl~BMI.power",  "cl",    "cov_BMI_power_cl",  "log(BMI / 26)",   FALSE
)

# ---- candidate UI builder (ODE), structural thetas optionally frozen -------
# freeze_struct = TRUE  -> wrap all structural thetas + omega + prop in fix(),
#                          leaving ONLY the covariate theta free (the frozen
#                          1-D profile used for the warm start).
make_cand_ui <- function(target, theta_name, term, cov_init,
                         base_theta, base_omega, bound = OPTS$cov_bound,
                         freeze_struct = FALSE) {
  gv <- function(nm, d) { v <- suppressWarnings(as.numeric(base_theta[nm]))
    if (length(v) == 0L || is.na(v)) d else v }
  lTVCL <- gv("lTVCL", log(0.6)); lTVQ <- gv("lTVQ", log(1.8))
  lTVVc <- gv("lTVVc", log(20));  lTVVp <- gv("lTVVp", log(80))
  prop  <- gv("prop.err", 0.1)
  o11 <- base_omega["eta.cl", "eta.cl"]; o22 <- base_omega["eta.vc", "eta.vc"]
  o21 <- base_omega["eta.vc", "eta.cl"]
  term_cl <- if (target == "cl") sprintf(" + %s * %s", theta_name, term) else ""
  term_vc <- if (target == "vc") sprintf(" + %s * %s", theta_name, term) else ""
  w <- if (freeze_struct) function(x) sprintf("fix(%s)", x) else function(x) x
  ini_lines <- c(
    sprintf("lTVCL <- %s", w(sprintf("%.10g", lTVCL))),
    sprintf("lTVQ  <- %s", w(sprintf("%.10g", lTVQ))),
    sprintf("lTVVc <- %s", w(sprintf("%.10g", lTVVc))),
    sprintf("lTVVp <- %s", w(sprintf("%.10g", lTVVp))),
    "lTVKA <- fix(log(0.7))",
    sprintf("%s <- c(%.6g, %.10g, %.6g)", theta_name, bound[1], cov_init, bound[2]),
    if (freeze_struct)
      sprintf("eta.cl + eta.vc ~ fix(c(%.8g, %.8g, %.8g))", o11, o21, o22)
    else
      sprintf("eta.cl + eta.vc ~ c(%.8g, %.8g, %.8g)", o11, o21, o22),
    sprintf("prop.err <- %s", w(sprintf("%.10g", prop)))
  )
  model_lines <- c(
    sprintf("cl <- exp(lTVCL%s + eta.cl)", term_cl),
    sprintf("vc <- exp(lTVVc%s + eta.vc)", term_vc),
    "q  <- exp(lTVQ)", "vp <- exp(lTVVp)", "ka <- exp(lTVKA)",
    "d/dt(depot)   <- -ka * depot",
    "d/dt(central) <-  ka * depot - (cl/vc)*central - (q/vc)*central + (q/vp)*periph",
    "d/dt(periph)  <-  (q/vc)*central - (q/vp)*periph",
    "cp <- central/vc",
    "cp ~ prop(prop.err)"
  )
  fn_text <- paste0(
    "function() {\n  ini({\n    ", paste(ini_lines, collapse = "\n    "),
    "\n  })\n  model({\n    ", paste(model_lines, collapse = "\n    "),
    "\n  })\n}"
  )
  eval(parse(text = fn_text), envir = globalenv())
}

# ---- fit summariser --------------------------------------------------------
summ <- function(fit, base_objf, theta_name) {
  if (is.null(fit) || isTRUE(fit$.failed))
    return(list(objf = NA_real_, dOFV = NA_real_, theta_hat = NA_real_,
                conv = NA_integer_, msg = "FIT ERROR", cn = NA_real_,
                switch = NA))
  d <- diagnose_fit(fit)
  objf <- as.numeric(fit$objf)
  th <- suppressWarnings(as.numeric(fit$theta[theta_name]))
  ri <- tryCatch(fit$runInfo, error = function(e) NULL)
  list(objf = objf, dOFV = base_objf - objf,
       theta_hat = if (length(th)) th else NA_real_,
       conv = d$convergence_code, msg = d$message %||% "",
       cn = d$cond_num_cor,
       switch = if (is.null(ri)) NA else any(grepl("could not be solved", ri)))
}

# ---- base fit: irlsfoceif / lbfgsb3c ---------------------------------------
message("===== base fit: irlsfoceif / lbfgsb3c =====")
base_ctrl <- make_est_control("irlsfoceif", "lbfgsb3c", "screen", fix1 = TRUE)$ctrl
base_fit <- nlmixr2(base_2cmt_oral_ode, ds_i,
                    est = nlmixr_est_name("irlsfoceif"), control = base_ctrl)
base_objf  <- as.numeric(base_fit$objf)
base_theta <- base_fit$theta
base_omega <- base_fit$omega
message(sprintf("  base OFV = %.3f  (conv=%s)",
                base_objf, diagnose_fit(base_fit)$convergence_code))

# control used for the frozen bobyqa 1-D warm-start profile (== runSCM)
prof_ctrl <- make_est_control("focei", "bobyqa", "screen")$ctrl

rows <- list()
for (k in seq_len(nrow(candidates))) {
  cc <- candidates[k, ]
  message(sprintf("\n---- candidate: %s (truth=%s) ----", cc$key, cc$truth))

  # step A: profileInit-style frozen 1-D bobyqa warm start ------------------
  prof_ui <- make_cand_ui(cc$target, cc$theta_name, cc$term, cov_init = 0,
                          base_theta, base_omega, freeze_struct = TRUE)
  pf <- tryCatch(suppressWarnings(nlmixr2(prof_ui, ds_i, est = "focei",
                                          control = prof_ctrl)),
                 error = function(e) NULL)
  warm_init <- tryCatch(unname(pf$theta[[cc$theta_name]]),
                        error = function(e) NA_real_)
  if (is.null(warm_init) || !is.finite(warm_init)) warm_init <- 0.001
  message(sprintf("  profileInit (focei/bobyqa 1-D) -> %s = %.5g",
                  cc$theta_name, warm_init))

  # step B: irlsfoceif / lbfgsb3c from the warm start, Fix 1 OFF then ON -----
  for (fx in c(FALSE, TRUE)) {
    ctrl <- make_est_control("irlsfoceif", "lbfgsb3c", "screen", fix1 = fx)$ctrl
    ui <- make_cand_ui(cc$target, cc$theta_name, cc$term, cov_init = warm_init,
                       base_theta, base_omega, freeze_struct = FALSE)
    t0 <- proc.time()[["elapsed"]]
    f <- tryCatch(suppressWarnings(nlmixr2(ui, ds_i,
                                           est = nlmixr_est_name("irlsfoceif"),
                                           control = ctrl)),
                  error = function(e) list(.failed = TRUE))
    dt <- proc.time()[["elapsed"]] - t0
    s <- summ(f, base_objf, cc$theta_name)
    moved <- is.finite(s$theta_hat) && abs(s$theta_hat - warm_init) > 1e-4
    rows[[length(rows) + 1L]] <- tibble::tibble(
      candidate = cc$key, truth = cc$truth, fix1 = fx,
      warm_init = warm_init, theta_hat = s$theta_hat,
      moved = moved, dOFV = s$dOFV, accept = is.finite(s$dOFV) && s$dOFV > 3.84,
      conv = s$conv, fd_switch = s$switch, cn = s$cn,
      secs = dt, msg = s$msg
    )
    message(sprintf("  [fix1=%-5s] init=%.4g -> theta=%s (%s)  dOFV=%8.2f  conv=%s  fd_switch=%s  %.1fs  '%s'",
                    fx, warm_init, formatC(s$theta_hat, format = "g", digits = 4),
                    if (moved) "MOVED" else "FROZEN", s$dOFV, s$conv,
                    s$switch, dt, s$msg))
  }
}

res <- dplyr::bind_rows(rows)

# ---- report ----------------------------------------------------------------
cat("\n\n============ FIX 1 A/B : irlsfoceif / lbfgsb3c (WARM start) ============\n")
cat(sprintf("N=%d  scenario=%d  dataset=%d  structure=%s  base_OFV=%.3f\n\n",
            OPTS$N, OPTS$scenario, OPTS$dataset, OPTS$structure, base_objf))

cat("---- dOFV (LRT chi2_1 @0.05 = 3.84) : Fix1 OFF vs ON ----\n")
print(as.data.frame(
  res |>
    dplyr::mutate(fix1 = ifelse(fix1, "ON", "OFF")) |>
    dplyr::select(candidate, truth, fix1, dOFV) |>
    tidyr::pivot_wider(names_from = fix1, values_from = dOFV)
), digits = 4)

cat("\n---- analytic->FD fallback (fd_switch) : Fix1 OFF vs ON ----\n")
print(as.data.frame(
  res |>
    dplyr::mutate(fix1 = ifelse(fix1, "ON", "OFF")) |>
    dplyr::select(candidate, fix1, fd_switch) |>
    tidyr::pivot_wider(names_from = fix1, values_from = fd_switch)
), right = FALSE)

cat("\n---- theta / MOVED / conv / runtime ----\n")
print(as.data.frame(
  res |> dplyr::mutate(fix1 = ifelse(fix1, "ON", "OFF")) |>
    dplyr::select(candidate, fix1, warm_init, theta_hat, moved, conv, secs, cn)
), right = FALSE, digits = 4)

dir.create(dirname(OPTS$out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(res = res, opts = OPTS, base_objf = base_objf), OPTS$out_rds)
cat(sprintf("\nSaved -> %s\n", OPTS$out_rds))

cat("\n---- VERDICT ----\n")
verdict <- res |>
  dplyr::group_by(candidate, truth) |>
  dplyr::summarise(
    off_switch = fd_switch[!fix1], on_switch = fd_switch[fix1],
    off_dOFV = dOFV[!fix1], on_dOFV = dOFV[fix1],
    rescued = isTRUE(off_switch) && identical(on_switch, FALSE),
    .groups = "drop"
  )
print(as.data.frame(verdict), right = FALSE, digits = 4)
cat("\nrescued = Fix1 flipped fd_switch TRUE->FALSE for that candidate.\n")
