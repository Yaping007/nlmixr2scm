# ============================================================================
# diagnose_scm_gradient_trace.R
# ----------------------------------------------------------------------------
# One-dataset mechanism trace: WHY do gradient optimizers (nlminb, lbfgsb3c)
# collapse in SCM while bobyqa succeeds, even though ALL of them recover the
# true model fine in the refit benchmark?
#
# Hypothesis to confirm/refute:
#   The failure is NOT that a gradient optimizer cannot fit a covariate model.
#   It is that each SCM candidate is started COLD (new covariate theta = 0,
#   inner EBEs recomputed from scratch). On the gravelly FOCEi surface the
#   finite-difference gradient at theta=0 drowns in inner-ETA ripple, so the
#   optimizer stops early ("false convergence (8)") at a near-zero slope ->
#   the candidate's dOFV is ~0 -> LRT p-value ~1 -> the TRUE covariate is
#   rejected -> Power collapses.
#
# The decisive test (this script):
#   For a handful of SINGLE-covariate candidates added to the base model,
#   fit each candidate three ways and compare the LRT-relevant dOFV:
#     (1) bobyqa   , cold  (cov theta init = 0)              <- reference
#     (2) nlminb   , cold  (cov theta init = 0)              <- the suspect
#     (3) nlminb   , WARM  (cov theta init = bobyqa's est)   <- proximity fix
#   plus an optional (4) nlminb WARM + etaMat from base fit  <- inner-eta fix.
#
#   Prediction if the hypothesis is correct:
#     - TRUE covariates: bobyqa dOFV >> 3.84 with a sensible theta;
#       nlminb-cold dOFV ~ 0 (theta stuck ~0, msg "false convergence (8)");
#       nlminb-WARM dOFV recovers to ~ bobyqa's value.
#     This isolates PROXIMITY-AT-START, not optimizer capability, as the cause.
#
# Runtime: ~3-5 min (≈16 short fits). Scenario 16 / N80 / dataset 1.
#
# Usage:
#   Rscript script/diagnose_scm_gradient_trace.R
#   # or source() interactively in an R session started at repo root.
# ============================================================================

suppressPackageStartupMessages({
  library(nlmixr2)
  library(nlmixr2est)
  library(rxode2)
  library(dplyr)
  library(tibble)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Config ----------------------------------------------------------------
OPTS <- list(
  N            = 80L,
  scenario     = 16L,
  dataset      = 1L,
  input_root   = "Inputdataset",
  # (estimator, optimizer) cells to trace. focei/bobyqa is the reference and
  # also supplies the warm starting theta values reused by the gradient cells.
  cells        = list(
    c("focei",      "bobyqa"),
    c("focei",      "nlminb"),
    c("focei",      "lbfgsb3c"),
    c("foceif",     "nlminb"),
    c("foceif",     "lbfgsb3c"),
    c("irlsfoceif", "lbfgsb3c")
  ),
  warm_src     = "focei/bobyqa",          # cell whose thetas seed the warm arm
  cov_bound    = c(-5, 5),                # runSCM default candidate bounds
  # structure: "ode" uses the explicit d/dt() base model (unlocks the Almquist
  # analytic outer gradient -> the fd_switch mechanism is testable); "linCmt"
  # uses the analytic linCmt() base (fast is auto-disabled, switch is always
  # NA/FALSE). Override with env DIAG_STRUCTURE=ode|linCmt.
  structure    = Sys.getenv("DIAG_STRUCTURE", "ode"),
  out_rds      = "output/diag/scm_gradient_trace_N80_scn16_ds01.rds"
)

.script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(fa)) return(normalizePath(dirname(fa[1]), mustWork = FALSE))
  "script"
})()
source(file.path(.script_dir, "scm_bench_helpers.R"))   # base model + to_nm_dataset
source(file.path(.script_dir, "estimator_factory.R"))   # make_est_control

# to_nm_dataset() reads the fixed study dose from the global DOSE_MG constant
# (defined in the bench drivers, not the helper). Mirror it here so the
# diagnostic can build the NM dataset standalone.
if (!exists("DOSE_MG")) DOSE_MG <- 100

# ---- Load the exact dataset the bench uses --------------------------------
sim_path <- file.path(OPTS$input_root, sprintf("sim_obs_N%d", OPTS$N),
                      sprintf("sim_obs_scenario_%02d.rds", OPTS$scenario))
stopifnot(file.exists(sim_path))
sim_all <- readRDS(sim_path)

ds_i <- to_nm_dataset(sim_all) |>
  dplyr::filter(DATASET == OPTS$dataset) |>
  dplyr::select(-SCENARIO, -DATASET) |>
  dplyr::mutate(ID = as.integer(ID),
                SEX = as.integer(SEX),
                RACE = as.integer(RACE))

# ---- Candidate catalogue --------------------------------------------------
# Each entry adds ONE covariate term to the base model. `truth` flags the
# covariates that are genuinely in scenario 16 (should be selected).
#   cl ~ BW power , cl ~ CrCL power , vc ~ BW power , vc ~ SEX  (TRUE)
#   cl ~ BMI power                                             (NULL/decoy)
candidates <- tibble::tribble(
  ~key,            ~target, ~theta_name,        ~term,                     ~truth,
  "cl~BW.power",   "cl",    "cov_BW_power_cl",   "log(BW / 70)",            TRUE,
  "cl~CrCL.power", "cl",    "cov_CrCL_power_cl", "log(CrCL / 95)",          TRUE,
  "vc~BW.power",   "vc",    "cov_BW_power_vc",   "log(BW / 70)",            TRUE,
  "vc~SEX.lin",    "vc",    "cov_SEX_lin_vc",    "SEX",                     TRUE,
  "cl~BMI.power",  "cl",    "cov_BMI_power_cl",  "log(BMI / 26)",           FALSE
)

# ---- Build one candidate UI, warm-starting structural theta from base -----
# base_theta : named numeric from base_fit$theta (log-scale structural)
# base_omega : 2x2 omega matrix (eta.cl, eta.vc)
# cov_init   : starting value for the new covariate theta (0 = cold)
make_candidate_ui <- function(target, theta_name, term, cov_init,
                              base_theta, base_omega, bound = OPTS$cov_bound) {
  gv <- function(nm, default) {
    v <- suppressWarnings(as.numeric(base_theta[nm]))
    if (length(v) == 0L || is.na(v)) default else v
  }
  lTVCL <- gv("lTVCL", log(0.6)); lTVQ <- gv("lTVQ", log(1.8))
  lTVVc <- gv("lTVVc", log(20));  lTVVp <- gv("lTVVp", log(80))
  prop  <- gv("prop.err", 0.1)
  o11 <- base_omega["eta.cl", "eta.cl"]; o22 <- base_omega["eta.vc", "eta.vc"]
  o21 <- base_omega["eta.vc", "eta.cl"]

  term_cl <- if (target == "cl") sprintf(" + %s * %s", theta_name, term) else ""
  term_vc <- if (target == "vc") sprintf(" + %s * %s", theta_name, term) else ""

  ini_lines <- c(
    sprintf("lTVCL <- %.10g", lTVCL),
    sprintf("lTVQ  <- %.10g", lTVQ),
    sprintf("lTVVc <- %.10g", lTVVc),
    sprintf("lTVVp <- %.10g", lTVVp),
    "lTVKA <- fix(log(0.7))",
    sprintf("%s <- c(%.6g, %.10g, %.6g)", theta_name, bound[1], cov_init, bound[2]),
    sprintf("eta.cl + eta.vc ~ c(%.8g, %.8g, %.8g)", o11, o21, o22),
    sprintf("prop.err <- %.10g", prop)
  )
  model_lines <- if (identical(OPTS$structure, "ode")) {
    c(
      sprintf("cl <- exp(lTVCL%s + eta.cl)", term_cl),
      sprintf("vc <- exp(lTVVc%s + eta.vc)", term_vc),
      "q  <- exp(lTVQ)",
      "vp <- exp(lTVVp)",
      "ka <- exp(lTVKA)",
      "d/dt(depot)   <- -ka * depot",
      "d/dt(central) <-  ka * depot - (cl/vc)*central - (q/vc)*central + (q/vp)*periph",
      "d/dt(periph)  <-  (q/vc)*central - (q/vp)*periph",
      "cp <- central/vc",
      "cp ~ prop(prop.err)"
    )
  } else {
    c(
      sprintf("cl <- exp(lTVCL%s + eta.cl)", term_cl),
      sprintf("vc <- exp(lTVVc%s + eta.vc)", term_vc),
      "q  <- exp(lTVQ)",
      "vp <- exp(lTVVp)",
      "ka <- exp(lTVKA)",
      "cp <- linCmt()",
      "cp ~ prop(prop.err)"
    )
  }
  fn_text <- paste0(
    "function() {\n  ini({\n    ",
    paste(ini_lines, collapse = "\n    "),
    "\n  })\n  model({\n    ",
    paste(model_lines, collapse = "\n    "),
    "\n  })\n}"
  )
  eval(parse(text = fn_text), envir = globalenv())
}

# ---- Fit summariser --------------------------------------------------------
summarise_fit <- function(fit, base_objf, theta_name) {
  if (is.null(fit) || isTRUE(fit$.failed)) {
    return(list(objf = NA_real_, dOFV = NA_real_, pchisqr = NA_real_,
                theta_hat = NA_real_, conv = NA_integer_, msg = "FIT ERROR",
                cn = NA_real_, cov_ok = NA, switch = NA))
  }
  d <- diagnose_fit(fit)                       # from scm_bench_helpers.R
  objf <- as.numeric(fit$objf)
  dOFV <- base_objf - objf                     # forward-step convention
  pch  <- if (is.finite(dOFV) && dOFV > 0) 1 - stats::pchisq(dOFV, df = 1) else 1
  th   <- suppressWarnings(as.numeric(fit$theta[theta_name]))
  # Did the Almquist analytic OUTER gradient fall back to finite differences
  # for any iteration? This is the mechanism-A signal: fast=TRUE estimators
  # (foceif/irlsfoceif) log "could not be solved at this point" in $runInfo
  # when fallbackFD fires. focei (fast=FALSE) never has this -> NA.
  ri   <- tryCatch(fit$runInfo, error = function(e) NULL)
  swtch <- if (is.null(ri)) NA else any(grepl("could not be solved", ri))
  list(objf = objf, dOFV = dOFV, pchisqr = pch,
       theta_hat = if (length(th)) th else NA_real_,
       conv = d$convergence_code, msg = d$message %||% "",
       cn = d$cond_num_cor, cov_ok = d$cov_ok, switch = swtch)
}

# ---- Main sweep ------------------------------------------------------------
rows <- list()
base_fits <- list()   # keyed by "est/opt"

for (cell in OPTS$cells) {
  est <- cell[1]; opt <- cell[2]; cell_key <- paste(est, opt, sep = "/")
  message(sprintf("\n===== base fit: %s =====", cell_key))
  ctrl <- make_est_control(est, opt, "screen")$ctrl
  base_model <- if (identical(OPTS$structure, "ode"))
    base_2cmt_oral_ode else base_2cmt_oral_linCmt
  bf <- tryCatch(nlmixr2(base_model, ds_i,
                         est = nlmixr_est_name(est), control = ctrl),
                 error = function(e) { message("base failed: ", conditionMessage(e)); NULL })
  if (is.null(bf)) next
  base_fits[[cell_key]] <- bf
  base_objf <- as.numeric(bf$objf)
  base_theta <- bf$theta
  base_omega <- bf$omega
  message(sprintf("  base OFV = %.2f  (conv=%s)", base_objf,
                  diagnose_fit(bf)$convergence_code))

  for (k in seq_len(nrow(candidates))) {
    cc <- candidates[k, ]
    ui <- make_candidate_ui(cc$target, cc$theta_name, cc$term, cov_init = 0,
                            base_theta = base_theta, base_omega = base_omega)
    f <- tryCatch(suppressWarnings(nlmixr2(ui, ds_i,
                                           est = nlmixr_est_name(est), control = ctrl)),
                  error = function(e) list(.failed = TRUE))
    s <- summarise_fit(f, base_objf, cc$theta_name)
    rows[[length(rows) + 1L]] <- tibble::tibble(
      arm = "cold", estimator = est, optimizer = opt, cell = cell_key,
      candidate = cc$key, truth = cc$truth,
      base_objf = base_objf, cand_objf = s$objf, dOFV = s$dOFV,
      pchisqr = s$pchisqr, theta_hat = s$theta_hat,
      conv = s$conv, msg = s$msg, cn = s$cn, switch = s$switch
    )
    message(sprintf("  [cold|%-18s] %-14s dOFV=%9.2f  p=%.3g  theta=%s  conv=%s  fd_switch=%s  '%s'",
                    cell_key, cc$key, s$dOFV, s$pchisqr,
                    formatC(s$theta_hat, format = "g", digits = 3),
                    s$conv, s$switch, s$msg))
  }
}

# ---- Warm-start arm --------------------------------------------------------
# Seed EVERY gradient cell's covariate theta with the reference cell's fitted
# value (focei/bobyqa). The decisive question for each optimizer:
#   Does it MOVE from the warm start and exit cleanly (rescue), or does it stay
#   FROZEN at the init (like foceif/lbfgsb3c did in the runSCM 0.5-init run)?
warm_cells <- Filter(function(c) c[2] %in% c("nlminb", "lbfgsb3c"), OPTS$cells)
ref_theta <- dplyr::bind_rows(rows) |>
  dplyr::filter(arm == "cold", cell == OPTS$warm_src) |>
  dplyr::select(candidate, ref_theta = theta_hat)

for (cell in warm_cells) {
  est <- cell[1]; opt <- cell[2]; cell_key <- paste(est, opt, sep = "/")
  bf <- base_fits[[cell_key]]; if (is.null(bf)) next
  base_objf <- as.numeric(bf$objf)
  ctrl <- make_est_control(est, opt, "screen")$ctrl
  for (k in seq_len(nrow(candidates))) {
    cc <- candidates[k, ]
    if (!cc$truth) next
    warm_init <- ref_theta$ref_theta[ref_theta$candidate == cc$key]
    if (!length(warm_init) || !is.finite(warm_init) || abs(warm_init) < 1e-6)
      warm_init <- 0.5   # fall back to 0.5 when the reference declined to move
    ui <- make_candidate_ui(cc$target, cc$theta_name, cc$term,
                            cov_init = warm_init,
                            base_theta = bf$theta, base_omega = bf$omega)
    f <- tryCatch(suppressWarnings(nlmixr2(ui, ds_i,
                                           est = nlmixr_est_name(est), control = ctrl)),
                  error = function(e) list(.failed = TRUE))
    s <- summarise_fit(f, base_objf, cc$theta_name)
    moved <- is.finite(s$theta_hat) && abs(s$theta_hat - warm_init) > 1e-4
    rows[[length(rows) + 1L]] <- tibble::tibble(
      arm = "warm_theta", estimator = est, optimizer = opt, cell = cell_key,
      candidate = cc$key, truth = cc$truth,
      base_objf = base_objf, cand_objf = s$objf, dOFV = s$dOFV,
      pchisqr = s$pchisqr, theta_hat = s$theta_hat,
      conv = s$conv, msg = s$msg, cn = s$cn, switch = s$switch
    )
    message(sprintf("  [warm|%-18s] %-14s init=%.3f -> theta=%s (%s)  dOFV=%9.2f  conv=%s  fd_switch=%s  '%s'",
                    cell_key, cc$key, warm_init,
                    formatC(s$theta_hat, format = "g", digits = 3),
                    if (moved) "MOVED" else "FROZEN",
                    s$dOFV, s$conv, s$switch, s$msg))
  }
}

trace_tbl <- dplyr::bind_rows(rows)

# ---- Report ----------------------------------------------------------------
cat("\n\n================ SCM GRADIENT MECHANISM TRACE ================\n")
cat(sprintf("N=%d  scenario=%d  dataset=%d  structure=%s\n\n",
            OPTS$N, OPTS$scenario, OPTS$dataset, OPTS$structure))

cat("---- dOFV per candidate (LRT threshold chi2_1 @0.05 = 3.84) ----\n")
wide <- trace_tbl |>
  dplyr::mutate(col = paste(arm, cell, sep = " | ")) |>
  dplyr::select(candidate, truth, col, dOFV) |>
  tidyr::pivot_wider(names_from = col, values_from = dOFV)
print(as.data.frame(wide), digits = 4)

cat("\n---- selection outcome (would this candidate be accepted? p<=0.05) ----\n")
sel <- trace_tbl |>
  dplyr::mutate(accept = pchisqr <= 0.05,
                col = paste(arm, cell, sep = " | ")) |>
  dplyr::select(candidate, truth, col, accept) |>
  tidyr::pivot_wider(names_from = col, values_from = accept)
print(as.data.frame(sel))

cat("\n---- warm arm: did the optimizer MOVE off the warm init? ----\n")
warm_diag <- trace_tbl |>
  dplyr::filter(arm == "warm_theta") |>
  dplyr::mutate(status = ifelse(is.finite(theta_hat) &
                                  abs(dOFV) < 1e4 & conv %in% c(0L, 4L),
                                "rescued", "still-broken")) |>
  dplyr::select(cell, candidate, theta_hat, dOFV, conv, msg, status) |>
  dplyr::arrange(cell, candidate)
print(as.data.frame(warm_diag), right = FALSE)

cat("\n---- analytic outer gradient -> FD fallback (mechanism A) ----\n")
cat("TRUE  = the Almquist analytic gradient 'could not be solved' for >=1\n",
    "        iteration and fallbackFD swapped in an FD gradient (poisons the\n",
    "        quasi-Newton curvature). NA = fast=FALSE estimator (focei).\n", sep = "")
switch_diag <- trace_tbl |>
  dplyr::mutate(col = paste(arm, cell, sep = " | ")) |>
  dplyr::select(candidate, truth, col, switch) |>
  tidyr::pivot_wider(names_from = col, values_from = switch)
print(as.data.frame(switch_diag), right = FALSE)

cat("\n---- convergence codes / messages ----\n")
print(as.data.frame(
  trace_tbl |>
    dplyr::select(candidate, arm, cell, conv, msg, theta_hat) |>
    dplyr::arrange(candidate, cell, arm)
), right = FALSE)

# ---- Persist ---------------------------------------------------------------
dir.create(dirname(OPTS$out_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(trace = trace_tbl, opts = OPTS,
             base_objf = vapply(base_fits, function(x) as.numeric(x$objf), numeric(1))),
        OPTS$out_rds)
cat(sprintf("\nSaved trace -> %s\n", OPTS$out_rds))

cat("\n---- INTERPRETATION GUIDE ----\n")
cat("* If TRUE-covariate dOFV is large under cold/bobyqa but ~0 under\n",
    "  cold/nlminb (theta stuck ~0, msg 'false convergence (8)'), the LRT\n",
    "  never sees the signal -> covariate rejected -> Power collapses.\n",
    "* If warm_theta/nlminb dOFV recovers to ~bobyqa level, the cause is\n",
    "  COLD START (proximity), not the optimizer -> warm-starting candidates\n",
    "  is the correct fix, and gradient optimizers CAN stand on their own.\n",
    "* If warm_theta/nlminb dOFV stays ~0 even when started at the answer,\n",
    "  the surface noise itself defeats the FD gradient -> derivative-free\n",
    "  screening is required (a legitimate scientific finding).\n",
    "* If irlsfoceif/lbfgsb3c shows switch=TRUE on the candidates it botches,\n",
    "  the corruption is mechanism A (analytic->FD gradient swapping poisons\n",
    "  the quasi-Newton Hessian). Fix 1 (sensitivity-matched atolSens/rtolSens\n",
    "  in estimator_factory.R) should flip those to switch=FALSE.\n", sep = "")
