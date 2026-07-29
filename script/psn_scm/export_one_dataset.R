# ==============================================================================
# export_one_dataset.R  --  PsN/NONMEM SCM benchmark: single-cell pilot
# ------------------------------------------------------------------------------
# Exports ONE (N, scenario, dataset) cell from the shared simulated populations
# into a NONMEM-ready CSV + control stream + PsN scm config, so it can be run on
# DaVinci HPCE with:
#
#     module load PsN/5.5.0-GCCcore-12.3.0
#     bsub scm run.scm -nm_version=7.5.1-gfortran --lsf_options="-n 12 -W 10000"
#
# This mirrors to_nm_dataset() from script/scm_bench_helpers.R exactly:
#   * single 100 mg oral dose into depot (CMT=1) at TIME=0, EVID=1
#   * observations of cp_obs in central (CMT=2), EVID=0
#   * covariates BW BMI CrCL SEX RACE carried on every row
#
# Structural model = 2-cmt oral (ADVAN4 TRANS4), KA fixed, seeded at the true
# values (TVCL=0.6, TVQ=1.8, TVVc=20, TVVp=80, KA=0.7; Omega 0.1/0.02/0.1;
# proportional error 0.1).  These are the same inits as base_2cmt_oral_ode().
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ---- cell selection (CLI-overridable so this scales to the full grid) -------
# Usage:
#   Rscript script/psn_scm/export_one_dataset.R [--N 300] [--scenario 16]
#           [--dataset 1] [--input_root Inputdataset] [--out_root output/psn_scm/runs]
.argv <- commandArgs(trailingOnly = TRUE)
.opt <- function(flag, default) {
  i <- which(.argv == flag)
  if (length(i) && i < length(.argv)) .argv[i + 1L] else default
}
N          <- as.integer(.opt("--N", "300"))
SCEN       <- as.integer(.opt("--scenario", "16"))
DATASET    <- as.integer(.opt("--dataset", "1"))
INPUT_ROOT <- .opt("--input_root", "Inputdataset")
OUT_ROOT   <- .opt("--out_root", file.path("output", "psn_scm", "runs"))
DOSE_MG    <- 100

# ---- resolve paths ---------------------------------------------------------
here     <- normalizePath(file.path(INPUT_ROOT, sprintf("sim_obs_N%d", N)),
                          mustWork = TRUE)
sim_rds  <- file.path(here, sprintf("sim_obs_scenario_%02d.rds", SCEN))
# runs/ layout: heavy, disposable, git-ignored.  The parser later emits the
# light record into output/psn_scm/records/N*/scn*/ds*/.  The exporter writes
# the input files directly into the cell dir; submit_scm.sh then creates
# base_run/ scm_dir/ refit/ logs/ alongside them.
cell_dir <- file.path(OUT_ROOT, sprintf("N%d", N),
                      sprintf("scn%02d", SCEN), sprintf("ds%03d", DATASET))
dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
out_dir  <- cell_dir

# ---- build the NONMEM dataset for this single replicate --------------------
sim <- readRDS(sim_rds) |>
  dplyr::filter(DATASET == !!DATASET)

stopifnot(nrow(sim) > 0)

obs_rows <- sim |>
  dplyr::filter(!(time == 0)) |>            # drop t=0 obs (DV=0 breaks prop error); parity with nlmixr2
  dplyr::transmute(
    ID   = as.integer(SUBJECT),
    TIME = time,
    EVID = 0L,
    MDV  = 0L,
    AMT  = 0,
    CMT  = 2L,               # central
    DV   = cp_obs,
    BW, BMI, CRCL = CrCL, SEX = as.integer(SEX), RACE = as.integer(RACE)
  )

dose_rows <- sim |>
  dplyr::distinct(SUBJECT, BW, BMI, CrCL, SEX, RACE) |>
  dplyr::transmute(
    ID   = as.integer(SUBJECT),
    TIME = 0,
    EVID = 1L,
    MDV  = 1L,
    AMT  = DOSE_MG,
    CMT  = 1L,               # depot
    DV   = 0,                # NONMEM: DV on dose row ignored (MDV=1)
    BW, BMI, CRCL = CrCL, SEX = as.integer(SEX), RACE = as.integer(RACE)
  )

nm <- dplyr::bind_rows(dose_rows, obs_rows) |>
  dplyr::arrange(ID, TIME, dplyr::desc(EVID))

csv_path <- file.path(out_dir, "data.csv")
# NONMEM prefers no quoting; write plain comma-separated with a header.
utils::write.csv(nm, csv_path, row.names = FALSE, quote = FALSE, na = ".")
cat("wrote", nrow(nm), "rows ->", csv_path, "\n")

# ---- base NONMEM control stream (2-cmt oral, ADVAN4 TRANS4) -----------------
# TRANS4 params: CL, V2 (central), Q, V3 (peripheral), KA.
# KA fixed at 0.7 (log space in nlmixr2; here on natural scale).
base_mod <- '$PROBLEM 2cmt oral base (scn16 N300 ds001) -- SCM benchmark
$INPUT ID TIME EVID MDV AMT CMT DV BW BMI CRCL SEX RACE
$DATA data.csv IGNORE=@

$SUBROUTINE ADVAN4 TRANS4

$PK
  TVCL = THETA(1)
  TVV2 = THETA(2)
  TVQ  = THETA(3)
  TVV3 = THETA(4)
  TVKA = THETA(5)

  CL = TVCL * EXP(ETA(1))
  V2 = TVV2 * EXP(ETA(2))
  Q  = TVQ
  V3 = TVV3
  KA = TVKA
  S2 = V2

$ERROR
  IPRED = F
  Y = IPRED * (1 + EPS(1))

; --- typical values seeded at the data-generating truth ---------------------
$THETA
  (0, 0.6)      ; TVCL
  (0, 20)       ; TVV2  (central volume)
  (0, 1.8)      ; TVQ
  (0, 80)       ; TVV3  (peripheral volume)
  0.7 FIX       ; TVKA  (unidentifiable from sparse design)

$OMEGA BLOCK(2)
  0.1           ; var(eta.CL)
  0.02  0.1     ; cov, var(eta.V2)

$SIGMA
  0.01          ; proportional error variance (~0.1 CV)

; FOCEI-INTER, SIGDIG=4 mirrors nlmixr2 focei+bobyqa (sigdig=4).
; NOTE: NO $COVARIANCE here -- SCM screening runs cov-step OFF, exactly like
; nlmixr2 make_est_control(tier="screen", covMethod="").  The covariance step
; is run ONCE on the final selected model (see final_refit.template below).
$ESTIMATION METHOD=1 INTER MAXEVAL=9999 SIGDIG=4 PRINT=5 NOABORT
'
writeLines(base_mod, file.path(out_dir, "base.mod"))
cat("wrote base.mod\n")

# ---- PsN scm config --------------------------------------------------------
# Search space = the SAME 16 candidates nlmixr2 tests:
#   3 continuous {BW,CRCL,BMI} x 2 params {CL,V2} x 2 shapes {exp,power} = 12
# + 2 categorical {SEX,RACE}   x 2 params {CL,V2}                        =  4  => 16
#
# valid_states (PsN standard numbering): 1=not-incl, 2=linear, 3=hockey-stick,
#   4=exponential, 5=power.  nlmixr2's shapes are log-additive (inside exp()):
#     "power" = log(cov/med)  -> TVP*(cov/med)^theta   == PsN state 5 (power)
#     "lin"   = (cov - med)   -> TVP*exp(theta*(cov-med)) == PsN state 4 (exponential)
#   so continuous=1,4,5 reproduces nlmixr2's {power, lin} shape set exactly.
#
# CATEGORICAL (SEX, RACE): nlmixr2 uses the EXPONENTIAL form  TVP*exp(theta_nlmixr*I),
# whereas PsN's native categorical (state 2) is LINEAR-multiplicative
# TVP*(1 + theta_psn*I).  PsN cannot emit exp() categorical without fragile
# custom-code injection, but the two are the SAME one-parameter fold-change for a
# binary level, so there is an EXACT algebraic bridge applied in the parser:
#     exp(theta_nlmixr) = 1 + theta_psn
#     => theta_nlmixr = log(1 + theta_psn)
#     => SE_nlmixr     = SE_psn / (1 + theta_psn)   (delta method)
# Reference level: for a binary covariate PsN's most-frequent reference == the
# nlmixr2 most-frequent reference, so the fold-change (and hence |beta|) matches;
# any reference flip only inverts the sign, which the aggregator already handles.
# Selection/detection scoring is unaffected (same relation, reparameterized).
# We therefore keep categorical=1,2 (native) and back-transform theta afterward.
#
# Centering: the STRUCTURAL typical value (TVCL/TVV2) is the parameter value AT
# the covariate center.  The truth is defined at fixed physiological references
# (BW=70, CrCL=95), so we pin PsN's center THERE via the [code] section -- PsN
# then estimates the intercept directly at 70/95 and NO post-hoc re-centering is
# needed.  theta is center-invariant (it only reanchors the intercept), so
# covariate SELECTION + coefficient comparison are unaffected:
#     exp(theta*(cov - c)) = exp(theta*cov) * exp(-theta*c)   [exp(-theta*c) -> TVCL]
#     (cov / c)^theta      = cov^theta      * c^(-theta)      [c^(-theta)    -> TVCL]
# The [code] section redefines only the exp (state 4) and power (state 5) forms
# for BW and CRCL, substituting the fixed anchor for PsN's default `median`.
# Both forms stay strictly positive, so PsN's default theta bounds
# (-1000000, 0.001, 1000000) keep the covariate function positive -- no
# [lower_bounds]/[upper_bounds] override required.  BMI (a distractor with no
# physiological reference) keeps PsN's median default.  Wildcard `*:` applies to
# both CL and V2.  (An earlier [reference_values] block was an INVALID PsN 5.5
# section -- "Found invalid section: reference_values" -- which aborted scm; the
# [code] section is the correct, documented mechanism.)
#
# Forward p=0.05, backward p=0.01.  linearize=0 => full re-estimation each step
# (LRT parity with nlmixr2 runSCM).
scm_cfg <- '
model=base.mod
search_direction=both
p_forward=0.05
p_backward=0.01
linearize=0

continuous_covariates=BW,CRCL,BMI
categorical_covariates=SEX,RACE

[test_relations]
CL=BW,CRCL,BMI,SEX,RACE
V2=BW,CRCL,BMI,SEX,RACE

[valid_states]
continuous=1,4,5
categorical=1,2

[code]
;fix centering at fixed references so TVCL/TVV2 are estimated at BW=70, CrCL=95
;state 4 = exponential, state 5 = power (PsN default numbering)
*:BW-4=PARCOV=EXP(THETA(1)*(COV-70))
*:BW-5=PARCOV=((COV/70)**THETA(1))
*:CRCL-4=PARCOV=EXP(THETA(1)*(COV-95))
*:CRCL-5=PARCOV=((COV/95)**THETA(1))

[lower_bounds]
*:BW-4=-5
*:BW-5=-5
*:CRCL-4=-5
*:CRCL-5=-5

[upper_bounds]
*:BW-4=5
*:BW-5=5
*:CRCL-4=5
*:CRCL-5=5
'
writeLines(scm_cfg, file.path(out_dir, "run.scm"))
cat("wrote run.scm\n")

# ---- final-refit template (covariance step ON, mirrors nlmixr2 tier=final) --
# After scm finishes, its winning model has NO $COVARIANCE (screening tier).
# Reproduce nlmixr2's single tight-tol R,S covariance refit on the winner:
# submit_scm.sh takes the winning candidate model, appends this line, and runs
# `execute` on it.
#   MATRIX=RSR   -> robust R,S sandwich  R^-1 S R^-1  (== nlmixr2 covMethod="r,s")
#                   NOT MATRIX=R (that would be Hessian-only R^-1, non-robust).
#   UNCONDITIONAL-> SEs for ALL params even near a boundary, so the correlation
#                   matrix (and its condition number) is always complete/comparable.
#   PRINT=E      -> print eigenvalues of the correlation matrix -> cond_num_cor.
final_refit_cov <- "$COVARIANCE UNCONDITIONAL MATRIX=RSR PRINT=E"
writeLines(final_refit_cov, file.path(out_dir, "final_refit_cov.txt"))
cat("wrote final_refit_cov.txt\n")

# ---- copy the timed submission wrapper into the cell dir -------------------
# submit_scm.sh records wall_sec (= runSCM_traced elapsed_s) into timing.json;
# the parser reconstructs cpu_sec + hog from per-subrun .lst files.
submit_src <- file.path("script", "psn_scm", "submit_scm.sh")
if (file.exists(submit_src)) {
  file.copy(submit_src, file.path(out_dir, "submit_scm.sh"), overwrite = TRUE)
  cat("copied submit_scm.sh\n")
}

# seed_refit.R = update_inits replacement (seeds refit initials from winner .ext
# so the covariance refit OFV matches the SCM winner exactly).  submit_scm.sh
# runs it from the cell dir, so ship a copy alongside submit_scm.sh.
seed_src <- file.path("script", "psn_scm", "seed_refit.R")
if (file.exists(seed_src)) {
  file.copy(seed_src, file.path(out_dir, "seed_refit.R"), overwrite = TRUE)
  cat("copied seed_refit.R\n")
}

cat("\n--- HPCE commands (run from", cell_dir, ") ---\n")
cat("bsub -n 4 -W 10000 \"bash submit_scm.sh\"        # scm(base+search) + mandatory cov refit\n")
cat("\n--- parse (from repo root) ---\n")
cat(sprintf("Rscript script/psn_scm/parse_psn_scm.R --cell %s --N %d --scenario %d --dataset %d\n",
            cell_dir, N, SCEN, DATASET))
