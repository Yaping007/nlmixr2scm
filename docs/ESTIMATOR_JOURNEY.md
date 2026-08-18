# The Estimator Expansion Journey — A PMx Beginner's Field Guide

*Cohort: `nlmixr2scm` benchmarking, N=40/80/300, 16 scenarios, Khandelwal 2019 replication*
*Last update: 2026-07-20*

---

## Table of Contents

1. [Estimators & Optimizers in Plain English](#1-estimators--optimizers-in-plain-english)
2. [Two Benchmarks, Two Purposes](#2-two-benchmarks-two-purposes)
3. [The Pitfalls We Fell Into & The Fixes](#3-the-pitfalls-we-fell-into--the-fixes)
4. [R Codebook — Inspecting Results](#4-r-codebook--inspecting-results)
5. [HPCE Codebook — LSF Submission & Monitoring](#5-hpce-codebook--lsf-submission--monitoring)
6. [Current Method Grid](#6-current-method-grid)

---

## 1. Estimators & Optimizers in Plain English

A **PopPK fit** has two nested tasks:

- **Outer task** — find the best *typical* parameters (θ) and variability (Ω) for the population. Runs an **outer optimizer** (bobyqa, nlminb, lbfgsb3c…) over the marginal likelihood.
- **Inner task** — for each subject, find the best individual "random-effect" parameters (η) given the current θ, Ω. Runs an **inner method** (Laplace, importance sampling, SAEM chain, VAE encoder…).

The **estimator** is the *combination* of an inner method + an approximation to the marginal likelihood.

### 1.1 The estimator family

| Estimator | Plain English | Inner method | Best for |
|---|---|---|---|
| **FOCE-I** (`focei`) | "First-Order Conditional Estimation with Interaction". The industry workhorse since the 1980s. Approximates the subject-level likelihood by a Taylor expansion around each subject's best η, correcting for the interaction between η and residual error. | Laplace + inner ETA optimization | General PopPK, especially rich sampling. |
| **FOCE-I with analytic gradient** (`foceif`) | Same math as FOCE-I but with an analytic gradient (Almquist 2015) instead of finite differences. Faster and less noisy when it works — but requires a compatible model. | Laplace + inner ETA + analytic ∂L/∂θ | Same as FOCE-I; faster when compatible. |
| **IRLS-FOCE-I** (`irlsfoceif`) | Iteratively-Reweighted Least Squares reformulation of FOCE-I. Should be more stable numerically. **Dropped from our grid** — its `cov` matrix is unpopulated in nlmixr2est's current class, so it fails our strict convergence check. | Same as foceif | Retry when nlmixr2est fixes cov step. |
| **SAEM** (`saem`) | "Stochastic Approximation Expectation-Maximization". Instead of optimizing η, it *samples* η from its posterior (an MCMC chain), then updates θ using expected values. Slower but robust to weird likelihood surfaces. | Metropolis-Hastings MCMC on η | Sparse data, complex random effects, poor-starting-value cases. |
| **VAE** (`vae`) | "Variational Autoencoder". A neural net learns to approximate the posterior of η. Cutting-edge; **dropped from our grid** because `runSCM` can't currently locate `cl/vc` symbols in the VAE-produced UI. | Neural encoder / decoder | Very large datasets; experimental. |

### 1.2 The outer optimizer family

| Optimizer | Type | How it works | Strengths | Weaknesses |
|---|---|---|---|---|
| **bobyqa** | **Derivative-free** ("trust-region") | Plants ~10 exploratory points around the current guess, fits a smooth quadratic bowl through them, walks to the bottom of that bowl, repeats. | *Deaf to noise* — averages tiny FOCE-I inner-ETA ripples out. Robust with bad starting values. Only optimizer that reliably starts SCM candidates at θ=0. | Slower per iteration when the surface is genuinely smooth. |
| **nlminb** | **Gradient (quasi-Newton)** | Measures the slope by taking two sideways steps (finite differences), walks steepest downhill, builds a Hessian approximation. | Fast when the surface is smooth. Widely used in NONMEM-era workflows. | Very sensitive to FD noise. Reports **"false convergence (8)"** when the local slope drowns in FOCE-I ripples — even when the fit is fine. |
| **lbfgsb3c** | **Gradient (limited-memory Newton, with bounds)** | Same idea as nlminb but keeps a *memory* of past gradients to guess the Hessian without storing it, and handles box constraints by projection. | Excellent when bounds are active and surface is smooth. | Same FD-noise sensitivity as nlminb. Its projected-gradient stopping rule can also mis-report code 8. |

### 1.3 Why the optimizer choice matters more than beginners expect

**FOCE-I's likelihood surface is *gravelly*.** The inner ETA re-optimization introduces tiny discontinuities (~0.01–0.1 OFV units) that look like a smooth surface plus noise. This is a **surface property, not a fit problem**. Different optimizers "hear" this noise differently:

- **bobyqa** builds a smooth bowl through many noisy measurements → deaf to ripples.
- **nlminb / lbfgsb3c** measure the local slope by tiny sideways steps → the ripples become the entire signal near the optimum.

The consequence: gradient optimizers frequently exit with `convergence = 8` ("false convergence") at points that are actually excellent fits. This is the **root cause of ~80% of the pitfalls** we hit.

---

## 2. Two Benchmarks, Two Purposes

### 2.1 `bench_refit_estimators` — "given the truth, can this estimator recover it?"

- **Model**: `make_true_model(scenario, boundary)` — the true structural model *with* the true covariates already present. Only continuous covariate slopes are bounded (±5 tight / ±10 narrow / ±10^5 wide).
- **Question answered**: parameter recovery accuracy (RMRSE, bias, CI coverage) *independent of model selection*.
- **Layout**:
  ```
  output/bench_refit_estimators/N{NN}/scn{SS}/{est}_{opt}/res_ds{DDD}.rds
  ```
- **Aggregator**: `script/aggregate_bench_refit_results.R`

### 2.2 `scm_bench` — "given only the base model, can this estimator find the true covariates?"

- **Model**: `base_2cmt_oral_linCmt()` — 2-cmt oral, **no covariates**. `runSCM()` then searches over `{CL, Vc} × {BW, CrCL, BMI, SEX, RACE} × {power, lin, cat}`.
- **Question answered**: operating characteristics of the SCM workflow (Power, false-positive rate, structural recovery).
- **Layout**:
  ```
  output/scm_bench/N{NN}/scn{SS}/{est}_{opt}/res_ds{DDD}.rds
  ```
- **Aggregator**: `script/aggregate_Bench_estimator_result.R`

### 2.3 The critical difference

- `bench_refit`: covariates always in the model. RMRSE is a **parameter estimation problem**.
- `bench_scm`: covariates must be **discovered by LRT**. Requires the estimator + optimizer to reliably distinguish OFV differences of ~3.84 units (χ²₁ at α=0.05). Any optimizer that adds >1 unit of OFV noise per candidate fit poisons the LRT and destroys Power.

**This is why `bench_refit` may look clean while `bench_scm` collapses for the same estimator.**

---

## 3. The Pitfalls We Fell Into & The Fixes

Presented chronologically with symptom → root cause → fix. Use as a checklist before running new arrays.

### 3.1 CRLF line endings on HPCE

- **Symptom**: bash scripts fail on HPCE with "bad interpreter" or R scripts fail with cryptic parse errors.
- **Root cause**: Windows editors save `\r\n` line endings; Linux tools expect `\n`.
- **Fix**: One-shot `sed -i 's/\r$//' script/**/*.{sh,lsf,R}` on HPCE after every `git pull`.

### 3.2 Stale HPCE code masquerading as new results

- **Symptom**: rerun on HPCE produces identical (unfixed) results.
- **Root cause**: forgot to `git push` from Windows, or forgot to `git pull` on HPCE.
- **Fix**: Always **push from Windows → pull on HPCE** before submitting arrays. Check with `git log -1 --oneline` on both sides.

### 3.3 lbfgsb3c "false convergence (8)" flagged as failure

- **Symptom**: 0% convergence for lbfgsb3c in `bench_refit`, despite excellent fits (finite OFV, cov OK, low CN).
- **Root cause**: old `diagnose_fit()` gated on `fit$convergence == 0`. lbfgsb3c returns code 8 for good fits.
- **Fix**: `diagnose_fit()` now uses **numerical evidence** (finite OFV + finite cov + finite CN). The raw code is kept as `convergence_code` for reference.

### 3.4 `saemControl(covMethod = "")` crashes

- **Symptom**: `Error in match.arg(covMethod): 'arg' should be one of "linFim","fim","sa","r,s","r","s"`. All 20 saem tasks exit in 3 seconds.
- **Root cause**: `foceiControl` accepts `""` (skip cov step); `saemControl` does not.
- **Fix**: hardcode `covMethod = "linFim"` for both saem tiers (cheap linearized FIM). Post-SAEM foceif refit still delivers `r,s` cov for the winning model.

### 3.5 `irlsfoceif` never satisfies strict convergence

- **Symptom**: 0% convergence forever. Model fits look fine.
- **Root cause**: `nlmixr2est`'s IRLS class does not populate `fit$conditionNumberCor` and the `cov2cor + kappa` fallback also fails.
- **Fix**: Dropped from `valid_combos()` and both HPCE submit defaults. Re-add once upstream fix ships.

### 3.6 PMx-strict convergence definition

- **Symptom**: mixing "converged" heuristics across scripts → inconsistent Power numbers.
- **Root cause**: no single agreed criterion.
- **Fix**: Adopted **PMx-strict flag** in both aggregators:
  ```
  converged_strict := (converged %in% TRUE)
                    & (cn_below_cutoff %in% TRUE)   # CN ≤ 1000
                    & !(est_bnd %in% TRUE)          # no bounded θ at boundary
  ```
  `%in% TRUE` idiom converts `NA → FALSE`, making aggregation NA-safe.

### 3.7 `nlminb` / `lbfgsb3c` collapse in SCM

- **Symptom**: `MedNSelected = 1–2`, Power ≈ 0. `bench_refit` for the same estimator looks fine.
- **Root cause**: `runSCM` starts every candidate covariate slope at **θ = 0**. At θ = 0 the covariate contributes nothing → the FD gradient of OFV w.r.t. that θ is essentially zero → gradient optimizer says "no signal" and quits with code 8 → LRT sees no improvement → real covariate rejected.
- **Fix (A)**: pass `inits = list(power = list(est = 0.5, lower = -5, upper = 5), lin = ..., cat = list(est = log(1.5), ...))` to `runSCM`. Non-zero starts give a real slope to detect.
- **Fix (B)**: per-optimizer `derivEps` — `nlminb → 1e-3` (was 5e-3), `lbfgsb3c → 5e-3`, `bobyqa → n/a`. Smaller FD step for nlminb sharpens the signal without dropping into ODE noise.
- **Non-fix**: `maxRetries` doesn't help. Same starting point → same noise → same wrong answer.

### 3.8 SAEM selects wrong covariate *shape* (lin instead of power)

- **Symptom**: saem_NA Power = 0 despite fits looking numerically fine.
- **Root cause**: screening tier used 300 burn / 400 EM. OFV noise ~5–10 units, larger than the χ²₁ threshold of 3.84 → shape-vs-shape LRT is a coin flip.
- **Fix**: unified SAEM chain to **500 burn / 1000 EM** across both tiers. Matches the paper's chain length.

### 3.9 Aggregator crashes on missing/NA fields

- **Symptom**: `NA/NaN argument` or dimension mismatch during aggregation.
- **Root cause**: mixing lenient converged (may be NA) with strict flags without guards.
- **Fix**: `%in% TRUE` idiom + explicit `is.finite()` guards in both aggregators.

### 3.10 SCM Power collapsed ~0.9 → ~0.6 (the forward-vs-backward winner bug)

- **Symptom**: the new `scm_bench` `focei_bobyqa` linCmt Power figure started at
  ~0.6 in easy/null scenarios where the ORIGINAL run
  (`output/full_scm_focei_bobyqa`) started at ~0.9. In the null scenario 1
  (`n_true = 0`, so `Power = 1 − false-positive rate`) the new run's FP-rate was
  ~0.36 and *rose* with N — the wrong direction for a correct α = 0.05 SCM.

- **The long false trail** (each *disproven* by a controlled A/B — record them as
  ruled-out, they are seductive but wrong):
  - ❌ **fork-over-OpenMP corruption** — 96% ODE convergence, sane CN/OFV.
  - ❌ **different Power-scoring machinery** — `compute_scm_power()` (new) and
    `compute_power_block()` (original) are byte-identical exact-match
    `Power = n_exact / N`.
  - ❌ **screening precision** (`sigdig` 3 vs 4, ODE `atol/rtol` 1e-6/1e-4 vs
    1e-8/1e-6) — a 10-dataset HPCE A/B (`SCREEN_SIGDIG=3` vs `4`) gave
    **byte-identical** results.
  - ❌ **inner-Hessian warm-start** (`foceiControl(warm=)`: nlmixr2est 6.2.0
    changed the default from classic `"save"` to `"calc"`) — a `WARM=save` vs
    `WARM=calc` A/B gave **byte-identical** results.
  - ❌ **base model / fit quality** — a per-dataset diff showed identical base
    OFVs and identical forward-step ΔOFVs / p-value counts across both runs.

- **The decisive evidence**: diff the SCM `step_hist` for a null dataset where
  the two runs disagreed (ds1/N300). BOTH runs' *backward* step dropped
  `cov_CrCL_power_vc` (ΔOFV 5.23 vs 5.22, p 0.0222 vs 0.0223,
  `included = "dropped"`), yet the ORIGINAL `selected` had **0 rows** while the
  NEW `selected` had **1 row** — it kept a covariate backward had explicitly
  dropped. The SCM search was the same; the divergence was purely in how
  `selected` was extracted.

- **Root cause**: `package_scm_schema21()` in `script/scm_bench_helpers.R`
  picked the SCM winner with **reversed precedence** —
  `winner <- .pickFit(scm_res$resFwd)` first, falling back to `resBck`. For a
  bidirectional `"scm"` search the FINAL model is the *backward-eliminated* one
  (`resBck`); taking `resFwd` retained every forward false-positive that
  backward correctly removed. The original `package_scm_result()` had the
  correct order (`resBck` first). In the null scenario this inflated the
  FP-rate and collapsed Power.

- **Fix**: reverse the precedence to match the original —
  ```r
  winner <- .pickFit(scm_res$resBck)
  if (is.null(winner)) winner <- .pickFit(scm_res$resFwd)  # forward-only search
  ```
  Verified by a `--force_repackage` pilot (10 datasets, cached SCM objects,
  no re-fitting): null scn-1 Power jumped `0.4/0.6/0.6 → 0.7/0.9/0.9` across
  N = 40/80/300. Committed in `7007b17`.

- **Incidental damage**: the same file carried merge corruption — a mangled
  `rec$wall_total_sec <- wall_totalchmark.` splice, an orphaned top-level
  fragment (which crashed `source()` with
  `object of type 'closure' is not subsettable` because `identity` resolved to
  base R's `identity()`), a missing `base_2cmt_oral_linCmt` definition, and a
  duplicated `# ...existing code...` tail block. All repaired in the same
  commit.

- **Lesson**: when the *metric* code and the *estimator* code both test clean,
  suspect the **packaging/extraction layer** between the search and the score.
  The step-history is ground truth — trust `included == "dropped"` over the
  mere presence of a `cov_*` theta in the final fit.

---

## 4. R Codebook — Inspecting Results

Store these snippets in an `.R` file, or copy-paste into the console.

### 4.1 Setup (once per session)

```r
library(dplyr); library(tidyr); library(tibble); library(purrr); library(readr)
`%||%` <- function(a, b) if (is.null(a)) b else a
setwd("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
```

### 4.2 Per-dataset peek — SCM run

```r
# Peek at res_ds01.rds for every cell in scn16/N80
cells <- c("focei_bobyqa","focei_nlminb","focei_lbfgsb3c",
           "foceif_nlminb","foceif_lbfgsb3c","saem_NA")

peek_scm <- function(cell, root = "output/scm_bench", N = 80, scn = 16, ds = 1) {
  fp <- file.path(root, sprintf("N%d", N), sprintf("scn%02d", scn), cell,
                  sprintf("res_ds%02d.rds", ds))
  if (!file.exists(fp)) { cat("MISSING:", fp, "\n"); return(invisible()) }
  r <- readRDS(fp)
  cat("=====", cell, "=====\n")
  cat(sprintf("t_total=%.1fs  t_scm=%.1fs  t_refit=%.1fs\n",
              r$t_total_sec %||% NA, r$t_scm_sec %||% NA, r$t_refit_sec %||% NA))
  cat("diag: converged=", r$test$diag$converged,
      " objf=", round(r$test$diag$objf, 1),
      " cn=", round(r$test$diag$cond_num_cor, 2),
      " cov_ok=", r$test$diag$cov_ok,
      " msg='", r$test$diag$message, "'\n", sep="")
  cat("selected covariates:\n"); print(r$test$selected)
  invisible(r)
}
invisible(lapply(cells, peek_scm))
```

### 4.3 Per-dataset peek — refit run

```r
peek_refit <- function(cell, root = "output/bench_refit_estimators",
                       N = 80, scn = 16, ds = 1) {
  fp <- file.path(root, sprintf("N%d", N), sprintf("scn%02d", scn), cell,
                  sprintf("res_ds%03d.rds", ds))
  if (!file.exists(fp)) { cat("MISSING:", fp, "\n"); return(invisible()) }
  r <- readRDS(fp)
  cat("=====", cell, "=====\n")
  cat(sprintf("t_fit=%.1fs  refit=%.1fs  boundary=%s status=%s\n",
              r$fit_runtime_sec %||% NA, r$refit_runtime_sec %||% NA,
              r$boundary, r$status))
  cat("diag:\n"); print(r$diag)
  cat("rel_err (head):\n"); print(head(r$rel_err, 12))
  invisible(r)
}
```

### 4.4 Aggregate the SCM benchmark

```r
source("script/aggregate_Bench_estimator_result.R")
tp  <- readRDS("Inputdataset/true_params_long.rds")

# root = "."  when files are in ./scm_bench (default download location)
# root = "output" when they're in output/scm_bench (aggregator's default)
res <- aggregate_bench_run(root = "output", true_params = tp,
                           write_outputs = FALSE, verbose = FALSE)

# The three key tables
res$diag_rates    # convergence % + CN + timing per cell
res$power         # Power / PowerCN / PowerMinSuc per cell
res$rmse_success  # RMRSE / bias per (cell, parameter), converged fits only
res$rmse_all      # same, all fits
res$relpower      # fraction recovering >=k of the n_true covariates




source("script/aggregate_Bench_estimator_result.R")
tp <- readRDS("Inputdataset/true_params_long.rds")
res <- aggregate_bench_run(root = "output", true_params = tp, write_outputs = FALSE, verbose = FALSE)
cat("cells discovered:", nrow(res$diag_rates), "\n\n===== diag_rates =====\n")
print(res$diag_rates |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 1))) |>
  dplyr::select(estimator, outer_opt, n_total, Converged_pct, ConvergedStrict_pct,
                CovStep_pct, CNBelowCutoff_pct, MedCN, MaxCN, MedObjF,
                MedNSelected, MedTotal_min), n = Inf, width = 200)
cat("\n===== power =====\n")
print(res$power |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 3))) |>
  dplyr::select(estimator, outer_opt, N, n_exact, Power, PowerCN, PowerMinSuc),
  n = Inf, width = 200)
cat("\n===== rmse_success (covariate coefs) =====\n")
print(res$rmse_success |>
  dplyr::filter(parameter %in% c("CLBW","CLcrCL","VcBW","VcSEX")) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(x) round(x, 1))) |>
  dplyr::select(estimator, outer_opt, parameter, N = n_used, MedRE_pct, RMRSE_pct, MARE_pct) |>
  dplyr::arrange(parameter, estimator, outer_opt), n = Inf)
```

### 4.5 Aggregate the refit benchmark

```r
source("script/aggregate_bench_refit_results.R")
res <- aggregate_bench_refit_run(root = "output", write_outputs = FALSE)

res$diag_rates    # convergence % + CN + timing
res$rmse_all      # RMRSE per (cell, parameter)
res$rmse_success  # restricted to converged_strict
```

### 4.6 Pretty-print the summary tables

```r
# Standard "print all rows" incantation
print_tbl <- function(x, digits = 1) {
  x |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), \(v) round(v, digits))) |>
    print(n = Inf, width = 220)
}

print_tbl(res$diag_rates)
print_tbl(res$power, digits = 3)
print_tbl(dplyr::filter(res$rmse_success,
                        parameter %in% c("CLBW","CLcrCL","VcBW","VcSEX")))
```

### 4.7 Quickly count files by cell

```r
# PowerShell equivalent baked into R
count_by_cell <- function(root = "scm_bench/N80/scn16") {
  list.files(root, pattern = "^res_ds\\d+\\.rds$",
             recursive = TRUE, full.names = FALSE) |>
    dirname() |>
    table()
}
count_by_cell()
```

### 4.8 Diff two estimator cells on the same dataset

```r
diff_cells <- function(cell_a, cell_b, ds = 1, scn = 16, N = 80) {
  ra <- peek_scm(cell_a, N = N, scn = scn, ds = ds)
  rb <- peek_scm(cell_b, N = N, scn = scn, ds = ds)
  tibble::tibble(
    metric   = c("objf", "cn", "n_selected", "power"),
    cell_a   = c(ra$test$diag$objf, ra$test$diag$cond_num_cor,
                 nrow(ra$test$selected %||% data.frame()), NA),
    cell_b   = c(rb$test$diag$objf, rb$test$diag$cond_num_cor,
                 nrow(rb$test$selected %||% data.frame()), NA)
  )
}
diff_cells("focei_bobyqa", "focei_nlminb")
```

---

## 5. HPCE Codebook — LSF Submission & Monitoring

**Queuing system**: IBM Spectrum LSF. Evidenced by `#BSUB` directives, `bsub` command, `LSB_JOBINDEX` env var, `%J`/`%I` log placeholders.

### 5.1 The daily hygiene sequence

```bash
# On HPCE, at start of each work session
cd ~/nlmixr2scm/nlmixr2scm
git pull
sed -i 's/\r$//' script/**/*.{sh,lsf,R} 2>/dev/null   # kill any CRLF
bjobs -A                                                # what's still running?
```

### 5.2 Submit ALL estimator × optimizer arrays for one (N, scenario)

```bash
# Default grid = focei/{bobyqa,nlminb,lbfgsb3c}, foceif/{nlminb,lbfgsb3c}, saem/NA
# Positional args: n_datasets [maxpar]
NS=80 SCENARIOS=16 bash script/hpce_bench/submit_all_arrays.sh 20 20

# Multiple scenarios in one shot
NS=80 SCENARIOS="1 4 8 16" bash script/hpce_bench/submit_all_arrays.sh 20 20

# Multiple cohorts (loops N × scn × est × opt)
NS="40 80 300" SCENARIOS="1 8 16" \
  bash script/hpce_bench/submit_all_arrays.sh 20 20
```

### 5.3 Submit ONE specific (est, opt) cell

```bash
# Positional: N scn est opt n_ds [maxpar] [ds_start]
bash script/hpce_bench/submit_one_array.sh 80 16 focei bobyqa 20 20
bash script/hpce_bench/submit_one_array.sh 80 16 saem   NA     20 20

# Resume from ds=5 through ds=20  (7th arg = ds_start)
bash script/hpce_bench/submit_one_array.sh 80 16 focei bobyqa 16 20 5
```

### 5.4 Restrict which estimators to submit

```bash
# Only saem  (useful after saem-only bugfix)
NS=80 SCENARIOS=16 ESTIMATORS="saem" \
  bash script/hpce_bench/submit_all_arrays.sh 20 20

# Only two gradient cells
NS=80 SCENARIOS=16 ESTIMATORS="focei foceif" \
  bash script/hpce_bench/submit_all_arrays.sh 20 20
```

### 5.5 Bench_refit submission (same structure, different top-level dir)

```bash
NS=80 SCENARIOS=16 BOUNDARY=tight \
  bash script/hpce_refit_estimators/submit_all_arrays.sh 20 20
```

### 5.6 Monitor running jobs

```bash
bjobs -A                                    # summary: NJOBS PEND DONE RUN EXIT SSUSP
bjobs -J 'bench_N8*'                        # only bench_N8* jobs
bjobs -l 240051 | head -60                  # everything about job 240051
bhist -a -J 'bench_N8*' | head -40          # historical: runtime / exit reason
```

**Reading `bjobs -A` column codes**:

| column | meaning | good/bad? |
|---|---|---|
| DONE | successfully finished | ✅ target: NJOBS == DONE |
| RUN | currently executing | 🟡 wait |
| PEND | queued, waiting for resources | 🟡 wait |
| EXIT | crashed | ❌ → check `.err` log |
| SSUSP | System-SUSPended (load throttling) | 🟡 auto-resumes; can `bresume` |

### 5.7 Inspect failed tasks

```bash
# Find latest error log
ls -lt logs/bench_N80/bench_N80_scn16_saem_NA.*.1.err | head

# Read the actual error (files are small, use cat not tail-glob)
cat logs/bench_N80/bench_N80_scn16_saem_NA.240053.1.err
echo "===== OUT ====="
cat logs/bench_N80/bench_N80_scn16_saem_NA.240053.1.out

# Post-mortem for a specific job
bhist -l 240051 2>/dev/null | tail -40
```

### 5.8 Kill jobs

```bash
bkill 240051                          # one job
bkill -J 'bench_N8*'                  # all matching name
bkill 0                                # ALL your jobs (use with caution)
```

### 5.9 Resubmit only the failed dataset indices

```bash
# Find missing res_ds*.rds under a cell
CELL=output/scm_bench/N80/scn16/focei_bobyqa
seq -w 1 20 | while read i; do
  [ -f "$CELL/res_ds${i}.rds" ] || echo "$i"
done
# then feed those indices back to submit_one_array.sh with a custom start/end
```

### 5.10 Rsync results back to Windows

```bash
# From HPCE (using rsync + a mounted share); or from Windows PowerShell:
# scp -r liuya8j@glchbs-sp221595:~/nlmixr2scm/nlmixr2scm/output/scm_bench .
```

---

## 6. Current Method Grid

**Included** (6 cells, as of 2026-07-13):

| Estimator | Outer optimizer | Notes |
|---|---|---|
| focei | bobyqa | The reliable reference. Highest Power in SCM. |
| focei | nlminb | Gradient method. Needs non-zero SCM inits + `derivEps=1e-3`. |
| focei | lbfgsb3c | Gradient method + bounds. Needs `derivEps=5e-3`. |
| foceif | nlminb | Analytic gradient version of nlminb. |
| foceif | lbfgsb3c | Analytic gradient version of lbfgsb3c. |
| saem | (none) | Chain 500 burn / 1000 EM. Cov via post-SAEM foceif+bobyqa refit. |

**Excluded (temporarily)**:

- `irlsfoceif × lbfgsb3c` — cov step unavailable in current `nlmixr2est`.
- `vae × *` — `runSCM` can't parse VAE-produced UI.

**Convergence rule** (both benchmarks):
```
converged        := is.finite(objf) & cov_ok & is.finite(cond_num_cor)
converged_strict := converged & (cond_num_cor <= 1000) & !any_theta_at_boundary
```

**Runtime budget (N=80, scn=16, per dataset)**:

| cell | typical wall-clock |
|---|---|
| focei_bobyqa | 8–10 min |
| focei_nlminb | 5–6 min |
| focei_lbfgsb3c | 7–8 min |
| foceif_* | 5–8 min |
| saem_NA | 15–20 min |

LSF template requests **4 CPU, 6 GB RAM, 3 h walltime** per task. `-R "span[hosts=1]"` keeps a task on a single node.

---

## 7. What's Next

- [ ] Verify (A) + (B) fix nlminb / lbfgsb3c SCM Power at scn16 × N80 → aim for Power ≥ 0.3 (was 0).
- [ ] Extend to all 16 scenarios × 3 cohorts (N=40, 80, 300).
- [ ] Fig 3 renderer: 4 cov params × 3 boundaries per cohort × scenario, Khandelwal 2019 layout.
- [ ] Consider optional "SCM screening with bobyqa, final polish with cell's optimizer" arm — the pragmatic PMx workflow.

---

*File: `docs/ESTIMATOR_JOURNEY.md`. Update it after every substantive fix — it will save your future self hours.*
