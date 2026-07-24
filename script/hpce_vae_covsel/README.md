# HPCE VAE covariate-selection pilot

One-run VAE covariate selection (`nlmixr2(..., est = "vae")`,
`covariateSelection = TRUE`) across the full pilot grid on the LSF cluster.
**`runSCM` is not used** — the stepwise search is internal to the VAE training
loop (BICc-ELBO penalty).

## Grid

| axis | values | count |
|------|--------|-------|
| cohort `N` | 40, 80, 300 | 3 |
| scenario | 1 … 16 | 16 |
| structure | linCmt, ode | 2 |
| datasets / cell (pilot) | 1 … `NDS` (default **5**) | 5 |

- One **LSF array per `(N, scenario, structure)` cell**; the array index
  (`$LSB_JOBINDEX`) **is the dataset id**.
- Full pilot = 3 × 16 × 2 = **96 arrays** × 5 datasets = **480 tasks**.

## Files

### Submission (this directory)

| file | role |
|------|------|
| `vae_covsel_array.lsf`   | LSF array template; reads `SAMPLE_N/SCN/STRUCTURE` from env, calls `script/vae_covsel_driver.R --dataset $LSB_JOBINDEX`. |
| `submit_one_array.sh`    | Submit ONE cell: `bsub -J "vaecov_N<N>_scn<SS>_<struct>[1-NDS]%MAXPAR"`. |
| `submit_all_arrays.sh`   | Triple loop over `NS × SCENARIOS × STRUCTURES`. |

### R scripts (under `script/`, required at runtime)

| file | role |
|------|------|
| `vae_covsel_driver.R`    | **Entry point** run by each array task. Loads the cell's dataset, fits VAE with `covariateSelection = TRUE`, derives the per-scenario `true_set`, and writes the schema-2.1 sidecars. |
| `refit_helpers.R`        | `to_nm_dataset()`, `PsN_scenarios` (scenario→covariate indicators). |
| `true_model_factory.R`   | Loaded for parity with the smoke file (model-attr helpers). |
| `scm_bench_helpers.R`    | Base models `base_2cmt_oral_linCmt` / `base_2cmt_oral_ode`. |
| `output_schema.R`        | `assemble_common()` (schema 2.1) and `write_fit_sidecar()`. |
| `estimator_factory.R`    | `seed_for_dataset()`, `seed_all()` for reproducible seeding. |

The driver `source()`s the five helpers automatically via its `.script_dir`
resolver, so only `vae_covsel_driver.R` is named on the command line.

### R package prerequisites

`nlmixr2 >= 3.0`, `nlmixr2est >= 7.0.1` (provides `est = "vae"` / `vaeControl`),
`rxode2` (dev build exporting `getIndCmt`), `nlmixr2scm`, `dplyr`, `tibble`,
`tidyr`.

**No `torch` / libtorch.** `est = "vae"` is a native C++/Armadillo engine.
The only optional extra is **`L0Learn`** (CRAN), used solely for large
covariate searches (`covSelectMethod = "l0learn"`/`"auto"`). The four-covariate
grid here uses the exact branch-and-bound, so `L0Learn` is not required — but it
is cheap to install and future-proofs wider grids.


## Output layout (schema 2.1)

```
output/vae_covsel_pilot/N<N>/scn<SS>_<structure>/covsel/
    res_ds<DDD>.rds        # full record (rel_err, diag, cov, $covsel block)
    res_ds<DDD>.fit.rds    # raw nlmixr2 fit
    res_ds<DDD>.meta.json  # greppable scalar manifest
    res_ds<DDD>_ERROR.txt  # only if the task failed (missing input, fit error)
```

The `$covsel` block carries `selected` (VAE's promoted `beta_*` terms mapped to
`var/covar`) and `true_set` (the scenario's true covariates, derived from the
`PsN_scenarios` indicators) — the basis for per-cell TP/FN/FP once aggregated.

## Engine prerequisite (no torch)

`est = "vae"` is a **native C++/Armadillo** LSTM encoder with an analytic
backward pass; the decoder is your ordinary `rxode2` model. There is **no
`torch` / libtorch dependency** — nothing to download on the login node.

What each node does need is a working **`nlmixr2est >= 7.0.1`** install (which
links against `rxode2`, `RcppParallel`/`tbb`, `stringfish`). If a node's build
is broken you will see a `LoadLibrary`/`getIndCmt` error at `library()` time,
not a torch error. Optionally install `L0Learn` for large covariate searches:

```r
install.packages("L0Learn")   # optional; only for covSelectMethod="l0learn"/"auto"
```

## Recommended workflow: probe first, then full sweep

VAE resource use on the cluster is unknown, so submit **one array** and read
its actual usage before launching all 96.

### 1. Single probe (known-good cell: N80 × scn16 × linCmt × 5 ds)

```bash
cd ~/nlmixr2scm && git pull
# normalize line endings + clear any stale ERROR sidecars first
sed -i 's/\r$//' script/hpce_vae_covsel/*.sh script/hpce_vae_covsel/*.lsf script/*.R
rm -f output/vae_covsel_pilot/N80/scn16_linCmt/covsel/*_ERROR.txt #Remove them from the probe cell:
find output/vae_covsel_pilot -name '*_ERROR.txt' -delete # clear stale ERROR files across the whole pilot tree

NS=80 SCENARIOS=16 STRUCTURES=linCmt \
  bash script/hpce_vae_covsel/submit_all_arrays.sh 5 20
# equivalently:
# bash script/hpce_vae_covsel/submit_one_array.sh 80 16 linCmt 5 20

FULL sweep: all cohorts x all 16 scenarios x both structures,
# datasets 6..250 (245 per cell), 40 concurrent tasks per array
bash script/hpce_vae_covsel/submit_all_arrays.sh 245 40 6
#submit_all_arrays.sh takes [n_ds] [maxpar] [ds_start], and DS_END = DS_START + NDS - 1, so ds_start=6, n_ds=245 gives indices 6–250.

NS="40" bash script/hpce_vae_covsel/submit_all_arrays.sh 245 40 6   
#Submitted 32 array(s), 7840 task(s) total. 2 model parameterization *16 scenarios- each scenario have 245 tasks. 
#The testing queue at novartis has 25 hosts/nodes; each node has 64-core. ~1000 jobs = 20-host limits 
#Occupy at most 20 nodes at once. With mostly 64-core nodes → ~1,280 cores → matches your ~1,066 running tasks. 
#So you're already near that ceiling; the rest pend until tasks finish.

bgadd -L 300 /liuya8j/vaecov 2>/dev/null   # # never run more than 150 tasks at once
NS="80" bash script/hpce_vae_covsel/submit_all_arrays.sh 245 10 6  
NS="300" bash script/hpce_vae_covsel/submit_all_arrays.sh 245 10 6
```

Monitor and inspect resource usage:

```bash
bjobs -J 'vaecov_*'                 # running/pending
bjobs -l <jobid>                    # live per-task detail
bacct -l <jobid>                    # peak MEM / CPU / walltime AFTER completion
ls output/vae_covsel_pilot/N80/scn16_linCmt/covsel/   # res_ds00{1..5}.*
```

Check for task failures:

```bash
find output/vae_covsel_pilot -name '*_ERROR.txt' -exec sed -n '1,6p' {} +
```

### 2. Tune resources

From `bacct -l`, read **MAX MEM** and **CPU TIME**; adjust `#BSUB -M` and
`#BSUB -W` in `vae_covsel_array.lsf` if the 6000 MB / 3 h probe values are off.

### 3. Full pilot

```bash
bash script/hpce_vae_covsel/submit_all_arrays.sh 5 20
bjobs -J 'vaecov_*'
```

## Aggregation

Once (some or all of) the sweep has landed on disk, roll the per-dataset
`res_ds*.rds` records up into operating characteristics with
`script/aggregate_vae_covsel.R`. It scans
`output/<sub>/N<N>/scn<SS>_<structure>/covsel/res_ds*.rds` and writes 11 CSVs +
1 bundled `.rds` to `output/vae_covsel_aggregated/`.

**Key argument — `--sub` / `sub`** selects which run tree to scan. The full
sweep lives in `output/vae_covsel_full/`, so pass `--sub vae_covsel_full`
(the default is `vae_covsel_pilot`). `root` stays `output` — it is the common
parent of both the input tree (`output/<sub>`) and the output tree
(`output/vae_covsel_aggregated`).

### Approach A — command line (HPCE-friendly)

```bash
Rscript script/aggregate_vae_covsel.R --root output --sub vae_covsel_full
# optional custom output dir:
# Rscript script/aggregate_vae_covsel.R --root output --sub vae_covsel_full --out_dir output/vae_covsel_full_aggregated
```

### Approach B — interactive R

```r
source("script/aggregate_vae_covsel.R")
res <- aggregate_vae_covsel_run(root = "output", sub = "vae_covsel_full")

res$power                                                 # Power / PowerCN / PowerMinSuc per cell
res$diag_rates                                            # %converged, CN, runtime
subset(res$estim_all, param_class == "covariate_beta")    # clean covariate-beta rel-err
res$covsel_by_cov                                         # per-covariate detection rate
res$relpower                                              # per-k fraction_at_least_k
```

Both approaches produce identical files in `output/vae_covsel_aggregated/`:

| file | contents |
|------|----------|
| `vae_file_index.csv`      | one row per discovered RDS (+ `has_error` flag) |
| `vae_diag_long.csv`       | per fit: convergence, CN, timing, selection tally |
| `vae_rse_long.csv`        | per (fit, parameter): `rel_err` + `param_class` |
| `vae_covsel_long.csv`     | per (fit, var, covar): `in_true`/`in_vae`/`verdict` (TP/FN/FP) |
| `vae_diag_rates.csv`      | per cell: %Converged(Strict), CN, MedObjF, runtime |
| `vae_estim_all.csv`       | per (cell, parameter): MedRE / MARE / RMRSE (all fits) |
| `vae_estim_success.csv`   | same, strict-converged fits only |
| `vae_estim_cond.csv`      | same, exact-match fits only |
| `vae_power.csv`           | per cell: Power / PowerCN / PowerMinSuc |
| `vae_relpower.csv`        | per (cell, k): fraction recovering ≥ k true covariates |
| `vae_covsel_by_covar.csv` | per (cell, var, covar): detection rate + TP/FP/FN |

It aggregates **whatever is on disk at call time** — partially-complete cells
just contribute fewer datasets (smaller `N` in `vae_power`/`vae_diag_rates`), so
re-run it after the sweep finishes for final numbers.

**Reference caveat** (baked into `param_class`): the DGP centres continuous
covariates at fixed references (`BW/70`, `CrCL/95`), while both VAE and runSCM
median-centre. Covariate-β coefficients are centring-invariant, so
`param_class == "covariate_beta"` rel-err is the clean, comparable estimation
metric. `structural_intercept` (TVCL/TVVc) rel-err is reference-dependent —
compare with care.

## Notes

- **`assemble_common` is untouched** — the VAE-specific `beta_*` → true-label
  mapping and `rel_err` backfill live only in the driver, so the FOCEi HPCE
  pipeline is unaffected.
- The **CRCL uppercase alias** (`CRCL = CrCL`) is added in the driver because
  VAE uppercases covariate names when building `beta_lTVCL_CRCL`.

## ⚠️ Known regression: `nlmixr2est ≥ 7.0.0` and `fit$theta`

**Symptom.** After a manager-side `nlmixr2est` update (≥ 7.0.0), a fresh full
re-run scored **`Power = 0` / all-FN** for every cell — even scenarios where the
raw fit clearly promoted the correct covariates (verified by hand for
`N300 × scn16 × ode`).

**Root cause.** The updated `nlmixr2est` makes **`fit$theta` return `NULL`**.
The old driver derived `covsel$selected` by grepping `names(fit$theta)` for
`^beta_<PARAM>_<COV>`, so it silently produced `selected = NULL` — no covariate
was ever counted as selected. A `tryCatch` masked it as "nothing selected"
rather than an error, so **no `_ERROR.txt` was written** (a loud failure became
a silent one). The promoted coefficients were never lost — they still live in
the record's **`parFixed`** table (rownames `beta_lTVCL_BW`, `beta_lTVCL_CRCL`,
`beta_lTVVc_BW`, `beta_lTVVc_SEX`; column `Estimate`).

**Fixes (both applied).**

1. **Driver (`vae_covsel_driver.R`) — prevents recurrence.**
   `selected` is now read via a version-robust accessor
   `.named_theta_estimates()` that prefers `parFixed` (→ `fit$parFixedDf` →
   `fit$parFixed` → legacy `fit$theta`). It also emits a **loud `warning()`**
   when 0 `beta_*` terms are extracted for a scenario that *has* true
   covariates, so a future API change cannot silently zero the results again.

2. **Aggregation (`aggregate_vae_covsel.R`) — salvages existing runs.**
   `.ensure_selected()` rebuilds `covsel$selected` from `parFixed` whenever the
   stored value is `NULL`/empty, and `.ensure_relerr_backfill()` then fills the
   NA covariate-β rows of `rel_err` (`CLBW`, `CLcrCL`, `VcBW`, `VcSEX`) from the
   same recovered estimates — so **power, FP/FN, detection, AND per-β relative
   error / RMRSE / MARE all come back**. **No re-fit needed** — just
   re-aggregate the affected tree, e.g.:

   ```bash
   Rscript script/aggregate_vae_covsel.R --root output --sub vae_covsel_full0722
   ```

   Validation on the hand-checked cell: `N300 × scn16 × ode` went from
   `Power = 0` to **`Power = 0.90`** (27/30 exact), all four true effects at
   `detection_rate = 1.0`, and covariate-β `rel_err` fully populated
   (e.g. `CLBW` −17 %, `CLcrCL` +10 %, `VcBW` −27 %, `VcSEX` −26 %).

**Scope.** The `diag` / `diag_t3` blocks and all scalars (convergence, CN,
`cov_ok`, `min_suc`, objf, runtime) were **unaffected** by the update, so
convergence / CN gating and the estimation-metric denominators never broke.
The only two casualties — `selected` and the covariate-β `rel_err` — are both
recovered from `parFixed` at aggregation time. **The VAE aggregation is now
fully compatible with `nlmixr2est ≥ 7.0.0`; every output (Power/PowerCN/
PowerMinSuc, FP/FN/detection, MedRE/MARE/RMRSE, relpower) populates correctly.**
Only the distractor-β rows (`CLBMI`, `VcBMI`, `VcCrCL`, `VcRACE`) remain `NA` —
by design, since they have no `true_value` to compare against.

> **Note.** This regression and its recovery are specific to the **VAE**
> covariate-selection workflow. The FOCEi / `runSCM` benchmark is a separate
> pipeline and is not addressed here.
