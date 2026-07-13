# SCM Estimator x Optimizer Benchmark

This directory contains the HPCE infrastructure for the SCM operating-
characteristics benchmark comparing estimation methods and outer optimizers
on the Khandelwal 2019 PopPK model.

## Files

| File | Purpose |
|---|---|
| `bench_array.lsf` | LSF job-array template (1 task per dataset in a cell) |
| `submit_one_array.sh` | Submit one (N, scenario, estimator, outer_opt) cell |
| `submit_all_arrays.sh` | Enumerate + submit the full grid |
| `README.md` | This file |

Companion R files (in `script/`):

| File | Purpose |
|---|---|
| `estimator_factory.R` | `make_est_control(est, outer_opt, tier)` dispatch + `valid_combos()` |
| `scm_bench_helpers.R` | Shared helpers (to_nm_dataset, diagnose_fit, package_scm_result, ...) |
| `PerformanceEvaluation_scm_bench.R` | Per-dataset driver (CLI) |
| `migrate_scm_output.R` | One-time relocation of legacy `output/output_N*/` tree |

## Output layout

All bench cells write to a single unified tree:

```
output/scm_bench/N{NN}/scn{SS}/{est}_{opt}/res_ds{DDD}.rds
```

Examples:

```
output/scm_bench/N80/scn16/focei_bobyqa/res_ds01.rds       (migrated from legacy)
output/scm_bench/N80/scn16/focei_nlminb/res_ds01.rds       (new)
output/scm_bench/N80/scn16/focei_lbfgsb3c/res_ds01.rds     (new)
output/scm_bench/N80/scn16/foceif_nlminb/res_ds01.rds
output/scm_bench/N80/scn16/foceif_lbfgsb3c/res_ds01.rds
output/scm_bench/N80/scn16/irlsfoceif_lbfgsb3c/res_ds01.rds
output/scm_bench/N80/scn16/saem_NA/res_ds01.rds
```

Each `res_ds{DDD}.rds` contains a named list with 5 coordinate keys
(`sample_N`, `scenario_id`, `dataset_id`, `estimator`, `outer_opt`),
provenance (`seed`, `timestamp`, `resumed`), timings (`t_base_sec`,
`t_scm_sec`, `t_refit_sec`, `t_total_sec`), and the packaged SCM result
in `$test`.

## Method grid

| Estimator | Outer optimizers | Notes |
|---|---|---|
| focei | bobyqa (migrated), nlminb, lbfgsb3c | Baseline; bobyqa preserved from legacy runs |
| foceif | nlminb, lbfgsb3c | Almquist analytic outer gradient (`fast=TRUE`) |
| irlsfoceif | lbfgsb3c | IRLS inner + fast outer (mu-referenced models only) |
| saem | NA | Post-SAEM `foceif` refit for covariance (`refit_estimator="foceif"`) |

`valid_combos()` in `estimator_factory.R` is the source of truth.

**VAE dropped**: `runSCM()` cannot locate `cl`/`vc` symbols in the ui returned
by `nlmixr2est::vae()` (labels are rewritten to `lTVCL`/`lTVVc`). All 16
step-1 candidates errored with *"'cl' has not been found in the model ui"*.
See `output/scm_bench/N80/scn16/vae_NA/ds01_ERROR.txt` for the reference
failure. VAE remains available for standalone base fits.

## Cell count

Full new-cell grid (per (N, scenario)): **6 cells**:

```
focei_nlminb, focei_lbfgsb3c,
foceif_nlminb, foceif_lbfgsb3c,
irlsfoceif_lbfgsb3c, saem_NA
```

Plus `focei_bobyqa` from migration for analysis (not re-run) = **7 cells
total** in downstream aggregation.

## Prerequisites (one-time)

### HPCE R environment

Install dev-branch nlmixr2 stack in the correct dependency order (see main
project README for details):

```r
remotes::install_github("nlmixr2/lotri",       upgrade = "never")
remotes::install_github("nlmixr2/rxode2ll",    upgrade = "never")
remotes::install_github("nlmixr2/rxode2",      upgrade = "never")
remotes::install_github("nlmixr2/nlmixr2data", upgrade = "never")
remotes::install_github("nlmixr2/nlmixr2est",  upgrade = "never")
remotes::install_github("nlmixr2/nlmixr2",     upgrade = "never")

```

Verify (should include foceif, irlsfoceif, mufocei):

```r
m <- as.character(utils::methods("nlmixr2Est"))
sort(sub("^nlmixr2Est\\.", "", m))
```

### Migrate legacy focei_bobyqa outputs

Runs once, non-destructive (copy, not move). Backfills coordinate keys.

```bash
cd ~/nlmixr2scm
Rscript script/migrate_scm_output.R
ls output/scm_bench/N80/scn16/focei_bobyqa/ | wc -l   # sanity check
```

## Pilot: scn16, N=80, 10 datasets, 6 new cells

6 cells x 10 datasets = 60 fits.

```bash
cd ~/nlmixr2scm
git pull

NS=80 SCENARIOS=16 bash script/hpce_bench/submit_all_arrays.sh 10 20

# Track
bjobs -J 'bench_*'
```

Runtime budget on HPCE (rough, from N=80 ds01 smoke test):
- focei_bobyqa (legacy reference): ~19 min / fit
- focei_{nlminb,lbfgsb3c}: ~20-30 min / fit
- foceif_{nlminb,lbfgsb3c}: ~20-25 min / fit
- irlsfoceif_lbfgsb3c: ~10-15 min / fit (fastest)
- saem_NA (+ foceif refit): ~25-30 min / fit

With 20 tasks in parallel, the pilot completes in ~1-2 hours.

## Scale-up (full grid)

Full grid = 3 N x 16 scenarios x 6 new cells = 288 arrays, 100 datasets
each = 28,800 new fits (plus the ~4,800 migrated focei_bobyqa results for
reference). Do NOT launch until the pilot has validated:

- All 6 cells produce parseable `res_ds*.rds`
- Runtimes fit inside the 3-hour LSF walltime

Then:

```bash
bash script/hpce_bench/submit_all_arrays.sh 100 50
```

## Resume / re-run behaviour

The driver has a 2-tier cache (mirrors `refit_one_dataset.R`):

- **Tier 1**: `res_ds{DDD}.rds` exists -> skip, return cached
- **Tier 2**: `scm_ds{DDD}.rds` exists but `res_*` missing (or `--force_repackage`)
  -> re-run only `package_scm_result` (~30-60 s for the refit)
- **Tier 3**: neither cache present -> full pipeline

Flags:
- `--force_rerun` bypasses both caches
- `--force_repackage` ignores stale `res_*` but reuses cached `scm_*`

## Aggregation (next step)

To be drafted as `script/aggregate_bench_results.R`. Will discover files
under `output/scm_bench/N*/scn*/*_*/res_ds*.rds`, parse the 4 path
coordinates (backstop for missing in-RDS keys), and produce per-cell
tables of Power, PowerCN, PowerMinSuc, RMRSE, and runtime.

## Manual smoke test (single fit, local Windows)

```r
setwd("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
source("script/PerformanceEvaluation_scm_bench.R")
run_bench_cell(list(
  N = 80, scenario = 16, dataset = 1,
  estimator = "foceif", outer_opt = "lbfgsb3c",
  force_rerun = FALSE, force_repackage = FALSE,
  out_root = "output/scm_bench",
  input_root = "Inputdataset",
  true_params_path = "Inputdataset/true_params_long.rds"
))
```
