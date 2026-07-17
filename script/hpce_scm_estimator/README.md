# HPCE SCM estimator × optimizer benchmark

Stepwise covariate modelling (SCM) operating-characteristics benchmark
comparing **estimation methods × outer optimizers × model structures** on the
Khandelwal 2019 PopPK model. Unlike the VAE pilot, the search here **is**
`runSCM()` (classical forward/backward LRT), with a two-tier control to keep it
affordable:

- **screen tier** — cheap control for the base fit and every SCM candidate LRT
  (no covariance, no residual tables);
- **final tier** — one tight-tolerance `r,s` covariance refit of the SCM winner.

All cells emit **schema-2.1** records via `package_scm_schema21()`.

## Scope / grid

| axis | values | count |
|------|--------|-------|
| cohort `N` | 40, 80, 300 | 3 |
| scenario | 1 … 16 | 16 |
| structure | linCmt, ode | 2 |
| estimator × outer_opt | see grid below | 6 |
| datasets / cell | 1 … `NDS` | tunable |

- One **LSF array per `(N, scenario, structure, estimator, outer_opt)` cell**;
  the array index (`$LSB_JOBINDEX`) **is the dataset id**.
- Full grid = 3 × 16 × 2 × 6 = **576 arrays** × `NDS` datasets.

### Method grid (`valid_combos()` is the source of truth)

| Estimator | Outer optimizers | Notes |
|---|---|---|
| focei | bobyqa, nlminb, lbfgsb3c | Baseline; `focei_bobyqa` re-run under schema 2.1 |
| foceif | nlminb, lbfgsb3c | Almquist analytic outer gradient (`fast=TRUE`) |
| irlsfoceif | lbfgsb3c | IRLS inner + fast outer (mu-referenced models only) |

**saem dropped** (2026-07-17): removed from `valid_combos()` and the job grid.

**VAE dropped**: `runSCM()` cannot locate `cl`/`vc` symbols in the ui returned
by `nlmixr2est::vae()` (labels rewritten to `lTVCL`/`lTVVc`). VAE covariate
selection has its own driver (`script/vae_covsel_driver.R` + `hpce_vae_covsel/`).

### Structure note

`foceif`/`irlsfoceif` need an ODE to interaction-linearise. On `linCmt` (analytic
solution) they **degrade to `focei`** internally; they are still run there to
document *"prefer linCmt when possible"*. The `ode` structure unlocks their
analytic gradients.

### Thread budget & timing (fair speed comparison)

Every cell uses the **same, fixed** parallelism policy so wall-times are
comparable across estimators and every run is bitwise reproducible. The driver
pins it once at the top of `run_bench_cell()`:

| knob | value | effect |
|------|-------|--------|
| SCM candidate workers | `workers = 3L` (hardcoded) | `runSCM()` **forks** 3 child processes for the candidate LRT fits |
| rxode2 ODE threads    | `rxThreads = 1` (pinned) | `rxode2::setRxThreads(1L)` + `OMP_NUM_THREADS=1` — OpenMP **off** |

**Why rxThreads is pinned to 1.** `runSCM(workers=3)` parallelises the SCM by
`fork()`ing child processes. rxode2's ODE solver uses OpenMP. Calling `fork()`
over a live OpenMP thread pool is **undefined behaviour** — it silently
corrupts the ODE solve in the children, producing non-deterministic, *wrong*
fits (frozen-at-init OFVs, junk LRT deltas, unstable covariate selection). We
therefore keep OpenMP off (`rxThreads=1`) so the fork is clean. Parallelism
comes from the fork (`workers`) and, across datasets, from the LSF array width.

> Never combine `workers > 1` **and** `rxThreads > 1`. Safe configurations are
> `workers=3, rxThreads=1` (fork only — what we use) or `workers=1,
> rxThreads>1` (OpenMP only). The forbidden combination is what caused earlier
> non-reproducible focei results.

The pinned values are recorded in every record and `meta.json` as
`scm_workers` (3) and `rx_threads` (1).

**Both wall-time and CPU-time are recorded.** CPU is measured with
`proc.time()` including the `user.child`/`sys.child` columns, so the forked
workers' CPU (invisible to wall-clock) is captured. The ratio
**`hog_factor = cpu_total / wall_total`** reflects the fork parallelism during
the candidate-fit phase (CPU can exceed wall there); serial phases (outer
optimizer, covariance, packaging) pull the overall factor back toward ~1.

## Files

### Submission (this directory)

| file | role |
|------|------|
| `bench_array.lsf`      | LSF array template; reads `SAMPLE_N/SCN/STRUCTURE/EST/OPT` from env, calls `PerformanceEvaluation_scm_bench.R --dataset $LSB_JOBINDEX`. |
| `submit_one_array.sh`  | Submit ONE cell: `bsub -J "bench_N<N>_scn<SS>_<struct>_<est>_<opt>[1-NDS]%MAXPAR"`. |
| `submit_all_arrays.sh` | Nested loop over `STRUCTURES × NS × SCENARIOS × valid(EST,OPT)`. |
| `README.md`            | This file. |

`submit_one_array.sh` positional args:

```
submit_one_array.sh <N> <SCN> <EST> <OPT> <STRUCTURE> <NDS> [MAXPAR] [DS_START]
```

### R scripts (under `script/`, required at runtime)

| file | role |
|------|------|
| `PerformanceEvaluation_scm_bench.R` | **Entry point** run by each array task. Loads the cell dataset, runs the two-tier SCM, packages the schema-2.1 record. |
| `estimator_factory.R` | `make_est_control(est, outer_opt, tier)` dispatch, `valid_combos()`, `nlmixr_est_name()`, `seed_*()`. |
| `scm_bench_helpers.R`  | Base models `base_2cmt_oral_linCmt` / `base_2cmt_oral_ode` (each `attr(,"structure")`-stamped), `runSCM_traced()`, `package_scm_schema21()`. |
| `refit_helpers.R`      | `to_nm_dataset()`, `PsN_scenarios`, tight-tol refit control. |
| `output_schema.R`      | `assemble_common()` (schema 2.1) + `write_fit_sidecar()`. |

The driver `source()`s the helpers automatically, so only
`PerformanceEvaluation_scm_bench.R` is named on the command line.

### R package prerequisites

Dev-branch nlmixr2 stack (install in dependency order):

```r
remotes::install_github("nlmixr2/lotri",       upgrade = "never")
remotes::install_github("nlmixr2/rxode2ll",    upgrade = "never")
remotes::install_github("nlmixr2/rxode2",      upgrade = "never")
remotes::install_github("nlmixr2/nlmixr2data", upgrade = "never")
remotes::install_github("nlmixr2/nlmixr2est",  upgrade = "never")
remotes::install_github("nlmixr2/nlmixr2",     upgrade = "never")
```

Verify the extra estimator methods are present:

```r
m <- as.character(utils::methods("nlmixr2Est"))
sort(sub("^nlmixr2Est\\.", "", m))   # expect foceif, irlsfoceif, mufocei
```

## Output layout (schema 2.1)

```
output/scm_bench/N<N>/scn<SS>_<structure>/<est>_<opt>/
    res_ds<DDD>.rds        # full record (rel_err, diag, diag_t3, cov, $scm block)
    res_ds<DDD>.fit.rds    # raw refit fit
    res_ds<DDD>.meta.json  # greppable scalar manifest
    scm_ds<DDD>.rds        # cached SCM search result (tier-2 resume)
    res_ds<DDD>_ERROR.txt  # only if the task failed
```

Examples:

```
output/scm_bench/N80/scn16_linCmt/focei_bobyqa/res_ds001.rds
output/scm_bench/N80/scn16_ode/foceif_lbfgsb3c/res_ds001.rds
output/scm_bench/N300/scn16_ode/irlsfoceif_lbfgsb3c/res_ds001.rds
```

Each record carries the coordinate keys (`sample_N`, `scenario_id`,
`dataset_id`, `estimator`, `outer_opt`), `model_type`/`structure`, the SCM
`$scm` block (`selected`, `step_hist`, `cov_done`), `rel_err`, the two
diagnostic blocks (`diag`, `diag_t3`), and **two parallel timing blocks**:

- `$runtime` — WALL-clock seconds per phase: `base_sec`, `scm_sec`,
  `refit_sec`, `total_sec`.
- `$cpu` — CPU seconds (self + forked-child) per phase plus
  `hog_factor = cpu$total_sec / runtime$total_sec` (the parallelism benefit).

The flat scalars `wall_total_sec`, `cpu_total_sec`, `hog_factor` are also
mirrored into `res_ds<DDD>.meta.json` so timing/benefit is greppable without
opening the `.rds`:

```bash
grep -h hog_factor output/scm_bench/**/res_ds*.meta.json
```

> **Legacy tree** `output/full_scm_focei_bobyqa/` (old schema) is left
> untouched. New schema-2.1 `focei_bobyqa` records land in the new
> `scn<SS>_<structure>` layout.

## Recommended workflow: probe first, then full sweep

### 1. Single probe (N80 × scn16 × ode × focei_bobyqa × few ds)

```bash
cd ~/nlmixr2scm && git pull
# normalize line endings first (edited on Windows)
sed -i 's/\r$//' script/hpce_scm_estimator/*.sh script/hpce_scm_estimator/*.lsf script/*.R

STRUCTURES=ode NS=80 SCENARIOS=16 ESTIMATORS=focei FOCEI_OPTS=bobyqa \
  bash script/hpce_scm_estimator/submit_all_arrays.sh 5 20 --force_rerun
# equivalently:
# bash script/hpce_scm_estimator/submit_one_array.sh 80 16 focei bobyqa ode 5 20
```

Monitor and read actual resource usage:

```bash
bjobs -J 'bench_*'                  # running/pending
bjobs -l <jobid>                    # live per-task detail
bacct -l <jobid>                    # peak MEM / CPU / walltime AFTER completion
ls output/scm_bench/N80/scn16_ode/focei_bobyqa/     # res_ds00{1..5}.*
```

Check for failures:

```bash
find output/scm_bench -name '*_ERROR.txt' -exec sed -n '1,6p' {} +
```

### 2. Tune resources

From `bacct -l`, read **MAX MEM** and **CPU TIME**; adjust `#BSUB -M` and
`#BSUB -W` in `bench_array.lsf` if the probe values are off. ODE + SCM is the
slow corner — size the walltime off an `ode × focei` cell.

### 3. Full sweep

```bash
# both structures, all cohorts, all 16 scenarios, all valid cells, 100 ds
bash script/hpce_scm_estimator/submit_all_arrays.sh 100 50

# or slice it — env filters are space-separated lists:
STRUCTURES=linCmt NS="80 300" SCENARIOS=16 \
  bash script/hpce_scm_estimator/submit_all_arrays.sh 100 50
```

`submit_all_arrays.sh` env knobs (all optional):

| var | default | meaning |
|-----|---------|---------|
| `STRUCTURES` | `linCmt ode` | structure sweep |
| `NS`         | `40 80 300`  | cohort sizes |
| `SCENARIOS`  | `1 … 16`     | scenario ids |
| `ESTIMATORS` | `focei foceif irlsfoceif` | estimators |
| `FOCEI_OPTS` | `bobyqa nlminb lbfgsb3c`  | focei outer optimizers |
| `DS_START`   | `1`          | first dataset index |

Positional args: `submit_all_arrays.sh <n_datasets> [maxpar]`.

Runtime budget (rough, from N=80 ds01 smoke test on linCmt):
- `focei_bobyqa`: ~19 min / fit
- `focei_{nlminb,lbfgsb3c}`: ~20–30 min / fit
- `foceif_{nlminb,lbfgsb3c}`: ~20–25 min / fit
- `irlsfoceif_lbfgsb3c`: ~10–15 min / fit (fastest)

ODE cells run longer than the linCmt equivalents.

## Resume / re-run behaviour

The driver has a 2-tier cache:

- **Tier 1**: `res_ds{DDD}.rds` exists → skip, return cached.
- **Tier 2**: `scm_ds{DDD}.rds` exists but `res_*` missing (or `--force_repackage`)
  → re-run only the final-tier packaging (covariance refit, ~30–60 s).
- **Tier 3**: neither cache present → full pipeline.

Flags:
- `--force_rerun` bypasses both caches;
- `--force_repackage` ignores stale `res_*` but reuses cached `scm_*`.

## Aggregation (next step — deferred)

To be drafted (mirrors `script/aggregate_vae_covsel.R`). It will discover files
under `output/scm_bench/N*/scn*_*/*_*/res_ds*.rds`, parse the path coordinates
(`scn<SS>_<structure>` → scenario + structure) as a backstop for missing in-RDS
keys, and roll `rec$scm$*` / `rec$rel_err` / `rec$runtime$*` / `rec$cpu$*` up
into per-cell tables of Power, PowerCN, PowerMinSuc, RMRSE, wall/CPU runtime,
and the `hog_factor` parallelism benefit.

## Manual smoke test (single fit, local Windows)

```r
setwd("C:/Users/LIUYA8J/OneDrive - Novartis Pharma AG/Internship/GithubRepo/nlmixr2scm")
source("script/PerformanceEvaluation_scm_bench.R")
run_bench_cell(list(
  N = 80, scenario = 16, dataset = 1,
  estimator = "foceif", outer_opt = "lbfgsb3c",
  structure = "ode",
  force_rerun = FALSE, force_repackage = FALSE,
  out_root = "output/scm_bench",
  input_root = "Inputdataset",
  true_params_path = "Inputdataset/true_params_long.rds"
))
```

Or via the CLI:

```bash
Rscript script/PerformanceEvaluation_scm_bench.R \
  --N 80 --scenario 16 --dataset 1 \
  --estimator focei --outer_opt bobyqa --structure ode
```
