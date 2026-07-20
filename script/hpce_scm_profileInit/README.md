# HPCE profileInit A/B — scn16

Runs the `profileInit` warm-start A/B experiment (profileOff vs profileOn) for
one SCM cell on the cluster, one LSF array task per `(cell, arm, structure)`.

---

## Files to upload to HPCE

The cluster job needs the **package (reinstalled)** plus a specific set of
`script/` sources and two data inputs. Simplest: sync the whole repo and
reinstall. If cherry-picking, upload exactly the following.

### 1. The package itself — MUST be reinstalled

The profileInit warm-start lives in `R/scm.R` (`.profileCovInit`,
`.freezeUiForProfile`), which is **part of the installed package namespace** —
sourcing scripts is not enough. `nlmixr2scm` is pure R (no C++), so installing
needs no compiler / `remotes` / dependency chain. After syncing the repo:

```bash
# from the repo root on HPCE, AFTER the nlmixr2est 6.2.0 dev stack is in place
cd ~/nlmixr2scm/nlmixr2scm
R CMD INSTALL --with-keep.source .

# GATE: confirm profileInit is a real argument BEFORE submitting jobs
Rscript -e 'cat("profileInit present:", "profileInit" %in% names(formals(nlmixr2scm::runSCM)), "\n")'
```

Submit only when the gate prints `TRUE`. Skipping the reinstall makes HPCE run
the OLD runSCM (error: `unused argument (profileInit = FALSE)`).

### 2. Scripts sourced by the driver (all under `script/`)

Traced from the `source()` block in `_ab_profileInit_scn16.R`:

| File | Provides |
|------|----------|
| `script/_ab_profileInit_scn16.R` | the driver (env-var aware) |
| `script/refit_helpers.R`         | `diagnose_fit_table3`, `rel_err_one` |
| `script/true_model_factory.R`    | `make_true_model(structure=)` |
| `script/bench_refit_estimators.R`| `to_nm_dataset`, `DOSE_MG` |
| `script/output_schema.R`         | `assemble_common`, `write_fit_sidecar` (schema 2.1) |
| `script/estimator_factory.R`     | `make_est_control`, `nlmixr_est_name` |
| `script/scm_bench_helpers.R`     | `base_2cmt_oral_ode/linCmt`, `runSCM_traced`, `package_scm_schema21` |

### 3. HPCE job files (under `script/hpce_scm_profileInit/`)

- `profileInit_array.lsf`
- `submit_one.sh`
- `submit_grid.sh`
- `README.md` (this file)

### 4. Data inputs

- `Inputdataset/sim_obs_N80/sim_obs_scenario_16.rds` — read at the exact path
- `Inputdataset/true_params_long.rds` — upload the cached `.rds`; the driver
  prefers an in-session `true_params`, else falls back to `true_params_long()`.
  Having the `.rds` avoids sourcing the 2800-line `PerformanceEvaluation04062026.R`.

### Software stack (already on HPCE)

- `nlmixr2est` **6.2.0** dev (required for `ifoceif` / analytic gradient)
- matching `rxode2` / `nlmixr2` dev stack
- Reinstall `nlmixr2scm` (step 1) AFTER the 6.2.0 stack so it links correctly.

---

## Environment variables (per array task)

## What runs

Driver: `script/_ab_profileInit_scn16.R`. It now reads its scope + parallelism
from environment variables (session var > env > default):

| env var        | meaning                                   | values / default |
|----------------|-------------------------------------------|------------------|
| `AB_CELL`      | which estimator cell                      | `irlsfoceif_lbfgsb3c` \| `focei_nlminb` \| `focei_lbfgsb3c` |
| `AB_ARM`       | warm-start off/on                         | `profileOff` \| `profileOn` |
| `AB_STRUCTURE` | model parameterization                    | `linCmt` \| `ode` (default `linCmt`) |
| `AB_WORKERS`   | future multisession workers (candidate fan-out) | default `3` |
| `AB_RXTHREADS` | rxode2 ODE threads per worker             | default `2` |

Keep `AB_WORKERS * AB_RXTHREADS <= -n cores` to avoid oversubscription.

## Output (schema 2.1)

```
output/ab_profileInit/scn16_<structure>/<cell>/<arm>/res_ds001.rds       # record
output/ab_profileInit/scn16_<structure>/<cell>/<arm>/res_ds001.fit.rds   # refit winner
output/ab_profileInit/scn16_<structure>/<cell>/<arm>/res_ds001.meta.json # manifest
```

The `linCmt` and `ode` runs write to separate roots, so they never overwrite.

## Submit

```bash
# ONE task (the heavy case that was too slow locally):
bash script/hpce_scm_profileInit/submit_one.sh irlsfoceif_lbfgsb3c profileOff ode
bash script/hpce_scm_profileInit/submit_one.sh irlsfoceif_lbfgsb3c profileOn  ode

# request more muscle: 6 workers x 2 threads = 12 cores, 10h wall:
bash script/hpce_scm_profileInit/submit_one.sh irlsfoceif_lbfgsb3c profileOn ode 6 2 12 10:00

# WHOLE grid for one structure (3 cells x 2 arms = 6 concurrent jobs):
bash script/hpce_scm_profileInit/submit_grid.sh ode
bash script/hpce_scm_profileInit/submit_grid.sh linCmt
```

`submit_one.sh <cell> <arm> <structure> [workers] [rxthreads] [cores] [walltime]`

## Monitor / logs

```bash
bjobs -w
tail -f logs/profileInit/ab_scn16_ode_irlsfoceif_lbfgsb3c_profileOn.<jobid>.out
```

## After it finishes

Read a record straight from R:

```r
on  <- readRDS("output/ab_profileInit/scn16_ode/irlsfoceif_lbfgsb3c/profileOn/res_ds001.rds")
off <- readRDS("output/ab_profileInit/scn16_ode/irlsfoceif_lbfgsb3c/profileOff/res_ds001.rds")
on$scm$selected; on$objf; on$cond_num_cor
```

## Resource notes

- `profileInit_array.lsf` defaults: `-n 6`, `-M 8000`, `-W 06:00`.
- ODE `irlsfoceif` is the heavy case (numerical integration + analytic-gradient
  compile per candidate); `linCmt` is much faster.
- To go bigger, pass `[workers] [rxthreads] [cores] [walltime]` -- e.g. a 16-core
  host: `... 8 2 16 12:00`.
