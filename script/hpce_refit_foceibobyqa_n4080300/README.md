# HPCE job-array templates for the refit-true-model study

This directory contains LSF submission templates for the Figure 3 / Table 3
replication of Khandelwal et al. 2019.  Each submission fits the "true" model
(known covariate structure per scenario) to one simulated dataset under one
of three boundary configurations.

## Study design

|                        | count | comment |
|------------------------|:-----:|---------|
| Cohorts (sample sizes) | 3     | `N40`, `N80`, `N300` |
| Scenarios              | 15    | 1 = no covariates (skipped); 2 = SEX only (single boundary); 3-16 = at least one continuous covariate (three boundaries each) |
| Boundary configs       | 3     | `none` (bare init), `wide` = c(-1e5, 1e5), `narrow` = c(-10, 10) |
| Datasets per combo     | 10 or 250 | 10 for pilot, 250 for full run |

- **Combos per cohort**: `1 + 14 * 3 = 43`
- **Total arrays**:      `3 * 43 = 129`
- **Pilot fits**:        `129 * 10  = 1,290`
- **Full-run fits**:     `129 * 250 = 32,250`

## Files

| file | purpose |
|------|---------|
| `refit_array.lsf`     | LSF job-array template (one task = one fit).  Not submitted directly; consumed by `submit_one_array.sh`. |
| `submit_one_array.sh` | Submit one array for a specific (cohort, scenario, boundary). |
| `submit_all_arrays.sh`| Loop over 129 (cohort, scenario, boundary) triples and submit one array per triple. |
| `pilot_smoke.sh`      | Convenience wrapper: `submit_all_arrays.sh 10 50` (all 3 cohorts, 10 datasets, 50 parallel). |

## Typical workflow

In the workspace these files live under `scripts/hpce/`.  You can rename or
relocate the whole directory on your cluster (for example `script/hpce_refit/`)
— the scripts auto-detect their own location.  All commands below use
`<hpce_dir>/` as a placeholder for wherever the directory lives on your
cluster.

If the scripts are NOT exactly 2 levels below your repo root, export
`REPO_ROOT` explicitly, e.g.:

```bash
export REPO_ROOT=/home/liuya8j/nlmixr2scm/nlmixr2scm
```

```bash
# 0. From the workspace root, make everything executable
chmod +x <hpce_dir>/*.sh <hpce_dir>/*.lsf

# 1. Pilot: 3 cohorts x 43 combos x 10 datasets = 1,290 fits (~90 min wall)
bash <hpce_dir>/pilot_smoke.sh

# 2. Track queue state
bjobs -J 'refit_*'
# ... or per cohort:
bjobs -J 'refit_N300_*'

# 3. When the pilot is green, launch the full run
bash <hpce_dir>/submit_all_arrays.sh 250 100

# 4. Cancel everything if needed
bkill -J 'refit_*'
```

## Output layout

```
outputs/refit_true_N40 /scn02_none  /res_ds001.rds ... res_dsNNN.rds
outputs/refit_true_N40 /scn03_none  /res_ds001.rds ...
outputs/refit_true_N40 /scn03_wide  /res_ds001.rds ...
outputs/refit_true_N40 /scn03_narrow/res_ds001.rds ...
                       ...
outputs/refit_true_N80 /scn02_none  /res_ds001.rds ...
                       ...
outputs/refit_true_N300/scn16_narrow/res_ds001.rds ... res_ds250.rds
```

Each `res_dsDDD.rds` is a slim named list with 15 fields (see
`scripts/refit_one_dataset.R` for full schema):

- `cohort`, `scenario_id`, `dataset_id`, `boundary`, `status`
- `final_est`  - 13-row tibble of estimates
- `rel_err`    - 13-row tibble with `estimate`, `true_value`, `rel_err`, ...
- `diag`       - Table-3 diagnostics (`min_suc`, `est_bnd`, `rnd_err`,
                 `zero_grad`, `cov_step`, **`phys_bnd`**, `cond_num_cor`,
                 `bound_hits`, `phys_hits`, `objf`, `message`)
- `parFixed`   - nlmixr2 `parFixedDf` (SE / %RSE / CI)
- `bounds_spec`, `cont_params`, `cat_params`, `fn_text`
- `runtime_sec`, `timestamp`

## Log files

Per-task stdout / stderr are written to:

```
logs/refit_<cohort>/refit_<cohort>_scn<NN>_<bnd>.<jobid>.<idx>.{out,err}
```

## Failure sentinels

If a fit throws inside `tryCatch()`, the driver writes a sentinel file
`res_dsDDD_ERROR.txt` alongside the RDS (which still contains
`status = "error"` and `error_msg`).  Aggregation should treat presence of
the sentinel as a failure regardless of the RDS.

## Resource assumptions (defaults)

| resource        | value    | rationale |
|-----------------|----------|-----------|
| cores per task  | 1        | nlmixr2 focei is serial per fit; parallelism = array width |
| memory per task | 4000 MB  | peak observed ~1.8 GB for N=300 |
| walltime        | 1 h      | typical N=300 fit ~5 min; margin for pathological trajectories in `none`/`wide` |
| max parallel    | 50 (pilot), 100 (full) | pass as 5th arg to `submit_one_array.sh` |

Adjust the `#BSUB` lines in `refit_array.lsf` if your queue policies differ.

## Pathology diagnostics

Because nlmixr2 focei has no equivalent of NONMEM's `SIGL` significant-digit
floor, unconstrained fits (`boundary in {none, wide}`) can lock onto a
degenerate MLE where a covariate theta escapes to O(100) so the proportional-
error sigma^2 collapses to zero and `-log(sigma^2) -> +Inf`.  This is exactly
the failure mode that the paper's `narrow` boundary is designed to prevent.

The driver flags these runs via `diag$phys_bnd = TRUE` (any covariate |theta|
> 10) with the offending parameters listed in `diag$phys_hits`.  Table 3 will
report **PhysBnd %** as an additional diagnostic column alongside MinSuc,
EstBnd, RndErr, ZeroGrad, and CovStep.
