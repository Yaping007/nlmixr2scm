# HPCE: profileInit warm-start A/B (scn16)

Offload the `runSCM(profileInit=...)` A/B experiment to the cluster, where each
`(cell, arm, structure)` runs as an independent LSF job (concurrent, and each
task can request more cores than a laptop has).

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
bash script/hpce_profileInit/submit_one.sh irlsfoceif_lbfgsb3c profileOff ode
bash script/hpce_profileInit/submit_one.sh irlsfoceif_lbfgsb3c profileOn  ode

# request more muscle: 6 workers x 2 threads = 12 cores, 10h wall:
bash script/hpce_profileInit/submit_one.sh irlsfoceif_lbfgsb3c profileOn ode 6 2 12 10:00

# WHOLE grid for one structure (3 cells x 2 arms = 6 concurrent jobs):
bash script/hpce_profileInit/submit_grid.sh ode
bash script/hpce_profileInit/submit_grid.sh linCmt
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
