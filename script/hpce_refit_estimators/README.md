# hpce_refit_estimators

LSF job-array infrastructure for the **true-model refit** benchmark
(`script/bench_refit_estimators.R`) — fits the true model per
(N, scenario, estimator, outer_opt, boundary, dataset) with **no SCM search**.

Mirrors the layout of `script/hpce_bench/` but shorter/lighter per-task
resources (no SCM step, no post-refit chain).

## Files

- `refit_estimators_array.lsf` — LSF task template. `-n 2 -M 4000 -W 01:00`.
- `submit_one_array.sh` — submit one cell as an array over datasets.
- `submit_all_arrays.sh` — enumerate cells and submit all arrays.

## Output layout

```
output/bench_refit_estimators/N<NN>/scn<SS>/<est>_<opt>/res_ds<DDD>.rds
logs/refitE_N<NN>/refitE_N<NN>_scn<SS>_<est>_<opt>_<bnd>.<jobid>.<idx>.{out,err}
```

## Quick start

```bash
# scn16 N=80 tight, ds 1..10, 20 parallel, all 7 cells
NS=80 SCENARIOS=16 BOUNDARIES=tight \
  bash script/hpce_refit_estimators/submit_all_arrays.sh 10 20

# One cell only (foceif + lbfgsb3c, ds 1..10)
bash script/hpce_refit_estimators/submit_one_array.sh 80 16 foceif lbfgsb3c tight 10 20

# Extend to ds 11..50 without redoing ds 1..10
DS_START=11 NS=80 SCENARIOS=16 BOUNDARIES=tight \
  bash script/hpce_refit_estimators/submit_all_arrays.sh 40 20

# Track / kill
bjobs -J 'refitE_*'
bkill  -J 'refitE_*'
```

## Env vars

| Var | Default | Notes |
|---|---|---|
| `NS` | `40 80 300` | sample-size cohorts |
| `SCENARIOS` | `1..16` | scenario indices |
| `BOUNDARIES` | `tight` | any of `none wide narrow tight` |
| `ESTIMATORS` | `focei foceif irlsfoceif saem` | |
| `FOCEI_OPTS` | `nlminb lbfgsb3c bobyqa` | outer opts for focei row |
| `DS_START` | `1` | shift LSB_JOBINDEX for restart / extend runs |
