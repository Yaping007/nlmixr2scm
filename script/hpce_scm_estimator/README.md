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
| focei | **bobyqa** (active), nlminb, lbfgsb3c | Baseline; `focei_bobyqa` re-run under schema 2.1 |
| irlsfocei | **bobyqa** (active) | IRLS mu-referenced FOCEi, **non-fast** `est="ifocei"` — derivative-free outer, no analytic-gradient build. The IRLS counterpart to `focei_bobyqa` |
| irlsfoceif | lbfgsb3c (parked) | Fast IRLS mu-referenced FOCEi (`est="ifoceif"`, analytic outer gradient) |
| foceif | nlminb, lbfgsb3c (parked) | Almquist analytic outer gradient (`fast=TRUE`) |

**Active comparison (2026-07-23):** the default sweep runs only `focei_bobyqa`
vs `irlsfocei_bobyqa` to isolate whether the IRLS mu-profiling speeds up the
SCM search under the same derivative-free outer optimiser. All other combos are
valid but parked (enable via the `ESTIMATORS`/`*_OPTS` env vars).

**saem dropped** (2026-07-17): removed from `valid_combos()` and the job grid.

**VAE dropped**: `runSCM()` cannot locate `cl`/`vc` symbols in the ui returned
by `nlmixr2est::vae()` (labels rewritten to `lTVCL`/`lTVVc`). VAE covariate
selection has its own driver (`script/vae_covsel_driver.R` + `hpce_vae_covsel/`).

### Structure note

The analytic outer gradient only matters for the **parked** analytic-gradient
cells (`foceif`, `irlsfoceif_lbfgsb3c`): those need an ODE to
interaction-linearise and **degrade to a finite-difference gradient** on
`linCmt`. The **active** `irlsfocei_bobyqa` cell dispatches non-fast
`est="ifocei"` and uses no outer gradient, so structure only changes its base
model, not its optimiser.

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

### Measuring the parallel benefit correctly (read this)

Every record stores both a WALL-clock and a CPU timing block, but **the
in-record CPU numbers do NOT prove the fork speedup** — for two reasons:

1. **`proc.time()` cannot see the fork workers.** `runSCM(workers=3)` forks
   child processes; `proc.time()`'s `user.child`/`sys.child` columns are only
   populated for children reaped via `wait()`, which the fork pool does not do.
   So the recorded `cpu$scm_sec` misses almost all of the parallel work.
2. **The driver never passes `base_cpu`/`scm_cpu`** to `package_scm_schema21()`,
   so `cpu$base_sec`/`cpu$scm_sec` are `NA` and `cpu_total_sec` is essentially
   just the refit CPU.

Consequently the in-record `hog_factor = cpu_total / wall_total` comes out
**≈ 0.03–0.25** — an artefact of unmeasured CPU, *not* a slowdown. **Do not use
it to claim (or deny) a parallelism benefit.**

**The authoritative CPU measure is LSF.** LSF accounts CPU for the whole
process tree (parent + all fork workers + OpenMP threads) and writes it into
every task's `-o .out` log (`Resource usage summary → CPU time` / `Run time`);
`bacct -l <jobid>` reports the same. Parse it with
`script/aggregate_lsf_cpu.R`.

**Why even the LSF whole-job `hog` is < 1.** An SCM task is
`base fit (serial) → SCM search (3 forks, parallel) → refit (serial) →
covariance (serial)`. The fork only accelerates the SCM-search phase, which is
a minority of total wall time, so the *whole-job* CPU/wall averages below 1.
The fork genuinely works — it is simply diluted by the serial phases and by
non-compute wall time (model compilation, disk I/O, host contention).

**The clean demonstration is a phase-isolated A/B:** run identical datasets with
`workers=1` vs `workers=3` (both `rxThreads=1`) on **exclusive nodes** and
compare the SCM-phase wall time directly:

```
SCM-phase speedup = scm_sec(workers=1) / scm_sec(workers=3)
```

Wall time is only trustworthy when cores are not shared, so run the A/B with
`EXCLUSIVE=1` (and `CHAIN=1` / small `MAXPAR`) — see *Contention-free timing*
below.

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
module load R          # same module the jobs use
Rscript -e 'install.packages(".", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))' #reinstall nlmixr2scm


```

Verify the extra estimator methods are present:

```r
m <- as.character(utils::methods("nlmixr2Est"))
  sort(sub("^nlmixr2Est\\.", "", m))   # expect foceif, ifocei, ifoceif, mufocei
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
  `refit_sec`, `total_sec`. **These are the usable timing numbers**; `scm_sec`
  is the phase to compare for the fork speedup.
- `$cpu` — CPU seconds per phase plus
  `hog_factor = cpu$total_sec / runtime$total_sec`. **Caveat:** `proc.time()`
  does not capture the fork workers' CPU and the driver omits `base_cpu`/
  `scm_cpu`, so `cpu$base_sec`/`cpu$scm_sec` are `NA` and `hog_factor` (≈0.03–
  0.25) is a measurement artefact — use LSF CPU (`aggregate_lsf_cpu.R`) for the
  real figure, not this.

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
  bash script/hpce_scm_estimator/submit_all_arrays.sh 5 20
# equivalently:
# bash script/hpce_scm_estimator/submit_one_array.sh 80 16 focei bobyqa ode 5 20
bkill 261294 261295 261296 261297 261298 261299 261300 261301 261302 261303 261304 261305 261306 261307 261308 261309

bash script/hpce_scm_estimator/submit_one_array.sh 40  2 focei bobyqa linCmt 1 1 212 #retun failed runs
bash script/hpce_scm_estimator/submit_one_array.sh 80  9 focei bobyqa linCmt 1 1  13


bash script/hpce_scm_estimator/submit_one_array.sh 300  16 focei bobyqa ode 1 1 

OUT_ROOT=output/scm_bench_rescue FORCE_RERUN=1 JOBTAG=rescue   bash script/hpce_scm_estimator/submit_one_array.sh 300 16 focei bobyqa ode 1 1


R CMD INSTALL --no-multiarch --with-keep.source .
sed -i 's/\r$//' script/hpce_scm_estimator/*.sh script/hpce_scm_estimator/*.lsf

FORCE_RERUN=1 \
OUT_ROOT=output/scm_bench_rescue_winner \
STRUCTURES=ode \
NS=300 \
SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16" \
  bash script/hpce_scm_estimator/submit_all_arrays.sh 5 30
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


FORCE_RERUN=1 \
OUT_ROOT=output/scm_bench \
STRUCTURES=ode \
NS="40 80 300" \
SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16" \
ESTIMATORS=focei \
FOCEI_OPTS="bobyqa lbfgsb3c" \
  bash script/hpce_scm_estimator/submit_all_arrays.sh 100 50


OUT_ROOT=output/scm_bench \
STRUCTURES=ode \
NS=300 \
SCENARIOS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16" \
ESTIMATORS=focei \
FOCEI_OPTS=bobyqa \
DS_START=6 \
  bash script/hpce_scm_estimator/submit_all_arrays.sh 100 50
```

`submit_all_arrays.sh` env knobs (all optional):

| var | default | meaning |
|-----|---------|---------|
| `STRUCTURES` | `linCmt ode` | structure sweep |
| `NS`         | `40 80 300`  | cohort sizes |
| `SCENARIOS`  | `1 … 16`     | scenario ids |
| `ESTIMATORS` | `focei irlsfocei` | estimators (active default) |
| `FOCEI_OPTS` | `bobyqa`     | focei outer optimizers |
| `IRLSFOCEI_OPTS`| `bobyqa`  | irlsfocei outer optimizers |
| `IRLS_OPTS`  | `lbfgsb3c`   | parked fast-IRLS (irlsfoceif) outer optimizers |
| `FOCEIF_OPTS`| `nlminb lbfgsb3c` | parked foceif outer optimizers |
| `DS_START`   | `1`          | first dataset index |

Positional args: `submit_all_arrays.sh <n_datasets> [maxpar]`.

Runtime budget (rough, from N=80 ds01 smoke test on linCmt):
- `focei_bobyqa`: ~19 min / fit
- `irlsfocei_bobyqa` (non-fast IRLS): comparable to `focei_bobyqa`; the IRLS mu-profiling is the variable under test
- parked cells: `focei_{nlminb,lbfgsb3c}` ~20–30 min, `foceif_{nlminb,lbfgsb3c}` ~20–25 min, `irlsfoceif_lbfgsb3c` ~10–15 min / fit
bash script/hpce_scm_estimator/submit_one_array.sh 300  16 focei nlminb ode 5 5 
bash script/hpce_scm_estimator/submit_one_array.sh 300  16 focei lbfgsb3c ode 5 5 

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

## Known bugs & fixes

### SCM Power collapse ~0.9 → ~0.6 — FIXED 2026-07-20 (`7007b17`)

`package_scm_schema21()` picked the SCM winner from the **forward** model
(`resFwd` first) instead of the **backward-eliminated** model (`resBck` first).
For a bidirectional `"scm"` search the final model is the backward one, so any
forward false-positive that backward elimination correctly *dropped* was still
counted in `$scm$selected`. In null/easy scenarios this inflated the
false-positive rate and collapsed Power from ~0.9 to ~0.6.

The search itself was never wrong — the `step_hist` showed the covariate as
`included = "dropped"` in both the old and new runs; only the winner-extraction
differed. Fix: `resBck` first, `resFwd` fallback (forward-only searches), matching
the original `package_scm_result()`. See
`docs/ESTIMATOR_JOURNEY.md` §3.10 for the full diagnosis (and the four
hypotheses that were disproven by controlled A/B first: fork corruption, screen
precision, `warm`, base-fit quality).

**Action required**: any `res_ds*.rds` produced before this commit has a wrong
`$scm$selected`. Re-derive them with `--force_repackage` (reuses cached
`scm_ds*.rds`, no re-fitting), then re-aggregate:

```bash
# repackage every cell that has a cached SCM object
for f in output/scm_bench/N*/scn*_*/*_*/scm_ds*.rds; do
  d=$(basename "$f"); ds=$(echo "$d" | sed 's/[^0-9]//g' | sed 's/^0*//')
  cell=$(dirname "$f"); IFS='_' read -r est opt <<< "$(basename "$cell")"
  scn_struct=$(basename "$(dirname "$cell")"); scn=$(echo "$scn_struct" | sed 's/scn0*\([0-9]*\)_.*/\1/'); struct=${scn_struct#*_}
  N=$(basename "$(dirname "$(dirname "$cell")")" | tr -d 'N')
  Rscript script/PerformanceEvaluation_scm_bench.R \
    --N "$N" --scenario "$scn" --estimator "$est" --outer_opt "$opt" \
    --structure "$struct" --out_root output/scm_bench --dataset "$ds" --force_repackage
done
Rscript script/aggregate_scm_estimator2.1.R --root output --sub scm_bench
```

### `irlsfoceif`/`foceif` analytic outer gradient never engages on ODE (`fast=TRUE` is inert) — INVESTIGATED 2026-07-23

> **RESOLVED for the active sweep (2026-07-23).** Adopted option 2 below for the
> IRLS cell: the active comparison uses a distinct `irlsfocei` grid label that
> dispatches the **non-fast** `est="ifocei"` (see `nlmixr_est_name()` in
> `estimator_factory.R`), which skips the wasted analytic augmented-model build
> entirely. Because the outer optimiser is derivative-free `bobyqa`, the
> analytic gradient was never going to be used anyway, so this is a strict win
> (no ~100 s setup, identical IRLS estimation). The investigation below is
> retained as the rationale and still applies to the **parked**
> analytic-gradient cells (`foceif`, `irlsfoceif_lbfgsb3c`).

**Symptom.** Every `foceif`/`irlsfoceif` ODE fit prints `grad: fd` and carries
the run-info note *"fast=TRUE: the analytic outer gradient could not be solved
at this point; using the finite-difference gradient for the affected
iteration(s)"*. The Almquist analytic outer gradient (the whole point of
`fast=TRUE`) therefore provides **no speedup** — `foceif`/`irlsfoceif` run as
finite-difference `focei` (with the IRLS inner step for `irlsfoceif`) plus a
wasted ~100 s symbolic augmented-gradient build.

**Not a corruption.** The selection is *correct*: on N300/scn16/ds002
`irlsfoceif_lbfgsb3c` recovers all four true covariates
(`CrCL_power_cl`, `SEX_vc`, `BW_power_cl`, `BW_power_vc`) with Cond#(Cor)=11,
while `focei_nlminb` on the same dataset mis-selects (`BW_lin`, `CrCL_lin`, the
`BMI` decoy), Cond#(Cor)=2218, `false convergence (8)`. The FD fallback is a
benign, *informational* note; the fit is sound.

**Root cause (traced in nlmixr2est 7.0.0 source, live).** The analytic core
`.foceiAnalyticGradCore()` returns `NULL` (→ FD fallback) on **every**
iteration. Instrumenting the `return(NULL)` gates pinned it to the `canVanish`
check:

```r
if (isTRUE(ef$canVanish)) {              # TRUE for a pure proportional error model
  .fa <- abs(E$f)
  if (any(!is.finite(.fa)) || min(.fa) < 1e-06 * max(.fa))
    return(NULL)                          # rejects -> finite-difference gradient
}
```

The 2-cmt **oral** model predicts concentration ≈ 0 during early absorption, and
every subject has a `TIME=0, DV=0` observation where the central compartment is
structurally 0 → `min|f| = 0`. With proportional error the residual variance at
`f=0` is `(prop·0)² = 0` (degenerate, infinite weight), so nlmixr2est
**refuses** the analytic gradient. This is a property of the *design*
(oral + proportional error + a t=0 record), not a tolerance or `sensMethod`
issue.

**What does NOT fix it (all tested live):**
- `sigdig` 5 → 6, or matching `atolSens/rtolSens` to the tight state tols
  ("Fix 1"): 0 % analytic either way; `fd_switch` stays `TRUE`. Fix 1's only
  proven benefit is on *warm, fully-specified* refits (≈5× faster, no switch),
  not on the ODE cells here.
- `sensMethod = forward|adjoint`, `covSolveTol`: the rejection happens in the
  `canVanish` gate *before* the sensitivity solve, so these knobs are inert.
- **Dropping the `TIME=0, DV=0` rows**: a short 15-iteration probe showed no
  gate firing, but a *full* run still ends `grad: fd` (the gate re-fires at
  later θ where oral troughs dip near zero), and deleting 300/1800 observations
  shifts the OBJF (−16868 → −6067) and the likelihood/LRT — so it changes the
  science and is **not** a valid fix.

**Conclusion / action.** On this proportional-error oral-ODE design the analytic
outer gradient cannot engage in nlmixr2est 7.0.0. Practical options:
1. **Report as a limitation** — `foceif`/`irlsfoceif` give no gradient speedup
   here and are numerically equivalent to FD `focei` (IRLS inner for
   `irlsfoceif`). No code change; honest and defensible.
2. **Set `fast=FALSE`** explicitly for the ODE cells to skip the wasted
   augmented-model build (byte-identical result, faster, honest `grad: fd`
   header).
3. **Upstream reprex** to nlmixr2est: the `canVanish` gate rejects on every
   iteration for oral + proportional-error models because early-absorption /
   t=0 predictions vanish — the only route to a genuine fix.

Diagnostic harnesses: `script/_test_fix1_irlsfoceif_scn16.R` (Fix 1 A/B) and
`script/diagnose_scm_gradient_trace.R` (candidate-level gradient trace).

## Diagnostic A/B knobs (screen precision & warm-start)

These were added while hunting the Power collapse and are kept as diagnostics.
**Defaults reproduce prior behaviour**, so they are inert unless set.

| driver flag | env (submit_one_array) | default | effect |
|-------------|------------------------|---------|--------|
| `--screen_sigdig` | `SCREEN_SIGDIG` | `NA` (→ 4) | screen-tier `foceiControl(sigdig=)` override |
| `--screen_atol`   | `SCREEN_ATOL`   | `NA` | screen-tier ODE `atol` override |
| `--screen_rtol`   | `SCREEN_RTOL`   | `NA` | screen-tier ODE `rtol` override |
| `--warm`          | `WARM`          | `calc` | `foceiControl(warm=)` — `save` = classic self-initialized inner Hessian |

Example A/B (original coarse screening vs current):

```bash
SCREEN_SIGDIG=3 SCREEN_ATOL=1e-6 SCREEN_RTOL=1e-4 \
  OUT_ROOT=output/scm_pilot_sd3 JOBTAG=sd3 \
  bash script/hpce_scm_estimator/submit_one_array.sh 300 1 focei bobyqa linCmt 10 10
```

## Contention-free timing (fork-parallelism benefit)

The bulk sweep runs many tasks per node (`%MAXPAR` shared hosts), so its wall
times are **contaminated by host contention** and understate the fork speedup.
To measure the benefit cleanly, run a small A/B on isolated nodes.

`submit_all_arrays.sh` / `submit_one_array.sh` expose two isolation knobs:

| env | effect |
|-----|--------|
| `EXCLUSIVE=1` | adds `-x` \u2014 each task owns its whole node, so the 3 fork workers are never starved for cores (trustworthy wall time). |
| `CHAIN=1`     | submits arrays with an LSF `ended()` dependency chain \u2014 only one array runs at a time (no cross-array contention). |

Strictest, cleanest configuration (single task per node, arrays serialized):

```bash
CHAIN=1 EXCLUSIVE=1 \
  STRUCTURES='linCmt ode' NS='40 80 300' SCENARIOS=16 \
  ESTIMATORS=focei FOCEI_OPTS=bobyqa \
  bash script/hpce_scm_estimator/submit_all_arrays.sh 20 1
```

(`MAXPAR=1` \u2192 no intra-array concurrency either.) Use a small dataset count
(20\u201330) \u2014 exclusive nodes are scarce and medians stabilise quickly.

Then extract the honest numbers:

```bash
# LSF whole-tree CPU vs wall (authoritative CPU; parses the .out summaries)
Rscript script/aggregate_lsf_cpu.R --logs logs --out_dir output/scm_timing
# in-record WALL phase split (base/scm/refit) for the scm_sec speedup ratio
Rscript script/aggregate_scm_timing.R --root output/scm_bench --out_dir output/scm_timing
```

### The workers=1 vs workers=3 A/B

The driver takes `--workers` (default 3) and `--rx_threads` (default 1, and
forced to 1 whenever `workers>1` for fork-safety). `submit_one_array.sh` /
`submit_all_arrays.sh` forward these via `WORKERS` / `RX_THREADS`, and route
each arm to its own `OUT_ROOT` + `JOBTAG` so tier-1 caches and LSF logs never
collide. The wrapper `submit_timing_ab.sh` submits both arms (serial, exclusive):

```bash
# focei_bobyqa, scn16, all N & structures, 20 datasets, both arms
bash script/hpce_scm_estimator/submit_timing_ab.sh 20
#   arm w1 -> output/scm_timing_ab/w1  (WORKERS=1, serial SCM)
#   arm w3 -> output/scm_timing_ab/w3  (WORKERS=3, forked SCM)

# per-arm WALL phase split, then compare:
Rscript script/aggregate_scm_timing.R --root output/scm_timing_ab/w1 --out_dir output/scm_timing_ab/agg_w1
Rscript script/aggregate_scm_timing.R --root output/scm_timing_ab/w3 --out_dir output/scm_timing_ab/agg_w3
```

The publishable metric is the **SCM-phase wall-time ratio**, per `(N, structure)`:

```
speedup = wall_scm_med(w1) / wall_scm_med(w3)
```

from the two `scm_timing_by_cell.csv` files. Both arms ran one task per
exclusive node, so the wall times are uncontended and the ratio isolates the
fork benefit from the serial base/refit/covariance phases. `aggregate_lsf_cpu.R`
provides the corroborating LSF CPU cross-check.

## Aggregation (operating characteristics)

Roll the per-dataset `res_ds*.rds` records up into operating characteristics
with `script/aggregate_scm_estimator2.1.R` (sibling of
`aggregate_vae_covsel.R`; the "2.1" suffix marks the schema-2.1 record format,
distinct from the old-schema `aggregate_bench_refit_results.R`). It scans
`output/scm_bench/N*/scn*_*/*_*/res_ds*.rds` (skipping the `.fit.rds`
sidecars), parses the path coordinates (`scn<SS>_<structure>` plus the
`<est>_<opt>` dir) as a backstop for in-RDS keys, **derives** the true
covariate set per scenario from the `PsN_scenarios` indicators (SCM records do
not store `true_set`), recomputes `converged` from numerical evidence, and
writes 11 CSVs + 1 bundled `.rds` to `output/scm_bench_aggregated/`.

### Command line (HPCE)

```bash
Rscript "script/aggregate_scm_estimator2.1.R" --root output --sub scm_bench
# custom output dir:
# Rscript "script/aggregate_scm_estimator2.1.R" --root output --sub scm_bench --out_dir output/scm_bench_aggregated
```

### Interactive R

```r
source("script/aggregate_scm_estimator2.1.R")
res <- aggregate_scm_bench_run(root = "output", sub = "scm_bench")

res$power           # Power / PowerCN / PowerMinSuc per cell
res$diag_rates      # %converged, CN, WALL/CPU timing, hog per cell
subset(res$estim_all, param_class == "covariate_beta")   # clean covariate-beta rel-err
res$covsel_by_cov   # per-covariate detection rate
res$relpower        # per-k fraction recovering >= k true covariates
```

Both approaches write to `output/scm_bench_aggregated/`:

| file | contents |
|------|----------|
| `scm_file_index.csv`      | one row per discovered RDS (+ `has_error`) |
| `scm_diag_long.csv`       | per fit: convergence, CN, WALL/CPU timing, selection tally |
| `scm_rse_long.csv`        | per (fit, parameter): `rel_err` + `param_class` |
| `scm_covsel_long.csv`     | per (fit, var, covar): `in_true`/`in_scm`/`verdict` (TP/FN/FP) |
| `scm_diag_rates.csv`      | per cell: %Converged(Strict), CN, MedObjF, WALL phase timing, hog |
| `scm_estim_all.csv`       | per (cell, parameter): MedRE / MARE / RMRSE (all fits) |
| `scm_estim_success.csv`   | same, strict-converged fits only |
| `scm_estim_cond.csv`      | same, exact-match fits only |
| `scm_power.csv`           | per cell: Power / PowerCN / PowerMinSuc |
| `scm_relpower.csv`        | per (cell, k): fraction recovering >= k true covariates |
| `scm_covsel_by_covar.csv` | per (cell, var, covar): detection rate + TP/FP/FN |

Grouping cell = `(sample_N, scenario, structure, estimator, outer_opt)`. The
completed sweep is `focei_bobyqa` only, so estimator/outer_opt are effectively
constant, but the keys generalise to additional cells.

**Selection scoring** is on `(var, covar)` (not functional shape) for parity
with the VAE pilot and the DGP truth. **Timing note:** `scm_diag_rates` reports
the WALL phase split (`MedWallBase/Scm/Refit/Total_sec`) -- the usable timing.
The `MedCpuTotal_sec` / `MedHogFactor` columns are carried through but
UNDER-report (fork workers not captured by `proc.time()`); use
`script/aggregate_lsf_cpu.R` for the honest CPU / parallelism story. See
*Measuring the parallel benefit correctly* above.

**Reference caveat** (estimation metrics only): `covariate_beta` rel-err is
centring-invariant and clean; `structural_intercept` (TVCL/TVVc) is
reference-dependent -- interpret with care.

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
