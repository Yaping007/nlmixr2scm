# PsN / NONMEM SCM benchmark

A NONMEM/PsN `scm` parallel of the nlmixr2 `focei_bobyqa` SCM benchmark, built to compare **covariate-selection behaviour** (power, powerCN, powerMinsuc, selection pattern, covariate-effect accuracy, runtime) between the two NLME tools on the **same** simulated populations (`Inputdataset/sim_obs_N{N}/`).

------------------------------------------------------------------------

## 1. Covariate parameterisations (PsN vs nlmixr2)

The two tools use different default algebra for covariate effects. Continuous forms coincide on the estimation scale; the **categorical** form differs and is bridged in the parser.

### Continuous covariates

| shape | PsN state | PsN model form | nlmixr2 form | transform |
|---------------|---------------|---------------|---------------|---------------|
| power | `5` | `PARAM = TVP * (cov/ref)**θ` | `TVP * (cov/ref)^θ` | none |
| exp / lin | `4` | `PARAM = TVP * EXP(θ*(cov-ref))` | `TVP * exp(θ*(cov-ref))` | none |

Both tools center continuous covariates on **fixed physiological reference values** — **BW = 70 kg, CrCL = 95 mL/min**. PsN's `scm` centers on the covariate **median** by default (there is *no* `[reference_values]` config keyword — an earlier attempt to add one aborted every run with *"Found invalid section: reference_values"*). Instead the anchors are pinned inside NONMEM via `run.scm`'s **`[code]` section**, which overrides the default exp/power algebra so the fixed reference replaces `median`:

```         
[code]
*:BW-4=PARCOV=EXP(THETA(1)*(COV-70))
*:BW-5=PARCOV=((COV/70)**THETA(1))
*:CRCL-4=PARCOV=EXP(THETA(1)*(COV-95))
*:CRCL-5=PARCOV=((COV/95)**THETA(1))

[lower_bounds]
*:BW-4=-5
... (θ bounds -5, 5 for each)
[upper_bounds]
*:BW-5=5
...
```

θ is **center-invariant** (only the structural intercept TVCL/TVV2 shifts), so pinning the anchor makes the structural typical values directly comparable to the truth (defined at 70/95) with **no post-hoc re-centering in the parser**. BMI (a distractor) keeps PsN's median default. `continuous=1,4,5` in `run.scm` reproduces nlmixr2's `{power, lin}` shape set (`1` = not-included). `REFERENCE_VALUES` in `parse_psn_scm.R` is kept only to *validate* that the winner `.mod` actually carries the expected 70/95 anchor.

### Shape-handling difference (PsN staged vs nlmixr2 simultaneous)

The two tools **enumerate continuous shapes differently**, recorded here as a known, deliberate method difference (not a bug to equalize):

|   | nlmixr2 `runSCM` | PsN `scm` |
|------------------------|------------------------|------------------------|
| shapes offered per step | **both** `lin` **and** `power` per covariate, simultaneously | **one** shape at a time: enter at state 4 (exp/lin); the state-4→5 **power upgrade** is a *separate later step* (a free `dDF=0` OFV move) |
| entry statistic | $\max(\Delta\text{OFV}_\text{lin},\,\Delta\text{OFV}_\text{power})$ vs 1-df threshold | single-shape 1-df LRT, then a 0-df upgrade test |
| fits per step | \~2× continuous candidates | fewer per step, but **more steps** |

**The shapes are mathematically identical**, so cross-tool shape comparison is a 1:1 relabel. nlmixr2 parameters are on the log scale ($\text{CL}=\exp(\theta+\eta)$), so its **linear** term sits inside the exponent and equals PsN's **exponential** (state 4):

$$
\exp\!\big(\theta_{\text{CL}}+\eta+\theta_{\text{cov}}(\text{cov}-\text{ref})\big)
= \text{TVCL}\cdot\exp\!\big(\theta_{\text{cov}}(\text{cov}-\text{ref})\big),
\qquad
\Big(\tfrac{\text{cov}}{\text{ref}}\Big)^{\theta}=\exp\!\big(\theta\log\tfrac{\text{cov}}{\text{ref}}\big).
$$

Hence `nlmixr2 lin ≡ PsN exp (state 4)`, `nlmixr2 power ≡ PsN power (state 5)`, `cat ≡ state 2` — no lossy collapse.

**Impact on operating characteristics** (second-order, partly cancelling; we **measure** rather than remove them):

-   **Power.** nlmixr2's simultaneous `max(lin, power)` entry is slightly more *sensitive* to a strongly-nonlinear (power-truth) covariate; PsN offers exp first and may need the later upgrade → nlmixr2 marginally higher raw power.
-   **False positives.** nlmixr2 runs \~2× continuous tests per step → higher family-wise look count → marginally higher FP; PsN's staged search is more conservative.
-   **Path dependence.** Staged vs simultaneous greedy order can reach different final sets on borderline datasets — genuine method variance, averaged over the 100 datasets per cell.
-   **Model-fit count.** PsN trades more sequential fits for narrower per-step menus; the parser records `scm$n_models_fit` (+ `n_forward_fit`, `n_backward_fit`) so effort is reported **alongside** `cpu_sec`.

**Which is "better"?** No universal winner: nlmixr2's simultaneous test is the cleaner hypothesis test (higher sensitivity); PsN's staged single-shape climb is classical Jonsson–Karlsson SCM (higher specificity, controlled complexity). The benchmark's job is to **quantify** the trade-off — report `power` (any relation), `power_exact` (relation + correct shape), and false-positive rate side by side, stratified by N and effect strength, using both tools **at their defaults**.

### Categorical covariates (SEX, RACE)

Both tools use the **most-frequent level as reference** (verified in source: `R/scm.R:1187` `ref <- names(tbl)[[1L]]  # most frequent level = reference`; and in the PsN `.mod`: `IF(SEX.EQ.0) V2SEX = 1  ; Most common`). Same reference ⇒ same sign, no flip.

The functional forms differ:

$$
\text{PsN (state 2):}\quad \text{PARAM} = \text{TVP}\cdot\big(1 + \theta_{\text{psn}}\,I\big)
\qquad
\text{nlmixr2:}\quad \text{PARAM} = \text{TVP}\cdot\exp\big(\theta_{\text{nlmixr}}\,I\big)
$$

They describe the *same* fold-change for a binary level, so there is an exact algebraic bridge (applied in `parse_psn_scm.R`):

$$
\exp(\theta_{\text{nlmixr}}) = 1 + \theta_{\text{psn}}
\;\Longrightarrow\;
\theta_{\text{nlmixr}} = \log(1+\theta_{\text{psn}}),
\qquad
\text{SE}_{\text{nlmixr}} = \frac{\text{SE}_{\text{psn}}}{1+\theta_{\text{psn}}}
\;\;(\text{delta method})
$$

Example (pilot): PsN `θ_psn=0.450` → `log(1.450)=0.372`, comparable to truth `VcSEX = log(1.5) = 0.405`. We keep `categorical=1,2` (native) and transform θ afterward.

> **Reference caveat:** both tools pick most-common per resample, so they always agree with *each other*; the truth anchor `VcSEX` is defined on `SEX=1`. The parser reads the `; Most common` comment from the `.mod` and, if a resample ever makes `SEX=1` the majority, flips the effect back onto `SEX=1` before comparing.

------------------------------------------------------------------------

## 2. File organisation

Code lives in `script/psn_scm/` (version-controlled, constant across all cells). Per-cell artifacts live under `output/psn_scm/`, split into **heavy/disposable** (`runs/`, git-ignored) and **light/committed** (`records/`, `aggregated/`).

```         
script/psn_scm/
  export_one_dataset.R        # emits one cell's inputs (CLI: --N --scenario --dataset)
  submit_scm.sh               # timed scm(base+search)+refit on one node
  parse_psn_scm.R             # PsN output -> schema-2.1 record (theta, OFV, cond#, omega/sigma)
  parse_full_grid.sh          # batch parse manifest -> records/ (+ optional STAGE_DST res)
  aggregate_psn_scm.R         # PsN unpack overrides -> shared compute_scm_* -> CSVs
  run_pilot_scn16_N300.sh     # pilot driver: export+submit ds 1..5
  parse_pilot_scn16_N300.R    # pilot parser: 5 cells -> records/
  OUTPUT_SCHEMA.md            # PsN output tree + schema-2.1 field mapping
  README.md

output/psn_scm/
  runs/       N{n}/scn{s}/ds{d}/            # HEAVY, git-ignored (HPCE scratch)
    data.csv  base.mod  run.scm  final_refit_cov.txt  submit_scm.sh   # inputs
    refit_tol.txt                           # ODE only: MAXEVAL=0 refit TOL/ATOL
    scm_dir/  refit/                        # PsN/NONMEM churn
    logs/     timing.json  scmlog.txt  scm_console.log
              refit.lst  refit.ext  refit_console.log
  records/    N{n}/scn{s}/ds{d}/            # LIGHT, committed
    psn_scm_record.rds  psn_scm_record.meta.json
  aggregated/                               # cohort roll-ups, committed
    psn_scm_records_N{n}.rds  psn_vs_nlmixr2_scm.csv
```

**Keys** are zero-padded and identical everywhere: `N{n}/scn{s}/ds{d}` (`ds001`), so globbing `records/**/psn_scm_record.rds` and joining against the nlmixr2 side is trivial. The `.gitignore` re-includes `records/` and `aggregated/` while keeping `runs/` out.

`export_one_dataset.R` → `submit_scm.sh` → `parse_psn_scm.R` is a **generate → run → parse** chain; the light `records/` output is the only thing the `compute_scm_*` aggregators (shared with nlmixr2) consume.

### Dataset format (parity with `to_nm_dataset()`)

-   single 100 mg oral dose → depot (`CMT=1`, `EVID=1`, `MDV=1`) at `TIME=0`
-   observations of `cp_obs` in central (`CMT=2`, `EVID=0`); `t=0` obs dropped
-   covariates `BW BMI CRCL SEX RACE` on every row; `SEX`/`RACE` are 0/1
-   N=300 → 1800 rows (300 subj × \[1 dose + 5 obs\])

------------------------------------------------------------------------

## 3. Time tracking

nlmixr2's benchmark reports `t_total = t_base + t_scm + t_refit` (three separately-timed stages summed; `PerformanceEvaluation_scm_bench.R:280-338`). PsN's `scm` fits the base model as its own **node-0** and runs forward+backward selection in one call, so:

$$
t_{\text{total}}^{\text{nlmixr2}} = t_{\text{base}} + t_{\text{scm}} + t_{\text{refit}},
\qquad
t_{\text{total}}^{\text{PsN}} = \underbrace{t_{\text{scm}}}_{\text{base + search}} + t_{\text{refit}}
$$

`submit_scm.sh` measures the whole pipeline as **one block**:

| stage | command | produces |
|------------------------|------------------------|------------------------|
| 1 — scm | `scm -config_file=run.scm` | base fit (node-0) + selection; `scmlog.txt`, winner `.mod/.ext/.lst` |
| 2 — refit (**mandatory**) | seed winner inits → `MAXEVAL=0` + `$COVARIANCE` → `execute` | SE + condition number (`refit.lst/.ext`) |

### Locating the winning model for the covariance refit

PsN 5.5 on this build has a **missing `pharmpy` binary**, so `scm`'s own `final_models/final_backward.mod` (and `update_inits`) are written **corrupt** and are **never used**. The winner is instead reconstructed by **content matching**:

1.  **Derive the final selected set from `scmlog.txt`** (not from `final_models/`). Every accepted step is logged as `... chosen in this forward step: CL-CRCL-4` / `... backward step: CL-SEX-1`. Replaying them — forward adds/updates a `(param,covar)→state`, a **backward state-1 removes it** — yields the post-backward signature, e.g. `CLBW-5 CLCRCL-5 V2BW-5`. **Backward elimination is fully captured here**, purely from the log; if backward dropped or re-stated a relation, the signature reflects it.
2.  **Match that signature to a real candidate `.mod`.** PsN fits every candidate as a genuine NONMEM run and leaves it in the nested `forward_scm_dir*/` / `backward_scm_dir*/` trees. `submit_scm.sh` (`mod_sig()`) and `parse_psn_scm.R` (`find_final_model()`) compute each candidate's `tag-state` signature (power=5, exp=4, cat=2) and pick the `.mod` whose set **exactly equals** the selected set. Because the final model was itself fit during the search, its `.mod/.lst/.ext` exist somewhere in the tree.
3.  **Refit it with `MAXEVAL=0`.** `seed_refit.R` injects the winner `.ext` final estimates back as initials and forces `MAXEVAL=0`, so NONMEM evaluates
    -   computes the covariance **at** the winner. The refit OFV therefore equals the SCM winner OFV **by construction** — a validity guard rejects the refit if `|refit OFV − winner OFV| > 0.01`. Executing the winner `.mod` as-is would instead re-minimise from screening initials (θ ≈ 0.001) and drift to a wrong optimum.

> **No base/wrong-state fallback.** If no exact signature match is found, the refit is **skipped** (`cov_done=FALSE`) rather than refitting the wrong model — a mismatched cond#/SE would be worse than a missing one. The one exception is a **null-selection** cell (nothing selected, e.g. scenario 1): there the **base model** *is* the winner, and `find_base_model()` records its OFV/cond#/TVs.

### Field provenance (which file each record field comes from)

| field | source | why |
|------------------------|------------------------|------------------------|
| selected set, shapes | `scmlog.txt` | the authoritative selection log |
| **θ, OBJ, status** | winner **search** `.lst`/`.ext` | the run that actually *performed* the minimisation search (whether it ended `MINIMIZATION SUCCESSFUL` or `TERMINATED`) |
| **SE, cond#** | **refit** `.lst`/`.ext` | the search had no `$COV`; `MAXEVAL=0` refit adds them at the same point |
| **ω², ω-cov, σ (random effects)** | **refit** `.ext` (else search `.ext`) | `OMEGA(1,1)`→`var_CL`, `OMEGA(2,2)`→`var_Vc`, `OMEGA(2,1)`→`cov_VcCL`; `ResErr = sqrt(SIGMA(1,1))` (NONMEM stores the proportional **variance**; the DGP defines the residual on the **SD** scale) |

Status is read from the **search** `.lst` because the `MAXEVAL=0` refit performs *no* minimisation and prints no status line. A winner that ended `MINIMIZATION TERMINATED` (e.g. rounding errors, ERROR=134) is recorded as `converged=FALSE` — its OFV is still valid and the refit still reproduces it exactly.

-   **Wall** = `date` around the pipeline; **CPU** = sum of per-command `/usr/bin/time -v` (user+sys). Runs on one node with `-no-run_on_lsf` so every NONMEM child is local and its CPU is captured (no LSF-fanout blind spot).
-   `hog = cpu/wall` — same footing as the nlmixr2 `hog_factor`.
-   Output: `logs/timing.json` → `{wall_sec, cpu_sec, hog_factor, scope:"scm+refit"}`.

### Parallelism mapping (Option A: controlled single-node)

| axis | nlmixr2 `runSCM` | PsN `scm` |
|------------------------|------------------------|------------------------|
| candidates within a step | `workers=3` (`future`) | `-threads=3`, `-no-run_on_lsf` |
| within a single fit | `rxThreads=1` | serial (no `-nodes`/`-parafile`) |
| timing | wall + cpu (`proc.time`) | wall + cpu (`/usr/bin/time`, local) |

`linearize=0` forces full re-estimation each SCM step (LRT parity with nlmixr2).

### HPCE module note

`submit_scm.sh` is self-contained: it sources the module init script (needed in non-interactive batch shells) then `module purge && module load PsN/5.5.0-GCCcore-12.3.0` (the node defaults to `GCCcore/14.2.0`, which conflicts with PsN 5.5's `GCCcore/12.3.0`). NONMEM `7.5.1-gfortran-GCCcore-12.3.0` shares the toolchain, so `NM_VERSION=7.5.1-gfortran` is the default.

------------------------------------------------------------------------

## 4. Truth (scenario 16)

All four effects active — the exact set SCM should recover:

| relation   | shape                | truth θ         |
|------------|----------------------|-----------------|
| CL \~ BW   | power                | 0.75            |
| CL \~ CrCL | power                | 0.50            |
| V2 \~ BW   | power                | 1.00            |
| V2 \~ SEX  | categorical `exp(θ)` | log 1.5 = 0.405 |

Distractors that must **not** be selected: BMI on CL/V2, RACE on CL/V2, CrCL on V2, SEX on CL, and wrong-shape BW/CrCL.

------------------------------------------------------------------------

## 5. Running the pilot (scn16, N=300, datasets 1–5, with refit)

From the **repo root** on DaVinci HPCE:

``` bash
# export 5 cells + submit one single-node job each (scm + mandatory cov refit)
bash script/psn_scm/run_pilot_scn16_N300.sh

# watch
bjobs -w

# once all jobs finish, parse the 5 cells into output/psn_scm/records/
Rscript script/psn_scm/parse_pilot_scn16_N300.R
```

Manual single cell (equivalent to what the driver does per dataset):

``` bash
Rscript script/psn_scm/export_one_dataset.R --N 300 --scenario 16 --dataset 1
cd output/psn_scm/runs/N300/scn16/ds001
sed -i 's/\r$//' submit_scm.sh                 # strip CRLF if edited on Windows
bsub -n 4 -W 10000 "bash submit_scm.sh"
# after it finishes, from repo root:
Rscript script/psn_scm/parse_psn_scm.R --cell output/psn_scm/runs/N300/scn16/ds001 \
        --N 300 --scenario 16 --dataset 1
```

`parse_pilot_scn16_N300.R` prints a one-line summary per dataset (selected set, `objf`, `converged`, condition number, wall time) and writes each `records/.../psn_scm_record.{rds,meta.json}`.

------------------------------------------------------------------------

## 6. Full sweep — submit → parse → aggregate → visualize

Full grid = **N∈{40,80,300} × 16 scenarios × ds 1–100 = 4800 cells**, run in an **isolated tree** `output/psn_scm_full0727/{runs,records}` so it never touches the pilot under `output/psn_scm/` (override with `BENCH_ROOT=`). All drivers run from the **repo root** on DaVinci HPCE and resolve `Rscript` via `R_MODULE` (default `R/4.3.1-gomkl-2022a-0.1`; override `RSCRIPT=/path` or `R_MODULE=<module>`).

### Step 1 — Submit (single LSF job array, preferred)

`run_full_grid_array.sh` exports every cell, writes `array_index.tsv` (idx → cell) + `array_task.sh`, then submits **one** array `psnscm_full[1-4800]%THROTTLE`. LSF paces concurrency via the `%` throttle (no sleep loop); one `bkill <arrayid>` cancels the whole sweep.

``` bash
# dry-run: export + build index, submit nothing
DRY_RUN=1 bash script/psn_scm/run_full_grid_array.sh
# smoke one live cell (confirm scm_dir appears + winner carries (BW/70)**THETA)
N_LIST="300" SCEN_LIST="16" DS_MAX=1 bash script/psn_scm/run_full_grid_array.sh
# FULL launch
bash script/psn_scm/run_full_grid_array.sh

# watch
bjobs -A                       # array summary (NJOBS/DONE/RUN/EXIT)
bjobs -a -J 'psnscm_full' | head
```

Defaults: `WALL_MIN=60`, `MEM_MB=2000`, `THROTTLE=400`, `CORES=4`. Each element reads its `$LSB_JOBINDEX` → cell → `bash submit_scm.sh`, and a manifest row `jobid=<arrayid>[<idx>]` is appended so `parse_full_grid.sh` works unchanged. (The older `run_full_grid.sh` submits \~4800 *individual* `bsub` jobs — kept as a fallback; the array is preferred.)

> **Force rerun / clean relaunch.** `scm` resumes an existing `scm_dir`, so wipe the tree first: `rm -rf output/psn_scm_full0727` (deletes the *disposable* runs; `records/` is regenerated). Strip CRLF on any Windows-edited scripts: `sed -i 's/\r$//' script/psn_scm/*.sh script/psn_scm/*.R`.

### Step 2 — Parse (runs → records)

``` bash
# once bjobs is empty, parse every finished cell (has logs/timing.json)
bash script/psn_scm/parse_full_grid.sh
# re-parse existing records too (after a parser change):
FORCE=1 bash script/psn_scm/parse_full_grid.sh
```

`parse_full_grid.sh` walks the manifest and runs `parse_psn_scm.R` on each cell, writing the light `records/N*/scn*/ds*/psn_scm_record.{rds,meta.json}` (skips cells whose record already exists unless `FORCE=1`). Null-selection cells record the base model; failed exact-match cells record a reduced record.

Each record's `rel_err` now also carries the **random-effect** parameters — `var_CL`, `var_Vc`, `cov_VcCL` (from the refit `$OMEGA BLOCK(2)`) and `ResErr` (= `sqrt(SIGMA(1,1))`) — so the estimation-accuracy figure shows the BSV variances/covariance and residual error alongside the structural TVs and covariate betas.

**One-pass staging (skip the separate restage loop).** `parse_full_grid.sh` can also drop each record straight into the aggregator layout, so Step 3's `find | cp` loop is unnecessary:

``` bash
# record-first, then stage each fresh record -> ResforAggregation1/ in one pass
STAGE_DST=output/psn_scm_full0727/ResforAggregation1 \
FORCE=1 bash script/psn_scm/parse_full_grid.sh

# direct-to-res: write res_ds<D>.rds ONLY (no records/ at all) via parser --out_rds
STAGE_ONLY=1 STAGE_DST=output/psn_scm_full0727/ResforAggregation1 \
FORCE=1 bash script/psn_scm/parse_full_grid.sh
```

`STAGE_DST` copies each record to `${STAGE_DST}/N<N>/scn<SS>_<STRUCT>/<EST>_<OPT>/res_ds<D>.rds` (defaults `STAGE_STRUCT=advan4`, `STAGE_EST=nonmem_scm`, `STAGE_OPT=focei`). `STAGE_ONLY=1` skips `records/` entirely by pointing the parser's `--out_rds` at the res path. Either way the pre-existing `records/` and any earlier `ResforAggregation/` are untouched — new output lands in the `…Aggregation1` sibling.

> Staging is a microsecond-scale `cp`; the cost is the 4800 `parse_psn_scm.R` invocations. `STAGE_ONLY` does **not** meaningfully speed the run — for that, parallelise the loop (`xargs -P`) or batch the parse in one R session.

### Step 3 — Aggregate (records → CSVs)

Use the **PsN-specific** aggregator `script/psn_scm/aggregate_psn_scm.R`. It `source()`s the shared `compute_scm_*` math from `script/aggregate_scm_estimator2.1.R` but overrides the unpack layer for the PsN schema (top-level `objf`/`converged`/`cond_num_cor`/`min_success`/`cov_done`, and **no** median re-centering — PsN records are already at 70/95 via the `[code]` section). This keeps the shared aggregator untouched for the nlmixr2 `runSCM` workflow. It also folds `CrCL`≡`CRCL` case-insensitively so CrCL relations score correctly (a case mismatch previously zeroed power for the CrCL-active scenarios 5–8, 13–16).

If you staged in Step 2 (`STAGE_DST=…`), the records are already in the aggregator layout — just point `--sub` at that dir:

``` bash
# aggregate -> 11 CSVs + 1 RDS in output/psn_scm_full0727_aggregated1/
Rscript script/psn_scm/aggregate_psn_scm.R \
  --root output/psn_scm_full0727 --sub ResforAggregation1 \
  --out_dir output/psn_scm_full0727_aggregated1
```

If you parsed **without** staging, restage the records into the layout `output/<sub>/N<N>/scn<SS>_<struct>/<est>_<opt>/res_ds<D>.rds` first (content is already schema-2.1):

``` bash
# restage records/N40/scn01/ds005/psn_scm_record.rds
#   -> ResforAggregation/N40/scn01_advan4/nonmem_scm_focei/res_ds5.rds
SRC=output/psn_scm_full0727/records ; DST=output/psn_scm_full0727/ResforAggregation
find "$SRC" -name psn_scm_record.rds | while read -r f; do
  n=$(echo "$f"   | sed -E 's#.*/(N[0-9]+)/scn[0-9]+/ds[0-9]+/.*#\1#')
  scn=$(echo "$f" | sed -E 's#.*/scn0*([0-9]+)/ds[0-9]+/.*#\1#')
  ds=$(echo "$f"  | sed -E 's#.*/ds0*([0-9]+)/.*#\1#')
  out="$DST/${n}/scn$(printf '%02d' "$scn")_advan4/nonmem_scm_focei"
  mkdir -p "$out" ; cp "$f" "$out/res_ds${ds}.rds"
done
```

Key outputs consumed by the figures:

| CSV | figure |
|------------------------------------|------------------------------------|
| `scm_diag_rates.csv` | convergence status |
| `scm_power.csv`, `scm_relpower.csv` | power curves |
| `scm_covsel_by_covar.csv` | covariate selection heatmap |
| `scm_estim_{all,success,cond}.csv`, `scm_rse_long.csv` | estimation accuracy |
| `scm_diag_long.csv` | runtime benchmark |

### Step 4 — Visualize

The `fig_*` functions in `script/viz/` each take an `agg_dir` (or a `sources` list) plus a `csv_name`. The PsN CSVs are `scm_*`-prefixed and the structure is `advan4`, so pass those explicitly:

``` r
source("script/viz/fig_diag_rates.R");  source("script/viz/fig_power.R")
source("script/viz/fig_covsel_heatmap.R"); source("script/viz/fig_error_metrics.R")
source("script/viz/fig_runtime.R")
agg <- "output/psn_scm_full0727_aggregated1" ; out <- "output/figures/psn_scm"

p_conv <- fig_diag_rates(agg, csv_name = "scm_diag_rates.csv",
                         metric = c("Converged", "CNBelowCutoff", "ConvergedStrict"),
                         structure = "advan4", title_prefix = "PsN SCM",
                         save = TRUE, out_dir = out)
p_pow  <- fig_power(agg, csv_name = "scm_power.csv",
                    metric = c("Power", "PowerCN", "PowerMinSuc"),
                    structure = "advan4", title_prefix = "PsN SCM",
                    save = TRUE, out_dir = out)
p_sel  <- fig_covsel_heatmap_scm(agg_dir = agg, csv_name = "scm_covsel_by_covar.csv",
                    structure = "advan4", estimator = "nonmem_scm",
                    outer_opt = "focei", save = TRUE, out_dir = out)
p_est  <- fig_error_metrics(agg, scenario = 16, structure = "advan4",
                    estimator = "nonmem_scm", outer_opt = "focei",
                    metrics = c("RMRSE", "MARE"),
                    csv_all = "scm_estim_all.csv", csv_cond = "scm_estim_cond.csv",
                    save = TRUE, out_dir = out)
p_rt   <- fig_runtime(sources = list(
  `PsN SCM` = list(agg_dir = agg, csv = "scm_diag_long.csv",
                   estimator = "nonmem_scm", outer_opt = "focei")),
  save = TRUE, out_dir = out)
```

The viz layer was extended for PsN: the `.STRUCT_LAB` maps in `fig_power`/`fig_diag_rates`/`fig_covsel_heatmap` now include `advan4 → "NONMEM (ADVAN4)"`, and `fig_covsel_heatmap_scm` matches `var`/`covar`/`shape` case-insensitively and recognises the `exp` shape (PsN tests both `exp` and `power` per continuous covariate, so an `exp` selection is a legitimate shape-mismatch distractor).

**PsN vs nlmixr2.** Run the nlmixr2 `focei_bobyqa` records through `aggregate_scm_estimator2.1.R`, then pass **both** `agg_dir`s (or both `sources`) to each figure to draw the two tools side-by-side.

> **Centering caveat for accuracy.** The PsN aggregator (`aggregate_psn_scm.R`) deliberately **skips** `.backtransform_intercepts()` — PsN records are already at 70/95 via the `[code]` section, so re-centering would double-correct. The shared `aggregate_scm_estimator2.1.R` still applies the median→70/95 back-transform for the **nlmixr2** side (runSCM centers on the median). Use the PsN aggregator for PsN records and the shared one for nlmixr2.

------------------------------------------------------------------------

## 7. ODE structure — PsN `ADVAN13` vs runSCM `ode`

The default benchmark compares the NONMEM **analytic** `ADVAN4 TRANS4` (`structure="advan4"`) against runSCM's **analytic** `linCmt()` — both closed-form. A second, independent comparison benchmarks the **general-ODE** path on both tools: NONMEM `ADVAN13` + `$DES` (`structure="ode"`) against runSCM `make_true_model(structure="ode")`. NONMEM's ODE integrator is a closed-source commercial solver (LSODA-family), so this measures whether the two solvers reach the **same covariate-selection operating characteristics** on the identical ODE model — not just the same analytic answer.

Everything downstream is structure-agnostic: the aggregator groups by `(sample_N, scenario, structure, estimator, outer_opt)` and discovers all `scn<SS>_<struct>` dirs, so **adding `ode` records never requires re-aggregating `advan4`** — they simply appear as extra `structure="ode"` rows.

### Model encoding (same 2-cmt oral, written as an ODE)

`export_one_dataset.R --structure ode` emits the identical model as `ADVAN13`:

```
$SUBROUTINE ADVAN13 TOL=6 ATOL=8
$MODEL
  COMP=(DEPOT, DEFDOSE)      ; 1
  COMP=(CENTRAL, DEFOBS)     ; 2  (S2=V2)
  COMP=(PERIPH)              ; 3
$PK   ... CL,V2,Q,V3,KA,S2 (same as ADVAN4) ...
$DES
  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1) - (CL/V2)*A(2) - (Q/V2)*A(2) + (Q/V3)*A(3)
  DADT(3) =  (Q/V2)*A(2) - (Q/V3)*A(3)
$ERROR
  IPRED = F
  Y = IPRED*(1+EPS(1))
$ESTIMATION METHOD=1 INTER MAXEVAL=9999 SIGDIG=4 PRINT=5 NOABORT SADDLE_RESET=1
```

The `$PK`, `$THETA/$OMEGA/$SIGMA`, `run.scm` search space, and the `[code]` 70/95 centering are **byte-identical** to the `advan4` branch — only the structural block (analytic ↔ `$DES`) and the `$SUBROUTINE`/`$ESTIMATION` lines differ. `$DES` matches the runSCM ODE RHS term-for-term (`d/dt(depot|central|periph)`), with `S2=V2` so `F = A(2)/S2` equals `central/vc`.

### Tolerance matching (screening vs final refit)

NONMEM `TOL`/`ATOL` are **significant-digit** style: `TOL=n ⇔ rtol=1e-n`, `ATOL=n ⇔ atol=1e-n`. runSCM's `focei/bobyqa` cell solves at `rtol=1e-6, atol=1e-8` (identical for screening and final — `estimator_factory.R` shares the tight tols across tiers). So the ODE defaults are **`TOL=6 ATOL=8`**, matching runSCM exactly.

The screening search and the final MAXEVAL=0 covariance refit can carry **different** tols (mirrors runSCM's *optional* screen-coarsening knob, which defaults off):

| tier | where the tol lives | knob | default |
|---|---|---|---|
| screening (SCM search) | `$SUBROUTINE ADVAN13 TOL=/ATOL=` in `base.mod` | `--screen_tol` / `--screen_atol` | `TOL=6 ATOL=8` |
| final refit (MAXEVAL=0 cov) | `refit_tol.txt` → `seed_refit.R` rewrites the winner `.mod`'s `$SUBROUTINE` | `--refit_tol` / `--refit_atol` | = screen |

`seed_refit.R` rewrites `ATOL=` then `TOL=` (perl lookbehind so the `TOL` inside `ATOL` is not double-matched) only when `--tol/--atol` are supplied; `advan4` (no ODE tol) is a no-op. Defaults keep both tiers identical for strict runSCM parity; tighten the refit with e.g. `--refit_tol 7 --refit_atol 9` if desired.

### Stall / bad-init escape (parity with runSCM `profileInitOnStall`)

The benchmark runs runSCM with **`maxRetries = 0`** (no perturb-init retries) and **`profileInitOnStall = TRUE`** (a stalled forward candidate is reseeded via a 1-D Brent FOCEi profile). The honest NONMEM analogs on the ODE `$ESTIMATION` line:

| NONMEM option | role | runSCM analog |
|---|---|---|
| `NOABORT` | suppress abort on a transient non-computable objective/gradient (e.g. an ODE step that momentarily fails) and keep iterating from the last good point — a passive safety net, **no** init change | (nlmixr2 internal solver recovery) |
| `SADDLE_RESET=1` | on reaching an apparent minimum, **perturb the full θ vector and restart** the quasi-Newton search to slide off a saddle / flat false minimum | `profileInitOnStall` |
| *(omitted)* `RETRIES=n` | multi-start perturb-init retries | `maxRetries` — **kept at 0**, so `RETRIES` is deliberately **not** set (fair comparison) |

> **Caveat — not algorithmically identical.** `SADDLE_RESET` is a *blind full-vector* perturb-and-restart of the multivariate quasi-Newton (BFGS-style) optimizer; it is **not** Brent-based. runSCM's `profileInitOnStall` is a *targeted 1-D Brent* profile of the specific stalled covariate coefficient, followed by a full FOCEi refit. They share the goal (escape a flat/false optimum without a full multi-start budget) but differ in mechanism — a close parity, not an equivalence. `NOABORT` and `SADDLE_RESET=1` are applied to the **ODE branch only**; the `advan4` branch is left exactly as before so its earlier results stay reproducible.

### Running an ODE sweep (targeted, isolated from advan4)

Run ODE into its **own** `BENCH_ROOT` so the `advan4` tree is never touched:

``` bash
# submit the ODE grid (one LSF array), tol-matched TOL=6 ATOL=8 by default
STRUCT=ode BENCH_ROOT=output/psn_scm_ode0729 \
  bash script/psn_scm/run_full_grid_array.sh
# looser screening / tighter refit example:
STRUCT=ode BENCH_ROOT=output/psn_scm_ode0729 \
  STRUCT_TOL_ARGS="--screen_tol 5 --screen_atol 7 --refit_tol 6 --refit_atol 8" \
  bash script/psn_scm/run_full_grid_array.sh

# parse + stage ONLY the ode rows (STRUCT filter), tagged scn<SS>_ode/
STRUCT=ode STAGE_DST=output/psn_scm_ode0729/ResforAggregation FORCE=1 \
  bash script/psn_scm/parse_full_grid.sh
```

The manifest carries a 7th `structure` column; `parse_full_grid.sh` reads it per-row (so `advan4` and `ode` land in `scn<SS>_advan4/` vs `scn<SS>_ode/`), and `STRUCT=ode` skips non-matching rows. Older 6-column manifests default to `advan4`.

``` bash
# aggregate the ode tree (separate out_dir); advan4 untouched
Rscript script/psn_scm/aggregate_psn_scm.R \
  --root output/psn_scm_ode0729 --sub ResforAggregation \
  --out_dir output/psn_scm_ode0729_aggregated
```

### Visualising the ODE comparison

The four `fig_*` `.STRUCT_LAB` maps already carry `ode → "ODE"`. Because **both** tools use `structure="ode"`, the NONMEM-vs-runSCM distinction is the **estimator** facet (`nonmem_scm` vs the runSCM label), not the structure. Pass `structure = "ode"` and the estimator explicitly:

``` r
p_sel <- fig_covsel_heatmap_scm(agg_dir = "output/psn_scm_ode0729_aggregated",
           csv_name = "scm_covsel_by_covar.csv", structure = "ode",
           estimator = "nonmem_scm", outer_opt = "focei",
           save = TRUE, out_dir = "output/figures/psn_scm_ode")
```

For a direct **NONMEM-ODE vs runSCM-ODE** figure, aggregate the runSCM `ode` records with `aggregate_scm_estimator2.1.R` and pass both `agg_dir`s (or both `sources`) — the two share `structure="ode"` and separate on `estimator`.
