# PsN / NONMEM SCM benchmark

A NONMEM/PsN `scm` parallel of the nlmixr2 `focei_bobyqa` SCM benchmark, built to
compare **covariate-selection behaviour** (power, powerCN, powerMinsuc, selection
pattern, covariate-effect accuracy, runtime) between the two NLME tools on the
**same** simulated populations (`Inputdataset/sim_obs_N{N}/`).

---

## 1. Covariate parameterisations (PsN vs nlmixr2)

The two tools use different default algebra for covariate effects. Continuous
forms coincide on the estimation scale; the **categorical** form differs and is
bridged in the parser.

### Continuous covariates

| shape | PsN state | PsN model form | nlmixr2 form | transform |
|---|---|---|---|---|
| power | `5` | `PARAM = TVP * (cov/ref)**θ` | `TVP * (cov/ref)^θ` | none |
| exp / lin | `4` | `PARAM = TVP * EXP(θ*(cov-ref))` | `TVP * exp(θ*(cov-ref))` | none |

Both tools center continuous covariates on **fixed physiological reference
values** — **BW = 70 kg, CrCL = 95 mL/min** — set in `run.scm`'s
`[reference_values]` block (BMI stays at its default median). The same anchor is
used on the nlmixr2 side, so the continuous θ's are directly comparable —
**no back-transform needed**. `continuous=1,4,5` in `run.scm` reproduces
nlmixr2's `{power, lin}` shape set (`1` = not-included).

### Shape-handling difference (PsN staged vs nlmixr2 simultaneous)

The two tools **enumerate continuous shapes differently**, recorded here as a
known, deliberate method difference (not a bug to equalize):

| | nlmixr2 `runSCM` | PsN `scm` |
|---|---|---|
| shapes offered per step | **both** `lin` **and** `power` per covariate, simultaneously | **one** shape at a time: enter at state 4 (exp/lin); the state-4→5 **power upgrade** is a *separate later step* (a free `dDF=0` OFV move) |
| entry statistic | $\max(\Delta\text{OFV}_\text{lin},\,\Delta\text{OFV}_\text{power})$ vs 1-df threshold | single-shape 1-df LRT, then a 0-df upgrade test |
| fits per step | ~2× continuous candidates | fewer per step, but **more steps** |

**The shapes are mathematically identical**, so cross-tool shape comparison is a
1:1 relabel. nlmixr2 parameters are on the log scale ($\text{CL}=\exp(\theta+\eta)$),
so its **linear** term sits inside the exponent and equals PsN's **exponential**
(state 4):

$$
\exp\!\big(\theta_{\text{CL}}+\eta+\theta_{\text{cov}}(\text{cov}-\text{ref})\big)
= \text{TVCL}\cdot\exp\!\big(\theta_{\text{cov}}(\text{cov}-\text{ref})\big),
\qquad
\Big(\tfrac{\text{cov}}{\text{ref}}\Big)^{\theta}=\exp\!\big(\theta\log\tfrac{\text{cov}}{\text{ref}}\big).
$$

Hence `nlmixr2 lin ≡ PsN exp (state 4)`, `nlmixr2 power ≡ PsN power (state 5)`,
`cat ≡ state 2` — no lossy collapse.

**Impact on operating characteristics** (second-order, partly cancelling; we
**measure** rather than remove them):

- **Power.** nlmixr2's simultaneous `max(lin, power)` entry is slightly more
  *sensitive* to a strongly-nonlinear (power-truth) covariate; PsN offers exp
  first and may need the later upgrade → nlmixr2 marginally higher raw power.
- **False positives.** nlmixr2 runs ~2× continuous tests per step → higher
  family-wise look count → marginally higher FP; PsN's staged search is more
  conservative.
- **Path dependence.** Staged vs simultaneous greedy order can reach different
  final sets on borderline datasets — genuine method variance, averaged over the
  100 datasets per cell.
- **Model-fit count.** PsN trades more sequential fits for narrower per-step
  menus; the parser records `scm$n_models_fit` (+ `n_forward_fit`,
  `n_backward_fit`) so effort is reported **alongside** `cpu_sec`.

**Which is "better"?** No universal winner: nlmixr2's simultaneous test is the
cleaner hypothesis test (higher sensitivity); PsN's staged single-shape climb is
classical Jonsson–Karlsson SCM (higher specificity, controlled complexity). The
benchmark's job is to **quantify** the trade-off — report `power` (any relation),
`power_exact` (relation + correct shape), and false-positive rate side by side,
stratified by N and effect strength, using both tools **at their defaults**.

### Categorical covariates (SEX, RACE)

Both tools use the **most-frequent level as reference** (verified in source:
`R/scm.R:1187` `ref <- names(tbl)[[1L]]  # most frequent level = reference`; and
in the PsN `.mod`: `IF(SEX.EQ.0) V2SEX = 1  ; Most common`). Same reference ⇒
same sign, no flip.

The functional forms differ:

$$
\text{PsN (state 2):}\quad \text{PARAM} = \text{TVP}\cdot\big(1 + \theta_{\text{psn}}\,I\big)
\qquad
\text{nlmixr2:}\quad \text{PARAM} = \text{TVP}\cdot\exp\big(\theta_{\text{nlmixr}}\,I\big)
$$

They describe the *same* fold-change for a binary level, so there is an exact
algebraic bridge (applied in `parse_psn_scm.R`):

$$
\exp(\theta_{\text{nlmixr}}) = 1 + \theta_{\text{psn}}
\;\Longrightarrow\;
\theta_{\text{nlmixr}} = \log(1+\theta_{\text{psn}}),
\qquad
\text{SE}_{\text{nlmixr}} = \frac{\text{SE}_{\text{psn}}}{1+\theta_{\text{psn}}}
\;\;(\text{delta method})
$$

Example (pilot): PsN `θ_psn=0.450` → `log(1.450)=0.372`, comparable to truth
`VcSEX = log(1.5) = 0.405`. We keep `categorical=1,2` (native) and transform θ
afterward.

> **Reference caveat:** both tools pick most-common per resample, so they always
> agree with *each other*; the truth anchor `VcSEX` is defined on `SEX=1`. The
> parser reads the `; Most common` comment from the `.mod` and, if a resample
> ever makes `SEX=1` the majority, flips the effect back onto `SEX=1` before
> comparing.

---

## 2. File organisation

Code lives in `script/psn_scm/` (version-controlled, constant across all cells).
Per-cell artifacts live under `output/psn_scm/`, split into **heavy/disposable**
(`runs/`, git-ignored) and **light/committed** (`records/`, `aggregated/`).

```
script/psn_scm/
  export_one_dataset.R        # emits one cell's inputs (CLI: --N --scenario --dataset)
  submit_scm.sh               # timed scm(base+search)+refit on one node
  parse_psn_scm.R             # PsN output -> schema-2.1 record
  run_pilot_scn16_N300.sh     # pilot driver: export+submit ds 1..5
  parse_pilot_scn16_N300.R    # pilot parser: 5 cells -> records/
  OUTPUT_SCHEMA.md            # PsN output tree + schema-2.1 field mapping
  README.md

output/psn_scm/
  runs/       N{n}/scn{s}/ds{d}/            # HEAVY, git-ignored (HPCE scratch)
    data.csv  base.mod  run.scm  final_refit_cov.txt  submit_scm.sh   # inputs
    scm_dir/  refit/                        # PsN/NONMEM churn
    logs/     timing.json  scmlog.txt  scm_console.log
              refit.lst  refit.ext  refit_console.log
  records/    N{n}/scn{s}/ds{d}/            # LIGHT, committed
    psn_scm_record.rds  psn_scm_record.meta.json
  aggregated/                               # cohort roll-ups, committed
    psn_scm_records_N{n}.rds  psn_vs_nlmixr2_scm.csv
```

**Keys** are zero-padded and identical everywhere: `N{n}/scn{s}/ds{d}`
(`ds001`), so globbing `records/**/psn_scm_record.rds` and joining against the
nlmixr2 side is trivial. The `.gitignore` re-includes `records/` and
`aggregated/` while keeping `runs/` out.

`export_one_dataset.R` → `submit_scm.sh` → `parse_psn_scm.R` is a
**generate → run → parse** chain; the light `records/` output is the only thing
the `compute_scm_*` aggregators (shared with nlmixr2) consume.

### Dataset format (parity with `to_nm_dataset()`)
- single 100 mg oral dose → depot (`CMT=1`, `EVID=1`, `MDV=1`) at `TIME=0`
- observations of `cp_obs` in central (`CMT=2`, `EVID=0`); `t=0` obs dropped
- covariates `BW BMI CRCL SEX RACE` on every row; `SEX`/`RACE` are 0/1
- N=300 → 1800 rows (300 subj × [1 dose + 5 obs])

---

## 3. Time tracking

nlmixr2's benchmark reports `t_total = t_base + t_scm + t_refit` (three
separately-timed stages summed; `PerformanceEvaluation_scm_bench.R:280-338`).
PsN's `scm` fits the base model as its own **node-0** and runs forward+backward
selection in one call, so:

$$
t_{\text{total}}^{\text{nlmixr2}} = t_{\text{base}} + t_{\text{scm}} + t_{\text{refit}},
\qquad
t_{\text{total}}^{\text{PsN}} = \underbrace{t_{\text{scm}}}_{\text{base + search}} + t_{\text{refit}}
$$

`submit_scm.sh` measures the whole pipeline as **one block**:

| stage | command | produces |
|---|---|---|
| 1 — scm | `scm -config_file=run.scm` | base fit (node-0) + selection; `scmlog.txt`, winner `.mod/.ext/.lst` |
| 2 — refit (**mandatory**) | `update_inits` winner → append `$COVARIANCE` → `execute` | SE + condition number (`refit.lst/.ext`) |

- **Wall** = `date` around the pipeline; **CPU** = sum of per-command
  `/usr/bin/time -v` (user+sys). Runs on one node with `-no-run_on_lsf` so every
  NONMEM child is local and its CPU is captured (no LSF-fanout blind spot).
- `hog = cpu/wall` — same footing as the nlmixr2 `hog_factor`.
- Output: `logs/timing.json` → `{wall_sec, cpu_sec, hog_factor, scope:"scm+refit"}`.

### Parallelism mapping (Option A: controlled single-node)
| axis | nlmixr2 `runSCM` | PsN `scm` |
|---|---|---|
| candidates within a step | `workers=3` (`future`) | `-threads=3`, `-no-run_on_lsf` |
| within a single fit | `rxThreads=1` | serial (no `-nodes`/`-parafile`) |
| timing | wall + cpu (`proc.time`) | wall + cpu (`/usr/bin/time`, local) |

`linearize=0` forces full re-estimation each SCM step (LRT parity with nlmixr2).

### HPCE module note
`submit_scm.sh` is self-contained: it sources the module init script (needed in
non-interactive batch shells) then `module purge && module load
PsN/5.5.0-GCCcore-12.3.0` (the node defaults to `GCCcore/14.2.0`, which
conflicts with PsN 5.5's `GCCcore/12.3.0`). NONMEM `7.5.1-gfortran-GCCcore-12.3.0`
shares the toolchain, so `NM_VERSION=7.5.1-gfortran` is the default.

---

## 4. Truth (scenario 16)

All four effects active — the exact set SCM should recover:

| relation | shape | truth θ |
|---|---|---|
| CL ~ BW | power | 0.75 |
| CL ~ CrCL | power | 0.50 |
| V2 ~ BW | power | 1.00 |
| V2 ~ SEX | categorical `exp(θ)` | log 1.5 = 0.405 |

Distractors that must **not** be selected: BMI on CL/V2, RACE on CL/V2, CrCL on
V2, SEX on CL, and wrong-shape BW/CrCL.

---

## 5. Running the pilot (scn16, N=300, datasets 1–5, with refit)

From the **repo root** on DaVinci HPCE:

```bash
# export 5 cells + submit one single-node job each (scm + mandatory cov refit)
bash script/psn_scm/run_pilot_scn16_N300.sh

# watch
bjobs -w

# once all jobs finish, parse the 5 cells into output/psn_scm/records/
Rscript script/psn_scm/parse_pilot_scn16_N300.R
```

Manual single cell (equivalent to what the driver does per dataset):

```bash
Rscript script/psn_scm/export_one_dataset.R --N 300 --scenario 16 --dataset 1
cd output/psn_scm/runs/N300/scn16/ds001
sed -i 's/\r$//' submit_scm.sh                 # strip CRLF if edited on Windows
bsub -n 4 -W 10000 "bash submit_scm.sh"
# after it finishes, from repo root:
Rscript script/psn_scm/parse_psn_scm.R --cell output/psn_scm/runs/N300/scn16/ds001 \
        --N 300 --scenario 16 --dataset 1
```

`parse_pilot_scn16_N300.R` prints a one-line summary per dataset (selected set,
`objf`, `converged`, condition number, wall time) and writes each
`records/.../psn_scm_record.{rds,meta.json}`.

---

## Next steps (after the pilot)

1. Feed `records/**` through the existing `compute_scm_*` aggregators to build
   `aggregated/psn_vs_nlmixr2_scm.csv`.
2. **Full grid** (N∈{40,80,300} × 16 scenarios × ds 1–100 = 4800 cells).
   Submit and parse with the grid drivers (run from the repo root on HPCE):

   ```bash
   # submit every cell as its own single-node job (throttled, manifest-tracked)
   bash script/psn_scm/run_full_grid.sh
   # smoke a slice first / dry-run:
   DRY_RUN=1 N_LIST="300" SCEN_LIST="16" DS_MAX=5 bash script/psn_scm/run_full_grid.sh

   # watch
   bjobs -o 'stat' | sort | uniq -c

   # once bjobs is empty, parse the whole grid into records/
   bash script/psn_scm/parse_full_grid.sh
   ```

   `run_full_grid.sh` writes `output/psn_scm_full0727/manifest.csv` (one row per
   cell: keys, jobid, path) for tracking and locating failures;
   `parse_full_grid.sh` walks it and parses only finished cells
   (`logs/timing.json` present). The full launch uses an **isolated tree**
   `output/psn_scm_full0727/{runs,records}` so it never overwrites the working
   pilot under `output/psn_scm/` (override with `BENCH_ROOT=`).
   Both resolve `Rscript` via `R_MODULE` (default `R/4.3.1-gomkl-2022a-0.1`);
   override with `RSCRIPT=/path` or `R_MODULE=<module>`. Throttle concurrent
   jobs with `THROTTLE` (default 400).
