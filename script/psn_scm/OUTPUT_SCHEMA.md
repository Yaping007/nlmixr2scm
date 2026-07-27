# PsN SCM output structure & schema 2.1 mapping

Reference: pilot cell `output/psn_scm/N300/scn16/ds001/` (PsN 5.5.0, NONMEM
7.5.1-gfortran, scenario 16, N=300, dataset 1). This documents (a) the PsN
`scm` output tree, and (b) how each schema-2.1 field is parsed from it so PsN
records slot into the same aggregators as the nlmixr2 `focei_bobyqa` bench.

---

## 1. What the run produced

The cell directory after `bash submit_scm.sh`:

```
ds001/
├── base.mod              # the base NONMEM control stream we generated
├── data.csv              # NONMEM dataset (1800 rows)
├── run.scm               # PsN scm config
├── final_refit_cov.txt   # $COVARIANCE line for the final refit (tier-2)
├── scm_console.log       # tee'd PsN stdout — THE STEP TRACE (see below)
├── timing.json           # wall_sec / cpu_sec / hog_factor
├── 273953.out / .err     # LSF job logs
└── scm_dir/              # PsN working directory (-directory=scm_dir)
    ├── backward_scm_dir1/
    │   ├── m1/            # backward candidate models (one relation removed each)
    │   │   ├── CLBW4.{mod,lst,ext,phi}
    │   │   ├── CLCRCL4.{mod,lst,ext,phi}
    │   │   ├── CLSEX1.{mod,lst,ext,phi}   <- the WINNER (SEX-on-CL removed)
    │   │   ├── V2BW4.{mod,lst,ext,phi}
    │   │   └── V2SEX1.{mod,lst,ext,phi}
    │   ├── modelfit_dir1/NM_run*/         # raw NONMEM scratch (FDATA, psn.*)
    │   ├── data.csv
    │   └── modelfit1.log
    ├── PsN_scm_plots.R
    └── .Rprofile
```

> **`-clean=1` caveat.** It pruned `scmlog.txt`, `raw_results.csv`, and the
> *forward* `scm_dir` — only the backward tree survived. For the benchmark set
> **`-clean=0`** so `scm_dir/scmlog.txt` and `scm_dir/raw_results.csv` persist
> (the cleanest parse sources). Until then, the **`scm_console.log`** (captured
> via `tee`) is a complete substitute: it logs every `Adding …` / `Removing …`
> decision, which is exactly the selection path.

### The step trace (from `scm_console.log`)

```
Forward:  Adding CRCL on CL state 4  -> upgraded to state 5 (power)
          Adding SEX  on V2 state 2  (categorical, linear)
          Adding BW   on V2 state 4  -> upgraded to state 5 (power)
          Adding BW   on CL state 4  -> upgraded to state 5 (power)
          Adding SEX  on CL state 2  (categorical)          <- false positive
Backward: Removing SEX on CL state 1 (dropped)
```

**Final selected set** = {CL~BW (power), CL~CrCL (power), V2~BW (power),
V2~SEX (categorical)} — this is an **EXACT MATCH** to the scenario-16 truth.
The forward false-positive (SEX-on-CL) was correctly eliminated in backward.

---

## 2. Key file formats

### `.mod` — covariate parameterization (drives the back-transform)
```
;;; power (continuous), median-centered:
   CLBW  = ((BW/80)**THETA(6))        ; BW median   = 80
   CLCRCL= ((CRCL/98.83)**THETA(7))   ; CrCL median = 98.83
   V2BW  = ((BW/80)**THETA(8))
;;; categorical SEX, state 2 (linear-multiplicative), ref = most common:
   IF(SEX.EQ.0) V2SEX = 1
   IF(SEX.EQ.1) V2SEX = (1 + THETA(9))
```
- **Median centering** confirmed (80 / 98.83), matching nlmixr2 runSCM — NOT
  the 70/95 reference of the refit study.
- Theta index map is read from the `$THETA` comments: `CLBW→6, CLCRCL→7,
  V2BW→8, V2SEX→9` (structural TVCL..TVKA = θ1..θ5).

### `.ext` — parameter estimates (NM7 iteration codes)
Special negative ITERATION rows:
| code | meaning |
|---|---|
| `-1000000000` | **final estimates** (the row to read θ/Ω/Σ + OBJ) |
| `-1000000001` | standard errors — **absent here** (no `$COV` in screening) |
| `-1000000006` | fixed-parameter flags |
| `-1000000007` | parameter lower bounds |

Final row for the winner `CLSEX1.ext`:
`θ6=0.653, θ7=0.544, θ8=0.868, θ9=0.450, OBJ=-6088.60`.

### `.lst` — status + objective
```
0MINIMIZATION SUCCESSFUL
 OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT:   -6088.5965...
```
Condition number appears here **only when `$COV` ran** (final refit), as
`EIGENVALUES` / `COND. NR.`.

### `timing.json`
`wall_sec`, `cpu_sec`, `hog_factor` — measured locally (Option A single node).

---

## 3. Schema 2.1 field mapping (PsN → package_scm_schema21)

| schema-2.1 field | source in PsN output |
|---|---|
| `schema_version` | `"2.1"` (constant) |
| `model_type` | `"ode"`? No — `ADVAN4 TRANS4` is analytic → `"linCmt"`-equivalent; tag `"advan4"` |
| `objf` | final model `.lst` "OBJECTIVE FUNCTION VALUE WITHOUT CONSTANT" |
| `converged` | `.lst` `MINIMIZATION SUCCESSFUL` (final refit); numeric fallback = finite OBJ + finite cov |
| `cond_num_cor` | final-refit `.lst` `COND. NR.` (only after `$COV`); NA for screening-only |
| `scm$selected` | parsed from `scm_console.log` `Adding/Removing` trace (or `scmlog.txt` if `-clean=0`) → tibble(var, covar, shape, theta_name, estimate) |
| `scm$step_hist` | full `Adding/Removing` sequence with step OFV (from `scmlog.txt`/console) |
| `rel_err` | final `.ext` θ6..θ9 vs `true_params_long`; **categorical back-transform** applied (below) |
| `parFixed` | final-refit `.ext` estimates + SE row (`-1000000001`) after `$COV` |
| `cov` | final-refit `.cov` matrix |
| `runtime` | `timing.json$wall_sec` → `total_sec` |
| `cpu` | `timing.json$cpu_sec`, `hog_factor` |
| identity keys | `sample_N=300, scenario_id=16, dataset_id=1, estimator="nonmem_scm", outer_opt="focei", structure="advan4"` |

### Shape mapping (PsN state → nlmixr2 shape → transform)
| PsN state | `.mod` form | nlmixr2 shape | θ transform for `rel_err` |
|---|---|---|---|
| 5 | `(cov/med)**θ` | `power` | none (same scale) |
| 4 | `exp(θ*(cov-med))` | `lin` | none (same scale) |
| 2 (categorical) | `(1 + θ*I)` | `cat` (`exp(θ·I)`) | **`θ_nlmixr = log(1+θ_psn)`**, `SE = SE_psn/(1+θ_psn)` |

Worked example (V2~SEX): PsN `θ9 = 0.450` → `log(1.450) = 0.372`; truth
`log(1.5) = 0.405` → relative error `(0.372-0.405)/0.405 = -8.2%`.

### Selection scoring (feeds compute_scm_power etc.)
- True set (scn16): `CL~BW(power)`, `CL~CrCL(power)`, `V2~BW(power)`,
  `V2~SEX(cat)`.
- Selected set (this ds): identical → `n_true_hit=4, n_false_pos=0,
  exact_match=TRUE` → contributes to Power numerator.

---

## 4. Parser inputs checklist (per cell)

For the parser (`parse_psn_scm.R`, to be written) each cell needs:
1. `scm_dir/scmlog.txt` **or** `scm_console.log` — selection path & step OFV.
2. Winning backward/forward model `.lst` + `.ext` — final θ, OBJ, MIN status.
3. Final-refit `.lst`/`.cov`/`.ext` (after applying `final_refit_cov.txt`) —
   SE + condition number for `converged`/`cond_num_cor`/`parFixed`.
4. `timing.json` — wall/cpu/hog.

**Action item:** rerun the benchmark with `-clean=0` (edit `submit_scm.sh`) so
`scmlog.txt`, `raw_results.csv`, and the forward tree survive for robust parsing.
