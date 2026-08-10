# `script/viz/` — Result Visualization Guide

Reusable plotting layer for the schema-2.1 aggregated simulation outputs.
Each figure is a self-contained function that reads an aggregated CSV and
returns a `ggplot` object (and optionally writes PNG + PDF).

Design goals:

- **Source-agnostic.** Figures read the aggregated CSVs (not raw fits), keyed
  on `sample_N`, `scenario`, `structure`. Future non-VAE estimator
  aggregations (schema 2.1) drop in with an added `estimator` / `outer_opt`
  facet — no figure rewrites.
- **Knob-driven.** Every figure exposes filter arguments (`sample_N`,
  `structure`, metric selectors) that default to "show everything," so you can
  eyeball all combinations first and then focus for publication.
- **Return-first.** During the testing phase functions return the `ggplot`
  object; pass `save = TRUE` to also write files.

Shared conventions (currently defined locally in each script, promotable to a
`viz_common.R` later):

| Element | Value |
|---|---|
| Theme | `theme_scm()` = `theme_minimal(12)`, minor grid off, legend on top, bold strip/title |
| Sample-size palette | `40` grey `#7F7F7F`, `80` orange `#E8820C`, `300` blue `#1F77B4` |
| Structure labels | `linCmt` → "linCmt (analytic)", `ode` → "ODE" |
| Output folder | `output/figures/vae_covsel/` |
| Save format | PNG (`dpi = 150`) + PDF, `width = 9`, `height = 7` |

---

## `fig_power.R` — covariate-selection power curves

Power across the 16 scenarios, faceted by model parameterization
(`structure`: linCmt / ode), coloured by sample size
(`sample_N`: 40 / 80 / 300).

**Data source:** `output/vae_covsel_aggregated/vae_power.csv`
Columns used: `sample_N, scenario, structure, Power, PowerCN, PowerMinSuc`
(the `Power*` columns are **fractions 0–1**; the figure scales to %).

**Power definitions:**

| metric | meaning | denominator |
|---|---|---|
| `Power` | exact-match rate | all `N` datasets |
| `PowerCN` | exact-match among well-conditioned fits | fits with CN < cutoff |
| `PowerMinSuc` | exact-match among converged fits | fits that minimised successfully |

### Signature

```r
fig_power(agg_dir   = "output/vae_covsel_aggregated",
          metric    = c("Power", "PowerCN", "PowerMinSuc"),
          sample_N  = NULL,     # NULL = all (40, 80, 300); or subset e.g. c(80, 300)
          structure = NULL,     # NULL = all (linCmt, ode);  or subset e.g. "linCmt"
          save      = FALSE,    # TRUE also writes PNG + PDF
          out_dir   = "output/figures/vae_covsel")


source("script/viz/fig_power.R")

fig_power(
  agg_dir  = "output/vae_covsel_aggregated0722",   # note: no leading space
  metric   = c("Power", "PowerCN", "PowerMinSuc"),
  sample_N = NULL,     # all N
  structure = NULL,    # both structures
  save     = TRUE      # writes PNG + PDF
)

```

### Knobs

- **`metric`** — accepts **one or several** definitions.
  - One metric → colour = `sample_N`, single line per N.
  - Several metrics → colour = `sample_N`, **linetype = metric**. When the CN /
    convergence gates are non-binding the linetypes overlap exactly (a visual
    proof the three definitions are equal).
- **`sample_N`** — `NULL` shows all; pass e.g. `c(80, 300)` to focus.
- **`structure`** — `NULL` shows both; pass `"linCmt"` or `"ode"` to focus.
- **`save`** — returns the `ggplot` only by default; `TRUE` also writes
  `fig_power_<metric>.{png,pdf}` (or `fig_power_multi.*` in overlay mode).

### Usage

```r
source("script/viz/fig_power.R")

fig_power()                                          # all N, both structures, Power
fig_power(metric = "PowerCN")                        # CN-gated
fig_power(metric = "PowerMinSuc")                    # convergence-gated
fig_power(metric = c("Power","PowerCN","PowerMinSuc")) # overlay by linetype
fig_power(sample_N = c(80, 300))                     # focus larger N
fig_power(structure = "linCmt")                     # focus one parameterization
fig_power(save = TRUE)                               # write PNG + PDF

# save all three single-metric figures at once
for (m in c("Power","PowerCN","PowerMinSuc")) fig_power(metric = m, save = TRUE)
```

### Reading notes (current VAE run)

- N = 300 stays at/above the 80% threshold across nearly all scenarios; N = 80
  straddles it; N = 40 degrades sharply as scenario complexity rises.
- `linCmt` and `ode` curves are nearly identical.
- `Power`, `PowerCN`, `PowerMinSuc` **coincide exactly** — in this dataset
  `n_ok_cn = n_ok_min = N = 250` for every cell, i.e. the CN and convergence
  gates removed nothing. For publication, lead with plain `Power` and footnote
  that the gates were non-binding (the linetype-overlay figure is good
  supporting evidence).

---

## `fig_covsel_heatmap.R` — where FP / FN errors concentrate

Per-effect error-rate heatmap answering *"under each scenario, which covariate
effects tend to become false positive (FP) or false negative (FN)?"* Designed
around the **BW↔BMI collinearity story**: correlated covariates steal each
other's signal, so a true `~BW` effect leaks into `~BMI` (and vice versa).

**Data source:** `output/vae_covsel_aggregated/vae_covsel_by_covar.csv`.
Columns used: `sample_N, scenario, structure, var, covar, shape, is_true,
n_datasets, n_TP, n_FN, n_FP`.

### How FP and FN rates are calculated

The rate in each tile is built up over three stages — per fit → per cell →
per tile.

**Stage 1 — verdict per fit (`.unpack_covsel` in `aggregate_vae_covsel.R`).**
For a single fitted dataset, VAE's *selected* covariate terms are matched
against the scenario's *true* set on `(var, covar)` **and shape** (a true
`power` recovered as `lin` does **not** match). Every effect gets one verdict:

| verdict | meaning | condition |
|---|---|---|
| **TP** | true term recovered on the right shape | in true set **and** selected, shapes match |
| **FN** | true term missed (or recovered on wrong shape) | in true set, **not** matched by a selected term |
| **FP** | spurious term (wrong pair *or* wrong shape) | selected, **not** matching any true term |

This is written one row per `(fit, var, covar, shape)` to
`vae_covsel_long.csv`.

**Stage 2 — counts per cell (`compute_vae_covsel_by_covar`).** Rows are grouped
by `(sample_N, scenario, structure, var, covar, shape)` — a *cell* — and the
verdicts tallied into `n_TP`, `n_FN`, `n_FP`, with `n_datasets` = number of
fits contributing and `is_true = any(in_true)` flagging whether that effect is
truly active in the scenario. Because grouping includes `shape`, a true
`~BW.power` effect and its spurious `~BW.lin` shape-flip are **separate rows**:
the true row accrues FN when the shape is wrong, the distractor row accrues FP.

**Stage 3 — rate per tile (`fig_covsel_heatmap.R`).** Each effect is in exactly
one regime, set by `is_true`, and the tile shows that regime's rate:

| regime | when | rate shown |
|---|---|---|
| **FN** (false negative) | `is_true = TRUE` — a real effect | `n_FN / n_datasets` |
| **FP** (false positive) | `is_true = FALSE` — a null/distractor effect | `n_FP / n_datasets` |

The signed value `err_signed` = `-FN%` (blue) for true effects, `+FP%` (red) for
distractors, so one diverging colour scale carries both: **blue = missed truth,
red = spurious selection, white = no error.** A true effect recovered on the
wrong shape shows as a **blue FN tile stacked on a red FP tile** (the wrong-shape
row), so shape flips read vertically.

> **Denominator (`n_datasets`).** This is the cell's **fitted-dataset count**,
> the common denominator for both regimes, so `n_FN / n_datasets` is a per-effect
> **false-negative rate** (1 − sensitivity) and `n_FP / n_datasets` a per-effect
> **false-positive rate**. The heatmap additionally runs `.fix_covsel_denom()`,
> which broadcasts the authoritative per-scenario total (the max `n_datasets`
> over the cell's *true-effect* rows) to every row. This guards against any
> upstream CSV that stored `n_datasets` as the per-pair *appearance* count —
> where a distractor selected `k` times would get `n_datasets = k` and a
> spurious ≈100 % FP rate — and falls back to the cell-wide total for the null
> scenario (which has no true-effect row).

**Layout:**

- **x** = scenario (1–16); **y** = effect label `PARAM~COVAR`.
- y-order = CL block then VC block; **within each block `~BW` and `~BMI` are
  adjacent**, so a true `~BW` (bold outline = FN) sits directly above its
  `~BMI` thief (plain tile = FP) — leakage reads vertically.
- **fill** = signed error rate (blue = FN on true effects, red = FP on null
  effects, white = no error).
- **bold black outline** = TRUE effect (its shading is the FN rate); plain tiles
  are null effects (shading = FP rate).
- **facet** = by default `sample_N` (rows) × `structure` (columns) so `linCmt`
  and `ode` sit side by side; pass a single `structure=` to get the one-structure
  `sample_N`-facet layout instead.

### Signature

```r
fig_covsel_heatmap(agg_dir   = "output/vae_covsel_aggregated",
                   sample_N  = NULL,       # NULL = all (40, 80, 300); or subset
                   structure = "both",     # "both" = linCmt + ode side by side (default) | "linCmt" | "ode"
                   labels    = TRUE,       # print error % inside each tile
                   min_fp    = 0,          # blank FP cells below this rate (de-clutter)
                   save      = FALSE,      # TRUE also writes PNG + PDF
                   out_dir   = "output/figures/vae_covsel")

# ---- Covariate-selection error pattern (FP / FN heatmap) -------------------
source("script/viz/fig_covsel_heatmap.R")
fig_covsel_heatmap(
  agg_dir   = "output/vae_covsel_aggregated0722",
  sample_N  = NULL,        # all N
  structure = "linCmt",    # one structure per figure; use "ode" to swap
  save      = TRUE
)
fig_covsel_heatmap(agg_dir = "output/vae_covsel_aggregated0722",
                   structure = "ode", save = TRUE)
```

### Knobs

- **`sample_N`** — `NULL` shows all three N; pass e.g. `300` for the large cohort.
- **`structure`** — `"both"` (default) draws `linCmt` + `ode` as side-by-side
  facet columns; pass `"linCmt"` or `"ode"` for a single-structure figure.
- **`labels`** — `TRUE` prints the rounded error % in each tile; `FALSE` for a
  cleaner fill-only grid.
- **`min_fp`** — FP cells below this rate are greyed out to hide near-zero noise
  (e.g. `0.05`).
- **`save`** — writes `fig_covsel_heatmap_<Ntag>_<structure>.{png,pdf}`.

Tune the effect row order via the module constant `.COVAR_ORD` (keeps BW/BMI
adjacent by default).

### Usage

```r
source("script/viz/fig_covsel_heatmap.R")

fig_covsel_heatmap()                                 # linCmt, all N
fig_covsel_heatmap(structure = "ode")                # swap structure
fig_covsel_heatmap(sample_N = 300)                   # large-cohort only
fig_covsel_heatmap(labels = FALSE)                   # fill-only grid
fig_covsel_heatmap(min_fp = 0.05)                    # hide trivial FP noise
fig_covsel_heatmap(save = TRUE)                      # write PNG + PDF
```

### Reading notes (current VAE run)

- FP rates are modest (median ≈ 3 %, max ≈ 31 %), consistent with the ~90 %
  power seen in `fig_power.R` — an early version that divided by `n_datasets`
  wrongly showed ~100 % FP everywhere.
- The headline pattern is the **`~BW` / `~BMI` pair**: scenarios where `~BW` is
  the true effect show elevated `~BMI` FP directly below (and some `~BW` FN),
  i.e. the two continuous covariates trade the signal.

---

## `fig_error_metrics.R` — estimation accuracy, conditioned on correct selection

Per-parameter **RMRSE** / **MARE** dot-plot answering *"how accurately is each
population parameter and covariate effect estimated, and does restricting to
datasets where the covariate model was selected **correctly** change that
accuracy?"* Currently focused on **scenario 16** (all four covariate effects
active).

**Data source:** `output/vae_covsel_aggregated0722/vae_estim_all.csv`
(unconditioned) and `vae_estim_cond.csv` (conditioned on `exact_match == TRUE`,
i.e. the selected covariate set equals truth — no FP, no FN). Both share the
schema-2.1 columns: `sample_N, scenario, structure, parameter, param_class,
true_value, n_used, MedRE_pct, MeanRE_pct, MARE_pct, RMRSE_pct, P90AbsRE_pct`.

**Metrics** (both %, lower = better), from `rel_err = (θ̂ − θ_true)/θ_true`:

| metric | formula | character |
|---|---|---|
| `RMRSE` | `100·√mean(rel_err²)` | outlier-sensitive accuracy |
| `MARE`  | `100·median\|rel_err\|` | robust central accuracy (the repo's "MASE") |

> There is **no metric literally named "MASE"** in the codebase; `MARE_pct` is
> the robust median-relative-error metric and is what "MASE" maps to.

**Parameter groups (`param_class`):**

- `covariate_beta` — `CLBW, CLcrCL, VcBW, VcSEX` (headline; centering-invariant).
- `structural_intercept` — `TVCL, TVVc` (**now reference-aligned**: the
  aggregators back-transform each estimate onto the DGP's fixed `BW = 70` /
  `CrCL = 95` anchor before computing error, so the relative error is genuine
  bias rather than a centering artefact — see *Reference alignment* below).
- `structural_other` — `TVQ, TVVp, var_CL, var_Vc, cov_VcCL, ResErr` (clean
  baseline; fixed `TVKA` is dropped by default since its error is always 0).

### Reference alignment (structural intercepts)

The data-generating model defines `TVCL` / `TVVc` at **fixed** references
`BW = 70`, `CrCL = 95`. The estimators instead centre each *selected* continuous
covariate on a *sample* statistic, so a raw `exp(lTVCL)` is the typical value at
that sample centre — a reference artefact, not bias. Both aggregators now
back-transform the intercept onto the 70 / 95 anchor using the fit's **own**
reference-invariant power betas:

$$\text{TVCL}_{70/95} = \widehat{\text{TVCL}}\cdot\left(\tfrac{70}{c_{BW}}\right)^{\beta_{CLBW}}\!\!\cdot\left(\tfrac{95}{c_{CrCL}}\right)^{\beta_{CLcrCL}},\qquad \text{TVVc}_{70} = \widehat{\text{TVVc}}\cdot\left(\tfrac{70}{c_{BW}}\right)^{\beta_{VcBW}}$$

The centre $c$ differs by pipeline — this is the **only** difference between the
two implementations:

| aggregator | script | centre statistic $c$ | matches engine |
|---|---|---|---|
| VAE covariate selection | `script/aggregate_vae_covsel.R` | `mean` | VAE centres on the sample **mean** |
| runSCM estimator bench   | `script/aggregate_scm_estimator2.1.R` | `median` | `R/scm.R` centres on the subject-level **median** |

Both are wired as the first line of `.unpack_rse()`, so *every* per-parameter
metric (`RMRSE`, `MARE`, `MedRE`, …) is computed on the aligned estimate. The
refit-**true** model is already anchored on 70 / 95 and is skipped. Empirically
verified on scn16 / N = 300: the VAE mean back-transform recovers
`TVCL = 0.594` vs truth `0.6` (median would give `0.628`; uncorrected `0.698`).

```r
# Re-aggregate so the CSVs the figures read carry the aligned intercepts
Rscript script/aggregate_vae_covsel.R \
  --root output --sub vae_covsel_full0722 \
  --out_dir output/vae_covsel_aggregated0722

Rscript script/aggregate_scm_estimator2.1.R \
  --root output --sub scm_bench \
  --out_dir output/scm_bench_aggregated
```

**Layout:**

- **y** = parameter, grouped into `param_class` blocks down the rows.
- **x** = metric value (%); **colour** = conditioning (grey = *Unconditioned*,
  blue = *True selection*), dodged as a lollipop (linerange 0→value + point).
- **facet** = `param_class` (rows, free/space y) × `metric` (cols). Extra
  `sample_N` / `structure` facet columns are added automatically when more than
  one is present (the 0722 run is scn16 / N=300 / ODE only, so those collapse).

### Signature

```r
fig_error_metrics(agg_dir    = "output/vae_covsel_aggregated0722",
                  scenario   = 16,
                  sample_N   = NULL,                # NULL = all present
                  structure  = NULL,                # NULL = all present
                  metrics    = c("RMRSE", "MARE"),
                  drop_fixed = TRUE,                # drop TVKA (fixed → 0)
                  labels     = FALSE,               # print value at each point
                  save       = FALSE,               # TRUE also writes PNG + PDF
                  out_dir    = "output/figures/vae_covsel")
```

### Usage

```r
source("script/viz/fig_error_metrics.R")

fig_error_metrics()                         # scn16, both metrics
fig_error_metrics(metrics = "MARE")         # robust metric only
fig_error_metrics(labels = TRUE, save = TRUE)
```

Saves `fig_error_metrics_scn<scn>_<Ntag>_<Struct>_<metrics>.{png,pdf}`.

---

## Planned figures (to be added)

- `viz_common.R` — promote the shared theme / palette / labellers / `save_fig()`
  helper once ≥ 2 figures reuse them.

