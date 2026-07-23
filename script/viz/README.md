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

**Data source:** `output/vae_covsel_aggregated/vae_covsel_by_covar.csv`
(plus `vae_diag_rates.csv` for the per-cell dataset total `n_total`).
Columns used: `sample_N, scenario, structure, var, covar, is_true, n_FN, n_FP`.

**Two error regimes per (scenario × effect) cell:**

| regime | when | error shown | denominator |
|---|---|---|---|
| **FN** (false negative) | `is_true = TRUE` — a real effect | `n_FN / N_cell` | cell dataset total |
| **FP** (false positive) | `is_true = FALSE` — a null effect | `n_FP / N_cell` | cell dataset total |

> **Denominator gotcha.** `n_datasets` in the by-covar CSV is the number of
> datasets in which that `(var, covar)` pair *appeared*, **not** the cell size.
> A distractor only appears when it is falsely selected, so `n_datasets == n_FP`
> and `n_FP / n_datasets` is always ≈100 %. The figure therefore divides by
> `N_cell` (= `n_total` from `vae_diag_rates.csv`, equivalently `N` in
> `vae_power.csv`) for **both** regimes.

**Layout:**

- **x** = scenario (1–16); **y** = effect label `PARAM~COVAR`.
- y-order = CL block then VC block; **within each block `~BW` and `~BMI` are
  adjacent**, so a true `~BW` (bold outline = FN) sits directly above its
  `~BMI` thief (plain tile = FP) — leakage reads vertically.
- **fill** = error rate 0–100 % (white → dark red = worse).
- **bold black outline** = TRUE effect (its shading is the FN rate); plain tiles
  are null effects (shading = FP rate).
- **facet** = `sample_N` (columns). One `structure` per figure (linCmt ≈ ode).

### Signature

```r
fig_covsel_heatmap(agg_dir   = "output/vae_covsel_aggregated",
                   sample_N  = NULL,       # NULL = all (40, 80, 300); or subset
                   structure = "linCmt",   # ONE structure per figure ("linCmt" | "ode")
                   labels    = TRUE,       # print error % inside each tile
                   min_fp    = 0,          # blank FP cells below this rate (de-clutter)
                   save      = FALSE,      # TRUE also writes PNG + PDF
                   out_dir   = "output/figures/vae_covsel")
```

### Knobs

- **`sample_N`** — `NULL` shows all three N; pass e.g. `300` for the large cohort.
- **`structure`** — one structure per figure (default `"linCmt"`; `"ode"` to swap).
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
- `structural_intercept` — `TVCL, TVVc` (**caveat**: absorb the covariate
  centering shift, so their relative error is *not* pure bias — footnoted).
- `structural_other` — `TVQ, TVVp, var_CL, var_Vc, cov_VcCL, ResErr` (clean
  baseline; fixed `TVKA` is dropped by default since its error is always 0).

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

