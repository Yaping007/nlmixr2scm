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

## Planned figures (to be added)

- `fig_sc16_selection.R` — SC16 per-covariate TP / FP / FN
  (source: `vae_covsel_by_covar.csv`).
- `fig_error_metrics.R` — `RMRSE_pct` and `MARE_pct` by scenario / N / structure,
  faceted by `param_class` (source: `vae_estim_{all,success,cond}.csv`).
- `viz_common.R` — promote the shared theme / palette / labellers / `save_fig()`
  helper once ≥ 2 figures reuse them.
