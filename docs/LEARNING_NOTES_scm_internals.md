# Learning Notes — nlmixr2 SCM Internals

Personal study notes on the internals of `nlmixr2scm::runSCM` and the
surrounding estimation machinery. **Not for publication** — these are
working notes to record how the code actually behaves, verified against
`R/scm.R` and `script/estimator_factory.R`.

---

## 1. The two-tier (screening → final) candidate fit

**Question:** Is the two-tier scheme a guard against FOCE-I's noisy likelihood
surface?

**Answer:** No. It is **purely a run-time optimisation**.

The two tiers use **identical estimation settings** — the same `sigdig`,
`derivEps`, and ODE/sensitivity tolerances. The *only* differences are:

| Tier | `covMethod` | `calcTables` | Purpose |
|------|-------------|--------------|---------|
| **screen** | `""` (skip) | `FALSE` | Fast; only the OFV is needed for the LRT |
| **final** | `"r,s"` | `TRUE` | Adds covariance step + diagnostic tables for the winner |

Because the covariance step and output-table construction are post-hoc
reporting work (not part of *finding* the estimate), skipping them during
screening changes nothing about the estimate — it just saves time. Only the
single winning candidate at each step pays for the covariance/tables.

*Source:* `script/estimator_factory.R` lines 29–36; corroborated by
`Trials&errors/smoke_irlsfoceif_vae_scn16.R` lines 113–114.

---

## 2. Default initial estimate for a covariate effect

**Question:** What is the default initial estimate for a covariate coefficient?

**Answer:** In `forwardSearch()` / `backwardSearch()`, when no per-pair `init`
is supplied, the hard-coded default is:

```r
cov_init  = 0.1     # estimate (on the coefficient / log scale)
cov_lower = -5      # lower bound
cov_upper =  5      # upper bound
```

- Every candidate covariate coefficient starts at **0.1**, bounded to
  **[−5, 5]**, regardless of shape (power / linear / categorical).
- The `inits` argument can override this per shape. The docstring *suggests*
  `list(power = 0.1, lin = 0.01, cat = 0.01)`, but absent an override it is a
  flat **0.1** for all shapes.
- Strictly, this is the `runSCM` default rather than a deeper `nlmixr2est`
  default.

*Source:* `R/scm.R` ~line 1885.

---

## 3. What "BSV fixed" means in the 1-D frozen-base profile

**Question:** In the frozen-base 1-D profile warm-start, the between-subject
variability (BSV) was fixed to what?

**Answer:** BSV (omega) is **fixed at the parent model's current estimates —
present, but not re-estimated. It is NOT zeroed.**

From `.freezeUiForProfile()`:

```r
ini$fix[isTheta] <- TRUE                              # fix all structural thetas
ini$fix[isTheta & ini$name == freeTheta] <- FALSE     # except the new covariate
isOmega <- !is.na(ini$neta1)
if (any(isOmega)) ini$fix[isOmega] <- TRUE            # fix omega at parent values
```

**Why not zero it?** Zeroing omega misspecifies the model and can **flip the
sign** of the profiled coefficient. Keeping the random effects present (FOCE-I
still integrates over them) but frozen means the single free covariate theta is
profiled under the *correct* mixed-effects likelihood.

**When does this fire?** Only when a forward candidate *stalls* — its nested
model fails to improve on its parent (ΔOFV ≤ `stallTol`, default `0`), which at
a genuine optimum cannot happen and therefore signals the outer optimiser never
left the flat zero-effect init. The profiled value (found by a bounded scalar
search; boundary solutions are rejected) then warm-starts a rescue refit.

*Source:* `R/scm.R` lines 1643–1657 (`.freezeUiForProfile`), 1690+
(`.profileCovInit`).

---

## 4. Glossary of the shared control parameters

These are the settings that are held **identical** across the screening and
final tiers, so the likelihood surface the optimiser "sees" is the same in both.

| Parameter | What it controls | Intuition |
|-----------|------------------|-----------|
| **`sigdig`** | Significant digits of the **outer** optimiser's convergence target (translated internally into stopping tolerances on the OFV and parameters). | "Stop when the OFV / parameters are stable to ~N significant figures." Higher `sigdig` = tighter, more iterations, more precise. |
| **`derivEps`** | Step size(s) for **finite-difference gradients** of the objective w.r.t. the population parameters (the outer problem). Usually a length-2 vector (relative, absolute step). | How far the optimiser nudges each θ to numerically estimate ∂OFV/∂θ. Too large → biased gradient; too small → amplifies numerical noise. |
| **`atol` / `rtol`** (ODE tolerances) | Absolute and relative error tolerances for the **ODE solver** integrating the structural model over time for each subject. | How accurately the concentration–time curve itself is solved. `atol` guards small values near zero; `rtol` controls relative accuracy of larger values. Looser → faster but noisier likelihood. |
| **`atolSens` / `rtolSens`** (sensitivity tolerances) | Same absolute/relative tolerances but for the **sensitivity equations** — the ODEs propagating ∂(state)/∂(parameter) needed for **analytic gradients** (`foceif` and other analytic-gradient estimators). | How accurately the *derivatives* of the ODE solution w.r.t. parameters are computed. Often as tight as (or tighter than) the state tolerances, because gradient noise is what stalls gradient optimisers. |

**Tie-back to §1:** The estimate an optimiser converges to is determined by
`sigdig`/`derivEps` (outer optimisation precision) and
`atol`/`rtol`/`atolSens`/`rtolSens` (how cleanly the inner ODE + gradients are
solved). Since these are identical across tiers, the *estimation* is identical
— the screening tier only drops the post-hoc covariance step and output tables.
That is exactly why the two-tier design is a **speed** optimisation, not an
accuracy safeguard.

---

## 5. Related facts worth remembering

- **Selection-pattern heatmaps** plot per-covariate **error rates**:
  false-negative rate for the four true effects (missed, or recovered in the
  wrong shape) and false-positive rate for the BMI/RACE distractors —
  *not* TP/FP selection frequencies. (`n_FN / n_datasets`, `n_FP / n_datasets`.)
- **RMRSE / MARE** are, by default, **unconditional** — computed over *all* fits
  (finite estimate), following Khandelwal et al. (2019), not restricted to
  converged fits. Strict-converged and selection-conditional (δ = 1) variants
  are secondary sensitivity analyses.
- **`stallTol = 0`** by default: a nested model can never be genuinely worse
  than its parent at a true optimum, so ΔOFV ≤ 0 is a reliable "stalled" signal.

---

## 6. The 1-D profile warm-start: theory → implementation

**Question:** The `.profileCovInit()` warm-start is a 1-D profile. Is it a Brent
method — fix all other parameters at their last estimates and estimate only the
one parameter of interest?

**Answer:** The *structure* is a profile likelihood (fix everything else, free
one parameter), but the **optimiser is BOBYQA, not Brent**, and each evaluation
is a **full FOCEi fit** — not a fixed-effects-only 1-parameter regression.

### 6.1 Theory — what a profile warm-start is here

The goal is a good *starting value* for a new covariate coefficient before the
candidate's real fit runs. At the flat zero-effect init the outer objective has
(near) zero gradient, so gradient optimisers (`nlminb` / `lbfgsb3c`) stall. A
1-D profile of the outer objective — holding every other population parameter at
the parent estimate and moving only the new covariate coefficient — locates a
non-flat, gradient-informative value to seed the real fit.

### 6.2 Implementation — three deliberate choices

1. **Population thetas fixed at parent** (`.freezeUiForProfile`): all θ
   `fix = TRUE` except the one new covariate coefficient. This is the classic
   "profile" restriction.

2. **BOBYQA, not Brent, for the outer search**
   (`foceiControl(outerOpt = "bobyqa")`): BOBYQA is Powell's derivative-free
   trust-region method (bounded, quadratic interpolation model). It is chosen
   because it is *already exposed* as the FOCEi outer optimiser **and** because —
   being derivative-free — it does **not** stall at the flat zero-effect point
   the way a gradient method does. Brent would share that robustness, but BOBYQA
   is what is wired in. BOBYQA is only the *profiler*; it is never the estimator
   of record.

3. **BSV fixed but PRESENT, not zeroed** (see §3): the fit is a genuine FOCEi
   fit, so the **inner loop re-estimates each subject's EBEs (η) at every trial
   value** of the covariate coefficient. This profiles the covariate under the
   correct mixed-effects likelihood. Zeroing omega would misspecify the model and
   can flip the coefficient's sign.

### 6.3 Two EBE re-estimations, one number handed forward

There are **two distinct fits**, and EBEs are re-estimated in *both*:

| Stage | Outer free params | Inner EBEs | Purpose |
|-------|-------------------|-----------|---------|
| **1. 1-D profile** | only the new covariate θ (others fixed, ω fixed-present) | re-estimated at every outer trial | produce a starting value |
| **2. formal candidate fit** | all θ + ω free | re-estimated (normal FOCEi) | the actual estimate |

The profile does **not** reuse EBEs as the warm start — it passes forward **only
the single profiled coefficient** as the starting init. Stage 2 recomputes its
own EBEs from scratch.

### 6.4 Boundary rejection

A profiled value that is non-finite or pinned to (within `edge_tol` of) the
parameter bracket carries no gradient information, so it is rejected
(`NA_real_`) and the caller falls back to the default init. `edge_tol` scales
with the bracket width, so "essentially on the bound" (e.g. `1.9999997` for an
upper bound of `2`) is also rejected.

### 6.5 Where it fires

This warm-start is invoked by **Mechanism B (`profileInitOnStall`, default
`TRUE`)** — the on-stall rescue that refits from the profiled init and keeps the
result *only if it strictly improves ΔOFV*. The unconditional pre-fit variant,
**Mechanism A (`profileInit`, default `FALSE`)**, calls the *same*
`.profileCovInit()` for every forward candidate but is off by default. Both are
forward-only (`add == TRUE`); backward search never profiles.

*Source:* `R/scm.R` `.freezeUiForProfile()` (~L1690), `.profileCovInit()`
(~L1737–1793); Mechanism A ~L1996–2019; Mechanism B ~L2213–2295. See also
`presentation/_diagrams/profileInit_warmup_mechanismB.md`.
