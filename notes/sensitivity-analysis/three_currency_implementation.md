# 3-Currency Sensitivity Analysis — Implementation Plan

## Overview

Create a new notebook `sensitivity_analysis.ipynb` for the 3-currency (USD, EUR, GBP) sensitivity analysis. This is the version presented in the thesis. The previous 2-currency notebook has been renamed to `sensitivity_analysis_2ccy.ipynb` and stays in the repo as development/exploration work.

The Saltelli/Sobol infrastructure is extracted into `src/sobol.py` so it can be shared and is presentable as library code.

No changes to the existing `src/` model modules are needed — the Riccati solver, policy extraction, and model infrastructure already handle d=3 correctly.

## Naming convention

Since the 3-currency analysis is the one shown in the thesis, all functions use clean generic names:

- `build_modified_params()` (not `build_modified_params_3ccy`)
- `_revenue_rate()` (not `_revenue_rate_3ccy`)
- `evaluate_qois()` (not `evaluate_qois_3ccy`)

The 2-currency notebook keeps its own `build_modified_params()` etc. in its own scope. No name collision.

## File changes

### New files

| File | Purpose |
|------|---------|
| `src/sobol.py` | Saltelli sampling, evaluation, Sobol index computation |
| `sensitivity_analysis.ipynb` | The sensitivity analysis notebook (thesis version, 3-currency) |

### Modified files

| File | Change |
|------|--------|
| `src/__init__.py` | Add imports from `src/sobol` |

### Untouched

| File | Why |
|------|-----|
| `sensitivity_analysis_2ccy.ipynb` | Previous 2-currency development notebook, renamed, stays as-is |
| `src/model.py` | Already handles d=3, PairParams, correlations |
| `src/riccati.py` | Already handles d×d matrices, correlation in Sigma |
| `src/policy.py` | Already handles arbitrary d, direction vectors |
| `src/hamiltonian.py` | Pure scalar functions, no dimension dependence |

## src/sobol.py

Extract three functions currently in `sensitivity_analysis.ipynb` cell `670125a5`:

```
saltelli_sample(N, ranges, seed=42) → (A, B)
evaluate_saltelli_samples(model_func, A, B) → (f_A, f_B, f_C)
compute_sobol_indices(f_A, f_B, f_C) → (S, S_T)
```

These are generic (no model-specific code). ~50 lines total.

## Notebook structure: sensitivity_analysis.ipynb

### Cell 1 — Title + imports

```markdown
# Sensitivity Analysis of ODE Market-Making Model

Global sensitivity analysis (Sobol indices) and forward uncertainty quantification
for the 3-currency (USD, EUR, GBP) ODE approximation.

See `notes/sensitivity-analysis/three_currency_setup.md` for methodology.
```

Imports: numpy, matplotlib, scipy.stats.gaussian_kde, tqdm, time, src model classes, src.sobol functions.

### Cell 2 — Parameter and QoI definitions

20 parameters, 6 QoIs. Constants: PARAM_NAMES, PARAM_LABELS, NOMINAL, RANGES, QOI_NAMES, N_PARAMS, N_QOIS.

Parameter ordering:
```
 0  sigma_EUR     global
 1  sigma_GBP     global
 2  rho           global
 3  gamma         global
 4  lambda_scale  global
 5  alpha_1       EUR/USD
 6  alpha_2       EUR/USD
 7  beta_1        EUR/USD
 8  beta_2        EUR/USD
 9  eta           EUR/USD
10  alpha_1       GBP/USD
11  alpha_2       GBP/USD
12  beta_1        GBP/USD
13  beta_2        GBP/USD
14  eta           GBP/USD
15  alpha_1       EUR/GBP
16  alpha_2       EUR/GBP
17  beta_1        EUR/GBP
18  beta_2        EUR/GBP
19  eta           EUR/GBP
```

Also define plotting metadata:
```python
PARAM_GROUPS = ["Global"]*5 + ["EUR/USD"]*5 + ["GBP/USD"]*5 + ["EUR/GBP"]*5

GROUP_COLORS = {
    "Global":  "...",
    "EUR/USD": "...",
    "GBP/USD": "...",
    "EUR/GBP": "...",
}
```

### Cell 3 — Model evaluation functions

**build_modified_params(params):**

Maps 20-param vector → ModelParams with currencies=["USD","EUR","GBP"], 3 pairs.

- Cache the base 3-ccy model: `restrict_currencies(build_paper_example_params(), ["USD","EUR","GBP"])`
- Pull sizes_musd and psi from each pair's base PairParams (these are fixed, not varied)
- Build 3 PairParams with tiers, lambda_scale, eta from the param vector
- Set corr = {("EUR","GBP"): rho}
- Set sigma = {"USD": 0.0, "EUR": sigma_eur*BP, "GBP": sigma_gbp*BP}

**_revenue_rate(mp, res, y):**

Loop over all pairs in mp.pairs:
```python
revenue = 0.0
for (a, b), pp in mp.pairs.items():
    for ccy_pay, ccy_sell in [(a, b), (b, a)]:
        for t_idx, tier in enumerate(pp.tiers):
            for z, lam in zip(pp.sizes_musd, pp.lambdas_per_day):
                delta = res.markup(t_idx, ccy_pay, ccy_sell, z, y)
                f = logistic_f(delta, tier.alpha, tier.beta)
                revenue += lam * f * delta * z
return revenue * 1e6
```

**evaluate_qois(params):**

Inventory states (3D):
```python
y_flat     = np.zeros(3)           # [USD, EUR, GBP] = [0, 0, 0]
y_long_eur = np.array([0, 10, 0])  # 10 M$ long EUR
y_long_gbp = np.array([0, 0, 10])  # 10 M$ long GBP
```

6 QoIs:
```python
# Set A: quoting policy
# 1. Tier spread differential (EUR/USD, y=0, z=1 M$)
tier_diff = (res.markup(0, "EUR","USD", 1.0, y_flat)
           - res.markup(1, "EUR","USD", 1.0, y_flat)) / BP

# 2. Own-inventory skew (EUR/USD, tier 1, z=1)
own_skew = (res.markup(0, "EUR","USD", 1.0, y_long_eur)
          - res.markup(0, "EUR","USD", 1.0, y_flat)) / BP

# 3. Cross-inventory skew (EUR/USD quote, y_GBP=10)
cross_skew = (res.markup(0, "EUR","USD", 1.0, y_long_gbp)
            - res.markup(0, "EUR","USD", 1.0, y_flat)) / BP

# Set B: hedging and economics
# 4. Own-pair hedge rate (EUR/USD, y_EUR=10)
xi_own = res.hedge_rate("EUR", "USD", y_long_eur)

# 5. Cross-pair hedge rate (EUR/GBP, y_EUR=10)
xi_cross = res.hedge_rate("EUR", "GBP", y_long_eur)

# 6. Net revenue at y_EUR=10
R_long = _revenue_rate(mp, res, y_long_eur)
L_total = 0.0
for (a, b), pp in mp.pairs.items():
    xi = res.hedge_rate(a, b, y_long_eur)
    L_total += (pp.psi * abs(xi) + pp.eta * xi**2) * 1e6
net_revenue = R_long - L_total
```

### Cell 4 — Sanity check

- Evaluate at nominal, print all 6 QoIs
- Compare n_steps=500 vs 2000
- Time a single evaluation (expect ~12-15 ms)
- Print estimated Sobol time: t_per * 10000 * 22 / 60

### Cell 5 — Run Saltelli sampling

```python
N_SOBOL = 10_000
A, B = saltelli_sample(N_SOBOL, RANGES, seed=42)
f_A, f_B, f_C = evaluate_saltelli_samples(evaluate_qois, A, B)
```

Expected: 220,000 evaluations, ~45-55 min.

### Cell 6 — Compute and print Sobol indices

20 × 6 table. Print grouped by parameter category for readability.

### Cell 7 — Sobol bar plots

2×3 grid (Set A top row, Set B bottom row), 20 bars per panel.

Bars color-coded by parameter group. S_i and S_Ti shown side by side (lighter/darker shade of group color, or solid/hatched).

x-axis: short abbreviated labels, rotated. Within each group the names repeat (α₁, α₂, β₁, β₂, η), separated visually by color.

Save: `report/figures/fig_sa_sobol.pdf`

### Cell 8 — Forward UQ statistics

Reuse A∪B evaluations (20,000 samples). Print mean/std/5%/95%/nominal table.

### Cell 9 — Forward UQ KDE plots

2×3 grid. KDE density + nominal line + 90% CI shading.

Save: `report/figures/fig_sa_forward_uq.pdf`

## Plotting details

### 20-bar Sobol chart readability

With 20 bars the x-axis is crowded. The plan:

1. **Color-code by group** (4 colors): makes it immediately visible which parameter family dominates each QoI.
2. **Short labels**: σ_E, σ_G, ρ, γ, λ | α₁, α₂, β₁, β₂, η | α₁, α₂, β₁, β₂, η | α₁, α₂, β₁, β₂, η. The repeating names within each group are disambiguated by color.
3. **Wider figure**: figsize=(20, 10) or similar.
4. **Group legend**: one legend entry per color/group, not per bar.

If the vertical bars are still too crowded after the first render, switch to horizontal bar charts.

### Report figures

The sensitivity chapter will need:
- Parameter table (20 rows) — already in three_currency_setup.md
- Sobol indices: either one composite figure (2×3) or individual per-QoI bar plots + tables
- Forward UQ: one composite figure (2×3) + summary table
- All saved as PDF to `report/figures/`

Figure naming (no `3ccy` suffix since this is the only version in the report):
```
report/figures/fig_sa_sobol.pdf
report/figures/fig_sa_forward_uq.pdf
```

## Execution order

1. Create `src/sobol.py`, update `src/__init__.py`
2. Create notebook cells 1-4 (constants, model functions, sanity check)
3. Verify sanity check passes, check timing
4. Add cells 5-9 (Saltelli run, Sobol computation, plots, forward UQ)
5. Run notebook (~50 min computation)
6. Inspect results, save figures
7. Write results note: `notes/sensitivity-analysis/three_currency_results.md`
8. Update report chapter: `report/chapters/sensitivity.tex`

Steps 7-8 happen after we have results.
