# PDE vs ODE Comparison: Design and Methodology

Design notes for the PDE chapter of the thesis. The PDE comparison validates the
ODE approximation so that the reader trusts the ODE-based sensitivity analysis.

## Role in the thesis

The main contribution is the 3-currency sensitivity analysis (20 parameters,
6 QoIs, Sobol indices + forward UQ). The PDE chapter supports it by showing
that the ODE approximation is accurate enough for the SA conclusions to hold.

The PDE solver is limited to d=2 (2 currencies, 1 pair) because d=3 requires
301^3 = 27.3M grid points — infeasible. See `pde_dimension_discussion.md` for
the full argument. The comparison must therefore bridge from d=2 PDE validation
to d=3 SA credibility.

## What the previous comparison did

The first version (`pde_ode_comparison.ipynb`, April 2026) ran 9 PDE solves:
1 nominal + 4 parameters x 2 edges. It tested 3 QoIs (neutral spread, inventory
skew, hedge rate) and compared absolute values + OAT sensitivity magnitudes.

Findings:
- Markup: ~0.2% agreement across all configs
- Inventory skew: 3-15% relative error
- Hedge rate: 0-76% relative error, worst at lambda extremes

Problems with this version:
1. QoIs don't align with the SA — the SA reports tier spread differential and
   net revenue, the PDE comparison tests single-tier spread and omits revenue.
2. Only 4 of 8 available parameters tested. The selection was justified by the
   2-currency SA, but that SA doesn't appear in the report.
3. Only tests parameter edges (OAT), missing multi-parameter interactions.
4. Comparison metrics limited to point values and bar charts.
5. No attempt to bridge to d=3 beyond the theoretical argument.

## Improved comparison

### QoIs (5 quantities, aligned with the SA)

Match the 3-currency SA's QoIs as closely as possible in d=2:

| # | QoI | SA match | Notes |
|---|-----|----------|-------|
| 1 | Tier spread differential: delta*(tier1) - delta*(tier2) at y=0, z=1 M$ [bps] | SA QoI 1 | Direct match |
| 2 | Own-inventory skew: delta*(y_EUR=10) - delta*(y=0) on EUR/USD, tier 1, z=1 [bps] | SA QoI 2 | Direct match |
| 3 | Hedge rate: xi*(EUR/USD) at y=[0, 10] [M$/day] | SA QoI 4 | Direct match |
| 4 | Net revenue: R(y=[0,10]) - L(xi*) across tiers [$/day] | SA QoI 6 | Direct match |
| 5 | delta*(y=0) single tier [bps] | (diagnostic) | Useful as baseline; not in SA |

Cannot test in d=2 (structurally zero, explain in thesis):
- SA QoI 3: cross-inventory skew (A is diagonal in d=2, no cross-pair coupling)
- SA QoI 5: cross-pair hedge rate (only one pair exists)

This gives 4/6 SA QoIs validated. The 2 untestable QoIs are the ones that probe
cross-pair effects, which are exact in the Riccati algebra (see bridging section).

### Parameter selection: all 8, OAT

Test all 8 parameters available in d=2 (the same 8 as the 2-currency SA):

| # | Parameter | Nominal | Range | Width |
|---|-----------|---------|-------|-------|
| 1 | sigma_EUR | 80 bps | [56, 104] | +/-30% |
| 2 | gamma | 20 | [10, 40] | +/-50% |
| 3 | lambda_scale | 1.0 | [0.5, 1.5] | +/-50% |
| 4 | alpha_1 | -1.9 | [-2.5, -1.3] | +/-30% |
| 5 | alpha_2 | -0.3 | [-0.9, 0.3] | +/-100% |
| 6 | beta_1 | 11.0 1/bps | [7.7, 14.3] | +/-30% |
| 7 | beta_2 | 3.5 1/bps | [2.45, 4.55] | +/-30% |
| 8 | eta | 1e-5 bps | [0.5e-5, 2e-5] | +/-50% |

Total: 17 PDE solves (1 nominal + 2 edges x 8 parameters).
Estimated compute: ~3 hours.

Doing all 8 avoids any parameter selection bias and produces an additional
finding: the demand curve parameters (alpha, beta) contribute negligible ODE
error compared to the market condition parameters (sigma, gamma, lambda, eta).
This makes structural sense — the quadratic-H approximation quality depends on
how far |p| deviates from zero, which is controlled by the Riccati dynamics
(sigma, gamma) and the Hamiltonian's curvature (lambda, eta), not the demand
curve shape.

### Latin hypercube sample (multi-parameter corners)

OAT moves one parameter at a time and cannot detect whether ODE error compounds
when multiple parameters are simultaneously extreme. To test this, draw 15-20
points from a Latin hypercube over the 8-dimensional parameter space and solve
the PDE at each.

Total: 15-20 additional PDE solves.
Estimated compute: ~3 hours.

This tests whether the ODE error at (sigma_high, gamma_high, eta_low) — a
multi-parameter corner where the quadratic-H is under pressure from all
directions — is worse than any single-parameter edge.

If it turns out that the LHS errors are bounded by the worst OAT errors, the
OAT analysis is sufficient and the LHS confirms it. If the LHS finds worse
corners, those are important to report.

### Comparison metrics

#### 1. Summary error table

For each parameter config, report a table with:
- ODE value, PDE value, absolute error, relative error
- Across all 5 QoIs

Present both max-over-configs and mean-over-configs for each QoI.

#### 2. Sensitivity rankings

For each QoI, rank the 8 parameters by their OAT sensitivity magnitude
|QoI(high) - QoI(low)|. Do this for ODE and PDE separately. Report:
- The rankings side by side
- Spearman rank correlation rho between ODE and PDE rankings

If rankings match (rho close to 1), the SA conclusions hold even where
magnitudes differ. This is the single most important metric: a reader
who sees "76% hedge rate error" will worry, but "rankings agree with
rho = 0.95" is reassuring.

#### 3. Parity plots

For each QoI, scatter-plot ODE value vs PDE value across all configs
(OAT + LHS). The 45-degree line shows perfect agreement. Deviations
show systematic bias (e.g., ODE over-predicts hedging at low lambda).

Report both absolute and relative errors. The hedge rate at lambda_high
has 76% relative error but the absolute error (501 M$/day) is smaller
than at nominal (302 M$/day) — the percentage is inflated by the small
denominator.

#### 4. Value function comparison (theta profile)

Plot theta_ODE(0, y) vs theta_PDE(0, y) along the y_EUR axis (y_USD = 0)
at nominal parameters. This shows:
- Where on the inventory grid the approximation is good (near y=0)
  vs poor (extreme y)
- The error profile is expected to look cubic (quadratic-H error is
  O(p^3), and p is linear in y through the Riccati)

This localises the error to its source and gives the reader geometric
intuition for why downstream QoIs have different error levels.

#### 5. Hessian comparison (A matrix)

The ODE ansatz is theta = -y^T A y - y^T B - C, so A encodes the
value function's curvature. Extract the PDE's effective Hessian at y=0:

  A_PDE approx -(1/2) (d^2 theta / dy^2)|_{y=0}

via central finite differences on the PDE grid. Compare A_ODE vs A_PDE.
In d=2 this is a single number (the EUR-EUR entry). This tests the
approximation at its source — all downstream QoIs flow from A and B.

## Handling the hedge rate error

The hedge rate has the largest ODE errors (up to 76% at lambda_high).
The thesis should address this directly, not hide it.

### Why the error is large

The hedge rate is xi* = -(p - psi * sign(p)) / (2 eta) when |p| > psi.
The 1/(2 eta) amplification factor (~5 x 10^8) converts small errors
in p into large errors in xi*. And p = -(2Ay + B) dot d, so any error
in A propagates linearly into p and then gets amplified.

The quadratic-H approximation makes A slightly wrong. For the markup,
delta* depends on p through the Lambert W, which is sublinear — errors
in p are damped. For the hedge rate, the linear/amplified dependence on
p means errors are magnified.

### Why the error is less alarming than it looks

1. The worst case (76% at lambda_high) has absolute values -158 vs -659
   M$/day, both small compared to nominal (-1879). In practical terms:
   both say "barely hedge."

2. The ODE error is systematic, not random. It over-predicts hedging at
   low lambda and under-predicts at high lambda. A parity plot will show
   this clearly.

3. The SA samples parameters uniformly. The 76% error occurs at parameter
   edges that have low sampling weight (most samples are near the
   centre). The median ODE error on xi* across the full SA distribution
   is likely 10-15%, not 76%.

### What to report

- Acknowledge the large relative errors at extremes.
- Show that absolute errors tell a different story.
- Show that sensitivity rankings agree despite magnitude differences.
- Frame the SA's hedge rate conclusions as qualitative: which parameters
  matter and interact, not exact numbers.

## Bridging d=2 PDE to d=3 SA

The report presents a 3-currency SA but validates with a 2-currency PDE.
Three arguments bridge this gap.

### Argument 1: theoretical (mathematical)

The ODE approximation error has a single source: the quadratic Taylor
expansion of H(p) around p=0, made independently per pair/tier/size.
Cross-pair interactions enter through the Riccati matrix algebra (Eq. 5),
which is exact — no approximation added when going from 1 to 3 pairs.

The Riccati ODE has the form A_dot = F(A) where F involves:
- Sigma (covariance): exact, enters through -gamma/2 Sigma
- M, M_tilde (arrival rate matrices): exact, linear combinations of lambda
- alpha_2 coefficients (Hamiltonian curvature): approximate, per-pair

Each pair contributes its own alpha_2 terms to the diagonal blocks of F(A).
The off-diagonal blocks (A[EUR,GBP] coupling) depend only on Sigma and M,
which are exact. Therefore: going from d=2 to d=3 adds exact cross-pair
coupling but no new approximation error.

Write this argument explicitly in the thesis (1 paragraph).

### Argument 2: empirical, Riccati matrix inspection (cheap, no PDE needed)

Solve the 3-currency ODE. Extract A(0). Report:
- Diagonal: A[EUR,EUR], A[GBP,GBP]
- Off-diagonal: A[EUR,GBP]
- Ratio: |A[EUR,GBP]| / A[EUR,EUR]

If cross-pair coupling is 5-10% of diagonal, it's a small perturbation on
per-pair dynamics already validated by the d=2 PDE.

### Argument 3: empirical, 2-ccy vs 3-ccy ODE comparison (cheap, no PDE needed)

Compute EUR/USD QoIs from both the 2-currency and 3-currency ODE at nominal
parameters. The difference quantifies how much adding GBP perturbs EUR/USD
results. If small (< a few percent), the d=2 PDE validation carries over:
the per-pair error (validated) dominates the cross-pair perturbation (small).

Could extend this to a sweep: compute the 2-ccy vs 3-ccy ODE difference
across the same parameter configurations used in the PDE comparison. If the
cross-pair perturbation is everywhere smaller than the ODE-PDE error, the
d=2 comparison already captures the dominant effect.

## Compute plan

### What to run

1. OAT solves: 17 PDE solves (1 nominal + 2 edges x 8 params). ~3 hours.
2. LHS solves: 15-20 PDE solves at Latin hypercube points. ~3 hours.
3. ODE evaluations: milliseconds, run inline.
4. Riccati matrix / 2-vs-3 ODE: milliseconds, run inline.

Total PDE compute: ~6 hours, run once.

### What to cache (data/ folder)

Save to `data/pde_comparison/`:
- `oat_configs.npz`: parameter vectors for all 17 OAT configs
- `oat_theta0/config_name.npy`: theta(0, y) grid for each config (301x301 float64 ~ 700 KB each)
- `oat_qois.npz`: extracted QoI values for ODE and PDE, all configs
- `lhs_configs.npz`: parameter vectors for all LHS points
- `lhs_theta0/lhs_00.npy` through `lhs_19.npy`: theta grids
- `lhs_qois.npz`: QoI values
- `metadata.json`: grid parameters, PDE solver settings, timestamp

Total disk: ~25 MB.

### Notebook structure

The analysis notebook loads from `data/pde_comparison/` and produces:
1. Summary error table (all QoIs, all configs)
2. Sensitivity rankings + Spearman rho
3. Parity plots (ODE vs PDE, each QoI)
4. Theta profile along y_EUR
5. Hessian comparison (A_ODE vs A_PDE)
6. Bridging analysis (Riccati matrix, 2-ccy vs 3-ccy)
7. Discussion of hedge rate errors

No PDE solves in the notebook itself — it stays fast and re-runnable.

A separate compute script (`scripts/run_pde_comparison.py`) handles the PDE
solves and writes to `data/`. This script can be run once and never again
unless parameters change.

## Connection to other notes

- Solver implementation: `implicit/eta_continuation_and_monotone_scheme.md`
- Why d=2: `pde_dimension_discussion.md`
- Grid and domain: `grid_and_domain.md`
- SA setup: `../sensitivity-analysis/three_currency_setup.md`
- SA results (2-ccy): `../sensitivity-analysis/results.md`
