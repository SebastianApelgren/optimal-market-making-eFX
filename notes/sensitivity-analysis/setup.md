# Sensitivity Analysis & Uncertainty Quantification — Setup

## Goal

Study how uncertain model parameters affect the optimal market-making strategy produced by the ODE approximation (Riccati system). Determine which parameters the solution is most sensitive to, and quantify how parameter uncertainty propagates to the quantities a market maker cares about.

This follows the thesis structure discussed with Anders (2026-03-29): the sensitivity/UQ chapter sits between the ODE approximation and the PDE solver, and the PDE is later used to validate the ODE-based conclusions in low dimension.

## Currency setup
 
We work with the **2-currency case (USD, EUR)** — the simplest non-trivial model with one currency pair. This keeps the parameter count manageable, the results easy to visualize, and allows direct validation against the 2D PDE solver.

Extension to 3+ currencies is a natural follow-up: it introduces cross-pair correlations (rho) and cross-Hamiltonian coupling, adding new sensitivity dimensions that don't exist in the 2-currency case.

## Parameters under study

We study 8 parameters, using the paper's Table 1 values (EUR/USD pair) as nominal. Where the paper specifies per-size arrival rates (5 buckets), we group them into a single scale factor to reduce dimension.

| # | Parameter | Symbol | Nominal value | Range | Distribution |
|---|-----------|--------|---------------|-------|-------------|
| 1 | EUR volatility | sigma_EUR | 80 bps | [56, 104] bps | Uniform |
| 2 | Risk aversion | gamma | 20 | [10, 40] | Uniform |
| 3 | Arrival rate scale | lambda_scale | 1.0 | [0.5, 1.5] | Uniform |
| 4 | Logistic shift, tier 1 | alpha_1 | -1.9 | [-2.5, -1.3] | Uniform |
| 5 | Logistic shift, tier 2 | alpha_2 | -0.3 | [-0.9, 0.3] | Uniform |
| 6 | Logistic slope, tier 1 | beta_1 | 11 * 1e4 | [7.7, 14.3] * 1e4 | Uniform |
| 7 | Logistic slope, tier 2 | beta_2 | 3.5 * 1e4 | [2.45, 4.55] * 1e4 | Uniform |
| 8 | Execution cost (quadratic) | eta | 1e-5 bps | [0.5e-5, 2e-5] bps | Uniform |

### Rationale for parameter selection

- **sigma, lambda, alpha, beta**: These are estimated from market data and carry genuine estimation uncertainty. Volatility varies day-to-day; arrival rates depend on market conditions; the logistic demand parameters (alpha, beta) are fitted from historical fill/reject data via logistic regression and are never known exactly.
- **gamma**: A design choice by the bank (how risk-averse the strategy should be), but in practice it is tuned and there is no single "correct" value. Treating it as uncertain captures the sensitivity of the strategy to this subjective choice.
- **eta**: The quadratic execution cost parameter, related to market impact. Hard to estimate precisely in practice.

### Parameters not varied (and why)

- **psi (hedge dead zone)**: Set by the execution desk, relatively well known, and enters the model in a simple thresholding way. Expected to matter only for the hedge rate near the dead zone boundary.
- **kappa (terminal penalty)**: Set to 0 in the paper's example and acts as a boundary condition at T. Could be studied separately but is less interesting for the sensitivity of the running strategy.
- **T_days (time horizon)**: A design parameter, fixed at 0.05 days. Not uncertain in the same sense as the others.
- **mu (drift)**: Set to 0 in the paper. Could be included if drift estimation is relevant.
- **k (exogenous trade rate)**: Very small values in the paper (5e-3 to 8e-3 bps), negligible effect.

### Grouping choices

- **Arrival rates** (5 size buckets: 1, 5, 10, 20, 50 M$) are scaled by a single multiplier lambda_scale rather than varied independently. This reflects the realistic scenario that market activity rises or falls broadly, not that 1 M$ trades become more frequent while 50 M$ trades become less frequent. If needed, size-specific effects could be explored in a second pass.
- **Tier parameters** (alpha, beta) are varied independently per tier, since the two tiers represent fundamentally different client segments (aggressive vs. passive) with separate calibrations.

## Quantities of interest

We compute two sets of QoIs: policy-level quantities that probe the optimal strategy directly, and economic quantities that measure what the desk actually manages to.

### Set A — Policy QoIs

These probe the shape of the optimal policy at specific inventory levels.

1. **delta_star(y=0, tier 1) - delta_star(y=0, tier 2)** — The tier spread differential at zero inventory, trade size z = 1 M$. This measures how much wider the model quotes for passive clients (tier 2) vs. aggressive clients (tier 1). Banks care deeply about this: it is the basis of their client segmentation / tiered pricing. Note: tier 2 (passive) has lower alpha and beta, so tier 2 markup > tier 1 markup, making this quantity negative in sign. We report it as tier 1 minus tier 2 (a negative number reflecting the discount given to aggressive clients).

2. **delta_star(y=10) - delta_star(y=0)** — The inventory skew for tier 1, z = 1 M$. How much the markup changes when the market maker holds 10 M$ EUR. Measures the strategy's responsiveness to inventory risk. Parameters that affect the value function curvature (gamma, sigma) should dominate here.

3. **xi_star(y=10)** — The optimal hedge rate when holding 10 M$ EUR. Measures how aggressively the model hedges inventory. Directly influenced by eta (execution cost) and by the gradient of the value function.

### Set B — Economic QoIs

These aggregate across tiers and sizes to measure the economic outcomes the desk manages to.

4. **R(y=0)** — Instantaneous expected revenue rate at zero inventory:
   R = sum over tiers t, sizes j of: lambda_j * f(delta_star_{t,j}(y=0)) * delta_star_{t,j}(y=0) * z_j.
   This is the strategy's "earning power" when flat. It naturally weights each parameter by its economic importance (large-lambda tiers and large-z buckets contribute more), avoiding arbitrary tier/size choices.

5. **L(xi_star(y=10))** — Hedging cost rate at y = 10 M$ EUR: L = psi * |xi_star| + eta * xi_star^2.
   What the desk pays to unwind inventory through the interbank market.

6. **R(y=10) - L(xi_star(y=10))** — Net revenue rate at y = 10 M$ EUR. Revenue minus hedging cost under inventory pressure. The single number a desk head cares about most: "given that I'm carrying inventory, am I still making money after hedging costs?"

### Why these QoIs

- All six are **deterministic** given the ODE solution — no Monte Carlo simulation noise. This keeps the Sobol analysis clean.
- **Set A** probes the policy shape: client discrimination (QoI 1), inventory responsiveness (QoI 2), and hedging aggressiveness (QoI 3).
- **Set B** probes economic outcomes: earning power (QoI 4), hedging cost (QoI 5), and bottom line (QoI 6).
- A key expected result: the two sets may give different sensitivity rankings. A parameter could barely move the markup but strongly affect revenue (or vice versa), and that contrast is itself a useful practical insight.
- Both sets should appear in the main text. They can be presented compactly: one 2×3 subplot figure for Sobol indices (one column per QoI, one row per set), one 2×3 subplot figure for forward UQ densities. If space becomes a concern, Set A can be moved to an appendix — but the default should be to keep both in the main chapter.

## Methods

### Sobol indices (global sensitivity analysis)

We compute first-order (S_i) and total-effect (S_Ti) Sobol indices using Saltelli's Monte Carlo estimator (Sobol, 2001). The algorithm:

1. Draw two independent uniform sample matrices A, B of shape (N, k), k = 8 parameters.
2. For each parameter i, construct C_i = A with column i replaced from B.
3. Evaluate the model (ODE solve + QoI extraction) at all A, B, C_i samples: N(k+2) evaluations total.
4. Compute S_i and S_Ti from the sample covariances.

The first-order index S_i measures the fraction of output variance explained by parameter i alone. The total-effect index S_Ti includes all interactions involving parameter i. If S_Ti is close to zero, the parameter can be safely fixed.

Implementation is done from scratch (not using SALib), following the formulas in Sobol (2001) and Saltelli et al. (2008).

### Forward uncertainty quantification

After identifying the important parameters, we:

1. Sample parameters from their prior distributions (uniform on the ranges above).
2. Solve the ODE for each sample.
3. Compute QoIs for each sample.
4. Estimate the probability density functions of the QoIs using kernel density estimation (Gaussian kernel, Silverman's rule for bandwidth).
5. Report mean, standard deviation, and confidence intervals.

This shows the "spread" of possible market-making strategies given parameter uncertainty — the forward propagation from parameter space to strategy space.

## Key references

- Sobol', I.M. (2001). "Global sensitivity indices for nonlinear mathematical models and their Monte Carlo estimates." Mathematics and Computers in Simulation, 55(1-3), 271-280.
- Saltelli, A. et al. (2008). Global Sensitivity Analysis: The Primer. Wiley.
- Jakobsson, P.H. and Warnberg, A. (2023). "SIR-models and uncertainty quantification." KTH Bachelor thesis. [Same supervisor (Szepessy), same methodology applied to epidemiological ODE models.]
- Piazzola, C., Tamellini, L., and Tempone, R. (2021). "A note on tools for prediction under uncertainty and identifiability of SIR-like dynamical systems for epidemiology." Mathematical Biosciences, 332, 108514.
- Barzykin, A., Bergault, P., and Gueant, O. (2023). "Algorithmic market making in foreign exchange cash markets." arXiv:2207.04100v4. [The base model.]

## Connection to the full thesis

The sensitivity analysis serves two purposes:

1. **Practical insight**: Which parameters does a bank need to estimate carefully, and which can be set approximately without meaningfully affecting the strategy?

2. **Validation bridge**: The ODE-based sensitivity results can be spot-checked against the PDE solver in the 2-currency case. If both agree on which parameters matter, this validates the ODE approximation for sensitivity purposes and justifies using the cheap ODE for the full 5-currency UQ (where the PDE is computationally infeasible).
