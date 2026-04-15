# Sensitivity Analysis — 3-Currency Extension (USD, EUR, GBP)

## Motivation

The 2-currency analysis (setup.md, results.md) studied the simplest non-trivial model: one pair, no correlation, no cross-hedging. It identified which parameter *types* matter for which aspects of the strategy (demand curves for spreads, market conditions for inventory management, everything for hedging).

The 3-currency case introduces three qualitatively new features that cannot exist in 2 currencies:

1. **Cross-pair risk coupling via correlation.** The covariance matrix Sigma acquires off-diagonal entries: Sigma[EUR,GBP] = rho * sigma_EUR * sigma_GBP. This feeds into the Riccati ODE through the term -0.5 * gamma * Sigma, so the Riccati solution A develops non-zero off-diagonal entries A_{EUR,GBP}. Holding EUR now creates risk that depends on GBP volatility and the EUR/GBP correlation.

2. **Cross-inventory effects on quoting.** The markup formula p = ((2y + z*d) * A + B) * d sums over all currencies via the matrix product. When A_{EUR,GBP} != 0, the GBP inventory component y_GBP shifts the EUR/USD markup even though the trade itself involves only EUR and USD. In 2 currencies this cross-influence is structurally zero.

3. **Cross-hedging.** With 3 tradeable pairs (EUR/USD, GBP/USD, EUR/GBP), the market maker can hedge EUR exposure through the cross pair EUR/GBP instead of only through EUR/USD. The model's optimal hedge rates for all 3 pairs are determined simultaneously by the value function gradient 2*A*y + B. Whether this channel is economically significant, and what drives it, are open questions.

The 3-currency sensitivity analysis is self-contained: it includes all model parameters and lets the Sobol indices determine which matter. We do not pre-filter parameters based on the 2-currency results. The two analyses can then be compared side by side in the report, and any conclusions about which parameters are unimportant are demonstrated independently in each case.

## Currency setup

We use **USD, EUR, GBP** (3 currencies, 3 pairs). This is the smallest extension from 2 currencies that includes all three new features above. The paper's Table 1 provides nominal parameter values for all three pairs. EUR/USD and GBP/USD are the two most liquid G10 majors; EUR/GBP is a liquid cross pair with different characteristics (wider spreads, lower arrival rates, higher execution cost).

## Parameters under study

We study 20 parameters: 5 global and 5 per pair. All are varied uniformly around the paper's Table 1 nominal values.

### Global parameters

| # | Parameter | Symbol | Nominal | Range | Distribution |
|---|-----------|--------|---------|-------|-------------|
| 1 | EUR volatility | sigma_EUR | 80 bps | [56, 104] | Uniform (+-30%) |
| 2 | GBP volatility | sigma_GBP | 70 bps | [49, 91] | Uniform (+-30%) |
| 3 | EUR/GBP correlation | rho | 0.6 | [0.3, 0.9] | Uniform (+-50%) |
| 4 | Risk aversion | gamma | 20 | [10, 40] | Uniform (+-50%) |
| 5 | Arrival rate scale | lambda_scale | 1.0 | [0.5, 1.5] | Uniform (+-50%) |

### EUR/USD pair parameters

| # | Parameter | Symbol | Nominal | Range | Distribution |
|---|-----------|--------|---------|-------|-------------|
| 6 | Logistic shift, tier 1 | alpha_1^{EU} | -1.9 | [-2.5, -1.3] | Uniform (+-30%) |
| 7 | Logistic shift, tier 2 | alpha_2^{EU} | -0.3 | [-0.9, 0.3] | Uniform (+-100%) |
| 8 | Logistic slope, tier 1 | beta_1^{EU} | 11.0 (1/bps) | [7.7, 14.3] | Uniform (+-30%) |
| 9 | Logistic slope, tier 2 | beta_2^{EU} | 3.5 (1/bps) | [2.45, 4.55] | Uniform (+-30%) |
| 10 | Execution cost | eta^{EU} | 1e-5 bps | [0.5e-5, 2e-5] | Uniform (+-50%) |

### GBP/USD pair parameters

| # | Parameter | Symbol | Nominal | Range | Distribution |
|---|-----------|--------|---------|-------|-------------|
| 11 | Logistic shift, tier 1 | alpha_1^{GU} | -1.4 | [-2.0, -0.8] | Uniform (+-~43%) |
| 12 | Logistic shift, tier 2 | alpha_2^{GU} | 0.0 | [-0.6, 0.6] | Uniform (+-0.6 abs) |
| 13 | Logistic slope, tier 1 | beta_1^{GU} | 5.5 (1/bps) | [3.85, 7.15] | Uniform (+-30%) |
| 14 | Logistic slope, tier 2 | beta_2^{GU} | 2.0 (1/bps) | [1.4, 2.6] | Uniform (+-30%) |
| 15 | Execution cost | eta^{GU} | 1.5e-5 bps | [0.75e-5, 3e-5] | Uniform (+-50%) |

### EUR/GBP pair parameters

| # | Parameter | Symbol | Nominal | Range | Distribution |
|---|-----------|--------|---------|-------|-------------|
| 16 | Logistic shift, tier 1 | alpha_1^{EG} | -0.5 | [-1.1, 0.1] | Uniform (+-~120%) |
| 17 | Logistic shift, tier 2 | alpha_2^{EG} | 0.5 | [-0.1, 1.1] | Uniform (+-~120%) |
| 18 | Logistic slope, tier 1 | beta_1^{EG} | 3.5 (1/bps) | [2.45, 4.55] | Uniform (+-30%) |
| 19 | Logistic slope, tier 2 | beta_2^{EG} | 2.5 (1/bps) | [1.75, 3.25] | Uniform (+-30%) |
| 20 | Execution cost | eta^{EG} | 3e-5 bps | [1.5e-5, 6e-5] | Uniform (+-50%) |

### Rationale for parameter selection

Same categories as the 2-currency analysis (see setup.md), extended to cover all 3 pairs:

- **sigma_EUR, sigma_GBP**: Estimated from market data, genuine daily variation.
- **rho**: Estimated from return correlations, varies across market regimes. This is the parameter most specific to the multi-currency extension; it controls the strength of cross-pair coupling through the covariance matrix.
- **gamma**: Design choice by the bank. Same argument as 2-currency.
- **lambda_scale**: Shared multiplier on all arrival rates across all 3 pairs. Same grouping argument: market activity rises or falls broadly, not pair-by-pair in opposite directions. Applied to all three pairs because EUR/USD and GBP/USD are both G10 majors, and EUR/GBP flow tracks overall FX activity.
- **alpha, beta per pair per tier**: Each pair has its own client population with separately calibrated demand curves. Tier 1 (aggressive) and tier 2 (passive) represent fundamentally different client segments. All are fitted from fill/reject data and carry estimation uncertainty. We include both tiers for all pairs, making no assumptions about which tiers or pairs are negligible.
- **eta per pair**: Market impact / execution cost, varies by pair liquidity (EUR/USD: 1e-5, GBP/USD: 1.5e-5, EUR/GBP: 3e-5 in bps). Hard to estimate precisely in practice.

### Parameters not varied (and why)

Same exclusions as 2-currency, for the same structural reasons:

- **psi (per pair)**: Hedge dead zone. Set by the execution desk, well known, enters only as a simple threshold.
- **kappa (terminal penalty)**: Set to 0, boundary condition.
- **T_days**: Fixed design parameter (0.05 days).
- **mu (drift)**: Set to 0 in the paper's example.
- **k (exogenous trade rate)**: Negligibly small values (5e-3 to 8e-3 bps).

### Grouping choices

- **Arrival rates** use a single shared lambda_scale across all 3 pairs (same rationale as 2-currency; see setup.md). If pair-specific flow effects are of interest, independent scale factors could be added in a follow-up analysis.
- **Demand curve and execution cost parameters** are varied independently per pair. Each pair has a distinct client base and liquidity profile. No parameters are shared across pairs.

### Why 20 parameters?

The 2-currency analysis studied 8 parameters for 1 pair. The natural extension to 3 pairs with cross-currency effects gives 5 global + 5 per pair x 3 = 20 parameters. This is a self-contained analysis: all model inputs with genuine uncertainty are included, and the Sobol indices themselves determine which are important and which are not. No parameters are excluded based on 2-currency results.

Computational cost: N(k+2) = 10,000 x 22 = 220,000 ODE evaluations. With d=3, each Riccati solve involves 3x3 matrix operations (vs 2x2 for 2 currencies), estimated at ~12-15 ms per evaluation. Total: ~45-55 minutes. Feasible without any approximation shortcuts.

## Quantities of interest

### Design principles

We choose QoIs to satisfy four requirements:

1. **Direct comparability**: Some QoIs should be identical in definition to the 2-currency QoIs, so we can compare Sobol rankings side by side and see whether adding a third currency changes the sensitivity structure.
2. **Cross-pair probing**: Some QoIs should measure effects that are structurally absent in 2 currencies (cross-inventory influence, cross-hedging).
3. **Practical relevance**: Each QoI should correspond to a quantity a trading desk actually monitors or manages.
4. **Numerical cleanliness**: All QoIs are deterministic given the ODE solution (no Monte Carlo noise) and have well-defined scaling.

### Inventory evaluation states

All QoIs use single-currency perturbations from the flat state. This isolates each mechanism cleanly.

| Label | y = [y_USD, y_EUR, y_GBP] | Purpose |
|-------|---------------------------|---------|
| Flat | [0, 0, 0] | Baseline for all differential QoIs, tier spread, revenue |
| Long EUR | [0, 10, 0] | Own-inventory skew, hedge rates, net revenue |
| Long GBP | [0, 0, 10] | Cross-inventory skew |

Using 10 M$ as the inventory perturbation (same as 2-currency) keeps the magnitude comparable across analyses.

### Set A — Quoting policy (what the client sees)

**QoI 1: Tier spread differential** [bps]

delta_star(y=0, tier 1, EUR/USD) - delta_star(y=0, tier 2, EUR/USD), z = 1 M$.

How much wider the model quotes for passive clients vs. aggressive clients on EUR/USD at flat inventory. Identical definition to 2-currency QoI 1.

Mathematical note: at y=0 with B=0, the momentum input to the markup formula is p = z * d^T A_0 d, which depends on the 2x2 block of A involving USD and EUR only. But the entries A[EUR,EUR] and A[USD,EUR] are shaped by the full Riccati dynamics, including sigma_GBP and rho through the covariance coupling term (-0.5 * gamma * Sigma) and the nonlinear Riccati term (2*A*M_big*A). So cross-currency parameters could in principle affect even the flat-inventory tier differential, though this is expected to be a second-order effect.

**QoI 2: Own-inventory skew** [bps]

delta_star(y_EUR=10, tier 1, EUR/USD) - delta_star(y=0, tier 1, EUR/USD), z = 1 M$.

How much the EUR/USD markup shifts when the market maker holds 10 M$ EUR. Identical definition to 2-currency QoI 2.

The change in p is: Delta_p = 20 * (A_0[EUR,EUR] - A_0[EUR,USD]), which depends on the Riccati solution. In 3 currencies, A is affected by the cross-currency coupling, so rho and sigma_GBP may appear in the Sobol ranking (they could not in 2 currencies, where the Riccati system was decoupled).

**QoI 3: Cross-inventory skew** [bps]

delta_star(y_GBP=10, tier 1, EUR/USD) - delta_star(y=0, tier 1, EUR/USD), z = 1 M$.

How much GBP inventory shifts the EUR/USD markup. This QoI is structurally zero in 2 currencies (GBP doesn't exist). It measures the cross-currency coupling through the value function.

The change in p is: Delta_p = 20 * (A_0[GBP,EUR] - A_0[GBP,USD]). This depends on the off-diagonal entries of A that couple GBP to EUR/USD. These entries arise from (a) the correlation rho entering through -0.5*gamma*Sigma, and (b) the nonlinear Riccati coupling through M_big (which aggregates Hamiltonian coefficients from all 3 pairs, including EUR/GBP).

Practical meaning: "If I'm holding 10 M$ GBP, how does that change my EUR/USD pricing?" This is exactly the multi-desk coordination question that motivates the multi-currency model.

### Set B — Hedging and economics (what the desk manages)

**QoI 4: Own-pair hedge rate** [M$/day]

xi_star(EUR/USD) at y = [0, 10, 0].

How aggressively the model hedges EUR exposure through the EUR/USD interbank market. Identical definition to 2-currency QoI 3.

In 3 currencies, the model has the option to hedge EUR through EUR/GBP as well. If cross-hedging is significant, the own-pair hedge rate might be smaller (some load shifted to the cross pair). Whether the Sobol ranking changes compared to 2 currencies is an empirical question.

**QoI 5: Cross-hedge momentum** [decimal]

p_{EUR/GBP} = -(2*A_0*y + B_0) . d_{EUR/GBP} at y = [0, 10, 0].

The value function gradient projected onto the EUR/GBP hedging direction. This is the "desire to cross-hedge" — the momentum that would drive cross-pair hedging if execution costs were zero. It measures the strength of the cross-hedging incentive from the value function, independent of whether the dead zone and execution cost make it economical.

We use the momentum rather than the clipped hedge rate xi_star(EUR/GBP) because the EUR/GBP dead zone (psi_EG = 0.25 bps) is wide enough that xi_star is zero for ~99.7% of the parameter samples. The clipped rate has essentially no variance to decompose, making Sobol indices meaningless. The momentum is continuous and well-behaved.

The momentum depends on A_0[EUR,EUR] and A_0[EUR,GBP], which encode the cross-pair risk coupling. With y = [0, 10, 0] and d_{EUR/GBP} = [0, +1, -1]:

p = -20 * (A_0[EUR,EUR] - A_0[GBP,EUR])

Practical meaning: "How strongly does the value function push toward hedging EUR through the cross pair?" Its Sobol indices show what drives the cross-hedging incentive. The fact that the actual hedge rate is zero (blocked by the dead zone) is itself a finding about the EUR/GBP execution parameters.

**QoI 6: Net revenue** [$/day]

R(y=[0,10,0]) - L_total(y=[0,10,0]).

Revenue across all 3 pairs minus total hedging cost across all active hedge instruments, at the long-EUR inventory state. Revenue is: R = sum over pairs, directions, tiers, sizes of lambda_j * f(delta_star) * delta_star * z_j. Hedging cost is: L_total = sum over pairs of psi*|xi_star| + eta*xi_star^2.

This is the aggregated "bottom line under inventory stress." It integrates quoting distortions, reduced fill rates, and hedging costs into the single number a desk head tracks: "given that I'm carrying EUR inventory, am I making money after hedging costs?" The hedging cost now spans multiple instruments, and the revenue now includes contributions from all 3 pairs.

### Why these QoIs and not others

**Included:**
- QoIs 1, 2, 4 are identical to 2-currency QoIs 1, 2, 3. Side-by-side comparison shows whether the sensitivity structure changes when cross-currency effects are present.
- QoIs 3, 5 probe the two genuinely new mechanisms (cross-inventory influence on quoting, cross-hedging incentive). Neither can exist in 2 currencies.
- QoI 6 is the aggregated economic outcome, the most downstream KPI.

**Considered and excluded:**
- **R(y=0) (flat revenue)**: Measures earning power at zero inventory. Used in 2-currency Set B. Dropped here because at y=0 the cross-pair effects are weak (pairs are largely decoupled), and the slot is better used for the cross QoIs.
- **L_total separately**: Partially redundant with QoIs 4 and 5 (hedge rates already contain the information; L just adds eta and psi scaling). Captured within QoI 6.
- **GBP/USD inventory skew** (delta_star_{GBP/USD}(y_GBP=10) - delta_star_{GBP/USD}(y=0)): Same structure as QoI 2 but for a different pair. Interesting for completeness but doesn't test a new mechanism. Could be an appendix item.
- **A_0[EUR,GBP] (Riccati off-diagonal)**: Directly measures cross-pair coupling strength. Mathematically clean, but not a decision variable the desk manages. The cross-inventory skew (QoI 3) is its observable downstream consequence.
- **Hedge rate ratio** (xi_EURGBP / xi_EURUSD): Would show the cross-hedging fraction, but ratios are numerically unstable when the denominator is near zero.
- **Multi-currency inventory states** like y = [0, 10, 10] (long both) or y = [0, 10, -10] (spread position): Would probe the interaction of own and cross inventory. Interesting follow-up, but single-currency perturbations are cleaner for isolating mechanisms.

### All QoIs are deterministic

As in the 2-currency analysis, all 6 QoIs are deterministic functions of the 20 parameters (no Monte Carlo noise). The mapping is: parameter vector -> Riccati ODE solve -> A_0, B_0 -> policy extraction -> QoI values. This keeps the Sobol analysis clean.

## Methods

Same methodology as the 2-currency analysis: Saltelli Monte Carlo estimator for Sobol indices, forward UQ via KDE on the parameter samples. See setup.md for details.

The only change is the increased parameter dimension (k=20 vs k=8), which increases the Saltelli sample count from N(k+2) = 10N to 22N. With N=10,000, this gives 220,000 evaluations.

## What we expect and what would be interesting

### Expected results

- **QoI 1 (tier spread differential)**: Still dominated by EUR/USD demand parameters (alpha_1^EU, beta_1^EU, alpha_2^EU, beta_2^EU). The cross-pair coupling through A_0 is a second-order effect at y=0. If sigma_GBP or rho appear with nonzero S_Ti, that is itself a finding worth reporting.

- **QoI 2 (own-inventory skew)**: Driven by sigma_EUR, gamma, lambda_scale as in 2 currencies. The question is whether rho and sigma_GBP now appear with nonzero total-effect indices, since they affect A_0[EUR,EUR] through the Riccati coupling.

- **QoI 3 (cross-inventory skew)**: Expected drivers: rho (the correlation controls the off-diagonal of Sigma, which is the primary source of A_{EUR,GBP}), sigma_EUR and sigma_GBP (which scale the off-diagonal covariance), and gamma (which scales the entire Sigma contribution to the Riccati ODE). Demand parameters of EUR/GBP may appear because they contribute to M, M_tilde, P, which feed into the Riccati nonlinear coupling. Whether they are significant is unknown.

- **QoI 4 (own-pair hedge rate)**: Similar to 2-currency results (sigma, gamma, lambda, eta^EU), possibly with reduced magnitude if some hedging load shifts to the cross pair.

- **QoI 5 (cross-hedge momentum)**: Expected drivers: sigma_EUR, rho, gamma, lambda_scale, sigma_GBP. Since this is the raw momentum (not clipped by the dead zone or divided by eta), it should be continuous and well-behaved for Sobol analysis.

- **QoI 6 (net revenue)**: Integrates many effects. Lambda_scale likely dominates (it scales all revenue). The decomposition of sensitivity across the remaining parameters is the unknown.

### Key open questions the analysis will answer

1. **Does rho appear as a significant driver for any QoI?** If so, how does it rank against the parameters that dominate in 2 currencies?

2. **Do demand parameters of one pair affect the strategy for another pair?** The M_big matrix in the Riccati ODE aggregates Hamiltonian coefficients from all pairs. EUR/GBP's alpha and beta contribute to M, which affects A, which affects EUR/USD quoting. Is this coupling negligible or significant?

3. **Is cross-hedging economically significant?** QoI 5 measures the cross-hedging incentive from the value function. If the momentum is small, cross-hedging is a minor correction. If it is large but blocked by the dead zone, that tells us the value function wants to cross-hedge but execution costs prevent it.

4. **Are the 2-currency sensitivity rankings robust?** Comparing QoIs 1, 2, 4 with their 2-currency counterparts shows whether the main conclusions ("spreads depend on demand curves, inventory management depends on market conditions, hedging depends on everything") survive the addition of cross-currency effects.

5. **Which of the 20 parameters can be safely fixed?** The Sobol analysis will identify parameters with S_Ti near zero across all QoIs. These are candidates for exclusion in any follow-up analysis or in practical model deployment.

## Connection to the thesis

The sensitivity chapter presents the 3-currency analysis as the main (and only) sensitivity study. The 2-currency analysis was development work that established the methodology; it remains in the repo as `sensitivity_analysis_2ccy.ipynb` but does not appear in the report.

The chapter should cover:
1. Setup: parameters, QoIs, methodology (Sobol + forward UQ)
2. Results: Sobol indices per QoI, with interpretation
3. Forward UQ: output distributions, robustness assessment
4. Discussion: what the results mean for model calibration and deployment

The notebook is `sensitivity_analysis.ipynb`. Saved data is in `data/sensitivity_analysis/`. Figures are saved to `report/figures/fig_sa_sobol.pdf` and `report/figures/fig_sa_forward_uq.pdf`.
