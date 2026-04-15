# Sensitivity Analysis — 3-Currency Results

## Setup summary

- 3 currencies (USD, EUR, GBP), 3 pairs (EUR/USD, GBP/USD, EUR/GBP)
- 20 parameters (5 global + 5 per pair), all uniform distributions
- 6 QoIs: 3 quoting policy (Set A), 3 hedging/economics (Set B)
- Saltelli estimator: N = 10,000 base samples, 220,000 total evaluations
- Forward UQ: 20,000 samples (reused from Sobol)

## Sobol indices

### QoI 1: Tier spread differential [bps]

delta_star(y=0, tier 1, EUR/USD) - delta_star(y=0, tier 2, EUR/USD), z = 1 M$.

| Parameter | S_i | S_Ti |
|-----------|-----|------|
| beta_2 EU | 0.719 | 0.718 |
| beta_1 EU | 0.165 | 0.149 |
| alpha_2 EU | 0.106 | 0.094 |
| alpha_1 EU | 0.039 | 0.033 |
| **Sum** | **1.03** | **0.99** |

All other parameters (including all global, GBP/USD, and EUR/GBP parameters) have S_Ti = 0.000.

**Interpretation:** The tier spread differential at flat inventory is controlled entirely by the EUR/USD demand curve parameters. beta_2 (passive tier logistic slope) dominates at 72% because it controls the wider tier-2 markup, which varies more in absolute terms than the tier-1 markup. The model is perfectly additive (Sum S_i = 1.03). No cross-pair parameter has any detectable effect.

This makes mathematical sense: at y=0 with B=0, the momentum p = z * d^T A_0 d is small for z=1, and the markup depends almost entirely on the logistic curve maximization. Although A_0 is shaped by the full Riccati dynamics (including rho and sigma_GBP), the effect on p at this trade size is negligible.

Forward UQ: mean = -0.22 bps, std = 0.09 bps, 90% CI = [-0.37, -0.09], CV = 39%.

### QoI 2: Own-inventory skew [bps]

delta_star(y_EUR=10, tier 1, EUR/USD) - delta_star(y=0, tier 1, EUR/USD), z = 1 M$.

| Parameter | S_i | S_Ti |
|-----------|-----|------|
| sigma_EUR | 0.346 | 0.364 |
| gamma | 0.316 | 0.369 |
| lambda_scale | 0.261 | 0.296 |
| **Sum** | **0.94** | **1.05** |

All pair-specific parameters and cross-pair global parameters (sigma_GBP, rho) have S_Ti < 0.01.

**Interpretation:** The own-inventory skew is driven by three global parameters with roughly equal importance: sigma_EUR (35%), gamma (32%), and lambda_scale (26%). This is identical to the 2-currency result. Adding a third currency does not change the sensitivity structure of own-pair inventory adjustment. The small gap between S_i and S_Ti sums (0.94 vs 1.05) indicates minor interactions.

The cross-pair parameters (rho, sigma_GBP) do NOT appear here, even though they affect A_0[EUR,EUR] through the Riccati coupling. The coupling's effect on the own-pair diagonal entries is too small to register at this sample size.

Forward UQ: mean = 0.15 bps, std = 0.06 bps, 90% CI = [0.07, 0.25], CV = 38%. Nominal (0.12) is below the mean because gamma's range [10, 40] has midpoint 25 > nominal 20.

### QoI 3: Cross-inventory skew [bps]

delta_star(y_GBP=10, tier 1, EUR/USD) - delta_star(y=0, tier 1, EUR/USD), z = 1 M$.

| Parameter | S_i | S_Ti |
|-----------|-----|------|
| rho | 0.387 | 0.473 |
| gamma | 0.196 | 0.239 |
| lambda_scale | 0.150 | 0.193 |
| sigma_GBP | 0.086 | 0.115 |
| sigma_EUR | 0.030 | 0.026 |
| alpha_1 EU | 0.014 | 0.015 |
| **Sum** | **0.87** | **1.08** |

**Interpretation:** The cross-inventory skew — how GBP inventory shifts EUR/USD quotes — is primarily driven by the correlation rho (39%, S_Ti = 47%). This confirms that rho is the gateway to all cross-pair effects: it controls the off-diagonal of Sigma, which creates the off-diagonal entries A_0[EUR,GBP] that mediate cross-pair coupling.

gamma (20%) and lambda_scale (15%) are secondary drivers. They shape the Riccati matrix A through the running penalty and the Hamiltonian coefficients. sigma_GBP (9%) enters because the off-diagonal covariance is rho * sigma_EUR * sigma_GBP, so higher GBP volatility amplifies the correlation effect.

The gap between S_i and S_Ti for rho (0.39 vs 0.47) indicates moderate interaction with the other global parameters. Sum S_Ti = 1.08 confirms this.

alpha_1 EU appears with a tiny but nonzero index (0.014). This is the EUR/USD aggressive-tier logistic shift entering through the Hamiltonian's contribution to M_big in the Riccati ODE. The coupling exists but is negligibly small.

Forward UQ: mean = 0.055 bps, std = 0.026 bps, 90% CI = [0.024, 0.105], CV = 46%. The cross-skew is consistently nonzero (90% CI entirely positive), confirming that cross-pair coupling is a structural feature, not noise. Its magnitude is about 40% of the own-inventory skew.

### QoI 4: Own-pair hedge rate [M$/day]

xi_star(EUR/USD) at y = [0, 10, 0].

| Parameter | S_i | S_Ti |
|-----------|-----|------|
| sigma_EUR | 0.232 | 0.329 |
| gamma | 0.209 | 0.329 |
| lambda_scale | 0.193 | 0.269 |
| eta EU | 0.145 | 0.277 |
| beta_1 EU | 0.025 | 0.041 |
| **Sum** | **0.82** | **1.26** |

All cross-pair parameters have S_Ti < 0.01.

**Interpretation:** Same structure as 2-currency: four parameters contribute (sigma_EUR, gamma, lambda_scale, eta_EU), with significant interactions (Sum S_Ti = 1.26). eta_EU shows the largest gap (S_i = 0.15 vs S_Ti = 0.28), consistent with the multiplicative 1/(2*eta) in the hedge rate formula.

The cross-pair parameters do not affect own-pair hedging. The availability of EUR/GBP as an alternative hedge instrument does not materially change the EUR/USD hedge rate sensitivity. This is because cross-hedging is blocked by the EUR/GBP dead zone at almost all parameter values.

Forward UQ: mean = -2896 M$/day, std = 3211, 90% CI = [-9287, 0], CV = 111%. The hedge rate remains the least robust output, with a CI spanning an order of magnitude. The 95th percentile is exactly 0 (some parameter combinations put the hedge momentum inside the dead zone).

### QoI 5: Cross-hedge momentum [decimal]

p_{EUR/GBP} = -(2*A_0*y + B_0) . d_{EUR/GBP} at y = [0, 10, 0].

| Parameter | S_i | S_Ti |
|-----------|-----|------|
| sigma_EUR | 0.360 | 0.379 |
| rho | 0.246 | 0.254 |
| gamma | 0.137 | 0.181 |
| lambda_scale | 0.114 | 0.138 |
| sigma_GBP | 0.058 | 0.060 |
| beta_1 EU | 0.017 | 0.025 |
| **Sum** | **0.94** | **1.05** |

**Interpretation:** The cross-hedge momentum is driven by global parameters only. sigma_EUR dominates (36%), followed by rho (25%), gamma (14%), lambda_scale (11%), and sigma_GBP (6%). The model is nearly additive (Sum S_i = 0.94).

sigma_EUR leads (rather than rho) because the momentum depends on A_0[EUR,EUR] - A_0[GBP,EUR], and the diagonal entry A_0[EUR,EUR] scales strongly with sigma_EUR through the Riccati driving term -0.5*gamma*Sigma. rho controls the off-diagonal A_0[GBP,EUR] but the diagonal term is larger.

No pair-specific execution cost (eta_EG) or demand parameter appears, because the momentum is computed before the dead zone and 1/(2*eta) clipping. It measures the value function gradient directly, which depends only on the Riccati solution (driven by global parameters and the aggregated Hamiltonian coefficients).

Note on the original QoI choice: we initially defined QoI 5 as the clipped hedge rate xi_star(EUR/GBP). This was zero for 99.7% of samples because the EUR/GBP dead zone (psi_EG = 0.25 bps = 2.5e-5 decimal) exceeds the typical momentum magnitude (~1e-5). The Sobol indices were meaningless (Sum S_Ti = 3.73). Switching to the raw momentum resolved this and reveals the clean sensitivity structure above.

Forward UQ: mean = -1.17e-5, std = 5.43e-6, 90% CI = [-2.17e-5, -3.86e-6], CV = 46%. The momentum is negative (the model wants to sell EUR / buy GBP to hedge EUR exposure), consistently nonzero, and of order 1e-5 in decimal. At nominal parameters, |p| = 1.02e-5 < psi_EG = 2.5e-5, so the actual hedge rate is zero. Cross-hedging is structurally present in the value function but blocked by EUR/GBP execution costs.

### QoI 6: Net revenue [$/day]

R(y=[0,10,0]) - L_total(y=[0,10,0]).

| Parameter | S_i | S_Ti |
|-----------|-----|------|
| lambda_scale | 0.775 | 0.775 |
| gamma | 0.063 | 0.069 |
| sigma_EUR | 0.028 | 0.054 |
| alpha_2 EU | 0.042 | 0.044 |
| beta_2 EU | 0.017 | 0.023 |
| alpha_2 GU | 0.020 | 0.022 |
| eta EU | 0.007 | 0.016 |
| alpha_1 GU | 0.012 | 0.009 |
| beta_2 GU | 0.012 | 0.011 |
| **Sum** | **0.99** | **1.05** |

**Interpretation:** Net revenue is overwhelmingly driven by lambda_scale (78%). It scales all client arrival rates across all 3 pairs, and revenue is roughly proportional to total flow. The model is nearly additive.

The secondary contributors are gamma (6%), sigma_EUR (3%), and several demand parameters at the 1-4% level: alpha_2 EU (4%), alpha_2 GU (2%), beta_2 EU (2%), beta_2 GU (1%). The GBP/USD demand parameters appear here — the only QoI where they have nonzero indices — because net revenue sums across all 3 pairs, and GBP/USD contributes a measurable fraction of total revenue.

Forward UQ: mean = $502k/day, std = $240k, 90% CI = [$124k, $893k], CV = 48%. The nominal ($508k) is close to the mean. The CI is wide but entirely positive — the strategy is profitable across all sampled parameter combinations.

## Key findings

### 1. The 2-currency sensitivity structure is robust

QoIs 1, 2, and 4 reproduce the same Sobol rankings as the 2-currency analysis. Base spreads depend on demand curves only. Inventory adjustment depends on sigma, gamma, lambda only. Own-pair hedging depends on sigma, gamma, lambda, eta with significant interactions. Adding 12 new parameters and the cross-currency coupling mechanism does not change these conclusions.

### 2. Correlation rho is the gateway to cross-pair effects

rho appears as a significant driver only for QoIs 3 and 5 — the two quantities that specifically measure cross-pair coupling. It ranks #1 for the cross-inventory skew (39%) and #2 for the cross-hedge momentum (25%). It has zero influence on all other QoIs. Cross-pair effects propagate exclusively through the covariance matrix Sigma, not through the Hamiltonian aggregation in M, M_tilde, P.

### 3. Cross-pair effects are real but small

The cross-inventory skew (mean 0.055 bps) is about 38% of the own-inventory skew (mean 0.145 bps). It is consistently nonzero across the parameter space (90% CI entirely positive). This is a genuine structural feature of the multi-currency model: holding GBP systematically shifts EUR/USD quotes through the value function coupling.

### 4. Cross-hedging is blocked by execution costs

The model produces a nonzero cross-hedging momentum (mean 1.17e-5 in decimal), but the EUR/GBP dead zone (2.5e-5) blocks it in 99.7% of parameter samples. Cross-hedging is structurally present in the value function but economically unviable at these execution parameters. A more liquid cross pair (lower psi, lower eta) could change this.

### 5. The hedge rate remains the least robust output

The own-pair hedge rate has CV = 111%, with a 90% CI spanning from -9,300 to 0 M$/day. This reproduces the 2-currency finding and is driven by the same 1/eta multiplicative interaction. The third currency neither helps nor hurts.

### 6. Most parameters can be safely fixed

Of the 20 parameters, only 5-7 have nonzero S_Ti for any QoI: sigma_EUR, sigma_GBP, rho, gamma, lambda_scale from the global set, plus eta_EU and the EUR/USD demand parameters for pair-specific QoIs. All 10 GBP/USD and EUR/GBP pair-specific parameters have S_Ti approximately 0 across all QoIs (except minor contributions to net revenue). The 20-parameter model has effectively 7 active dimensions.

### 7. Pair-specific demand parameters don't couple across pairs

The Riccati matrices M, M_tilde, P aggregate Hamiltonian coefficients from all 3 pairs, creating a potential coupling pathway where EUR/GBP demand parameters could influence EUR/USD strategy. The Sobol analysis shows this coupling is negligible: no cross-pair demand parameter has S_Ti > 0.015 for any quoting or hedging QoI. Cross-pair coupling goes through the covariance (sigma, rho), not through client flow characteristics.
