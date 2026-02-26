# Figure 1 reproduction — investigation of discrepancies

Comparison of our ODE approximation output with the paper's Figure 1
(Barzykin, Bergault, Gueant 2023, p.7).

## Setup

Figure 1 plots top-of-book quotes (tier 1, z = 1 M$) for GBPUSD, EURUSD,
and EURGBP as functions of GBP inventory (other inventories = 0).
The paper's caption says the curves are `delta^{1,X,Y}` and `-delta^{1,Y,X}`
for each pair XY. Risk aversion gamma = 20.

## What we verified

### ODE convergence

A(0) is fully stationary — identical (to machine precision) for
T = 0.05, 0.1, 0.5, 1.0, 5.0 days. The Euler scheme also converges:
results are the same at n_steps = 2000 and 50000. Not the issue.

### Lambda interpretation (per tier vs total)

The paper says "intensity amplitudes lambda(z) are taken to be the same for
each tier." We tested halving lambda (treating the table values as totals
shared across 2 tiers). Result: negligible change (~1%). Not the issue.

### Minimum spreads at zero inventory

The zero-inventory spread is bounded below by the (alpha, beta) parameters
alone (via the supremum delta*(0)):

| Pair    | alpha | beta (1/decimal) | min delta*(0) (bps) |
|---------|-------|-------------------|---------------------|
| EURUSD  | -1.9  | 110 000           | 0.177               |
| GBPUSD  | -1.4  |  55 000           | 0.313               |
| EURGBP  | -0.5  |  35 000           | 0.401               |

These are hard floors: no A matrix can make the spread smaller. Our computed
values at inv = 0 (which add a small p from z * d^T A d) are:

| Pair    | Our delta*(0) (bps) |
|---------|---------------------|
| EURUSD  | 0.182               |
| GBPUSD  | 0.321               |
| EURGBP  | 0.411               |

These are consistent with the floors above.

### 5-currency vs 3-currency restriction

Restricting to {USD, EUR, GBP} changes the quotes by < 1%. The other
currencies (JPY, CHF) have minimal influence on the EUR-GBP-USD subsystem
at these parameters.

## Discrepancies found

### 1. Plotting convention (sign/slope mismatch)

In the paper's figure, ALL curves slope **downward** (left to right).
Our code (Convention 1) produces curves where:
- GBPUSD slopes **upward** (positive slope with GBP inventory)
- EURGBP slopes **downward**

This is only consistent with the paper if the figure actually plots

    -delta^{X,Y}   and   +delta^{Y,X}

i.e., bid and ask **prices** relative to mid, NOT the markup convention
stated in the caption (`delta^{X,Y}` and `-delta^{Y,X}`).

In this "price" convention (Convention 2):

| inv (M$) | GBPUSD -d^{G,U} | GBPUSD +d^{U,G} | EURGBP -d^{E,G} | EURGBP +d^{G,E} |
|----------|------------------|------------------|------------------|------------------|
|     -100 |          +0.192  |          +2.696  |          -2.050  |          -0.264  |
|      -50 |          +0.056  |          +1.445  |          -1.182  |          -0.038  |
|        0 |          -0.321  |          +0.321  |          -0.411  |          +0.411  |
|      +50 |          -1.445  |          -0.056  |          +0.038  |          +1.182  |
|     +100 |          -2.696  |          -0.192  |          +0.264  |          +2.050  |

All columns now have negative slope, matching the paper's visual. The paper's
caption likely has the convention swapped relative to the actual figure.

### 2. Magnitude discrepancy for direct pairs

Comparing Convention 2 extremes with the paper's figure:

| Pair    | Our extreme (bps) | Paper's extreme (bps, approx) | Ratio |
|---------|-------------------|-------------------------------|-------|
| EURGBP  | +/- 2.05          | +/- 2.0                       | ~1.0  |
| GBPUSD  | +/- 2.70          | +/- 0.9                       | ~3.0  |
| EURUSD  | +/- 0.85          | +/- 0.4                       | ~2.1  |

The cross pair EURGBP matches well. Direct pairs (GBPUSD, EURUSD) have
inventory sensitivity that is 2-3x larger than the paper's figure.

The inventory sensitivity is controlled by the A matrix diagonal:
- dp/d(inv_GBP) for GBPUSD ~ 2 * A[GBP,GBP] = 2.50e-6
- dp/d(inv_GBP) for EURGBP ~ 2 * (A[EUR,GBP] - A[GBP,GBP]) = -1.75e-6

Since A[USD, :] = A[:, USD] = 0 (reference currency has zero volatility),
the GBPUSD sensitivity is entirely determined by A[GBP,GBP]. A smaller
A[GBP,GBP] would reduce the GBPUSD sensitivity without affecting EURGBP
as strongly (since EURGBP depends on A[EUR,GBP] - A[GBP,GBP]).

## Possible explanations for the magnitude discrepancy

1. **Different parameters**: the paper may have generated the figure with
   slightly different parameters (or an earlier calibration) than what ended
   up in the published Table 1.

2. **Different ODE formulation detail**: a subtle factor (e.g., in how
   M_big is symmetrized, or in the Hadamard product term) could change the
   effective damping and hence A[GBP,GBP].

3. **Full PDE vs ODE approximation**: the paper validated the ODE against
   a full PDE grid solver for d = 2. The figure may have been generated
   from the PDE solution rather than the ODE approximation, which could
   produce different inventory sensitivities.

4. **Convention in M matrix**: if the paper's M matrix counts each
   directed pair (i,j) only once (upper triangle), while our code fills
   both (i,j) and (j,i), the effective M_big would differ by a factor.
   This would change the balance in the Riccati ODE and hence A.

None of these have been conclusively confirmed or ruled out. The cross pair
matching well while direct pairs are off suggests the issue is specifically
in how the reference currency (USD) row interacts with the rest of the
A matrix dynamics.
