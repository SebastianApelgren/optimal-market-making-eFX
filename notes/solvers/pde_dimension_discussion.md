# Why the PDE comparison stays at d=2

Discussion from 2026-04-14 about whether to extend the PDE vs ODE comparison
from 2 currencies (1 pair) to 3 currencies (3 pairs), matching the sensitivity
analysis which uses 3 currencies.

## Decision: keep PDE at d=2

## Computational argument

- d=2: 301² = 90,601 grid points. 9 PDE solves take ~90 minutes.
- d=3: 301³ = 27.3M points — completely impractical.
- Coarsening to dy=5 gives 31³–41³ ≈ 30k–69k points, but drops the z=1 M$
  trade size (1 is not a multiple of 5). This changes the model, making the
  comparison not apples-to-apples with the ODE.
- Even at comparable grid sizes, 3 pairs means ~3× more nonzeros per row in
  the sparse matrix (more quoting, hedging, and diffusion cross-terms).
  Estimated ~20–30 min per solve, 3–4 hours total.

## Theoretical argument (the important one)

The ODE approximation error has a single source: the quadratic Taylor expansion
of the per-pair Hamiltonian H(p) around p=0, replacing the exact H with
Ĥ = α₀ + α₁p + ½α₂p². This approximation is made independently for each
pair, tier, and trade size.

The cross-pair interactions (off-diagonal elements of A, coupling between e.g.
EUR and GBP inventories) are captured exactly by the Riccati matrix algebra in
Eq. 5. No additional approximation is introduced when going from 1 pair to 3.

Therefore the 2-currency PDE comparison already tests the dominant (and only)
source of ODE error. Extending to d=3 tests the same per-pair phenomenon on a
coarser grid, reducing the fidelity of the PDE baseline without probing a
qualitatively different error source.

## Suggested thesis text (adapt as needed)

The ODE approximation error arises from the quadratic expansion of the per-pair
Hamiltonian H(p). Cross-pair interactions enter through the Riccati matrix
algebra, which is exact. The 2-currency PDE comparison therefore characterises
the dominant approximation error. Extending to d=3 would test the same per-pair
phenomenon on a necessarily coarser grid, reducing the fidelity of the PDE
baseline without testing a qualitatively different source of error.

## Existing d=2 results (from pde_ode_comparison.ipynb)

- Markup δ*(y=0): ODE vs PDE agree to ~0.2% across all parameter configs
- Inventory skew: 3–8% relative difference
- Hedge rate: larger deviations at extreme parameters (14–76% rel diff),
  driven by the quadratic-H approximation breaking down at large |p|
- These numbers characterise the quadratic-H error fully.

## Bridging d=2 validation to d=3 SA credibility

Keeping PDE at d=2 means the thesis presents a 3-currency SA validated by a
2-currency PDE. Three arguments bridge this gap (details in
`pde_ode_comparison_design.md`):

1. **Theoretical**: the ODE error is per-pair (quadratic-H), cross-pair
   coupling is exact in the Riccati algebra. Going from d=2 to d=3 adds
   exact cross-pair terms but no new approximation error.

2. **Riccati matrix inspection**: solve the 3-ccy ODE, extract A(0), show
   that |A[EUR,GBP]| is small relative to A[EUR,EUR]. Cross-pair coupling
   is a small perturbation on per-pair dynamics already validated.

3. **2-ccy vs 3-ccy ODE comparison**: compute EUR/USD QoIs from both 2-ccy
   and 3-ccy ODE. If adding GBP barely changes the EUR/USD results, the
   d=2 PDE validation carries over.

Arguments 2 and 3 are cheap (milliseconds, ODE only) and produce concrete
numbers instead of a hand-wave.

## Supervisor context

Anders confirmed (2026-04-14 discussion) that the PDE comparison does not need
to be super thorough (fine grid, many time steps) since it is validation, not
the main contribution. This reinforces keeping d=2 with a clean, convincing
comparison rather than a noisy d=3 one.
