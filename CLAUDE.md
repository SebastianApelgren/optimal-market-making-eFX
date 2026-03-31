# CLAUDE.md

## Project

Master thesis: numerical methods for optimal multi-currency market making in eFX.
Based on Barzykin, Bergault, Guéant (2023), arXiv:2207.04100v4 (in `references/`).

Thesis structure:
1. Model and ODE approximation (complete)
2. Sensitivity analysis / UQ on ODE parameters (active — next major work)
3. PDE solver, PDE vs ODE comparison (complete)
4. Connect ODE sensitivity results to the full problem via PDE validation

## Environment

**Always use the `.venv` virtual environment. Never install packages globally.**

```bash
source .venv/bin/activate
pip install -r requirements.txt   # numpy, matplotlib, scipy
```

When running Python scripts or notebooks, use `.venv/bin/python` or ensure the venv is activated. Dependencies: numpy, matplotlib, scipy. No package manager beyond pip.

## Architecture

```
src/                        # Shared library — all notebooks import from here
  model.py                  # ModelParams, TierParams, PairParams, build_paper_example_params()
  hamiltonian.py            # Logistic intensity f(δ), optimal δ*, Hamiltonian H(p) via Lambert W
  riccati.py                # Σ, M, M̃, P matrix builders; Riccati ODE solver (Eq. 5)
  policy.py                 # Quote formula, hedge rate, MMResult, run_multicurrency_mm()
  pde.py                    # PDE spatial operators, explicit/semi-implicit/implicit solvers
  simulation.py             # Monte Carlo tau-leaping inventory simulator
  plotting.py               # Figures 1 & 2 reproduction helpers

ode_approximation_reproduction.ipynb   # Complete — reproduces [1] Table 1, Figures 1-2, MC sim
ode_solver_verification.ipynb          # Complete — manufactured solution tests for Riccati ODE
hjb_pde_solver.ipynb                   # Complete — explicit Euler PDE with η×1000 scaling
eta_continuation_solver.ipynb          # Active — implicit solver with original η, PDE vs ODE

notes/
  figure1_investigation.md             # ODE Figure 1 reproduction discrepancy
  spatial-operators/                   # Derivations for each PDE term (6 notes)
  solvers/
    summary.md                         # Overview of all solvers and what worked
    grid_and_domain.md                 # Shared: PDE, grid, boundary, validation
    explicit_euler.md                  # Explicit Euler + CFL + stiffness analysis
    semi-implicit/                     # IMEX scheme (tried, dropped)
    implicit/                          # Policy iteration + monotone scheme (working)
```

## Units and Conventions

Getting units right is critical. The paper mixes basis points and decimals.

| Quantity | Internal unit | Conversion |
|----------|--------------|------------|
| Inventory `Y` | M$ (millions USD) | — |
| Markup `δ` | decimal | × 10⁴ → bps |
| Volatility `σ` | decimal / √day | paper gives bps → multiply by `BP = 1e-4` |
| Risk aversion `γ` | 1/M$ | — |
| Hedge rate `ξ` | M$/day | — |
| Execution cost `L(ξ)` | decimal | `ψ|ξ| + ηξ²` |
| Arrival rate `λ` | per day | ÷ `DAY_SECONDS` → per second (simulation only) |
| `α` (logistic) | dimensionless | — |
| `β` (logistic) | 1/decimal | paper gives 1/bps → multiply by `1e4` |

Constants: `BP = 1e-4`, `DAY_SECONDS = 86400.0` (in `src/model.py`).

## Mathematical Correspondence

The code maps to the paper as follows:

- **Value function ansatz**: `θ ≈ -Y^T A(t) Y - Y^T B(t) - C(t)` (Eq. 4)
- **Riccati ODE**: `solve_AB_euler` solves Eq. 5 backward from `T` with `A(T) = κ`, `B(T) = 0`
- **Quadratic approximation**: `Ĥ = α₀ + α₁p + ½α₂p²` where `α₂ = H''(0)` (the ½ is explicit in the formula)
- **Quote formula** (p.6): `p = ((2y + z·d) · A + B) · d`, then `δ* = argmax_δ f(δ)·(δ - p)`
- **Hedge rate** (p.6): `ξ* = H'_L(-(2Ay + B) · d + market_impact)`
- **PDE solver**: fully implicit Euler with policy iteration (Howard's algorithm) and monotone (M-matrix) differencing for the hedging advection term

### Known paper typo (p.6)

The approximate hedging formula in the paper writes `(AY + B)` but the correct gradient of the ansatz gives `(2AY + B)`. The quote formula on the same page correctly has the factor 2. The code uses the correct `2AY + B` everywhere.

## PDE Solver Summary

Three solvers were implemented and tested:

1. **Explicit Euler** (`solve_hjb_explicit`): Simple, validates spatial operators. Requires η×1000 scaling — the hedging Hamiltonian's gain 1/(4η) ≈ 2.5×10⁸ makes it unconditionally unstable with original η.

2. **Semi-implicit IMEX** (`solve_hjb_semi_implicit`): Treats diffusion implicitly. Provides no improvement — diffusion was never the binding constraint. Dropped.

3. **Fully implicit + policy iteration** (`solve_hjb_implicit`): The working solver. Uses monotone (M-matrix) differencing for the hedging advection term — this is what the paper means by "monotone implicit Euler." PI converges in 2-7 iterations with tolerance 1e-4. Works with original η.

Also available: `solve_hjb_implicit_continuation` for η-continuation (solving with decreasing η scales). Tested but found unnecessary — the monotone scheme alone is sufficient.

## Code Conventions

- Frozen dataclasses for immutable parameter containers (`TierParams`, `PairParams`)
- Currency pairs stored as sorted tuples via `canon_pair(a, b) → tuple(sorted((a, b)))`
- Currency indexing: `mp.currencies` list defines the ordering; `d = len(mp.currencies)`
- Direction vectors: `dvec[i] = +1, dvec[j] = -1` for a trade paying currency `i`, selling `j`
- `MMResult` wraps a solved model — use `.markup()` and `.hedge_rate()` for policy queries
- `build_paper_example_params()` returns the full 5-currency setup from Table 1
- `restrict_currencies(mp, keep)` slices down to a subset (must include `ref_ccy`)
- Shared code goes in `src/`, not in notebook cells. Notebooks import via `from src import ...`
- `src/__init__.py` re-exports all public names — update it when adding new functions
- No test suite yet — verify correctness by running the reproduction notebook and comparing against paper values
- `B0 = 0` in the paper's example because parameters are direction-symmetric (`M̄` symmetric, `V = 0`, `μ = 0`)

## Git

- Never add `Co-Authored-By` lines to commits.

## Key Functions (Entry Points)

```python
# ODE approximation
from src import build_paper_example_params, restrict_currencies, run_multicurrency_mm

mp = build_paper_example_params()            # Full 5-currency model
mp2 = restrict_currencies(mp, ["USD","EUR"]) # 2-currency subset
res = run_multicurrency_mm(mp2)              # Solve Riccati ODE → MMResult

res.markup(tier_idx=0, ccy_pay="EUR", ccy_sell="USD", z_musd=1.0, y=np.zeros(2))
res.hedge_rate(ccy_buy="EUR", ccy_sell="USD", y=np.array([5.0, -5.0]))

# PDE solver (implicit, original η)
from src.pde import solve_hjb_implicit, build_pde_spec

y_grids = [np.arange(-150, 151, 1.0) for _ in range(2)]
result = solve_hjb_implicit(y_grids, mp2, n_steps=200, pi_tol=1e-4, pi_max_iter=50)
theta_0 = result['theta_0']  # θ(0, y) on the grid
```
