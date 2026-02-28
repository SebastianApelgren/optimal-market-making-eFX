# CLAUDE.md

## Project

Master thesis: numerical methods for optimal multi-currency market making in eFX.
Based on Barzykin, Bergault, Guéant (2023), arXiv:2207.04100v4 (in `references/`).

## Environment

**Always use the `.venv` virtual environment. Never install packages globally.**

```bash
source .venv/bin/activate
pip install -r requirements.txt   # numpy, matplotlib
```

When running Python scripts or notebooks, use `.venv/bin/python` or ensure the venv is activated. Dependencies: numpy, matplotlib. No package manager beyond pip.

## Architecture

```
src/                        # Shared library — all notebooks import from here
  model.py                  # ModelParams, TierParams, PairParams, build_paper_example_params()
  hamiltonian.py            # Logistic intensity f(δ), optimal δ*, Hamiltonian H(p)
  riccati.py                # Σ, M, M̃, P matrix builders; Riccati ODE solver (Eq. 5)
  policy.py                 # Quote formula, hedge rate, MMResult, run_multicurrency_mm()
  simulation.py             # Monte Carlo tau-leaping inventory simulator
  plotting.py               # Figures 1 & 2 reproduction helpers

ode_approximation_reproduction.ipynb   # Complete — reproduces [1] Table 1, Figures 1-2
hjb_pde_solver.ipynb                   # Stub — grid-based FD solver for full HJB PDE
pde_vs_ode_comparison.ipynb            # Stub — PDE vs ODE comparison
sensitivity_analysis.ipynb             # Stub — parameter sweeps, uncertainty quantification
pinns_solver.ipynb                     # Stub — physics-informed neural networks (extension)
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
- **Riccati ODE**: `solve_AB_euler` solves Eq. 5 backward from `T` with `A(T) = γΣ/2`, `B(T) = 0`
- **Quote formula** (p.6): `p = ((2y + z·d) · A + B) · d`, then `δ* = argmax_δ f(δ)·(δ - p)`
- **Hedge rate** (p.6): `ξ* = H'_L(-(2Ay + B) · d + market_impact)`

### Known paper typo (p.6)

The approximate hedging formula in the paper writes `(AY + B)` but the correct gradient of the ansatz gives `(2AY + B)`. The quote formula on the same page correctly has the factor 2. The code uses the correct `2AY + B` everywhere.

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
from src import build_paper_example_params, restrict_currencies, run_multicurrency_mm

mp = build_paper_example_params()            # Full 5-currency model
mp2 = restrict_currencies(mp, ["USD","EUR"]) # 2-currency subset
res = run_multicurrency_mm(mp2)              # Solve Riccati ODE → MMResult

res.markup(tier_idx=0, ccy_pay="EUR", ccy_sell="USD", z_musd=1.0, y=np.zeros(2))
res.hedge_rate(ccy_buy="EUR", ccy_sell="USD", y=np.array([5.0, -5.0]))
```
