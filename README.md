# Optimal Market Making in eFX

Code accompanying the master thesis

> **Understanding Optimal Multi-Currency FX Market Making: A Sensitivity and Numerical Study**
> Sebastian Apelgren, KTH Royal Institute of Technology, 2026.
> Degree project in Applied and Computational Mathematics, second cycle, 30 credits.

The thesis builds on the multi-currency market making model of Barzykin, Bergault, and Guéant (2023) [[arXiv:2207.04100](https://arxiv.org/abs/2207.04100)]. It implements (i) the ODE approximation of the Hamilton-Jacobi-Bellman (HJB) equation via a quadratic ansatz, (ii) a fully implicit backward Euler grid PDE solver with policy iteration and a monotone upwind discretization of the hedging gradient term, and (iii) a global Sobol sensitivity analysis and forward uncertainty quantification on the three-currency model.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
jupyter notebook
```

Notebooks import shared code from `src/`. The PDE comparison and sensitivity notebooks read precomputed Saltelli, Sobol, and PDE-vs-ODE samples from `data/`; cells that regenerate those samples are clearly labelled and can be re-run if you want to reproduce them from scratch.

## Repository layout

```
optimal-market-making-eFX/
├── README.md
├── LICENSE                                  # MIT
├── CITATION.cff
├── requirements.txt                         # numpy, scipy, matplotlib, tqdm
├── src/                                     # Shared library imported by all notebooks
│   ├── __init__.py
│   ├── model.py                             # ModelParams, parameter builders, Table 1 setup
│   ├── hamiltonian.py                       # Logistic intensity, optimal δ, Hamiltonian via Lambert W
│   ├── riccati.py                           # Σ, M, M̃, P matrices; Riccati ODE solver (Eq. 5)
│   ├── policy.py                            # Quote/hedge formulae, MMResult, run_multicurrency_mm
│   ├── pde.py                               # Spatial operators, explicit/semi-implicit/implicit HJB solvers
│   ├── pde_comparison.py                    # PDE vs ODE comparison utilities
│   ├── sensitivity.py                       # Parameter perturbation and quantity-of-interest helpers
│   ├── sobol.py                             # Saltelli sampling and Sobol index estimation
│   ├── simulation.py                        # Monte Carlo tau-leaping inventory simulator
│   └── plotting.py                          # Figure 1 and 2 reproduction helpers
├── data/                                    # Precomputed Saltelli / Sobol / PDE-vs-ODE samples
│   ├── pde_comparison/
│   └── sensitivity_analysis/
├── ode_approximation_reproduction.ipynb     # Reproduces Table 1, Figures 1-2, MC simulation
├── ode_solver_verification.ipynb            # Manufactured-solution tests for the Riccati ODE
├── hjb_pde_solver.ipynb                     # Explicit Euler PDE with η scaled, validates spatial operators
├── eta_continuation_solver.ipynb            # Implicit PDE solver with the paper's original η
├── pde_ode_comparison.ipynb                 # PDE vs ODE comparison across parameter regimes
└── sensitivity_analysis.ipynb               # Three-currency Sobol indices and forward UQ
```

## Notebooks

| Notebook | Description |
|----------|-------------|
| `ode_approximation_reproduction.ipynb` | Reproduces the ODE-based quadratic approximation of Barzykin et al. Includes Table 1 parameters, Riccati ODE solver, optimal quotes (Figure 1), hedging rates (Figure 2), and a Monte Carlo inventory simulation. |
| `ode_solver_verification.ipynb` | Verifies the Riccati ODE solver using manufactured solutions with known closed-form answers. Confirms first-order convergence in time. |
| `hjb_pde_solver.ipynb` | Explicit Euler PDE solver with η×1000 scaling. Validates each spatial operator (diffusion, quoting, hedging, drift) against the ODE solution. |
| `eta_continuation_solver.ipynb` | Fully implicit PDE solver with the paper's original parameters, monotone (M-matrix) upwinding for the hedging term, and policy iteration. |
| `pde_ode_comparison.ipynb` | Systematic PDE vs ODE comparison across one-at-a-time and Latin-hypercube parameter perturbations in the two-currency case. |
| `sensitivity_analysis.ipynb` | Sobol indices and forward UQ for the three-currency ODE model, the main sensitivity results of the thesis. |

## Mathematical correspondence

| Paper symbol | Code |
|--------------|------|
| Value function ansatz `θ ≈ -Yᵀ A(t) Y - Yᵀ B(t) - C(t)` | `solve_AB_euler` in `riccati.py` |
| Riccati system (Eq. 5) | `solve_AB_euler` (backward Euler from `T` with `A(T) = κ`, `B(T) = 0`) |
| Optimal client markup | `optimal_delta_logistic` in `hamiltonian.py`, exposed as `MMResult.markup(...)` |
| Optimal hedge rate | `optimal_hedge_rate` in `policy.py`, exposed as `MMResult.hedge_rate(...)` |
| Full HJB PDE | `solve_hjb_implicit` in `pde.py` (fully implicit Euler with policy iteration and monotone upwinding) |

### Note on a paper typo (p. 6)

The approximate hedging formula in Barzykin et al. (2023) is written with `(AY + B)`, but the gradient of the quadratic ansatz `θ̂ = -yᵀ A y - yᵀ B - C` is `-(2Ay + B)`. The quote formula on the same page already carries the correct factor of two. The code uses `2Ay + B` everywhere.

## Units and conventions

| Quantity | Internal unit | Conversion |
|----------|--------------|-----------|
| Inventory `Y` | M$ (millions USD) | — |
| Markup `δ` | decimal | × 10⁴ → bps |
| Volatility `σ` | decimal / √day | paper gives bps → multiply by `BP = 1e-4` |
| Risk aversion `γ` | 1/M$ | — |
| Hedge rate `ξ` | M$/day | — |
| Execution cost `L(ξ)` | decimal | `ψ\|ξ\| + ηξ²` |
| Arrival rate `λ` | per day | ÷ `DAY_SECONDS` → per second (simulation only) |
| `α` (logistic) | dimensionless | — |
| `β` (logistic) | 1/decimal | paper gives 1/bps → multiply by `1e4` |

Constants `BP = 1e-4` and `DAY_SECONDS = 86400.0` live in `src/model.py`.

## Citation

If you use this code, please cite the thesis (see `CITATION.cff`):

```
Apelgren, S. (2026). Understanding Optimal Multi-Currency FX Market Making:
A Sensitivity and Numerical Study. Master's thesis, KTH Royal Institute of Technology.
```

## License

MIT. See `LICENSE`.

## References

1. A. Barzykin, P. Bergault, O. Guéant. *Dealing with multi-currency inventory risk in FX cash markets*. arXiv:2207.04100v4, 2023.
