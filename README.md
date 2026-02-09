# Optimal Market Making in eFX

Master thesis project on numerical methods for optimal multi-currency market making in electronic foreign exchange (eFX) markets.

## Project Overview

This project is based on the work of Barzykin, Bergault, and Guéant [1], who formulate an optimal market making problem as a stochastic control problem, leading to a Hamilton-Jacobi-Bellman (HJB) partial differential equation (PDE) for the market maker's value function. The model optimizes bid and ask quotes across multiple currencies while accounting for inventory risk, exchange rate volatility, and cross-currency correlations.

### Thesis Objectives

1. **Grid-based numerical PDE solver for the HJB equation.** Due to the curse of dimensionality, the numerical analysis focuses on low-dimensional settings (2-3 currencies). The numerically computed HJB solution in the two-currency case will be compared to the ODE approximation from [1] to evaluate its accuracy and validate the PDE solver.

2. **Sensitivity analysis and uncertainty quantification.** Once reliable numerical solutions are obtained, study the robustness of optimal strategies with respect to model parameters (risk aversion, volatility, order flow intensities) under different market conditions.

3. **Physics-informed neural networks (PINNs) (extension).** If time permits, solve the HJB system using PINNs to handle higher dimensions without grid-based discretization, and compare results to the grid-based solver.

## Mathematical Framework

The market maker quotes bid and ask prices on multiple FX currency pairs. The model uses:

- **Logistic client-demand model:** Trade arrival intensity per tier, direction, and size is `Lambda(z, delta) = lambda(z) * f(z, delta)` where `f(z, delta) = 1 / (1 + exp(alpha + beta * delta))`.
- **Hamiltonian:** `H(z, p) = sup_delta f(z, delta) * (delta - p)`, where the optimal quote `delta_bar(z, p)` is the argmax.
- **Quadratic approximation** of the Hamiltonian around `p=0`: `H_check(z, p) = alpha_0(z) + alpha_1(z)*p + 0.5*alpha_2(z)*p^2`, yielding closed-form matrices `M`, `M_tilde`, `P`.
- **Matrix Riccati ODE** for `A(t)` and `B(t)` (Eq. (5) in [1]), solved backward from terminal time `T`.
- **Optimal quotes** and **hedging rates** are then computed from `A(t)`, `B(t)`, and the current inventory vector `Y`.

### Key Quantities

| Symbol | Description | Units |
|--------|-------------|-------|
| `Y` | Inventory vector (one entry per currency) | M$ |
| `delta` | Client markup / quote | decimal (multiply by 1e4 for bps) |
| `sigma_i` | Volatility of currency `i` vs reference | decimal / sqrt(day) |
| `gamma` | Risk aversion parameter | 1/M$ |
| `xi` | Hedging (externalization) rate | M$ / day |
| `A(t)` | Riccati solution matrix (d x d, symmetric) | - |
| `B(t)` | Riccati solution vector (d) | - |
| `L(xi)` | D2D execution cost: `psi*|xi| + eta*xi^2` | decimal |

## Repository Structure

```
optimal-market-making-eFX/
├── README.md
├── requirements.txt                              # numpy, matplotlib
├── src/                                          # Shared Python modules
│   ├── __init__.py
│   ├── model.py                                  # ModelParams, parameter builders
│   ├── hamiltonian.py                            # Logistic intensity, optimal quotes
│   ├── riccati.py                                # Sigma, M, Mtilde, P, Riccati ODE solver
│   ├── policy.py                                 # Quote/hedge policies, MMResult, run_multicurrency_mm
│   ├── simulation.py                             # Monte Carlo inventory simulation
│   └── plotting.py                               # Visualization helpers
├── ode_approximation_reproduction.ipynb           # ODE-based reproduction of [1]
├── hjb_pde_solver.ipynb                           # Grid-based HJB PDE solver (2-3 currencies)
├── pde_vs_ode_comparison.ipynb                    # Comparison: PDE solution vs ODE approximation
├── sensitivity_analysis.ipynb                     # Sensitivity analysis & uncertainty quantification
├── pinns_solver.ipynb                             # PINNs-based HJB solver (extension)
├── references/                                    # Reference papers
└── .gitignore
```

## Notebooks

| Notebook | Description |
|----------|-------------|
| `ode_approximation_reproduction.ipynb` | Reproduces the ODE-based quadratic approximation from [1]. Includes parameters from Table 1, Riccati ODE solver, optimal quotes, hedging rates, and Monte Carlo simulation. |
| `hjb_pde_solver.ipynb` | Grid-based finite difference solver for the full HJB PDE in 2-3 currency settings. |
| `pde_vs_ode_comparison.ipynb` | Compares the PDE solution against the ODE approximation: value functions, optimal controls, error analysis, and P&L under both policies. |
| `sensitivity_analysis.ipynb` | Parameter sweeps and uncertainty quantification for risk aversion, volatility, correlations, order flow, and hedging costs. |
| `pinns_solver.ipynb` | Physics-informed neural network solver for the HJB equation, targeting higher-dimensional cases beyond the grid-based approach. |

## Usage

```bash
pip install -r requirements.txt
jupyter notebook ode_approximation_reproduction.ipynb
```

## References

[1] A. Barzykin, P. Bergault, O. Guéant, "Dealing with multi-currency inventory risk in FX cash markets," arXiv:2207.04100v4, 2023.
