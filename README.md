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
├── requirements.txt                        # numpy, matplotlib
├── multicurrency_mm_reproduction.ipynb     # ODE-based reproduction of [1]
└── .gitignore
```

### `multicurrency_mm_reproduction.ipynb`

Jupyter notebook that reproduces the ODE-based approximation and key figures from [1]. Contains all logic in self-contained functions. Sections:

| Section | Description |
|---------|-------------|
| 1. Parameters / configuration | `ModelParams` dataclass with all model inputs. `build_paper_example_params()` reproduces Table 1 from [1] (5 currencies: USD, EUR, JPY, GBP, CHF; 2 tiers; sizes 1/5/10/20/50 M$). `restrict_currencies()` projects to a subset. |
| 2. Logistic intensity + optimal quotes | `logistic_f()`, `optimal_delta_logistic()` (bisection on `g'(delta)=0`), `H_logistic()` returning `(H, delta_bar, f_star)`. |
| 3. Quadratic approximation coefficients | `quadratic_coeffs_H_logistic()` computing `alpha_0, alpha_1, alpha_2` via envelope theorem + finite differences. |
| 4. Build Sigma, M, M_tilde, P and Riccati ODE | `build_Sigma()`, `build_M_tildeM_P()`, `solve_AB_euler()` (backward explicit Euler). |
| 5. Policy functions | `optimal_client_markup()` and `optimal_hedge_rate()` using `A(0), B(0)`. |
| 6. Main wrapper | `run_multicurrency_mm()` returning `MMResult` with helper methods `.markup()` and `.hedge_rate()`. |
| 7-8. Reproduction plots | Figures 1-2 from [1]: top-of-book quotes vs inventory, hedging rates vs inventory. |
| 9. Diagnostics | Eigenvalue checks for Sigma and A(0), verification that B(0) = 0 for symmetric parameters. |
| 10. Monte Carlo simulation skeleton | Tau-leaping Poisson simulator for inventory paths under the approximate policy. |

### Key Functions

| Function | Purpose |
|----------|---------|
| `build_paper_example_params()` | Constructs `ModelParams` matching Table 1 of [1] |
| `restrict_currencies(mp, keep)` | Projects model to a currency subset |
| `optimal_delta_logistic(p, alpha, beta)` | Computes optimal quote via bisection |
| `H_logistic(p, alpha, beta)` | Evaluates Hamiltonian and optimal quote |
| `quadratic_coeffs_H_logistic(alpha, beta)` | Computes quadratic approximation coefficients |
| `build_Sigma(mp)` | Constructs covariance matrix from vols and correlations |
| `build_M_tildeM_P(mp)` | Constructs aggregated matrices from client parameters |
| `solve_AB_euler(mp, M, Mtilde, P, Sigma)` | Solves Riccati ODE backward via explicit Euler |
| `run_multicurrency_mm(mp)` | End-to-end solver returning `MMResult` |
| `simulate_inventory_path_tau_leap(res, ...)` | Monte Carlo inventory simulation |

## Planned Development

The following components are planned but not yet implemented:

- **HJB PDE solver** (grid-based finite difference/element method for 2-3 currencies)
- **Comparison framework** between PDE solution and ODE approximation
- **Sensitivity analysis** (parameter sweeps over gamma, sigma, lambda, etc.)
- **Uncertainty quantification**
- **PINN-based solver** (extension, if time permits)

## Usage

```bash
pip install -r requirements.txt
jupyter notebook multicurrency_mm_reproduction.ipynb
```

To run the model programmatically (within the notebook or after extracting to `.py` files):

```python
mp = build_paper_example_params()
# Optionally restrict to fewer currencies:
# mp = restrict_currencies(mp, ["USD", "EUR", "GBP"])
res = run_multicurrency_mm(mp, n_steps=2000)

# Query optimal quotes
delta = res.markup(tier_idx=0, ccy_pay="EUR", ccy_sell="USD", z_musd=1.0, y=np.zeros(5))

# Query optimal hedging rate
xi = res.hedge_rate(ccy_buy="EUR", ccy_sell="USD", y=np.zeros(5))
```

## References

[1] A. Barzykin, P. Bergault, O. Guéant, "Dealing with multi-currency inventory risk in FX cash markets," arXiv:2207.04100v4, 2023.
