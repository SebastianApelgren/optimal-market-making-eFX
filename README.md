# Optimal Market Making in eFX

Master thesis project on numerical methods for optimal multi-currency market making in electronic foreign exchange (eFX) markets.

## Project Overview

This project is based on the work of Barzykin, Bergault, and Guéant [1], who formulate an optimal market making problem as a stochastic control problem, leading to a Hamilton-Jacobi-Bellman (HJB) partial differential equation (PDE) for the market maker's value function. The model optimizes bid and ask quotes across multiple currencies while accounting for inventory risk, exchange rate volatility, and cross-currency correlations.

### Thesis Structure

1. **Model and ODE approximation.** Present the multi-currency market making model and the quadratic ansatz that reduces the HJB PDE to a system of matrix Riccati ODEs. Reproduce the paper's results (Table 1, Figures 1-2).

2. **Sensitivity analysis on the ODE system.** Study how the optimal strategy depends on model parameters (risk aversion, volatility, correlations, order flow intensities, hedging costs). Identify which parameters the solution is most sensitive to.

3. **Grid-based PDE solver.** Solve the full HJB PDE numerically using finite differences with fully implicit Euler and policy iteration (Howard's algorithm). Compare against the ODE approximation to quantify the error of the quadratic ansatz.

4. **Connecting sensitivity to the full problem.** If PDE ≈ ODE across parameter variations, the ODE sensitivity analysis extends to the full nonlinear problem.

## Mathematical Framework

The market maker simultaneously quotes bid and ask prices on multiple FX currency pairs. Incoming client trades arrive according to a **logistic intensity model**: each trade size and tier has arrival rate `Lambda(z, delta) = lambda(z) * f(z, delta)` where `f(z, delta) = 1 / (1 + exp(alpha + beta * delta))`. The market maker also hedges inventory risk by trading on the dealer-to-dealer (D2D) market at a cost `L(xi) = psi*|xi| + eta*xi^2`.

The market maker's objective is to maximize expected terminal P&L minus a risk penalty on inventory, leading to a **Hamilton-Jacobi-Bellman (HJB) PDE** for the value function `theta(t, S, Y)`, where `S` is the vector of exchange rates and `Y` is the inventory vector.

### Solution approaches

1. **ODE approximation** [1]: Under a quadratic ansatz `theta ~ -Y^T A(t) Y - Y^T B(t) - C(t)`, the HJB PDE reduces to a system of matrix Riccati ODEs for `A(t)` and `B(t)`, solved backward from terminal time. This yields closed-form optimal quotes and hedging rates. Scales to many currencies but relies on the quadratic approximation of the Hamiltonian.

2. **Full HJB PDE** (grid-based finite differences): Solves the PDE directly on a discretized state space using fully implicit Euler with policy iteration. Uses the exact (non-quadratic) Hamiltonian. Limited to low dimensions (2-3 currencies) due to the curse of dimensionality, but captures the full nonlinear structure. Comparison with the ODE shows ~1% agreement for the paper's parameters.

### Key Quantities

| Symbol | Description | Units |
|--------|-------------|-------|
| `Y` | Inventory vector (one entry per currency) | M$ |
| `delta` | Client markup / quote | decimal (multiply by 1e4 for bps) |
| `sigma_i` | Volatility of currency `i` vs reference | decimal / sqrt(day) |
| `gamma` | Risk aversion parameter | 1/M$ |
| `xi` | Hedging (externalization) rate | M$ / day |
| `L(xi)` | D2D execution cost: `psi*|xi| + eta*xi^2` | decimal |

## Repository Structure

```
optimal-market-making-eFX/
├── README.md
├── CLAUDE.md                                        # Development instructions
├── requirements.txt                                 # numpy, matplotlib, scipy
├── src/                                             # Shared Python modules
│   ├── __init__.py
│   ├── model.py                                     # ModelParams, parameter builders
│   ├── hamiltonian.py                               # Logistic intensity, optimal quotes (Lambert W)
│   ├── riccati.py                                   # Sigma, M, Mtilde, P, Riccati ODE solver
│   ├── policy.py                                    # Quote/hedge policies, MMResult
│   ├── pde.py                                       # PDE spatial operators and solvers
│   ├── simulation.py                                # Monte Carlo inventory simulation
│   └── plotting.py                                  # Visualization helpers
├── ode_approximation_reproduction.ipynb              # ODE reproduction of [1]
├── ode_solver_verification.ipynb                     # Manufactured solution tests for ODE
├── hjb_pde_solver.ipynb                              # Explicit Euler PDE (scaled η)
├── eta_continuation_solver.ipynb                     # Implicit PDE solver (original η)
├── notes/                                           # Design notes and working logs
│   ├── solvers/                                     # Solver-specific notes
│   └── spatial-operators/                           # PDE term derivations
└── references/                                      # Reference papers
```

## Notebooks

| Notebook | Description |
|----------|-------------|
| `ode_approximation_reproduction.ipynb` | Reproduces the ODE-based quadratic approximation from [1]. Includes Table 1 parameters, Riccati ODE solver, optimal quotes (Figure 1), hedging rates (Figure 2), and Monte Carlo simulation. |
| `ode_solver_verification.ipynb` | Verifies the Riccati ODE solver using manufactured solutions with known closed-form answers. Confirms first-order convergence. |
| `hjb_pde_solver.ipynb` | Explicit Euler PDE solver with η×1000 scaling. Validates all spatial operators (diffusion, quoting, hedging, drift) against the ODE solution. Includes value function, quote, and hedge rate comparisons. |
| `eta_continuation_solver.ipynb` | Fully implicit PDE solver with the paper's original parameters. Uses monotone (M-matrix) differencing and policy iteration. Produces PDE vs ODE comparison showing ~1% agreement. |

## Usage

```bash
pip install -r requirements.txt
jupyter notebook ode_approximation_reproduction.ipynb
```

## References

[1] A. Barzykin, P. Bergault, O. Guéant, "Dealing with multi-currency inventory risk in FX cash markets," arXiv:2207.04100v4, 2023.
