# HJB PDE Solver — Summary

Overview of what the solver does, what was tried, and what worked.
For detailed derivations and experiments, see the notes in `spatial-operators/` and `solvers/`.


## What the solver does

We solve the HJB PDE (Eq. 1) backward from T to 0 using fully implicit Euler with policy iteration. At each time step, we solve a nonlinear equation for θ^{m+1}:

$$\theta^{m+1} - \Delta t \cdot F[\theta^{m+1}] = \theta^m$$

This is nonlinear because F contains the quoting and hedging Hamiltonians, which involve suprema over controls (optimal quotes and hedge rates).


## Policy iteration (Howard's algorithm)

Policy iteration linearises the HJB by alternating two steps:

1. **Fix controls.** Given current θ, compute the optimal quotes δ*(y) and hedge rates ξ*(y). These are known formulas: δ* via Lambert W (see `spatial-operators/hamiltonian_supremum_logistic.md`), ξ* via the piecewise formula from the Legendre-Fenchel transform of L(ξ) = ψ|ξ| + ηξ².

2. **Solve the linear system.** With controls frozen, the HJB becomes a linear equation:

$$(I - \Delta t \cdot A) \, \theta^{m+1} = \theta^m + \Delta t \cdot \text{source}$$

where A is a sparse matrix encoding how each grid point's θ depends on its neighbors (through finite differences), and `source` collects all terms that don't depend on θ^{m+1}.

Repeat until θ stops changing. This is equivalent to Newton's method on the HJB equation.


## The sparse matrix (I - dt·A)

The matrix is N×N where N = 301² = 90,601 (one row per grid point). A is the linearised spatial operator with contributions from:

- **Diffusion** (½σ²y² ∂²θ/∂y²) — second derivative stencil, 3 entries/row
- **Quoting** (λf*·(θ(y) - θ(y+zd))) — couples each point to shifted neighbors, up to 40 entries/row
- **Hedging advection** (ξ*·(1+ky)·∂θ/∂y) — first derivative stencil, 2 entries/row/axis
- **Drift** (μy·∂θ/∂y) — first derivative stencil (zero for the paper's μ=0)

The **source** vector collects known terms:
- Running penalty: -γ/2 · y^T Σ y
- Quoting profit: z·λ·f*·δ*
- Hedging cost: -ψ|ξ*| - η(ξ*)²
- Drift source: Σᵢ μᵢyᵢ (zero for the paper's example)


## What went wrong and how it was fixed

### Problem: central differences broke the matrix

With the paper's original η = 10⁻⁹, the optimal hedge rates reach ~50,000 M$/day. When the hedging advection v·∂θ/∂y is discretised with central differences, the off-diagonal matrix entries (~60) dwarf the diagonal (~1). The matrix loses diagonal dominance, `spsolve` produces garbage, and PI diverges.

### Fix: monotone (M-matrix) differencing

The hedging advection term uses directional differences instead of central differences:
- v > 0: forward diff (θ_{j+1} - θ_j)/dy
- v < 0: backward diff (θ_j - θ_{j-1})/dy

This guarantees (I - dt·A) is an **M-matrix**: non-positive off-diagonals and dominant diagonal. The M-matrix property ensures a discrete maximum principle, which is required for PI convergence. This is what the paper means by "monotone implicit Euler scheme."

Note: this is the **opposite** of the standard transport-equation upwind convention, because our advection term enters with a positive sign in the RHS (v·∂θ/∂y makes θ grow).

### Fix: overflow protection in Lambert W

For very negative quoting momentum p, the exp argument in H_logistic overflows float64. The fix clamps the argument to 700 and uses the asymptotic W(x) ≈ ln(x) for large x, giving physically correct values (H → arg/β, f* → 1, δ* → -α/β).


## What was tried but wasn't needed

### η-continuation (dropped)

At each time step, solve the implicit equation multiple times with decreasing η-scale factors (e.g., η×10, η×5, η×2, η×1), using each solution to warm-start PI for the next level. This was the supervisor's suggestion (continuation/homotopy method).

Testing showed that the monotone differencing alone is sufficient for stability — the solution with schedule (1,) matches the continuation results exactly. The η-continuation adds computational cost without improving convergence. Dropped.

### Semi-implicit IMEX scheme (dropped)

Treat diffusion implicitly, Hamiltonians explicitly. This was tried before the fully implicit solver. Testing showed it provides no improvement over explicit Euler — the binding stability constraint is the hedging Hamiltonian (1/(4η) ≈ 2.5×10⁸), not the diffusion. The IMEX scheme required the same η×1000 scaling as explicit Euler. See `solvers/semi-implicit/` for details. Dropped in favour of the fully implicit scheme.


## PI convergence behaviour

- PI converges the value function in ~5 iterations (rel_diff drops from 3.6e-2 to 2.8e-4)
- The residual then stalls at ~2.5e-5 due to first-order accuracy of the upwind scheme at boundary points
- `theta(0,0)` is unchanged from iteration 5 onward — the stall is at the grid edges, not in the region of interest
- PI tolerance set to 1e-4 (above the stall level): step 0 needs ~7 PI, subsequent steps need 2-4 PI
- 20-step run completes in ~9 minutes on 301×301 grid


## First result: PDE vs ODE comparison (20 steps, original η)

- PDE and ODE value functions agree to within ~1% over [-100, 100] M$
- Max |θ_PDE - θ_ODE| ≈ 1.1e-4 (relative to θ ≈ 0.01)
- The PDE value function is slightly higher than the ODE at large |y| (peaks at |y| ≈ 75 M$)
- Physical interpretation: the ODE's quadratic approximation of H underestimates the true Hamiltonian at large |p|, so the market maker extracts slightly more value from quoting at extreme inventories than the approximation predicts
- The ODE approximation is quite accurate for the paper's parameter set


## PDE vs ODE comparison methodology

See `pde_ode_comparison_design.md` for the full design of the improved
comparison. Key elements:

- 5 QoIs aligned with the 3-currency SA (tier spread differential, inventory
  skew, hedge rate, net revenue, plus single-tier spread as diagnostic)
- All 8 parameters tested (OAT, 17 PDE solves) + 15-20 Latin hypercube
  points to probe multi-parameter corners
- Comparison metrics: summary error table, sensitivity rankings with
  Spearman rho, parity plots, value function theta profile, Hessian A
  comparison
- Bridging d=2 to d=3: theoretical argument + Riccati matrix inspection +
  2-ccy vs 3-ccy ODE comparison (both cheap, ODE only)
- Hedge rate error (up to 76% at extremes): framed in absolute terms,
  shown to preserve sensitivity rankings, SA conclusions interpreted
  qualitatively for this QoI

Compute cached in `data/pde_comparison/`. Analysis notebook loads from
disk and stays fast.


## Notes index

| Topic | File |
|-------|------|
| **Shared** | |
| Grid, domain, boundary treatment | `grid_and_domain.md` |
| Figure 1 discrepancy investigation | `../figure1_investigation.md` |
| **Spatial operators** | |
| Quoting Hamiltonian integral | `../spatial-operators/pde_quoting_integral_design.md` |
| Hamiltonian supremum (Lambert W) | `../spatial-operators/hamiltonian_supremum_logistic.md` |
| Hedging Hamiltonian | `../spatial-operators/hedging_hamiltonian.md` |
| Diffusion term | `../spatial-operators/diffusion_term.md` |
| Drift term | `../spatial-operators/drift_term.md` |
| Terminal condition | `../spatial-operators/terminal_condition.md` |
| **Explicit Euler** | |
| Explicit solver design + CFL | `explicit_euler.md` |
| **Semi-implicit (IMEX)** | |
| IMEX design | `semi-implicit/semi_implicit_time_stepping.md` |
| IMEX results (didn't help) | `semi-implicit/semi_implicit_results_and_next_steps.md` |
| **Implicit + policy iteration** | |
| Implicit solver design (PI) | `implicit/implicit_euler_policy_iteration.md` |
| Implementation log | `implicit/implicit_solver_implementation_log.md` |
| η-continuation and monotone scheme | `implicit/eta_continuation_and_monotone_scheme.md` |
| **PDE vs ODE comparison** | |
| Comparison design and methodology | `pde_ode_comparison_design.md` |
| Why PDE stays at d=2 | `pde_dimension_discussion.md` |
