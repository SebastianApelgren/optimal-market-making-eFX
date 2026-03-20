# Semi-Implicit Solver: Results and Next Steps

## What was implemented (2026-03-13)

Added to `src/pde.py`:

- **`pde_rhs_nonlinear(theta, spec)`** — evaluates the explicit (nonlinear) part of the RHS: penalty, drift, quoting Hamiltonian, hedging Hamiltonian, and cross-axis diffusion. Excludes diagonal diffusion (handled implicitly).
- **`_diffusion_term_cross(...)`** — computes only the cross-derivative part of the diffusion (zero for d=2 with σ_USD=0; needed for d≥3).
- **`_thomas_factor(lower, diag, upper)`** — precomputes forward-sweep factors for the Thomas algorithm. Called once before the time loop.
- **`_thomas_solve_batch(lower, c_prime, denom_inv, rhs)`** — solves a tridiagonal system for multiple RHS vectors simultaneously (vectorized over the batch dimension, Python loop over the n=301 grid points).
- **`_ImplicitAxisOp`** — frozen dataclass storing Thomas factors for one grid axis.
- **`_build_implicit_ops(spec, dt)`** — builds the tridiagonal matrix $(I - \Delta t \mathcal{L}_k)$ for each axis $k$ with nonzero diagonal diffusion ($\Sigma_{kk} > 0$). Uses central differences for interior points, identity rows at boundaries (j=0 and j=n-1). Returns a tuple of `_ImplicitAxisOp`.
- **`_apply_implicit_step(theta, ops)`** — applies the implicit solve sequentially per axis. Uses `np.moveaxis` to generalize across dimensions.
- **`solve_hjb_semi_implicit(y_grids, mp, n_steps, snapshot_times)`** — main IMEX Euler time stepper. Same interface and return format as `solve_hjb_explicit`.

Updated `src/__init__.py` to export `pde_rhs_nonlinear` and `solve_hjb_semi_implicit`.

### Design decisions

- **Implicit operator = diagonal diffusion only.** No drift gradient (zero for the paper's μ=0), no cross-derivative diffusion (zero for d=2). For d≥3 or μ≠0, the drift gradient and cross terms stay explicit. The drift gradient could be moved implicit for better advection CFL; see `notes/semi_implicit_time_stepping.md`.
- **Boundary treatment = identity rows.** The first and last grid points (j=0, j=n-1) are not modified by the implicit step — they pass through from the explicit RHS. This is acceptable because they lie in the 50 M$ buffer zone ($|y| > 100$ M$ region of interest). The one-sided FD stencils used elsewhere would break the tridiagonal bandwidth.
- **Thomas algorithm in pure numpy** (no scipy dependency). Python loop over 301 grid points, vectorized over 301 batch slices. About 20 steps/s for the 301×301 grid.
- **Sequential axis application for d≥3.** When multiple axes have nonzero diffusion, the tridiagonal solves are applied sequentially (Lie splitting). This introduces an $O(\Delta t^2)$ splitting error, acceptable for our first-order-in-time scheme.


## What was tested

### 1. Step-by-step debugging (original η, 301×301 grid, 1000 steps)

Ran individual time steps to trace the blowup:

| Step | Q (quoting) max | H (hedging) max | θ max |
|------|----------------|-----------------|-------|
| 0 | 3.56e-01 | 1.06e+00 | 6.67e-04 |
| 1 | 3.54e-01 | 8.03e-01 | 1.35e-03 |
| 2 | 3.72e-01 | 6.43e-01 | 2.04e-03 |
| 3 | 4.36e-01 | 5.69e-01 | 2.74e-03 |
| 4 | 5.23e-01 | 6.87e-01 | 3.45e-03 |
| 5 | 6.25e-01 | 1.17e+00 | 4.17e-03 |
| 6 | 7.39e-01 | 3.21e+00 | 4.89e-03 |
| 7 | 9.33e-01 | 1.92e+01 | 5.61e-03 |
| 8 | 2.23e+00 | 4.94e+02 | 1.98e-02 |
| 9 | NaN | 2.62e+05 | NaN |

The hedging Hamiltonian grows exponentially: 0.6 → 3.2 → 19 → 494 → 262,000 → NaN. The quoting Hamiltonian also produces NaN at step 9, but it's a consequence of θ blowing up from H, not the root cause.

**Root cause:** The hedging Hamiltonian $H(p) = (\max(|p| - \psi, 0))^2 / (4\eta)$ with $\eta = 10^{-9}$ (decimal) amplifies gradients of θ by a factor $\sim 1/(4\eta) \approx 2.5 \times 10^8$. Even tiny gradients in θ produce large H values, which feed back into θ on the next step. This is a positive feedback loop (stiff nonlinearity). Implicit diffusion damps high-frequency modes but cannot counteract this $O(10^8)$ gain in the explicit Hamiltonian.


### 2. η scaling sweep (IMEX, 1000 steps, 301×301 grid)

Tested which η scaling factor makes the IMEX scheme stable:

| η scale | Stable? | θ(0,0) |
|---------|---------|--------|
| 1 (original) | No | NaN |
| 2 | No | NaN |
| 5 | No | NaN |
| 10 | No | NaN |
| 50 | No | NaN |
| 100 | No | NaN |
| 500 | No | NaN |
| 1000 | **Yes** | 0.010082 |

**Conclusion:** With 1000 time steps, the IMEX scheme needs the same η×1000 scaling as the fully explicit scheme. The implicit diffusion does not improve the stability threshold for the hedging Hamiltonian.


### 3. Correctness verification (η×1000, 1000 steps)

Compared IMEX vs explicit Euler with the same (scaled) parameters:

```
Max |θ_IMEX - θ_explicit| = 2.10e-06
Mean |θ_IMEX - θ_explicit| = 3.00e-10
```

The two schemes agree to ~2e-6 (expected: they are different numerical schemes with the same order of accuracy, applied to the same PDE). This confirms the IMEX infrastructure is correct.


### 4. Why implicit diffusion doesn't help

The CFL analysis (from `notes/hjb_pde_solver_design.md`) already showed that the diffusion CFL is very generous:

$$\Delta t_{\text{CFL, diffusion}} < \frac{\Delta y^2}{\sigma^2 Y_{\max}^2} = \frac{1}{0.008^2 \times 150^2} \approx 0.69 \text{ days}$$

With $\Delta t = 5 \times 10^{-5}$ days, the diffusion CFL number is $\sim 3.6 \times 10^{-5}$ — the diffusion was never the binding constraint. The stiffness is entirely from the hedging Hamiltonian's nonlinear gain $1/(4\eta)$.

Making diffusion implicit removes a non-binding CFL constraint. The binding constraint (hedging Hamiltonian stiffness) remains fully explicit.


## What's needed: policy iteration (Howard's algorithm)

To handle the original η, the hedging Hamiltonian must be treated implicitly. The standard approach for HJB equations is **policy iteration** (also called Howard's algorithm), which is what the paper describes as their "monotone implicit Euler scheme."

### The idea

At each time step, instead of evaluating the nonlinear Hamiltonians at $\theta^m$ (explicit), we:

1. **Fix the control policy** from $\theta^m$:
   - Quoting: for each (pair, tier, size), compute the optimal markup $\delta^*(y)$ from $\theta^m$
   - Hedging: for each pair, compute the optimal hedge rate $\xi^*(y)$ from $\nabla\theta^m$

2. **Solve the resulting LINEAR PDE** for $\theta^{m+1}$:
   - With fixed controls, the HJB PDE becomes a linear PDE (no supremum/infimum)
   - The quoting term becomes: $\sum \lambda f(\delta^*)\bigl(\delta^* - p(\theta^{m+1})\bigr)$, where $p(\theta) = (\theta(y) - \theta(y + z\,d_{ij}))/z$ is linear in $\theta$
   - The hedging term becomes: $\sum \bigl(\xi^* \cdot p(\nabla\theta^{m+1}) - L(\xi^*)\bigr)$, where $p$ depends linearly on $\nabla\theta^{m+1}$
   - This linear PDE can be solved implicitly via a sparse linear system

3. **Iterate**: recompute the policy from $\theta^{m+1}$, re-solve, until convergence (usually 2-5 iterations per time step).

### Why it works

- With fixed controls, the PDE is linear → the implicit system is a single matrix solve (no Newton iteration needed).
- The matrix changes at each iteration (because the policy changes), so we can't precompute LU factors, but the system is sparse and can be solved efficiently.
- Policy iteration converges superlinearly (quadratically in practice) because it's equivalent to Newton's method on the HJB equation.

### Implementation considerations

- The linear system with fixed controls couples all grid points (the quoting term shifts θ by $z \cdot d_{ij}$, connecting distant grid points). It's no longer tridiagonal.
- For d=2 with 301×301 = 90,601 unknowns, the sparse system can be solved with `scipy.sparse.linalg.spsolve` or iterative methods (GMRES, BiCGSTAB).
- This requires adding `scipy` as a dependency.
- The existing Thomas solver infrastructure can still be used for the diffusion-only part if we use operator splitting within the policy iteration.

### Simpler alternative: linearize hedging only

A middle ground between pure IMEX and full policy iteration:

1. At each time step, compute the hedging momentum $p^m$ from $\nabla\theta^m$
2. **Linearize** $H(p) \approx H(p^m) + H'(p^m)(p - p^m)$ around $p^m$
3. The linearized hedging term is $H'(p^m) \cdot p^{m+1}$, which is linear in $\nabla\theta^{m+1}$
4. Absorb this into the implicit system

**Problem:** For d=2, the hedging momentum is $p = \partial\theta/\partial y_{\text{EUR}} - \partial\theta/\partial y_{\text{USD}} + \ldots$, which involves gradients along both axes. This couples the axes and breaks the per-axis tridiagonal structure. A full sparse system solve would be needed.

### Recommendation

Go directly to policy iteration. The "linearize hedging only" approach has the same implementation complexity (need a sparse solver) but less generality. Policy iteration handles both Hamiltonians simultaneously and is the standard method for these problems. The paper used it, so we know it works.

The existing code (`pde_rhs_nonlinear`, Thomas solver, `_ImplicitAxisOp`, etc.) is still useful — the diffusion part of the implicit system remains tridiagonal and can be used as a preconditioner or within an operator-splitting variant of policy iteration.
