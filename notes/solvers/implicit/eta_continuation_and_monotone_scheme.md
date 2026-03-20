# η-Continuation and Monotone Discretisation for the Implicit Solver

Working notes from debugging the implicit Euler solver with the paper's original parameters.


## 1. Starting point: implicit solver blows up at original η

The implicit Euler solver with policy iteration (Howard's algorithm) was validated with η×1000 (scaled execution cost). With original η = 10⁻⁹ (decimal), it produced NaN within 3-4 time steps. The implementation log (`implicit_solver_implementation_log.md`) documented the progression of fixes (upwind, clamping, Lambert W overflow protection) but those were made by an AI agent without review and would be reverted for independent evaluation.


## 2. Diagnosing the root cause

Starting from the reverted code (no overflow protection, central differences, no clamping), we ran the implicit solver step by step and traced the failure chain:

1. The running penalty drives θ to ~-0.01 at the grid boundary after a few steps
2. This creates gradients ∇θ ~ 10⁻⁴ in the value function
3. The hedging momentum p ≈ ∇θ produces optimal hedge rates ξ* = (|p| - ψ)/(2η) ~ 50,000 M$/day (because 1/(2η) ≈ 5×10⁸)
4. These hedge rates enter the sparse matrix as advection coupling entries with magnitude dt·ξ*/dy ~ 60, while the diagonal is ~1
5. The matrix is not diagonally dominant → `spsolve` produces garbage
6. The garbage θ has huge gradients → quoting momentum p pushes `exp(-(1+α+βp))` past float64 overflow → NaN
7. NaN poisons the entire grid

The Lambert W overflow is a symptom. The root cause is the ill-conditioned sparse matrix from the hedging advection term.


## 3. Supervisor's advice: continuation method

The supervisor suggested:
1. Solve a simpler problem first (damped nonlinearity) — already done with η×1000
2. Use implicit + Newton for stability — we have policy iteration (= Newton for HJB)
3. **From damped nonlinearity, successively iterate to the desired nonlinearity** — the key new idea

This is a standard continuation (homotopy) method: at each time step, solve the implicit equation multiple times with decreasing η-scale factors. Each sub-solve's converged θ warm-starts PI for the next level. The implicit equation is always θ^{m+1} - dt·F_η[θ^{m+1}] = θ^m — only the spatial operator F_η changes; θ^m stays fixed.


## 4. Implementing and testing η-continuation

Added `_scale_hedging_eta(spec, scale)` helper and `solve_hjb_implicit_continuation(...)` solver to `src/pde.py`. The solver builds PDESpecs for each η level (sharing the expensive quoting spec), then at each time step loops over the schedule from most-damped to original.

### Test 1: schedule (1000, 100, 10, 1), no overflow protection

```
Step 0, η×1000: theta(0,0) = 5.749e-04, PI=14  ← OK
Step 0, η×100:  theta(0,0) = 5.751e-04, PI=9   ← OK, barely changed
Step 0, η×10:   theta(0,0) = 5.769e-04, PI=11  ← OK
Step 0, η×1:    NaN, PI=20, singular=True        ← blew up
```

The η×10 → η×1 jump is too large. The overflow in H_logistic produces NaN on the first bad PI iterate, killing everything.

### Test 2: same schedule, with overflow protection restored

Restored the Lambert W overflow protection (clamp exp argument to 700, use W(x) ≈ arg for large x). This prevents NaN but reveals the underlying divergence:

```
Step 0, η×10:   theta(0,0) = 5.769e-04     ← OK
Step 0, η×1:    theta(0,0) = 6.9e+13, PI=20 ← exploded (finite but garbage)
```

PI hits max iterations and θ reaches 10¹⁴. The overflow protection keeps values finite but doesn't fix convergence.

### Test 3: finer schedule (10, 5, 2, 1), with overflow protection

```
Step 0, η×10: theta(0,0) = 5.769e-04, PI=15  ← OK
Step 0, η×5:  theta(0,0) = 5.785e-04, PI=11  ← OK
Step 0, η×2:  theta(0,0) = 6.9e+18, PI=50    ← exploded
```

Even the η×5 → η×2 jump diverges. The finer schedule helps with η×10 → η×5 but the problem is not the schedule — it's the linear system itself.


## 5. The real problem: central differences break monotonicity

The continuation can't work because the underlying linear system at each PI step is ill-conditioned. The hedging gradient term is an advection operator v·∂θ/∂y with velocity v = ξ*·(1+k·y). With central differences:

- Off-diagonal entries at j±1: ±dt·v/(2dy), which can be **positive or negative**
- Diagonal: no contribution from the advection term
- With v ~ 50,000 and dt = 2.5e-3, dy = 1: |off-diag| ≈ 62, diagonal ≈ 1

The matrix is not diagonally dominant and not an M-matrix. This breaks policy iteration — **Howard's algorithm converges only for monotone schemes** (i.e., schemes where the discrete operator satisfies a comparison principle, which requires the matrix to be an M-matrix).

This is what the paper means by "monotone implicit Euler scheme." The word "monotone" refers to the spatial discretisation, not just the time stepping.


## 6. Fix: monotone (M-matrix) differencing for the hedging advection

For the spatial operator A = +v·∂/∂y entering the RHS F[θ], the implicit system is (I - dt·A)θ = rhs. For (I - dt·A) to be an M-matrix (non-positive off-diagonals, dominant diagonal):

- **v > 0: forward diff** (θ_{j+1} - θ_j)/dy
  → A[j,j] = -v/dy, A[j,j+1] = +v/dy
  → (I-dtA): diag = 1 + dt·v/dy > 1, off-diag = -dt·v/dy < 0  ✓

- **v < 0: backward diff** (θ_j - θ_{j-1})/dy
  → A[j,j] = v/dy, A[j,j-1] = -v/dy
  → (I-dtA): diag = 1 + dt·|v|/dy > 1, off-diag = dt·v/dy < 0  ✓

Note: this is the **opposite** of the standard transport-equation upwind convention (where v > 0 uses backward diff). The difference arises because our advection term has a positive sign in the RHS (+v·∂θ/∂y makes θ grow), whereas the transport equation has ∂u/∂t = -v·∂u/∂x.

Combined entries in (I - dt·A):
- diagonal:  += dt·|v|/dy          (always positive)
- j+1:       -= dt·max(v, 0)/dy   (≤ 0)
- j-1:       += dt·min(v, 0)/dy   (≤ 0)

This is first-order accurate in space (vs second-order for central differences), but it guarantees the M-matrix property that policy iteration requires for convergence.

Changed `assemble_implicit_system` in `src/pde.py` to use this monotone differencing for the hedging gradient terms.


## 7. Testing monotone scheme + continuation

### Test with schedule (10, 5, 2, 1), PI_TOL = 1e-6

The η×10 and η×5 levels converge fast (6-9 PI), but η×2 and η×1 hit 50 PI without converging to 1e-6. However, `theta(0,0)` is fully stable — it doesn't change between PI iteration 5 and 50. The remaining ~2.5e-5 rel_diff is at the grid boundary where the first-order upwind stencil is less accurate. The continuation warm-start is working (θ barely changes between η levels: 5.767e-4 → 5.782e-4 → 5.818e-4 → 5.858e-4), but PI stalls before reaching the tolerance.

### Test with schedule (1,) — no continuation, monotone scheme only

The monotone scheme alone is stable with original η. No blowup, no NaN, no singular matrices. Values match the continuation run exactly (`theta(0,0)` = 5.858e-4 at step 0 in both cases). The continuation was not needed — the M-matrix property was the key fix.

**Conclusion: the η-continuation is unnecessary for stability.** The monotone differencing is sufficient. The continuation is dropped from further tests.


## 8. PI convergence analysis

With PI_TOL = 1e-6, PI stalls at ~2.5e-5 rel_diff (oscillating, never reaching tolerance). Detailed iteration trace for step 0:

```
PI 1:  rel_diff = 3.6e-2     theta(0,0) = -1.047e-3
PI 2:  rel_diff = 1.6e-2     theta(0,0) =  3.907e-4
PI 3:  rel_diff = 6.3e-3     theta(0,0) =  4.981e-4
PI 4:  rel_diff = 1.7e-3     theta(0,0) =  5.659e-4
PI 5:  rel_diff = 2.8e-4     theta(0,0) =  5.843e-4
PI 10: rel_diff = 2.4e-5     theta(0,0) =  5.858e-4   ← converged in value
PI 50: rel_diff = 2.7e-5     theta(0,0) =  5.858e-4   ← rel_diff stalled
```

PI converges the value function in ~5 iterations. The residual stalls at ~2.5e-5 due to the first-order spatial accuracy of the upwind scheme at boundary/buffer points. `theta(0,0)` is identical from iteration 5 onward — the remaining rel_diff is at the grid edges and doesn't affect the solution in the region of interest.

**Fix: set PI_TOL = 1e-4.** This is above the stall level, so PI converges in:
- Step 0: 7 iterations (cold start from θ(T) = 0)
- Steps 1+: 2-4 iterations (warm-started from previous step)
- Some late steps (15-19): up to 16-50 iterations as the solution evolves further

20-step run with PI_TOL = 1e-4 completes in ~9 minutes on 301×301 grid.


## 9. First PDE vs ODE comparison (original η, 20 time steps)

With 20 steps (dt = 2.5e-3 days), PI_TOL = 1e-4, original η:

**Value function along EUR axis (y_USD = 0):**
- PDE and ODE curves are visually indistinguishable in [-100, 100]
- Max |θ_PDE - θ_ODE| = 1.09e-4 (after shifting ODE by C₀ to match at origin)
- θ(0,0) = 1.059e-2

**Difference pattern:**
- Smooth, symmetric, non-negative: PDE > ODE everywhere
- Peaks at |y_EUR| ≈ 75 M$ with magnitude ~1.1e-4
- Minimum near y = 0 (~0) and y = ±50 (~0)
- Relative difference: ~1.1e-4 / 1.06e-2 ≈ 1%

**Physical interpretation:** The ODE uses the quadratic Taylor approximation Ĥ(p) ≈ α₀ + α₁p + ½α₂p², which underestimates the true Hamiltonian H(p) at large |p| (large inventory). The exact H from the PDE captures the full nonlinear structure, giving a slightly higher value function — the market maker extracts slightly more value from quoting at extreme inventories than the quadratic approximation predicts. The ~1% correction confirms that the ODE approximation is quite accurate for the paper's parameter set.

**Remaining work:**
- Run with more time steps (200) for better time discretisation accuracy
- Compare optimal quotes and hedge rates (policy comparison), not just value functions
- Investigate the late-step PI slowdown (steps 15-19 need more iterations)
