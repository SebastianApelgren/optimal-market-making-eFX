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


## 7. Next: retest with monotone scheme + continuation

With the M-matrix discretisation, the sparse solve should produce well-conditioned results at each PI step. The η-continuation then provides the good initial guesses that PI needs to converge quickly. Both pieces are needed:
- Monotone scheme → each PI step is well-conditioned
- Continuation → PI starts close to the answer at each η level
