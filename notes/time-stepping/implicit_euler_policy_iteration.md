# Implicit Euler with Policy Iteration (Howard's Algorithm)

## Motivation

The semi-implicit (IMEX) Euler scheme treats diffusion implicitly and the Hamiltonians explicitly. As documented in `notes/semi_implicit_results_and_next_steps.md`, this does not cure the stiffness because the binding constraint is the hedging Hamiltonian's nonlinear gain $1/(4\eta) \approx 2.5 \times 10^8$, not the diffusion. The IMEX scheme requires the same $\eta \times 1000$ scaling as fully explicit Euler.

The paper (p.6) states they used a "monotone implicit Euler scheme" for the $d = 2$ PDE validation. This is a **fully implicit** scheme where the nonlinear HJB equation at each time step is solved via policy iteration (Howard's algorithm). Policy iteration linearizes the HJB by alternating between fixing the control policy and solving the resulting linear PDE. It converges quadratically (equivalent to Newton's method) and typically requires only 2-5 inner iterations per time step.


## The fully implicit Euler scheme

The PDE in backward time $\tau = T - t$:

$$\partial_\tau \theta = F[\theta]$$

where $F[\theta]$ is the full spatial RHS (diffusion + drift + penalty + quoting + hedging). The fully implicit Euler update is:

$$\frac{\theta^{m+1} - \theta^m}{\Delta t} = F[\theta^{m+1}]$$

Rearranged:

$$\theta^{m+1} - \Delta t \cdot F[\theta^{m+1}] = \theta^m$$

This is a **nonlinear** equation for $\theta^{m+1}$ because $F$ contains the quoting and hedging Hamiltonians, which involve suprema over controls.


## Policy iteration

At each time step, solve the nonlinear implicit equation by iterating:

### Step A: Fix the controls from $\theta^k$

Given current estimate $\theta^k$ (initialized to $\theta^m$ at the start of each time step):

**Quoting controls.** For each directed pair $(i, j)$, tier $n$, and size $z_k$, compute the momentum:

$$p_k(y) = \frac{\theta^k(y) - \theta^k(y + z_k d_{ij})}{z_k}$$

and the optimal markup via Lambert $W$:

$$\delta^*_k(y) = p_k(y) + \frac{W_0(e^{-(1 + \alpha_n + \beta_n p_k)}) + 1}{\beta_n}$$

Also record the fill probability at the optimum:

$$f^*_k(y) = f(\delta^*_k;\, \alpha_n, \beta_n) = \frac{w}{w + 1}, \qquad w = W_0(e^{-(1 + \alpha_n + \beta_n p_k)})$$

These are the same formulas in `hamiltonian.py:H_logistic` and `hamiltonian.py:optimal_delta_logistic`.

**Hedging controls.** For each canonical pair $(i, j)$ with $i < j$, compute the hedging momentum:

$$p^{i,j}(y) = \nabla\theta^k_i - \nabla\theta^k_j + k_i y_i (1 + \nabla\theta^k_i) - k_j y_j (1 + \nabla\theta^k_j)$$

and the optimal hedge rate:

$$\xi^{*,i,j}(y) = \operatorname{sign}(p) \cdot \frac{\max(|p| - \psi, 0)}{2\eta}$$

This is the derivative $H'(p)$ of the execution cost Hamiltonian $H(p) = (\max(|p| - \psi, 0))^2 / (4\eta)$.

At grid points where $|p| \leq \psi$ (inside the dead zone), $\xi^* = 0$ and the hedging term contributes nothing to the linearized operator — only to the source.


### Step B: Solve the linearized system

With $\delta^*$ and $\xi^*$ frozen, the HJB becomes linear in $\theta^{k+1}$. We need to identify what goes into the sparse matrix $A$ (linear dependence on $\theta^{k+1}$) and what goes into the source vector $b$ (known terms).

Write the implicit equation as:

$$\theta^{k+1} - \Delta t \cdot F_{\text{linear}}[\theta^{k+1}] = \theta^m + \Delta t \cdot \text{source}$$

i.e. $(I - \Delta t \cdot A)\,\theta^{k+1} = \theta^m + \Delta t \cdot s$, where $A$ is the linearized spatial operator and $s$ collects all terms that don't depend on $\theta^{k+1}$.


#### Term-by-term decomposition

**1. Running penalty:** $-\frac{\gamma}{2} y^\top \Sigma y$

No $\theta$ dependence. Goes entirely into the source:

$$s_{\text{penalty}}(y) = -\frac{\gamma}{2} y^\top \Sigma y$$

This is precomputed once (same as `running_penalty` in `pde.py`).


**2. Diffusion:** $\frac{1}{2} \operatorname{Tr}(\mathcal{D}(y) \Sigma \mathcal{D}(y)\, D^2_{yy} \theta)$

Entirely linear in $\theta$. Goes into $A$:

$$A_{\text{diff}}[\theta](y) = \sum_{i} \frac{1}{2} \Sigma_{ii} y_i^2 \frac{\partial^2 \theta}{\partial y_i^2} + \sum_{i < j} \Sigma_{ij} y_i y_j \frac{\partial^2 \theta}{\partial y_i \partial y_j}$$

For the diagonal terms, the FD stencil at interior point $\mathbf{j}$ along axis $i$ is:

$$\frac{\partial^2 \theta}{\partial y_i^2}\bigg|_{\mathbf{j}} \approx \frac{\theta_{\mathbf{j}+e_i} - 2\theta_{\mathbf{j}} + \theta_{\mathbf{j}-e_i}}{\Delta y_i^2}$$

This contributes three entries per row of $A$ (for each axis with $\Sigma_{ii} > 0$).

For the cross terms ($i \neq j$), the central FD stencil is:

$$\frac{\partial^2 \theta}{\partial y_i \partial y_j}\bigg|_{\mathbf{j}} \approx \frac{\theta_{\mathbf{j}+e_i+e_j} - \theta_{\mathbf{j}+e_i-e_j} - \theta_{\mathbf{j}-e_i+e_j} + \theta_{\mathbf{j}-e_i-e_j}}{4 \Delta y_i \Delta y_j}$$

This contributes four entries per row. For $d = 2$ with $\sigma_{\text{USD}} = 0$, the cross term vanishes and only the EUR diagonal term survives.


**3. Drift:** $y^\top \mu + y^\top \mathcal{D}(\mu) \nabla_y \theta = \sum_i \mu_i y_i (1 + \partial\theta/\partial y_i)$

Split into source and linear parts:

- Source: $s_{\text{drift}}(y) = \sum_i \mu_i y_i$
- Linear: $A_{\text{drift}}[\theta](y) = \sum_i \mu_i y_i \frac{\partial \theta}{\partial y_i}$

The gradient FD stencil (central differences at interior, one-sided at boundaries) contributes two entries per row per axis with $\mu_i \neq 0$.

For the paper's example, $\mu = 0$ everywhere, so the drift term is identically zero (both source and linear parts). The code should still handle nonzero $\mu$ for generality.


**4. Quoting Hamiltonian with fixed $\delta^*$:**

The original quoting term for a single contribution $(n, i, j, k)$ is:

$$z_k \lambda_k H^{n}(p_k) = z_k \lambda_k \sup_\delta f(\delta)(\delta - p_k)$$

With $\delta^*$ fixed, the supremum is removed and the contribution becomes:

$$z_k \lambda_k f^*_k \cdot (\delta^*_k - p_k)$$

where $p_k(y) = (\theta(y) - \theta(y + z_k d_{ij})) / z_k$. Expanding:

$$= z_k \lambda_k f^*_k \cdot \delta^*_k - \lambda_k f^*_k \cdot (\theta(y) - \theta(y + z_k d_{ij}))$$

The first term is a source (depends on $\theta^k$ through $\delta^*_k$ and $f^*_k$, but not on $\theta^{k+1}$):

$$s_{\text{quote},c}(y) = z_k \lambda_k f^*_k(y) \cdot \delta^*_k(y)$$

The second term is linear in $\theta^{k+1}$:

$$A_{\text{quote},c}[\theta](y) = -\lambda_k f^*_k(y) \cdot \bigl(\theta(y) - \theta(y + z_k d_{ij})\bigr)$$

In the sparse matrix, this contributes **two entries per row** for each contribution $c$: a diagonal entry $-\lambda_k f^*_k(y)$ at position $(y, y)$, and an off-diagonal entry $+\lambda_k f^*_k(y)$ at position $(y, y + z_k d_{ij})$.

**Boundary handling:** When $y + z_k d_{ij}$ falls outside the grid, the contribution is dropped entirely (no matrix entry, no source). This matches the existing explicit treatment.

**Note on the number of contributions:** For the $d = 2$ case (EURUSD only), there are 2 directions × 5 sizes × 2 tiers = 20 quoting contributions. Each adds 2 matrix entries per row, so the quoting operator adds up to 40 nonzeros per row (some may coincide if different contributions have the same shift).


**5. Hedging Hamiltonian with fixed $\xi^*$:**

The hedging Hamiltonian for canonical pair $(i, j)$ is:

$$H^{i,j}(p) = \sup_\xi \bigl[\xi \cdot p - \psi|\xi| - \eta\xi^2\bigr]$$

With $\xi^*$ fixed:

$$\xi^* \cdot p^{k+1} - \psi|\xi^*| - \eta(\xi^*)^2$$

where $p^{k+1} = \nabla\theta^{k+1}_i - \nabla\theta^{k+1}_j + k_i y_i (1 + \nabla\theta^{k+1}_i) - k_j y_j (1 + \nabla\theta^{k+1}_j)$.

The source part (known from Step A):

$$s_{\text{hedge}}^{i,j}(y) = -\psi|\xi^*| - \eta(\xi^*)^2 + \xi^* \cdot (k_i y_i - k_j y_j)$$

The linear part:

$$A_{\text{hedge}}^{i,j}[\theta](y) = \xi^*(y) \cdot \Bigl[(1 + k_i y_i) \frac{\partial\theta}{\partial y_i} - (1 + k_j y_j) \frac{\partial\theta}{\partial y_j}\Bigr]$$

Using central FD for the gradients, this contributes up to 4 entries per row per hedging pair (two per axis for $\partial\theta/\partial y_i$ and $\partial\theta/\partial y_j$). For $d = 2$ with $k_{\text{USD}} = 0$ and $k_{\text{EUR}} \approx 5 \times 10^{-7}$, the market impact terms are negligible but included for correctness.

**Dead zone:** Where $|\xi^*| = 0$ (i.e. $|p| \leq \psi$), the hedging contribution to both $A$ and $s$ is zero. No matrix entries needed at those grid points.


### Convergence criterion

Policy iteration terminates when:

$$\frac{\|\theta^{k+1} - \theta^k\|_\infty}{\max(1, \|\theta^{k+1}\|_\infty)} < \varepsilon_{\text{PI}}$$

with $\varepsilon_{\text{PI}} = 10^{-10}$ (well below the FD discretization error). A maximum of 20 inner iterations guards against non-convergence; in practice 2-5 are expected.


## Sparse matrix assembly

### Flattening convention

The $d$-dimensional grid array $\theta$ with shape $(n_0, n_1, \ldots, n_{d-1})$ is flattened to a 1D vector of length $N = \prod_k n_k$ using C-order (row-major) flattening, matching `numpy.ravel(order='C')`.

A multi-index $\mathbf{j} = (j_0, j_1, \ldots, j_{d-1})$ maps to flat index:

$$I(\mathbf{j}) = \sum_{k=0}^{d-1} j_k \cdot \text{stride}_k, \qquad \text{stride}_k = \prod_{l=k+1}^{d-1} n_l$$

Shifting $\mathbf{j}$ by $+1$ along axis $k$ adds $\text{stride}_k$ to the flat index. Shifting by $s$ steps along axis $k$ adds $s \cdot \text{stride}_k$.


### Assembly strategy

Build the matrix in COO format (lists of row, col, val), then convert to CSR for the solve:

```python
rows, cols, vals = [], [], []
source = np.zeros(N)

# 1. Diffusion: loop over axes with Sigma_kk > 0
#    For interior points along axis k:
#      row = I(j), col = I(j ± e_k), val = ±dt * 0.5 * Sigma_kk * y_k^2 / dy_k^2
#      row = I(j), col = I(j),        val = -dt * Sigma_kk * y_k^2 / dy_k^2  (diagonal)

# 2. Quoting: loop over contributions (n, i, j, k)
#    For each grid point y in the overlap region:
#      row = flat(y), col = flat(y),             val = +dt * lam * f_star
#      row = flat(y), col = flat(y + z*d_{ij}),  val = -dt * lam * f_star
#    source[flat(y)] += dt * z * lam * f_star * delta_star

# 3. Hedging: loop over canonical pairs
#    For each grid point y where |xi_star| > 0:
#      Gradient of theta along axis i: contributes (j-1, j+1) entries
#      row = flat(y), col = flat(y ± e_i), val = ±dt * xi_star * (1 + k_i*y_i) / (2*dy_i)
#      Similarly for axis j with opposite sign

# 4. Drift: same structure as hedging gradient terms

# Add identity diagonal
rows += range(N)
cols += range(N)
vals += [1.0] * N

A = scipy.sparse.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()
b = theta_old_flat + source
theta_new_flat = scipy.sparse.linalg.spsolve(A, b)
```

### Expected sparsity

For $d = 2$ with 301×301 grid ($N = 90{,}601$):

| Operator | Entries/row | Total nnz |
|----------|------------|-----------|
| Identity diagonal | 1 | 90k |
| Diffusion (EUR axis only) | 3 | 270k |
| Quoting (20 contributions × 2) | ≤ 40 | ≤ 3.6M |
| Hedging (1 pair, 4 gradient entries) | ≤ 4 | ≤ 360k |
| **Total** | **≤ 48** | **≤ 4.3M** |

With 8 bytes per entry (CSR float64), the matrix data occupies ~35 MB. The direct solver (SuperLU via `spsolve`) should handle this in ~0.05-0.2 seconds per solve.


## Boundary treatment

Boundary points ($j_k = 0$ or $j_k = n_k - 1$ for any axis $k$) are handled as follows:

**Option A (identity rows):** Set the matrix row to the identity ($A_{ii} = 1$, all other entries zero) and set $b_i = \theta^m_i + \Delta t \cdot s_i$ where $s_i$ includes only the source terms (penalty). This effectively applies explicit Euler at boundaries, which is acceptable since boundaries lie in the 50 M$ buffer zone.

**Option B (one-sided stencils):** Use forward/backward FD stencils at boundaries to assemble the same operators. This is more accurate but complicates the stencil logic.

**Recommendation:** Start with Option A (simpler, boundaries are in the buffer). If boundary artifacts appear, upgrade to Option B.


## Algorithm summary

```
Precompute (once):
    spec = build_pde_spec(y_grids, mp)
    strides = compute_strides(grid_shape)
    penalty_source = running_penalty(...)

Time loop (m = 0, 1, ..., n_steps - 1):
    theta_old = theta.copy()

    Policy iteration (k = 0, 1, ...):
        # Step A: extract controls from theta
        grad = compute_gradient(theta, dy_list)
        delta_star, f_star = extract_quoting_controls(theta, spec)
        xi_star = extract_hedging_controls(grad, y_grids, spec)

        # Step B: assemble and solve linear system
        A, source = assemble_linearized_system(
            delta_star, f_star, xi_star, spec, dt, strides)
        b = theta_old.ravel() + dt * source
        theta_new = spsolve(A, b).reshape(grid_shape)

        # Check convergence
        if max|theta_new - theta| / max(1, max|theta_new|) < eps:
            theta = theta_new
            break
        theta = theta_new
```


## What changes in the code

### New dependency

Add `scipy` to `requirements.txt`. Specifically: `scipy.sparse` (CSR matrix construction) and `scipy.sparse.linalg` (spsolve).

### New functions in `src/pde.py`

- **`extract_quoting_controls(theta, spec)`** — compute $(\delta^*, f^*)$ arrays for each quoting contribution. Returns a list of $(delta\_star, f\_star)$ arrays, one per contribution, defined on the overlap sub-grid.

- **`extract_hedging_controls(grad_theta, y_grids, spec)`** — compute $\xi^*$ for each canonical pair. Returns a list of $\xi^*$ arrays on the full grid.

- **`assemble_implicit_system(delta_star_list, f_star_list, xi_star_list, spec, dt, strides)`** — build the sparse CSR matrix $(I - \Delta t \cdot A)$ and source vector $s$. Returns `(A_csr, source_flat)`.

- **`solve_hjb_implicit(y_grids, mp, n_steps, snapshot_times)`** — main solver. Same interface as `solve_hjb_explicit` and `solve_hjb_semi_implicit`.

### Reused code

- `build_pde_spec`, `QuotingSpec`, `HedgingSpec` — unchanged
- `compute_gradient` — used in hedging control extraction
- `running_penalty`, `terminal_condition` — unchanged
- `H_logistic`, `optimal_delta_logistic` from `hamiltonian.py` — used in quoting control extraction

### Existing code kept but not used by the implicit solver

- `pde_rhs`, `pde_rhs_nonlinear` — still available for explicit/IMEX solvers
- Thomas solver (`_thomas_factor`, `_thomas_solve_batch`) — kept for potential use as preconditioner in d≥3
- `_build_implicit_ops`, `_apply_implicit_step` — IMEX infrastructure, kept for reference


## Extension to $d = 3$

For $d = 3$ with 301 points per axis, the system has $\sim 27$M unknowns. A direct sparse solve is infeasible.

**Approach:** Replace `spsolve` with an iterative solver (`scipy.sparse.linalg.gmres` or `bicgstab`) using the existing Thomas/ADI infrastructure as a preconditioner. The preconditioner applies the tridiagonal diffusion solve per axis (already implemented in `_apply_implicit_step`), which captures the dominant stiffness at low cost. The iterative solver then converges in a small number of iterations.

Alternatively, use operator splitting within the policy iteration: instead of assembling the full sparse matrix, split the implicit solve into per-axis tridiagonal sweeps. This avoids forming the sparse matrix entirely and reuses the existing Thomas solver. The splitting error is $O(\Delta t^2)$, acceptable for a first-order scheme.

For now, target $d = 2$ only. The $d = 3$ extension is a separate design step.


## Validation plan

1. **η × 1000 comparison:** Run the implicit solver with scaled η and compare against explicit/IMEX solutions. They should agree to within discretization error ($\sim 10^{-6}$), confirming the linearization and assembly are correct.

2. **Original η:** Run the implicit solver with the paper's original η. This is the main test — it should produce a finite, smooth solution where explicit/IMEX blow up.

3. **ODE comparison:** Compare $\theta_{\text{PDE}}(0, y)$ against the Riccati ODE ansatz along inventory slices. Should agree for moderate $|y|$, diverge at large $|y|$.
