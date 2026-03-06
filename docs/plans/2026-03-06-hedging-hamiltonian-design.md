# Hedging Hamiltonian — Implementation Design

## Goal

Implement the hedging Hamiltonian term from the HJB PDE (Eq. 1, third line) as a reusable function in `src/pde.py`. Given $\theta$ on a grid, compute:

$$\mathcal{H}_{\text{hedge}}(y) = \sum_{1 \leq i < j \leq d} \mathcal{H}^{i,j}\!\bigl(p^{i,j}(y)\bigr)$$

where $p^{i,j}(y) = \partial_{y^i}\theta - \partial_{y^j}\theta + k^i y^i(1 + \partial_{y^i}\theta) - k^j y^j(1 + \partial_{y^j}\theta)$ and $\mathcal{H}(p) = (\max(|p| - \psi, 0))^2 / (4\eta)$.


## Approach

Follow the same precomputed-spec pattern as `QuotingSpec` / `quoting_hamiltonian_integral`. The gradient of $\theta$ is computed separately and passed in, so it can be reused by other PDE terms (drift, diffusion) later.


## New functions in `src/pde.py`

### `H_execution_cost(p, psi, eta) -> np.ndarray`

Vectorized closed-form evaluation:

$$\mathcal{H}(p) = \frac{(\max(|p| - \psi, 0))^2}{4\eta}$$

Takes scalar or array `p`, returns same shape.

### `HedgingSpec` (frozen dataclass)

Precomputed specification, built once per grid/model. Fields:
- `contributions`: tuple of `(i, j, psi, eta, k_i, k_j)` for each D2D pair with $i < j$
- `d`: number of currencies

### `build_hedging_spec(y_grids, mp) -> HedgingSpec`

Loops over canonical pairs in `mp.pairs` where both currencies are in `mp.currencies`. For each pair $(i, j)$ with $i < j$, extracts `psi`, `eta`, `k_i = mp.k[ccy_i]`, `k_j = mp.k[ccy_j]` and stores them.

### `hedging_hamiltonian(grad_theta, y_grids, spec) -> np.ndarray`

Arguments:
- `grad_theta`: list of $d$ arrays (one per currency axis), each same shape as $\theta$. Entry `grad_theta[i]` is $\partial_{y^i}\theta$.
- `y_grids`: list of $d$ 1D arrays (the grid axes).
- `spec`: a `HedgingSpec`.

For each `(i, j, psi, eta, k_i, k_j)` in `spec.contributions`:
1. Broadcast `y_grids[i]` and `y_grids[j]` to the full grid shape.
2. Assemble $p = \nabla_i\theta - \nabla_j\theta + k^i y^i(1 + \nabla_i\theta) - k^j y^j(1 + \nabla_j\theta)$.
3. Accumulate `H_execution_cost(p, psi, eta)` into the result.

Returns array same shape as $\theta$.

### `compute_gradient(theta, dy_list) -> List[np.ndarray]`

General-purpose gradient on a uniform grid. Returns list of $d$ arrays.
- **Interior**: central differences, second-order: $(\ theta[k+1] - \theta[k-1]) / (2 \Delta y_i)$.
- **Boundaries**: one-sided (forward at left, backward at right), first-order.

This will be reused by drift and diffusion terms later.


## Exports

Add to `src/__init__.py`:
- `H_execution_cost`
- `HedgingSpec`, `build_hedging_spec`, `hedging_hamiltonian`
- `compute_gradient`
