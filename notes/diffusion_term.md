# Diffusion Term — Design Notes

## The term in the HJB PDE

The second-order (diffusion) term in Eq. 1 of the paper is:

$$\frac{1}{2} \operatorname{Tr}\!\Big(\mathcal{D}(y)\,\Sigma\,\mathcal{D}(y)\,D^2_{yy}\theta(t,y)\Big)$$

where $\mathcal{D}(y)$ is the diagonal matrix with $y_i$ on the diagonal, $\Sigma = (\rho^{i,j}\sigma^i\sigma^j)_{1 \leq i,j \leq d}$ is the covariance matrix, and $D^2_{yy}\theta$ is the Hessian of the value function.


## What it represents

The inventories $Y^i = q^i S^i$ are measured in reference currency. Since the exchange rates follow geometric Brownian motion ($dS^i = \mu^i S^i dt + \sigma^i S^i dW^i + \ldots$), the volatility of inventory $Y^i$ is proportional to $Y^i$ itself. The diffusion term captures how this stochastic volatility affects the expected evolution of the value function $\theta$ via Ito's formula.

The $y_i y_j$ prefactor means:
- **Zero inventory = zero diffusion**: if the market maker holds none of currency $i$, there is no volatility risk from $S^i$ fluctuations. This is the classic GBM feature.
- **Large inventory = large diffusion**: the further the market maker is from flat, the more the value function is affected by price uncertainty.


## Expanding the trace

Both $\Sigma$ and $D^2_{yy}\theta$ are symmetric, so:

$$\frac{1}{2} \sum_{i=1}^d \sum_{j=1}^d \Sigma_{ij}\,y_i\,y_j\,\frac{\partial^2\theta}{\partial y_i\,\partial y_j} \;=\; \frac{1}{2} \sum_{i=1}^d \Sigma_{ii}\,y_i^2\,\frac{\partial^2\theta}{\partial y_i^2} \;+\; \sum_{i < j} \Sigma_{ij}\,y_i\,y_j\,\frac{\partial^2\theta}{\partial y_i\,\partial y_j}$$


## Simplification from $\sigma_{\text{ref}} = 0$

The reference currency (USD) has $\sigma^1 = 0$, so $\Sigma_{1,j} = \Sigma_{i,1} = 0$ for all $i, j$. The entire first row and column of $\Sigma$ vanish. This means:

- No second derivatives involving the USD axis are ever needed.
- The number of nonzero terms depends only on the non-reference currencies and their correlations.

### Concrete cases

**$d = 2$ (USD, EUR):** Only one term survives — a single diagonal second derivative:

$$\frac{1}{2}\,\sigma_{\text{EUR}}^2\,y_{\text{EUR}}^2\,\frac{\partial^2\theta}{\partial y_{\text{EUR}}^2}$$

**$d = 3$ (USD, EUR, GBP):** Two diagonal terms plus one cross term:

$$\frac{1}{2}\,\sigma_{\text{EUR}}^2\,y_{\text{EUR}}^2\,\frac{\partial^2\theta}{\partial y_{\text{EUR}}^2} \;+\; \frac{1}{2}\,\sigma_{\text{GBP}}^2\,y_{\text{GBP}}^2\,\frac{\partial^2\theta}{\partial y_{\text{GBP}}^2} \;+\; \rho_{\text{EUR,GBP}}\,\sigma_{\text{EUR}}\,\sigma_{\text{GBP}}\,y_{\text{EUR}}\,y_{\text{GBP}}\,\frac{\partial^2\theta}{\partial y_{\text{EUR}}\,\partial y_{\text{GBP}}}$$


## Computing second derivatives on the grid

### Diagonal: $\partial^2\theta / \partial y_i^2$

Standard second-order central difference (3-point stencil):

$$\frac{\partial^2\theta}{\partial y_i^2}\bigg|_k \;\approx\; \frac{\theta_{...,k+1,...} - 2\,\theta_{...,k,...} + \theta_{...,k-1,...}}{\Delta y_i^2}$$

At boundaries, one-sided second-order stencils:

- Left ($k = 0$): $\;\frac{\theta_0 - 2\theta_1 + \theta_2}{\Delta y_i^2}$
- Right ($k = n-1$): $\;\frac{\theta_{n-3} - 2\theta_{n-2} + \theta_{n-1}}{\Delta y_i^2}$

These are the same 3-point formula shifted to use only interior-side points, preserving second-order accuracy.

### Cross-derivative: $\partial^2\theta / (\partial y_i\,\partial y_j)$

Standard 4-point central stencil:

$$\frac{\partial^2\theta}{\partial y_i\,\partial y_j}\bigg|_{k_i, k_j} \;\approx\; \frac{\theta_{k_i+1,k_j+1} - \theta_{k_i+1,k_j-1} - \theta_{k_i-1,k_j+1} + \theta_{k_i-1,k_j-1}}{4\,\Delta y_i\,\Delta y_j}$$

At boundaries where any of the four stencil points falls outside the grid, we use the same approach as for the gradient: shift the stencil inward (one-sided differences). For instance, at the left boundary of axis $i$ ($k_i = 0$), replace $k_i - 1$ with $k_i$ and $k_i + 1$ with $k_i + 1$, giving a forward-difference approximation in that direction. This is first-order in $\Delta y_i$ at the boundary, but we rely on $Y_{\max}$ being large enough that boundary quality is irrelevant.


## The running penalty $-\frac{\gamma}{2}\,y^\top \Sigma\,y$

This term also appears in the HJB PDE but has no $\theta$ dependence — it is a known function of $y$ on the grid:

$$-\frac{\gamma}{2} \sum_{i,j} \Sigma_{ij}\,y_i\,y_j$$

It uses the same $\Sigma_{ij}\,y_i\,y_j$ products as the diffusion term (without the Hessian factor). We implement it as a separate helper since it needs to be computed only once per grid (it does not change with $\theta$ during time-stepping).


## Implementation plan

The diffusion term is simpler than the quoting and hedging terms — no pair-specific market microstructure data, just $\Sigma$ and the grid. No need for a precomputed `Spec` dataclass.

```python
def diffusion_term(theta, y_grids, dy_list, Sigma):
    """Compute ½ Tr(D(y) Sigma D(y) D²_yy theta) at every grid point."""
```

1. Loop over unique $(i, j)$ with $i \leq j$ where $\Sigma_{ij} \neq 0$.
2. Compute the second derivative $\partial^2\theta / (\partial y_i\,\partial y_j)$ via finite differences.
3. Broadcast $y_i$ and $y_j$ along their respective axes.
4. Multiply: $\Sigma_{ij} \cdot y_i \cdot y_j \cdot \partial^2\theta/(\partial y_i\,\partial y_j)$.
5. Accumulate with the correct prefactor ($\frac{1}{2}$ for diagonal, $1$ for off-diagonal).

```python
def running_penalty(y_grids, Sigma, gamma):
    """Compute -gamma/2 * y^T Sigma y at every grid point."""
```

Same loop structure but without the Hessian — just $\Sigma_{ij} \cdot y_i \cdot y_j$ products.
