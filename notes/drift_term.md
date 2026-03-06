# Drift Term — Design Notes

## The term in the HJB PDE

The drift term in Eq. 1 of the paper is:

$$y^\top \mu(t) + y^\top \mathcal{D}(\mu(t))\,\nabla_y\theta(t,y)$$

where $\mu = (\mu^1, \ldots, \mu^d)^\top$ is the deterministic drift of the exchange rates and $\mathcal{D}(\mu)$ is the diagonal matrix with $\mu_i$ on the diagonal.


## What it represents

The exchange rates follow $dS^i = \mu^i S^i dt + \sigma^i S^i dW^i + \ldots$, so each inventory $Y^i = q^i S^i$ drifts at rate $\mu^i Y^i$ even without any trading. The drift term captures how this deterministic trend affects the expected evolution of the value function:

- $y^\top \mu$ — direct effect: the portfolio value changes at rate $\sum_i \mu_i y_i$ due to price drift.
- $y^\top \mathcal{D}(\mu)\,\nabla_y\theta$ — indirect effect: the drift moves the inventory state, which changes $\theta$ through its gradient.


## Combining into a single sum

$$y^\top \mu + y^\top \mathcal{D}(\mu)\,\nabla_y\theta \;=\; \sum_{i=1}^d \mu_i\,y_i + \sum_{i=1}^d \mu_i\,y_i\,\frac{\partial\theta}{\partial y_i} \;=\; \sum_{i=1}^d \mu_i\,y_i\!\left(1 + \frac{\partial\theta}{\partial y_i}\right)$$


## Simplification for the paper's example

In the paper's numerical example, $\mu_i = 0$ for all currencies, so the entire drift term vanishes. The implementation still handles the general case for completeness.


## Numerical method

The gradient $\partial\theta/\partial y_i$ is already computed by `compute_gradient` (shared with the hedging Hamiltonian). It uses:

- **Interior**: second-order central differences $\;\frac{\theta_{k+1} - \theta_{k-1}}{2\,\Delta y_i}$
- **Boundaries**: first-order one-sided differences

The $y_i$ factor is broadcast along axis $i$ of the grid, matching the pattern used in `diffusion_term` and `hedging_hamiltonian`.

No new finite-difference stencils are needed — the drift term reuses the existing gradient computation.


## Implementation

```python
def drift_term(grad_theta, y_grids, mu_vec):
    """Compute sum_i mu_i y_i (1 + dtheta/dy_i) at every grid point."""
```

Loops over dimensions $i$, skipping any where $\mu_i = 0$. For each nonzero $\mu_i$, broadcasts $y_i$ along axis $i$ and accumulates $\mu_i \cdot y_i \cdot (1 + \partial\theta/\partial y_i)$.
