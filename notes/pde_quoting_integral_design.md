# PDE Quoting Integral — Design Notes

## Computing $Q(y)$

The quoting term in the HJB PDE is

$$Q(y) = \int_{\mathbb{R}^*_+} \sum_{n=1}^{N} \sum_{1 \leq i \neq j \leq d} z \, H^{n,i,j}\!\left(z,\, \frac{\theta(y) - \theta(y + z e^i - z e^j)}{z}\right) \lambda^{n,i,j}(z) \, dz$$

where $d_{ij} = e_i - e_j$ and $H(p, \alpha, \beta) = \sup_{\delta} f(\delta;\, \alpha,\, \beta)\,(\delta - p)$ is the exact Hamiltonian, evaluated via Lambert $W$ (no quadratic Taylor approximation as in the ODE).

In the paper's model, the arrival measure $\Lambda_n(z)$ is discrete with support at $z_k \in \{1, 5, 10, 20, 50\}$ M\$, so the integral reduces to an exact finite sum:

$$Q(y) = \sum_{\substack{(i,j) \\ i \neq j}} \sum_{n} \sum_{k} z_k \, \lambda_{k} \, H\!\left(\frac{\theta(y) - \theta(y + z_k \, d_{ij})}{z_k},\, \alpha_n,\, \beta_n\right)$$

No numerical quadrature is needed.

We require $\Delta y \mid z_k$ (grid spacing divides every trade size) so that the shifted point $y + z_k d_{ij}$ always lands exactly on a grid node. With the paper's sizes, $\Delta y = 1$ M\$ is the natural choice. This means no interpolation of $\theta$ is ever needed.


## Boundary treatment

The grid covers a finite domain $y_i \in [-Y_{\max}, Y_{\max}]$ per axis. For points $y$ near the boundary, some shifts $y + z_k d_{ij}$ fall outside this domain. For those terms we set

$$H\bigl(p_k(y),\, \alpha_n,\, \beta_n\bigr) = 0 \qquad \text{if} \quad y + z_k \, d_{ij} \notin \Omega$$

i.e. the market maker does not quote that (pair, size) at that inventory.

This means $Q(y)$ near the boundary is a partial sum: large-$z$ contributions drop out first (they need larger shifts), small-$z$ contributions survive longer. At the very edge only the $z = 1$ terms may remain, and at the corners $Q = 0$ entirely.

The domain boundary $\partial \Omega$ therefore acts as an implicit inventory limit. As long as $Y_{\max}$ is large enough that the solution in the region of interest is not affected, this is fine. In practice, choosing $Y_{\max}$ well beyond the inventory range where optimal quotes are economically relevant is sufficient.


## Why exact $H$ matters

The ODE approximation expands

$$H(p) \approx \alpha_0 + \alpha_1 \, p + \alpha_2 \, p^2$$

around $p = 0$. This is accurate for small $|p|$, i.e. when inventory is moderate and $\theta$ is well-described by the quadratic ansatz $\theta \approx -y^\top A \, y - y^\top B$.

At large $|y|$, the true $\theta$ deviates from the quadratic form, $p$ becomes large, and the Taylor expansion breaks down. Evaluating the exact $H$ via Lambert $W$ removes this approximation, which is the whole point of solving the PDE instead of the Riccati ODE.
