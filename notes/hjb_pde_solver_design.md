# HJB PDE Solver — Design Notes

## The PDE

We solve PDE (1) from Barzykin, Bergault, Guéant (2023) directly on an inventory grid:

$$
0 = \partial_t \theta + y^\top \mu + y^\top \mathcal{D}(\mu) \nabla_y \theta
  - \tfrac{\gamma}{2} y^\top \Sigma y
  + \tfrac{1}{2} \operatorname{Tr}\!\bigl(\mathcal{D}(y)\Sigma\mathcal{D}(y)\, D^2_{yy}\theta\bigr)
  + Q(y)
  + \mathcal{H}_{\text{hedge}}(y),
$$

with terminal condition $\theta(T, y) = -y^\top \kappa\, y$ and:

- $Q(y)$ = quoting Hamiltonian (exact $H$ via Lambert $W$, not the quadratic Taylor used in the ODE),
- $\mathcal{H}_{\text{hedge}}(y)$ = hedging Hamiltonian via Legendre-Fenchel transform of $L(\xi) = \psi|\xi| + \eta\xi^2$.

The goal is to solve this without the quadratic ansatz $\hat\theta = -y^\top A y - y^\top B - C$ used in the Riccati ODE approximation (Eq. 4-5). This captures the true nonlinear structure of the value function, especially at large inventories where the Taylor expansion of $H$ around $p=0$ breaks down.


## Time stepping: explicit Euler

We march backward from $t = T$ to $t = 0$. Defining $\text{RHS}(\theta)$ as the sum of all spatial terms, the update is:

$$\theta^{m-1} = \theta^m + \Delta t \cdot \text{RHS}(\theta^m)$$

where $m$ indexes time steps from $M$ (terminal) down to $1$.

**Why explicit Euler?** It is the simplest scheme to implement and debug. Each time step is a single evaluation of the RHS — no linear systems, no iteration. This makes it easy to verify each spatial operator independently before worrying about solver convergence. The plan is to upgrade to semi-implicit or fully implicit schemes later if needed.

**Possible future upgrades:**
- Semi-implicit (operator splitting): treat diffusion implicitly, Hamiltonians explicitly. Removes the diffusion CFL constraint.
- Fully implicit + policy iteration (Howard's algorithm): the approach the paper used for their d=2 validation. Most robust but most complex.


## CFL stability constraint

For explicit Euler, the time step must satisfy the CFL condition. The most restrictive term is typically the diffusion:

$$\Delta t < \frac{\Delta y^2}{\sigma^2 \, Y_{\max}^2}$$

This arises because the diffusion coefficient is $\tfrac{1}{2} y_i^2 \sigma_i^2$, which is largest at the grid boundary $y_i = Y_{\max}$.

With the chosen parameters ($\sigma_{\text{EUR}} = 80$ bps $= 0.008$ per $\sqrt{\text{day}}$, $Y_{\max} = 150$ M\$, $\Delta y = 1$ M\$):

$$\Delta t < \frac{1}{0.008^2 \times 150^2} = \frac{1}{1.44} \approx 0.69 \text{ days}$$

This is very generous — the diffusion CFL is not the binding constraint here because inventory volatility in M\$ units is modest ($\sigma \cdot Y_{\max} \approx 1.2$ M\$/day$^{1/2}$).

The Hamiltonian terms may impose tighter (empirical) constraints, but with $\Delta t = 5 \times 10^{-5}$ days and 1000 time steps over $T = 0.05$ days, we have a large safety margin.


### Hedging Hamiltonian stiffness (discovered during implementation)

The initial CFL analysis above only considered the diffusion term. In practice, **the hedging Hamiltonian is the binding stability constraint** and makes explicit Euler unusable with the paper's parameters.

The hedging Hamiltonian is $H^{i,j}(p) = (\max(|p| - \psi, 0))^2 / (4\eta)$. For EURUSD, $\eta = 10^{-5} \text{ bps} \times 10^{-4} = 10^{-9}$ in decimal units. This means $H \sim p^2 / (4 \times 10^{-9}) \approx 2.5 \times 10^8 \, p^2$, which amplifies gradients of $\theta$ by a factor of $\sim 10^8$. The resulting positive feedback loop causes numerical blow-up regardless of how small $\Delta t$ is chosen.

This is not a bug — it reflects genuine stiffness in the PDE. The paper used an implicit scheme ("monotone implicit Euler") precisely because of this.

**Workaround for initial development:** We scale $\eta$ up by a factor (e.g., $\times 1000$) to make explicit Euler stable. This lets us validate the spatial discretization (quoting Hamiltonian, diffusion, drift) independently of the time-stepping scheme. The artificially large $\eta$ reduces the hedging Hamiltonian's gain, making explicit time stepping feasible.

**Future upgrade path:** Replace explicit Euler with a semi-implicit scheme that treats the hedging Hamiltonian implicitly while keeping the other terms explicit. This removes the stiffness constraint entirely and allows use of the paper's original parameters.


## Grid and domain

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| $Y_{\max}$ | 150 M\$ | 50 M\$ buffer beyond the ±100 M\$ region of interest |
| $\Delta y$ | 1 M\$ | Divides all trade sizes $\{1, 5, 10, 20, 50\}$ exactly |
| Grid points/axis | 301 | $[-150, 150]$ with spacing 1 |
| Region of interest | $\pm 100$ M\$ | Matches the paper's Figures 1-2 |
| $T$ | 0.05 days (72 min) | Paper's value; long enough for stationarity |
| $N_t$ | 1000 | $\Delta t = 5 \times 10^{-5}$ days $\approx$ 4.3 seconds |

**Why $\Delta y = 1$?** The paper's trade sizes are $\{1, 5, 10, 20, 50\}$ M\$. The quoting Hamiltonian evaluates $\theta(y + z_k d_{ij})$, which must land on a grid node. With $\Delta y = 1$, all sizes are exact multiples — no interpolation needed.

**Why $Y_{\max} = 150$ not $100$?** At $y = 100$ M\$, the largest trade ($z = 50$) shifts the lookup to $y = 150$, which is still on the grid. Without the buffer, the $z = 50$ contribution would be zeroed out for $|y| > 50$, artificially degrading the solution in the region of interest. The 50 M\$ buffer ensures all trade sizes contribute throughout $[-100, 100]$.


## Boundary treatment

**Quoting Hamiltonian:** No explicit boundary condition needed. When the shifted lookup $y + z_k d_{ij}$ falls outside $\Omega$, that contribution is set to $H = 0$ (the market maker doesn't quote that size at that inventory). This is an implicit inventory limit — see `notes/pde_quoting_integral_design.md` for details.

**Diffusion and drift terms:** Use one-sided finite difference stencils at boundary points (forward at left, backward at right), already implemented in `compute_gradient`, `_second_deriv_diagonal`, and `_second_deriv_cross`.

**No Dirichlet pinning:** All grid points (including boundaries) are evolved with the same explicit Euler update. No values are pinned to the ODE solution or any other prescribed function. This keeps the solver self-contained and avoids coupling to the Riccati approximation.

**Why this works:** If $Y_{\max}$ is large enough, the solution in the interior is insensitive to the boundary values. The 50 M\$ buffer was chosen to ensure this. If boundary artifacts are observed, the first remedy is to increase $Y_{\max}$; the second is to pin boundary values to the ODE ansatz.


## Dimensionality

**Primary target: $d = 2$** (e.g., USD-EUR). This gives a 1D inventory grid (301 points) since $Y_{\text{USD}}$ is determined by $Y_{\text{EUR}}$. This is the validation case — the paper mentions comparing ODE vs PDE for $d = 2$.

**Stretch goal: $d = 3$** (e.g., USD-EUR-GBP). This gives a 2D inventory grid ($301 \times 301 \approx 91$k points). This is the interesting case where cross-pair effects (EURGBP) appear in the PDE but are already captured by the ODE's off-diagonal elements in $A$. The comparison reveals whether the quadratic ansatz adequately captures cross-pair interactions.

The code is written dimension-generically (loops over axes, broadcasted arrays), so extending from $d = 2$ to $d = 3$ should require no changes to the solver — only the grid construction and parameter setup in the notebook.


## Validation strategy

**Step 1: Value function comparison.** Plot $\theta_{\text{PDE}}(0, y)$ against the ODE ansatz $-y^\top A(0)\, y - y^\top B(0) - C(0)$ along inventory slices. They should:
- Agree well for moderate $|y|$ (where the quadratic $\hat{H}$ approximation is accurate),
- Diverge at large $|y|$ (where the exact $H$ differs from the Taylor expansion).

This is the primary sanity check. If they disagree everywhere, there is a bug.

**Step 2: Policy comparison.** Extract optimal quotes $\delta^*$ and hedge rates $\xi^*$ from the PDE solution using Eqs. (2)-(3) with the exact $H$, and compare against the ODE-derived policy from the Riccati solution. This is the thesis-relevant result — it quantifies the approximation error of the ODE approach.
