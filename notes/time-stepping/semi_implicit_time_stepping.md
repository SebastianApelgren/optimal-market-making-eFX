# Semi-Implicit Time Stepping for the HJB PDE

## Motivation

The explicit Euler scheme for PDE (1) requires us to scale $\eta$ up by $\times 1000$ to avoid blow-up. The root cause is the hedging Hamiltonian $\mathcal{H}^{i,j}(p) = (\max(|p| - \psi, 0))^2 / (4\eta)$, which amplifies gradients of $\theta$ by $\sim 1/(4\eta) \approx 2.5 \times 10^8$ at the paper's $\eta = 10^{-9}$ (decimal). This positive feedback between $\theta$ and its own spatial derivatives makes the system stiff: explicit time stepping is unconditionally unstable (or requires impractically small $\Delta t$) regardless of the CFL condition from diffusion alone.

We want a scheme that:
1. Handles the paper's original $\eta$ values without rescaling.
2. Does not require solving a full nonlinear system at each time step (which a fully implicit scheme would).
3. Remains straightforward to implement.

The standard approach is **operator splitting**: treat the stiff linear part implicitly and the nonlinear Hamiltonian parts explicitly. This is sometimes called an IMEX (implicit-explicit) or semi-implicit scheme.


## Identifying the linear and nonlinear parts

Recall the full PDE (1). We write $\partial_t \theta + \mathcal{L}[\theta] + \mathcal{N}[\theta] = 0$, where we split the right-hand side into:

**Linear operator $\mathcal{L}[\theta]$** (depends linearly on $\theta$ and its derivatives):

$$\mathcal{L}[\theta] \;=\; \frac{1}{2}\operatorname{Tr}\!\Big(\mathcal{D}(y)\,\Sigma\,\mathcal{D}(y)\,D^2_{yy}\theta\Big) \;+\; y^\top \mathcal{D}(\mu)\,\nabla_y \theta$$

This is the diffusion and the $\nabla_y\theta$-dependent part of the drift, with the same signs as in the paper's PDE (1). On a grid, $\mathcal{L}$ becomes a sparse matrix acting on the vector of $\theta$ values.

**Nonlinear operator $\mathcal{N}[\theta]$** (nonlinear dependence on $\theta$):

$$\mathcal{N}[\theta] \;=\; y^\top \mu \;-\; \frac{\gamma}{2}\,y^\top \Sigma\,y \;+\; Q[\theta] \;+\; \mathcal{H}_{\text{hedge}}[\theta]$$

where:
- $y^\top \mu$ and $-\frac{\gamma}{2} y^\top \Sigma y$ are source terms (no $\theta$ dependence at all),
- $Q[\theta]$ is the quoting Hamiltonian (nonlinear: involves $H^{n,i,j}$ evaluated at finite differences of $\theta$),
- $\mathcal{H}_{\text{hedge}}[\theta]$ is the hedging Hamiltonian (nonlinear: involves $\mathcal{H}^{i,j}$ applied to $\nabla_y \theta$).

Note: the drift term $y^\top \mathcal{D}(\mu) \nabla_y \theta = \sum_i \mu_i y_i \partial_{y_i}\theta$ is linear in $\theta$ and goes into $\mathcal{L}$. The remaining part of the drift, $y^\top \mu = \sum_i \mu_i y_i$, is a pure source term in $\mathcal{N}$.


## The semi-implicit (IMEX Euler) scheme

We march backward from $t = T$ to $t = 0$, with time step $\Delta t = T / N_t$. The PDE in backward time $\tau = T - t$ reads:

$$\partial_\tau \theta \;=\; \mathcal{L}[\theta] + \mathcal{N}[\theta].$$

The IMEX Euler update from $\theta^m$ (at backward time $m \Delta t$) to $\theta^{m+1}$ is:

$$\frac{\theta^{m+1} - \theta^m}{\Delta t} \;=\; \mathcal{L}[\theta^{m+1}] \;+\; \mathcal{N}[\theta^m].$$

That is:
- The **linear part $\mathcal{L}$** is evaluated at the **new** (unknown) time level — **implicit**.
- The **nonlinear part $\mathcal{N}$** is evaluated at the **old** (known) time level — **explicit**.

Rearranging:

$$\boxed{\bigl(I - \Delta t\,\mathcal{L}\bigr)\,\theta^{m+1} \;=\; \theta^m \;+\; \Delta t\,\mathcal{N}[\theta^m]}$$

The right-hand side is known (it only involves $\theta^m$). The left-hand side is a **linear system** $A\,\theta^{m+1} = b$ where $A = I - \Delta t\,\mathcal{L}$ is a sparse matrix that depends only on the grid geometry and model parameters — it can be **precomputed once** before the time loop.


## Why this was expected to cure the stiffness (but didn't)

The reasoning below motivated the IMEX approach. It turned out to be **wrong in practice** — see `semi_implicit_results_and_next_steps.md` for the actual results. The argument is preserved here to document the thinking and where it broke down.

The idea was: implicit diffusion provides unconditional damping of high-frequency modes, which should counteract the amplification from the explicit Hamiltonian evaluation, stabilising the scheme for moderate $\Delta t$.

More precisely: the discrete diffusion operator $\mathcal{L}$ is built from second derivatives ($\partial_{yy}$), whose eigenvalues are negative: $\lambda_k < 0$. In an explicit forward-$\tau$ scheme the amplification factor is $|1 + \Delta t \lambda_k|$, which exceeds 1 when $\Delta t |\lambda_k| > 2$ — hence the CFL constraint. In the IMEX scheme, the implicit treatment gives amplification $|1/(1 - \Delta t \lambda_k)| = 1/(1 + \Delta t |\lambda_k|) < 1$ for all $\lambda_k < 0$ and all $\Delta t > 0$ — unconditional stability for the linear part. The original hope was that the explicit Hamiltonian contributions would be "bounded" in the sense of not having eigenvalues that grow with grid refinement.

**Why it failed:** The hedging Hamiltonian $\mathcal{H}(p) = (\max(|p| - \psi, 0))^2 / (4\eta)$ is *not* bounded as a functional of $\theta$. The argument $p$ depends on $\nabla\theta$, so $\mathcal{H} \sim |\nabla\theta|^2 / (4\eta)$. With $1/(4\eta) \approx 2.5 \times 10^8$, even small gradients in $\theta$ produce $O(10^8)$ contributions that feed back into $\theta$ on the next step. This is a nonlinear positive feedback loop — not a fixed-eigenvalue problem that implicit diffusion can damp. The diffusion CFL ($\Delta t < 0.69$ days) was never the binding constraint; the hedging stiffness dominates by many orders of magnitude.

The fix requires treating the hedging Hamiltonian implicitly, either via linearisation or full policy iteration — see `implicit_euler_policy_iteration.md`.


## Structure of the linear system

### For $d = 2$ (USD, EUR): effectively 1D

Since $\sigma_{\text{USD}} = 0$, the covariance matrix has $\Sigma_{1,j} = \Sigma_{j,1} = 0$. The diffusion operator reduces to a single second derivative along the EUR axis:

$$\mathcal{L}[\theta] = \frac{1}{2}\sigma_{\text{EUR}}^2\,y_{\text{EUR}}^2\,\frac{\partial^2 \theta}{\partial y_{\text{EUR}}^2} \;+\; \mu_{\text{EUR}}\,y_{\text{EUR}}\,\frac{\partial \theta}{\partial y_{\text{EUR}}}$$

(The USD drift vanishes too since $\mu_{\text{USD}} = 0$.)

On the grid, fix any value of $y_{\text{USD}}$. For each such row, $\mathcal{L}$ acts only along the EUR axis via a **tridiagonal** finite difference stencil. Therefore $A = I - \Delta t\,\mathcal{L}$ is tridiagonal along the EUR axis for each fixed $y_{\text{USD}}$. This can be solved in $O(N)$ per row via the Thomas algorithm. Sweeping over all $y_{\text{USD}}$ rows gives $O(N^2)$ total work per time step — the same as the explicit scheme.

### For $d = 3$ (USD, EUR, GBP): ADI splitting

With two non-reference currencies, the diffusion operator has diagonal terms in each axis plus a cross derivative:

$$\mathcal{L} = \frac{1}{2}\sigma_{\text{EUR}}^2 y_{\text{EUR}}^2 \partial_{y_{\text{EUR}}}^2 + \frac{1}{2}\sigma_{\text{GBP}}^2 y_{\text{GBP}}^2 \partial_{y_{\text{GBP}}}^2 + \rho_{\text{EUR,GBP}}\,\sigma_{\text{EUR}}\,\sigma_{\text{GBP}}\,y_{\text{EUR}}\,y_{\text{GBP}}\,\partial^2_{y_{\text{EUR}}\,y_{\text{GBP}}} + \text{drift terms}$$

Inverting $I - \Delta t\,\mathcal{L}$ directly would require solving a banded (but not tridiagonal) 2D system. Instead, we use **alternating direction implicit (ADI)** splitting — e.g. the Douglas-Rachford or Craig-Sneyd scheme — which factors each time step into a sequence of tridiagonal solves along each axis. The cross-derivative term is treated explicitly within the ADI substeps.

This is a well-studied approach for multi-dimensional diffusion PDEs in finance (Black-Scholes in multiple assets uses the same technique).


## Summary of the update algorithm

Each time step from $\theta^m$ to $\theta^{m+1}$:

1. **Evaluate the nonlinear terms at $\theta^m$** (explicit, same as current code):
   - Compute $\nabla_y \theta^m$ via finite differences.
   - Evaluate the quoting Hamiltonian $Q[\theta^m]$.
   - Evaluate the hedging Hamiltonian $\mathcal{H}_{\text{hedge}}[\theta^m]$.
   - Evaluate the source terms ($+y^\top\mu$, $-\frac{\gamma}{2}y^\top\Sigma y$).
   - Form the explicit right-hand side: $b = \theta^m + \Delta t \cdot \mathcal{N}[\theta^m]$.

2. **Solve the implicit diffusion step** (new, replaces the explicit Euler update):
   - Solve $(I - \Delta t\,\mathcal{L})\,\theta^{m+1} = b$.
   - For $d = 2$: tridiagonal solve along EUR axis for each $y_{\text{USD}}$ row.
   - For $d = 3$: ADI splitting into tridiagonal solves along each axis.

The computational cost per time step is essentially the same as explicit Euler (the tridiagonal solves are $O(N)$ per row), but the scheme is stable for much larger $\Delta t$.


## What changes in the code

- `pde_rhs` is split into two parts: `pde_rhs_nonlinear` (quoting, hedging, source terms) and a representation of the linear diffusion+drift operator.
- The diffusion+drift operator is precomputed as tridiagonal coefficient arrays (sub-diagonal, diagonal, super-diagonal) for each row of the grid.
- `solve_hjb_explicit` is replaced by `solve_hjb_semi_implicit` which, at each time step, evaluates the explicit RHS, then solves the tridiagonal system.
- For $d = 2$, this uses `numpy` or `scipy.linalg.solve_banded` (Thomas algorithm). No new dependencies beyond what is already available.


## Relation to the paper's approach

The paper (p.6) states they used a "monotone implicit Euler scheme" for the $d = 2$ PDE validation. That is a **fully implicit** scheme where both the linear and nonlinear parts are evaluated at $\theta^{m+1}$, requiring an iterative solve (policy iteration / Howard's algorithm) at each time step. Our semi-implicit scheme is a lighter-weight alternative that should suffice for moderate $\eta$ values. If it does not fully resolve the stiffness, the path forward is clear: add policy iteration on top of the implicit diffusion step, effectively upgrading to the paper's fully implicit scheme.
