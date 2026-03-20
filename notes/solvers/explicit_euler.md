# Explicit Euler Solver

The simplest time-stepping scheme — used for initial development and validation of the spatial operators.


## Time stepping

We march forward in backward time $\tau = T - t$, from $\tau = 0$ (terminal) to $\tau = T$ (initial). Defining $\text{RHS}(\theta)$ as the sum of all spatial terms, the PDE becomes $\partial_\tau \theta = \text{RHS}(\theta)$, and the explicit Euler update is:

$$\theta^{m+1} = \theta^m + \Delta t \cdot \text{RHS}(\theta^m)$$

where $m = 0, 1, \ldots, N_t - 1$ indexes time steps forward in $\tau$.

**Why explicit Euler?** It is the simplest scheme to implement and debug. Each time step is a single evaluation of the RHS — no linear systems, no iteration. This makes it easy to verify each spatial operator independently before worrying about solver convergence.


## CFL stability constraint

For explicit Euler, the time step must satisfy the CFL condition. The most restrictive term is typically the diffusion:

$$\Delta t < \frac{\Delta y^2}{\sigma^2 \, Y_{\max}^2}$$

This arises because the diffusion coefficient is $\tfrac{1}{2} y_i^2 \sigma_i^2$, which is largest at the grid boundary $y_i = Y_{\max}$.

With the chosen parameters ($\sigma_{\text{EUR}} = 80$ bps $= 0.008$ per $\sqrt{\text{day}}$, $Y_{\max} = 150$ M\$, $\Delta y = 1$ M\$):

$$\Delta t < \frac{1}{0.008^2 \times 150^2} = \frac{1}{1.44} \approx 0.69 \text{ days}$$

This is very generous — the diffusion CFL is not the binding constraint here because inventory volatility in M\$ units is modest ($\sigma \cdot Y_{\max} \approx 1.2$ M\$/day$^{1/2}$).

With $\Delta t = 5 \times 10^{-5}$ days and 1000 time steps over $T = 0.05$ days, the diffusion CFL has a large safety margin.


## Hedging Hamiltonian stiffness (discovered during implementation)

The initial CFL analysis above only considered the diffusion term. In practice, **the hedging Hamiltonian is the binding stability constraint** and makes explicit Euler unusable with the paper's parameters.

The hedging Hamiltonian is $H^{i,j}(p) = (\max(|p| - \psi, 0))^2 / (4\eta)$. For EURUSD, $\eta = 10^{-5} \text{ bps} \times 10^{-4} = 10^{-9}$ in decimal units. This means $H \sim p^2 / (4 \times 10^{-9}) \approx 2.5 \times 10^8 \, p^2$, which amplifies gradients of $\theta$ by a factor of $\sim 10^8$. The resulting positive feedback loop causes numerical blow-up regardless of how small $\Delta t$ is chosen.

This is not a bug — it reflects genuine stiffness in the PDE. The paper used an implicit scheme ("monotone implicit Euler") precisely because of this.

**Workaround for initial development:** We scale $\eta$ up by a factor (e.g., $\times 1000$) to make explicit Euler stable. This lets us validate the spatial discretization (quoting Hamiltonian, diffusion, drift) independently of the time-stepping scheme. The artificially large $\eta$ reduces the hedging Hamiltonian's gain, making explicit time stepping feasible.

**Upgrade path:** The semi-implicit and fully implicit solvers address this stiffness — see `semi-implicit/` and `implicit/`.
