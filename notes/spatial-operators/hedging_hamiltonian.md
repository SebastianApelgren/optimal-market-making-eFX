# Hedging Hamiltonian — Design Notes

## The term in the HJB PDE

The third line of the HJB equation (Eq. 1 in the paper) is the **hedging Hamiltonian**:

$$\sum_{1 \leq i < j \leq d} \mathcal{H}^{i,j}\!\Big(\partial_{y^i}\theta - \partial_{y^j}\theta + k^i y^i (1 + \partial_{y^i}\theta) - k^j y^j (1 + \partial_{y^j}\theta)\Big)$$

where

$$\mathcal{H}^{i,j}(p) = \sup_\xi \Big(p\,\xi - L^{i,j}(\xi)\Big) \quad \text{with} \quad L^{i,j}(\xi) = \psi^{i,j}|\xi| + \eta^{i,j}\xi^2.$$


## What it represents

The market maker can hedge inventory risk by trading on the D2D (dealer-to-dealer) market. For each currency pair $(i,j)$, the hedging rate $\xi$ is a continuous control: how many M\$/day to externalize.

At each instant the market maker chooses $\xi$ to maximize the net benefit:

$$p\,\xi - L(\xi)$$

- $p\,\xi$ is the **benefit** of hedging. The argument $p$ is constructed from $\nabla_y \theta$ and captures how sensitive the value function is to inventory changes along the $(i,j)$ direction — it is the marginal value of shifting inventory.
- $L(\xi) = \psi|\xi| + \eta\xi^2$ is the **cost** of trading on the D2D market:
  - $\psi|\xi|$ — a **proportional (spread) cost**. The market maker pays $\psi$ per unit traded, regardless of size. This models the bid-ask spread on D2D platforms.
  - $\eta\xi^2$ — a **quadratic market impact** cost. Trading faster is disproportionately more expensive.

The Hamiltonian $\mathcal{H}(p)$ is the optimal net value: benefit minus cost at the best choice of $\xi$.


## The argument $p$ — marginal value of hedging

For each D2D pair $(i,j)$ with $i < j$, the argument to $\mathcal{H}^{i,j}$ is:

$$p^{i,j}(y) = \partial_{y^i}\theta - \partial_{y^j}\theta + k^i y^i (1 + \partial_{y^i}\theta) - k^j y^j (1 + \partial_{y^j}\theta).$$

This is the **total marginal value of hedging one unit** along the $(i,j)$ direction. It collects every channel through which a small hedge $d\xi$ (buying currency $i$, selling currency $j$ on D2D) affects the market maker's expected value:

**1. Direct inventory shift.** Buying $d\xi$ of currency $i$ and selling $d\xi$ of currency $j$ directly changes the inventory vector. The effect on the value function is:

$$(\partial_{y^i}\theta - \partial_{y^j}\theta) \cdot d\xi$$

This is the chain rule: how much $\theta$ changes when $y^i$ goes up and $y^j$ goes down.

**2. Price impact on existing inventory (currency $i$).** From the price dynamics (p.3 of the paper), trading at rate $\xi^{i,j}$ permanently moves the exchange rate $S^i$ via the market impact parameter $k^i$:

$$dS^i \ni k^i \, S^i \, \xi^{i,j} \, dt$$

This price change affects the existing inventory $Y^i = q^i S^i$: it grows by $k^i y^i \, d\xi$. This inventory change affects both the P&L account directly (the position $y^i$ is now worth more, contributing $k^i y^i \cdot d\xi$) and the value function (contributing $k^i y^i \cdot \partial_{y^i}\theta \cdot d\xi$). Together:

$$k^i y^i (1 + \partial_{y^i}\theta) \cdot d\xi$$

**3. Same for currency $j$, with opposite sign.** Selling currency $j$ pushes $S^j$ down:

$$-k^j y^j (1 + \partial_{y^j}\theta) \cdot d\xi$$

The HJB collects all three channels into $p \cdot d\xi$ and then optimizes over $\xi$: $\sup_\xi (p\xi - L(\xi)) = \mathcal{H}(p)$.

**Simplification in the paper's base case:** The paper sets $k^i = 0$ for all currencies (no permanent market impact), so $p$ reduces to just $\partial_{y^i}\theta - \partial_{y^j}\theta$ — the pure inventory sensitivity of the value function.


## Closed-form derivation

We compute $\mathcal{H}(p) = \sup_\xi \bigl(p\,\xi - \psi|\xi| - \eta\xi^2\bigr)$ by case analysis. Define $g(\xi) = p\,\xi - \psi|\xi| - \eta\xi^2$.

### Case 1: $\xi > 0$

Here $|\xi| = \xi$, so

$$g(\xi) = (p - \psi)\,\xi - \eta\,\xi^2.$$

This is a downward parabola ($\eta > 0$). Setting $g'(\xi) = 0$:

$$p - \psi - 2\eta\xi = 0 \quad \Longrightarrow \quad \xi^* = \frac{p - \psi}{2\eta}.$$

Valid only when $\xi^* > 0$, i.e. $p > \psi$. Substituting back:

$$g(\xi^*) = (p - \psi) \cdot \frac{p - \psi}{2\eta} - \eta \cdot \frac{(p - \psi)^2}{4\eta^2} = \frac{(p - \psi)^2}{4\eta}.$$

### Case 2: $\xi < 0$

Here $|\xi| = -\xi$, so

$$g(\xi) = (p + \psi)\,\xi - \eta\,\xi^2.$$

Setting $g'(\xi) = 0$:

$$\xi^* = \frac{p + \psi}{2\eta}.$$

Valid ($\xi^* < 0$) only when $p < -\psi$. The value is $g(\xi^*) = \frac{(p + \psi)^2}{4\eta}$.

### Case 3: $|p| \leq \psi$ (dead zone)

When $|p| \leq \psi$, neither half-line yields a valid interior maximum:

- For $\xi > 0$: $g'(\xi) = (p - \psi) - 2\eta\xi$. Since $p \leq \psi$, the slope at $\xi = 0^+$ is $(p - \psi) \leq 0$, and $-2\eta\xi < 0$ makes it only more negative. So $g$ is strictly decreasing on $(0, \infty)$.
- For $\xi < 0$: $g'(\xi) = (p + \psi) - 2\eta\xi$. Since $p \geq -\psi$, the slope is $(p + \psi) \geq 0$, and $-2\eta\xi > 0$ for $\xi < 0$. So $g$ is strictly increasing on $(-\infty, 0)$.

Both sides slope toward $\xi = 0$, so the supremum is $g(0) = 0$. Economically: the marginal benefit $|p|$ does not cover the per-unit spread cost $\psi$, so any hedge — no matter how small — loses money. The market maker is better off doing nothing.

### Result

$$\boxed{\mathcal{H}(p) = \frac{\bigl(\max(|p| - \psi,\; 0)\bigr)^2}{4\eta}}$$

and the optimizer is

$$\xi^*(p) = \begin{cases} \frac{p - \psi}{2\eta} & \text{if } p > \psi, \\[4pt] \frac{p + \psi}{2\eta} & \text{if } p < -\psi, \\[4pt] 0 & \text{if } |p| \leq \psi. \end{cases}$$

Note: $\xi^*$ is exactly what `Hprime_execution_cost` in `policy.py` computes (since $\xi^* = \mathcal{H}'(p)$).


## PDE solver implementation plan

### Computing the argument $p(y)$

At each grid point $y$, for each D2D pair $(i,j)$ with $i < j$, we need:

$$p(y) = \partial_{y^i}\theta - \partial_{y^j}\theta + k^i y^i (1 + \partial_{y^i}\theta) - k^j y^j (1 + \partial_{y^j}\theta).$$

This requires the gradient $\nabla_y \theta$ on the grid. We use **central finite differences** in the interior:

$$\partial_{y^i}\theta \;\approx\; \frac{\theta[\ldots, k+1, \ldots] - \theta[\ldots, k-1, \ldots]}{2\,\Delta y_i}$$

and one-sided (forward/backward) differences at the boundaries.

### Evaluating $\mathcal{H}(p)$

Once $p$ is assembled at every grid point, the closed-form piecewise quadratic above is applied **pointwise and vectorized** — no optimization loop or iterative solve is needed.

### Summing over pairs

The full hedging contribution at grid point $y$ is:

$$\mathcal{H}_{\text{hedge}}(y) = \sum_{1 \leq i < j \leq d} \mathcal{H}^{i,j}\!\bigl(p^{i,j}(y)\bigr)$$

where each pair $(i,j)$ has its own $(\psi^{i,j}, \eta^{i,j})$ from `mp.pairs`.

### Boundary considerations

The gradient approximation degrades at the grid boundary (one-sided differences are first-order instead of second-order for central differences). However, this only affects the outermost layer of grid points. We rely on $Y_{\max}$ being large enough that the boundary does not influence the solution in the region of interest.
