# Terminal Condition

`terminal_condition(y_grids, kappa)` computes θ(T, y) = -ℓ(y) = -y^T κ y on the grid.

Same broadcasting pattern as `running_penalty`.

In the paper's example κ = 0, so θ(T, y) = 0 everywhere. The paper states (p.6): "both the drift vector μ and the terminal penalty κ are assumed to be 0." With the short time horizon T = 0.05 days (72 min), the solution converges to a stationary regime well before T, so the terminal condition has no influence on the quotes at t = 0. Setting κ = 0 avoids end-of-horizon effects that would distort the stationary solution.
