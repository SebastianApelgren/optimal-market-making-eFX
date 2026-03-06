# Terminal Condition

`terminal_condition(y_grids, kappa)` computes θ(T, y) = -ℓ(y) = -y^T κ y on the grid.

Same broadcasting pattern as `running_penalty`. In the paper's example κ = 0, so θ(T, y) = 0 everywhere.
