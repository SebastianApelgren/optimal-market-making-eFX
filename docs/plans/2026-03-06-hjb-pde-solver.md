# HJB PDE Solver Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Solve PDE (1) from Barzykin, Bergault, Guéant (2023) via explicit Euler on a finite inventory grid, and compare the solution against the Riccati ODE approximation.

**Architecture:** Add `PDESpec`, `build_pde_spec`, `pde_rhs`, and `solve_hjb_explicit` to `src/pde.py`. Add PDE-based policy extraction functions. Fill in `hjb_pde_solver.ipynb` with grid setup, solver execution, and validation plots (value function + policy comparison vs ODE).

**Tech Stack:** numpy, matplotlib. No new dependencies.

---

### Task 1: Add PDESpec and build_pde_spec

**Files:**
- Modify: `src/pde.py` (append after `terminal_condition`, ~line 597)

**Step 1: Implement PDESpec dataclass**

```python
@dataclass(frozen=True)
class PDESpec:
    """All precomputed data needed to evaluate the PDE RHS."""
    quoting: QuotingSpec
    hedging: HedgingSpec
    Sigma: np.ndarray       # (d, d) covariance matrix
    mu_vec: np.ndarray      # (d,) drift vector
    gamma: float
    dy_list: List[float]
    y_grids: List[np.ndarray]
    penalty: np.ndarray     # precomputed running_penalty (static, no theta dependence)
```

**Step 2: Implement build_pde_spec**

```python
def build_pde_spec(y_grids: List[np.ndarray], mp: ModelParams) -> PDESpec:
    """One-time setup: build all sub-specs and precompute static arrays."""
    from .riccati import build_Sigma

    dy_list = validate_pde_grid(y_grids, mp)
    Sigma = build_Sigma(mp)
    mu_vec = np.array([mp.mu.get(c, 0.0) for c in mp.currencies], dtype=float)

    quoting = build_quoting_spec(y_grids, mp)
    hedging = build_hedging_spec(y_grids, mp)
    penalty = running_penalty(y_grids, Sigma, mp.gamma)

    return PDESpec(
        quoting=quoting,
        hedging=hedging,
        Sigma=Sigma,
        mu_vec=mu_vec,
        gamma=mp.gamma,
        dy_list=dy_list,
        y_grids=y_grids,
        penalty=penalty,
    )
```

**Step 3: Verify in a scratch notebook cell or Python REPL**

```python
from src import build_paper_example_params, restrict_currencies
from src.pde import build_pde_spec
import numpy as np

mp = restrict_currencies(build_paper_example_params(), ["USD", "EUR"])
y_grids = [np.arange(-150, 151, 1.0) for _ in range(len(mp.currencies))]
spec = build_pde_spec(y_grids, mp)
print(f"QuotingSpec contributions: {len(spec.quoting.contributions)}")
print(f"HedgingSpec contributions: {len(spec.hedging.contributions)}")
print(f"Penalty shape: {spec.penalty.shape}")
print(f"Sigma:\n{spec.Sigma}")
```

Expected: no errors, shapes consistent with d=2 and 301-point grids.

**Step 4: Commit**

```bash
git add src/pde.py
git commit -m "add PDESpec dataclass and build_pde_spec factory"
```

---

### Task 2: Add pde_rhs

**Files:**
- Modify: `src/pde.py` (append after `build_pde_spec`)

**Step 1: Implement pde_rhs**

This assembles the full right-hand side of PDE (1): all five spatial terms summed together. The sign convention is: `dθ/dt + RHS = 0`, so we return `RHS` and the time stepper does `θ^{m-1} = θ^m + Δt * RHS(θ^m)`.

```python
def pde_rhs(theta: np.ndarray, spec: PDESpec) -> np.ndarray:
    """Evaluate the spatial RHS of PDE (1) at every grid point.

    Returns RHS such that dθ/dt = -RHS, i.e. the backward Euler update is
    θ^{m-1} = θ^m + Δt * RHS.
    """
    grad = compute_gradient(theta, spec.dy_list)

    rhs = spec.penalty.copy()                                         # -γ/2 y^T Σ y
    rhs += diffusion_term(theta, spec.y_grids, spec.dy_list, spec.Sigma)  # ½ Tr(...)
    rhs += drift_term(grad, spec.y_grids, spec.mu_vec)                    # y^T μ (1+∂θ)
    rhs += quoting_hamiltonian_integral(theta, spec.quoting)              # Q(y)
    rhs += hedging_hamiltonian(grad, spec.y_grids, spec.hedging)          # H_hedge(y)

    return rhs
```

**Step 2: Verify with a quick smoke test**

```python
from src.pde import build_pde_spec, pde_rhs, terminal_condition
import numpy as np

# Use d=2 (USD-EUR) with a small grid for speed
mp = restrict_currencies(build_paper_example_params(), ["USD", "EUR"])
y_grids = [np.arange(-10, 11, 1.0) for _ in range(len(mp.currencies))]
spec = build_pde_spec(y_grids, mp)
kappa = mp.kappa if mp.kappa is not None else np.zeros((len(mp.currencies), len(mp.currencies)))
theta = terminal_condition(y_grids, kappa)

rhs = pde_rhs(theta, spec)
print(f"RHS shape: {rhs.shape}, range: [{rhs.min():.4f}, {rhs.max():.4f}]")
assert rhs.shape == theta.shape
assert np.all(np.isfinite(rhs))
```

Expected: finite values, same shape as theta.

**Step 3: Commit**

```bash
git add src/pde.py
git commit -m "add pde_rhs assembling all spatial operators for PDE (1)"
```

---

### Task 3: Add solve_hjb_explicit

**Files:**
- Modify: `src/pde.py` (append after `pde_rhs`)

**Step 1: Implement the explicit Euler time stepper**

```python
def solve_hjb_explicit(
    y_grids: List[np.ndarray],
    mp: ModelParams,
    n_steps: int = 1000,
    snapshot_times: List[float] | None = None,
) -> dict:
    """Solve PDE (1) backward from T to 0 via explicit Euler.

    Parameters
    ----------
    y_grids : list of d 1D arrays (grid axes).
    mp : ModelParams with T_days, gamma, kappa, etc.
    n_steps : number of time steps.
    snapshot_times : optional list of times (in days) at which to save θ.
        Always saves t=0.

    Returns
    -------
    dict with keys:
        'theta_0'   : θ(0, y), d-dimensional array on the grid.
        'snapshots' : dict mapping time -> θ(t, y) for requested snapshot_times.
        'spec'      : the PDESpec used.
        'dt'        : time step size.
    """
    d = len(mp.currencies)
    kappa = mp.kappa if mp.kappa is not None else np.zeros((d, d))
    T = mp.T_days
    dt = T / n_steps

    spec = build_pde_spec(y_grids, mp)
    theta = terminal_condition(y_grids, kappa)

    # Prepare snapshot collection
    if snapshot_times is None:
        snapshot_times = []
    snap_steps = {}
    for t_snap in snapshot_times:
        step_idx = round((T - t_snap) / dt)
        step_idx = max(0, min(n_steps, step_idx))
        snap_steps[step_idx] = t_snap
    snapshots = {}

    # Backward march: step m goes from n_steps down to 1
    for m in range(n_steps):
        rhs = pde_rhs(theta, spec)
        theta = theta + dt * rhs

        if (m + 1) in snap_steps:
            t_snap = snap_steps[m + 1]
            snapshots[t_snap] = theta.copy()

    return {
        'theta_0': theta,
        'snapshots': snapshots,
        'spec': spec,
        'dt': dt,
    }
```

**Step 2: Verify on a tiny grid that the solver runs without errors**

```python
from src import build_paper_example_params, restrict_currencies
from src.pde import solve_hjb_explicit
import numpy as np

mp = restrict_currencies(build_paper_example_params(), ["USD", "EUR"])
y_grids = [np.arange(-10, 11, 1.0) for _ in range(len(mp.currencies))]
result = solve_hjb_explicit(y_grids, mp, n_steps=100)
theta_0 = result['theta_0']
print(f"theta_0 shape: {theta_0.shape}")
print(f"theta_0 range: [{theta_0.min():.6f}, {theta_0.max():.6f}]")
print(f"theta_0 at origin: {theta_0[10, 10]:.6f}")
assert np.all(np.isfinite(theta_0))
```

Expected: finite values, θ(0,0) should be positive (value of having zero inventory).

**Step 3: Commit**

```bash
git add src/pde.py
git commit -m "add solve_hjb_explicit backward Euler time stepper"
```

---

### Task 4: Update src/__init__.py with new exports

**Files:**
- Modify: `src/__init__.py` (add new names to the pde import block, ~lines 24-38)

**Step 1: Add new exports**

Add `PDESpec`, `build_pde_spec`, `pde_rhs`, `solve_hjb_explicit` to the import block from `.pde`.

**Step 2: Verify imports work**

```python
from src import PDESpec, build_pde_spec, pde_rhs, solve_hjb_explicit
print("All imports OK")
```

**Step 3: Commit**

```bash
git add src/__init__.py
git commit -m "export PDE solver functions from src/__init__.py"
```

---

### Task 5: Notebook — d=2 grid setup and solver run

**Files:**
- Modify: `hjb_pde_solver.ipynb` (fill in the stub notebook)

**Step 1: Problem setup cell**

```python
import numpy as np
import matplotlib.pyplot as plt
from src import (
    build_paper_example_params, restrict_currencies, run_multicurrency_mm,
    solve_hjb_explicit,
)
from src.riccati import build_Sigma

# d=2: USD-EUR
mp_full = build_paper_example_params()
mp = restrict_currencies(mp_full, ["USD", "EUR"])
d = len(mp.currencies)
print(f"Currencies: {mp.currencies}, d={d}")
```

**Step 2: Grid construction cell**

```python
Y_MAX = 150.0
DY = 1.0
y_axis = np.arange(-Y_MAX, Y_MAX + DY/2, DY)  # [-150, -149, ..., 150]
y_grids = [y_axis for _ in range(d)]
print(f"Grid: {len(y_axis)} points per axis, range [{y_axis[0]}, {y_axis[-1]}]")
print(f"Total grid points: {len(y_axis)**d}")
```

**Step 3: Run solver cell**

```python
N_STEPS = 1000
result = solve_hjb_explicit(y_grids, mp, n_steps=N_STEPS,
                            snapshot_times=[mp.T_days * 0.5, mp.T_days * 0.25])
theta_0 = result['theta_0']
print(f"Solver complete. dt = {result['dt']:.2e} days")
print(f"theta_0 shape: {theta_0.shape}")
print(f"theta_0 at origin: {theta_0[150, 150]:.6f}")
```

**Step 4: Commit**

```bash
git add hjb_pde_solver.ipynb
git commit -m "add d=2 grid setup and PDE solver execution to HJB notebook"
```

---

### Task 6: Notebook — validation A: value function comparison

**Files:**
- Modify: `hjb_pde_solver.ipynb` (add cells after Task 5)

**Step 1: Compute ODE reference solution**

```python
ode_result = run_multicurrency_mm(mp, n_steps=2000)
A0 = ode_result.A0
B0 = ode_result.B0
print(f"A0:\n{A0}")
print(f"B0: {B0}")
```

**Step 2: Plot θ_PDE(0, y) vs ODE ansatz along the EUR inventory axis**

The ODE ansatz is `θ_hat(y) = -y^T A0 y - y^T B0 - C0`. Since `B0 = 0` and `C0` is an additive constant we don't track, compare shapes by shifting to match at y=0.

```python
# Slice along EUR axis (axis 1) with USD inventory = 0
usd_idx = 150  # y_USD = 0
y_eur = y_axis

theta_pde_slice = theta_0[usd_idx, :]

# ODE ansatz: theta_hat = -y^T A0 y - y^T B0 (ignore C0, shift to match at origin)
theta_ode_slice = np.zeros_like(y_eur)
for idx, ye in enumerate(y_eur):
    y_vec = np.array([0.0, ye])  # [y_USD=0, y_EUR=ye]
    theta_ode_slice[idx] = -y_vec @ A0 @ y_vec - B0 @ y_vec

# Shift ODE to match PDE at origin for fair comparison (absorb unknown C0)
shift = theta_pde_slice[150] - theta_ode_slice[150]
theta_ode_slice += shift

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Left: value functions
roi = np.abs(y_eur) <= 100  # region of interest
ax1.plot(y_eur[roi], theta_pde_slice[roi], label='PDE (exact H)', linewidth=2)
ax1.plot(y_eur[roi], theta_ode_slice[roi], '--', label='ODE (quadratic H)', linewidth=2)
ax1.set_xlabel('EUR Inventory (M$)')
ax1.set_ylabel('θ(0, y)')
ax1.set_title('Value Function: PDE vs ODE')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Right: difference
diff = theta_pde_slice[roi] - theta_ode_slice[roi]
ax2.plot(y_eur[roi], diff, 'r-', linewidth=2)
ax2.set_xlabel('EUR Inventory (M$)')
ax2.set_ylabel('θ_PDE - θ_ODE')
ax2.set_title('Value Function Difference')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

**Step 3: Commit**

```bash
git add hjb_pde_solver.ipynb
git commit -m "add value function comparison PDE vs ODE in HJB notebook"
```

---

### Task 7: Notebook — validation B: policy comparison

**Files:**
- Modify: `hjb_pde_solver.ipynb` (add cells after Task 6)

**Step 1: Extract optimal quotes from PDE solution**

The PDE-based optimal quote uses Eq. (2): `δ*(z, p) = δ_bar(z, p)` where `p = (θ(y) - θ(y + z*d_ij)) / z`. We read this directly from the θ grid.

```python
from src.hamiltonian import optimal_delta_logistic, H_logistic
from src.model import BP, canon_pair

def pde_optimal_quote(theta_0, y_grids, mp, tier_idx, ccy_pay, ccy_sell, z_musd, y_slice_axis, y_slice_val):
    """Extract optimal quote from PDE θ grid along a 1D inventory slice.

    Parameters
    ----------
    y_slice_axis : int, which axis to vary (0=USD, 1=EUR, ...)
    y_slice_val : float, value of the fixed axis

    Returns (y_values, delta_values_bps) along the slice.
    """
    ccy = mp.currencies
    i = ccy.index(ccy_pay)
    j = ccy.index(ccy_sell)
    d = len(ccy)
    key = canon_pair(ccy_pay, ccy_sell)
    pp = mp.pairs[key]
    tier = pp.tiers[tier_idx]

    dy = y_grids[y_slice_axis][1] - y_grids[y_slice_axis][0]
    shift_i = round(z_musd / dy) if i == y_slice_axis else 0
    shift_j = round(-z_musd / dy) if j == y_slice_axis else 0

    # For the fixed axis, find the index
    fixed_axis = 1 - y_slice_axis  # works for d=2
    fixed_idx = int(round((y_slice_val - y_grids[fixed_axis][0]) / dy))

    n = len(y_grids[y_slice_axis])
    y_vals = []
    delta_vals = []

    for k in range(n):
        # Build index for theta(y)
        idx_here = [0, 0]
        idx_here[y_slice_axis] = k
        idx_here[fixed_axis] = fixed_idx

        # Build index for theta(y + z * d_ij)
        idx_shifted = list(idx_here)
        idx_shifted[i] += round(z_musd / dy)
        idx_shifted[j] += round(-z_musd / dy)

        # Check bounds
        if any(s < 0 or s >= len(y_grids[ax]) for ax, s in enumerate(idx_shifted)):
            continue

        theta_here = theta_0[tuple(idx_here)]
        theta_shifted = theta_0[tuple(idx_shifted)]
        p = (theta_here - theta_shifted) / z_musd

        delta_star = optimal_delta_logistic(p, tier.alpha, tier.beta)
        y_vals.append(y_grids[y_slice_axis][k])
        delta_vals.append(delta_star / BP)  # convert to bps

    return np.array(y_vals), np.array(delta_vals)
```

**Step 2: Compare PDE vs ODE quotes (reproducing Figure 1 style)**

```python
# Optimal EURUSD top-of-book quote (tier 1, z=1 M$) as function of EUR inventory
y_pde, delta_pde = pde_optimal_quote(
    theta_0, y_grids, mp, tier_idx=0,
    ccy_pay="EUR", ccy_sell="USD", z_musd=1.0,
    y_slice_axis=1, y_slice_val=0.0,
)

# ODE reference: use the existing optimal_client_markup
from src.policy import optimal_client_markup
y_ode = np.arange(-100, 101, 1.0)
delta_ode = np.array([
    optimal_client_markup(mp, 0, "EUR", "USD", 1.0,
                          np.array([0.0, ye]), A0, B0) / BP
    for ye in y_ode
])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

roi_pde = np.abs(y_pde) <= 100
ax1.plot(y_pde[roi_pde], delta_pde[roi_pde], label='PDE', linewidth=2)
ax1.plot(y_ode, delta_ode, '--', label='ODE', linewidth=2)
ax1.set_xlabel('EUR Inventory (M$)')
ax1.set_ylabel('δ* (bps)')
ax1.set_title('EURUSD Optimal Quote: PDE vs ODE (tier 1, z=1 M$)')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Difference
common_y = y_ode
common_pde = np.interp(common_y, y_pde[roi_pde], delta_pde[roi_pde])
ax2.plot(common_y, common_pde - delta_ode, 'r-', linewidth=2)
ax2.set_xlabel('EUR Inventory (M$)')
ax2.set_ylabel('δ*_PDE - δ*_ODE (bps)')
ax2.set_title('Quote Difference')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

**Step 3: Compare PDE vs ODE hedge rates**

```python
from src.pde import compute_gradient, H_execution_cost

def pde_optimal_hedge_rate(theta_0, y_grids, dy_list, mp, ccy_buy, ccy_sell, y_slice_axis, y_slice_val):
    """Extract hedge rate from PDE θ grid along a 1D slice."""
    ccy = mp.currencies
    i = ccy.index(ccy_buy)
    j = ccy.index(ccy_sell)
    key = canon_pair(ccy_buy, ccy_sell)
    pp = mp.pairs[key]
    d = len(ccy)

    grad = compute_gradient(theta_0, dy_list)

    dy = dy_list[y_slice_axis]
    fixed_axis = 1 - y_slice_axis
    fixed_idx = int(round((y_slice_val - y_grids[fixed_axis][0]) / dy_list[fixed_axis]))

    k_i = mp.k.get(ccy_buy, 0.0)
    k_j = mp.k.get(ccy_sell, 0.0)

    y_vals = []
    xi_vals = []

    for k in range(len(y_grids[y_slice_axis])):
        idx = [0, 0]
        idx[y_slice_axis] = k
        idx[fixed_axis] = fixed_idx
        idx = tuple(idx)

        gi = grad[i][idx]
        gj = grad[j][idx]
        yi = y_grids[i][idx[i]] if i < d else 0.0
        yj = y_grids[j][idx[j]] if j < d else 0.0

        # Eq. (3): p = ∂θ/∂y_i - ∂θ/∂y_j + k_i y_i (1+∂θ/∂y_i) - k_j y_j (1+∂θ/∂y_j)
        p = (gi - gj) + k_i * yi * (1.0 + gi) - k_j * yj * (1.0 + gj)

        # ξ* = H'(p) = sign(p) * max(|p| - ψ, 0) / (2η)
        if p > pp.psi:
            xi = (p - pp.psi) / (2.0 * pp.eta)
        elif p < -pp.psi:
            xi = (p + pp.psi) / (2.0 * pp.eta)
        else:
            xi = 0.0

        y_vals.append(y_grids[y_slice_axis][k])
        xi_vals.append(xi)

    return np.array(y_vals), np.array(xi_vals)


# Compare
from src.policy import optimal_hedge_rate

y_pde_h, xi_pde = pde_optimal_hedge_rate(
    theta_0, y_grids, result['spec'].dy_list, mp,
    "EUR", "USD", y_slice_axis=1, y_slice_val=0.0,
)

y_ode_h = np.arange(-100, 101, 1.0)
xi_ode = np.array([
    optimal_hedge_rate(mp, "EUR", "USD", np.array([0.0, ye]), A0, B0)
    for ye in y_ode_h
])

fig, ax = plt.subplots(figsize=(8, 5))
roi_h = np.abs(y_pde_h) <= 100
ax.plot(y_pde_h[roi_h], xi_pde[roi_h], label='PDE', linewidth=2)
ax.plot(y_ode_h, xi_ode, '--', label='ODE', linewidth=2)
ax.set_xlabel('EUR Inventory (M$)')
ax.set_ylabel('ξ* (M$/day)')
ax.set_title('EURUSD Optimal Hedge Rate: PDE vs ODE')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

**Step 4: Commit**

```bash
git add hjb_pde_solver.ipynb
git commit -m "add policy comparison (quotes + hedge rates) PDE vs ODE"
```

---

### Task 8: Stability check and stationarity verification

**Files:**
- Modify: `hjb_pde_solver.ipynb` (add cells after Task 7)

**Step 1: Check that θ has reached stationarity**

Run the solver with snapshots at several times and verify the solution has converged.

```python
# Re-run with more snapshots to check convergence
snap_times = [mp.T_days * f for f in [0.75, 0.5, 0.25, 0.1, 0.0]]
result_snap = solve_hjb_explicit(y_grids, mp, n_steps=N_STEPS, snapshot_times=snap_times)

fig, ax = plt.subplots(figsize=(8, 5))
usd_idx = 150
roi = slice(50, 251)  # y in [-100, 100]

for t_snap in sorted(result_snap['snapshots'].keys(), reverse=True):
    theta_snap = result_snap['snapshots'][t_snap]
    label = f't = {t_snap:.4f} days ({t_snap*24*60:.1f} min)'
    ax.plot(y_axis[roi], theta_snap[usd_idx, roi], label=label, alpha=0.7)

ax.plot(y_axis[roi], result_snap['theta_0'][usd_idx, roi],
        'k-', linewidth=2, label='t = 0')
ax.set_xlabel('EUR Inventory (M$)')
ax.set_ylabel('θ(t, y)')
ax.set_title('Convergence to Stationarity')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
```

**Step 2: Check for numerical stability (no NaN/Inf, no blow-up)**

```python
print(f"theta_0 finite: {np.all(np.isfinite(result_snap['theta_0']))}")
print(f"theta_0 range: [{result_snap['theta_0'].min():.4f}, {result_snap['theta_0'].max():.4f}]")

# Check max absolute difference between last two snapshots
times_sorted = sorted(result_snap['snapshots'].keys())
if len(times_sorted) >= 2:
    t1, t2 = times_sorted[0], times_sorted[1]
    max_diff = np.max(np.abs(result_snap['snapshots'][t1] - result_snap['snapshots'][t2]))
    print(f"Max |θ(t={t1:.4f}) - θ(t={t2:.4f})| = {max_diff:.2e}")
```

**Step 3: Commit**

```bash
git add hjb_pde_solver.ipynb
git commit -m "add stationarity and stability checks to HJB notebook"
```

---

## Summary of tasks

| Task | What | Files |
|------|------|-------|
| 1 | PDESpec + build_pde_spec | `src/pde.py` |
| 2 | pde_rhs (assemble all operators) | `src/pde.py` |
| 3 | solve_hjb_explicit (time stepper) | `src/pde.py` |
| 4 | Update __init__.py exports | `src/__init__.py` |
| 5 | Notebook: grid setup + solver run | `hjb_pde_solver.ipynb` |
| 6 | Notebook: value function comparison | `hjb_pde_solver.ipynb` |
| 7 | Notebook: policy comparison | `hjb_pde_solver.ipynb` |
| 8 | Notebook: stationarity + stability | `hjb_pde_solver.ipynb` |

Tasks 1-4 are pure library code. Tasks 5-8 are notebook cells. Each task builds on the previous one.
