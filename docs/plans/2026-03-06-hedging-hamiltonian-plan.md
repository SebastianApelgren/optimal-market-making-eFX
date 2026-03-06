# Hedging Hamiltonian Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement the hedging Hamiltonian term from HJB PDE Eq. 1 as reusable functions in `src/pde.py`.

**Architecture:** Precomputed-spec pattern (matching existing `QuotingSpec`). Gradient of theta computed separately via central FD and passed in. Closed-form piecewise quadratic H(p) evaluated vectorized.

**Tech Stack:** numpy

---

### Task 1: `H_execution_cost`

**Files:**
- Modify: `src/pde.py` (append after `quoting_hamiltonian_integral`)
- Create: `tests/test_pde_hedging.py`

**Step 1: Write tests**

Create `tests/test_pde_hedging.py`:

```python
"""Tests for the hedging Hamiltonian functions in src/pde.py."""
from __future__ import annotations

import numpy as np
import pytest

from src.pde import H_execution_cost


class TestHExecutionCost:
    """Tests for H(p) = (max(|p| - psi, 0))^2 / (4*eta)."""

    def test_dead_zone_zero(self):
        """H(p) = 0 when |p| <= psi."""
        psi, eta = 0.5, 1.0
        assert H_execution_cost(0.0, psi, eta) == 0.0
        assert H_execution_cost(0.5, psi, eta) == 0.0
        assert H_execution_cost(-0.5, psi, eta) == 0.0
        assert H_execution_cost(0.3, psi, eta) == 0.0

    def test_positive_p_above_psi(self):
        """H(p) = (p - psi)^2 / (4*eta) when p > psi."""
        psi, eta = 0.5, 2.0
        p = 1.5
        expected = (1.5 - 0.5) ** 2 / (4.0 * 2.0)  # 1.0 / 8.0 = 0.125
        assert H_execution_cost(p, psi, eta) == pytest.approx(expected)

    def test_negative_p_below_neg_psi(self):
        """H(p) = (|p| - psi)^2 / (4*eta) when p < -psi. Symmetric."""
        psi, eta = 0.5, 2.0
        assert H_execution_cost(-1.5, psi, eta) == pytest.approx(
            H_execution_cost(1.5, psi, eta)
        )

    def test_vectorized(self):
        """Works on numpy arrays."""
        psi, eta = 1.0, 1.0
        p = np.array([-3.0, -0.5, 0.0, 0.5, 3.0])
        result = H_execution_cost(p, psi, eta)
        expected = np.array([
            (3.0 - 1.0) ** 2 / 4.0,  # 1.0
            0.0,
            0.0,
            0.0,
            (3.0 - 1.0) ** 2 / 4.0,  # 1.0
        ])
        np.testing.assert_allclose(result, expected)

    def test_continuity_at_psi(self):
        """H is continuous at |p| = psi (value is 0 from both sides)."""
        psi, eta = 0.5, 1.0
        eps = 1e-10
        assert H_execution_cost(psi + eps, psi, eta) == pytest.approx(0.0, abs=1e-15)
        assert H_execution_cost(psi, psi, eta) == 0.0
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_pde_hedging.py -v`
Expected: ImportError — `H_execution_cost` does not exist yet.

**Step 3: Implement `H_execution_cost` in `src/pde.py`**

Append after the `quoting_hamiltonian_integral` function:

```python
def H_execution_cost(p, psi: float, eta: float):
    """Hedging Hamiltonian H(p) = (max(|p| - psi, 0))^2 / (4*eta).

    Legendre-Fenchel transform of L(xi) = psi*|xi| + eta*xi^2.
    Returns 0 in the dead zone |p| <= psi.
    """
    p = np.asarray(p, dtype=float)
    excess = np.maximum(np.abs(p) - psi, 0.0)
    return excess * excess / (4.0 * eta)
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_pde_hedging.py::TestHExecutionCost -v`
Expected: All 5 tests PASS.

**Step 5: Commit**

```bash
git add src/pde.py tests/test_pde_hedging.py
git commit -m "add H_execution_cost: closed-form hedging Hamiltonian"
```

---

### Task 2: `compute_gradient`

**Files:**
- Modify: `src/pde.py` (append after `H_execution_cost`)
- Modify: `tests/test_pde_hedging.py` (add new test class)

**Step 1: Write tests**

Append to `tests/test_pde_hedging.py`:

```python
from src.pde import compute_gradient


class TestComputeGradient:
    """Tests for central-FD gradient on uniform grids."""

    def test_linear_1d_exact(self):
        """Gradient of theta = 3*y is exactly 3 everywhere."""
        y = np.linspace(-5.0, 5.0, 11)  # dy = 1.0
        theta = 3.0 * y
        grad = compute_gradient(theta, [1.0])
        np.testing.assert_allclose(grad[0], 3.0, atol=1e-12)

    def test_quadratic_1d_interior(self):
        """Gradient of theta = y^2 is 2*y. Central FD is exact for quadratics."""
        y = np.linspace(-5.0, 5.0, 101)  # dy = 0.1
        dy = y[1] - y[0]
        theta = y ** 2
        grad = compute_gradient(theta, [dy])
        # Interior: central differences are exact for quadratics
        np.testing.assert_allclose(grad[0][1:-1], 2.0 * y[1:-1], atol=1e-10)

    def test_2d_grid(self):
        """Gradient on a 2D grid: theta = 2*y0 + 3*y1."""
        y0 = np.linspace(-5.0, 5.0, 11)
        y1 = np.linspace(-5.0, 5.0, 11)
        dy0, dy1 = y0[1] - y0[0], y1[1] - y1[0]
        Y0, Y1 = np.meshgrid(y0, y1, indexing="ij")
        theta = 2.0 * Y0 + 3.0 * Y1
        grad = compute_gradient(theta, [dy0, dy1])
        np.testing.assert_allclose(grad[0], 2.0, atol=1e-10)
        np.testing.assert_allclose(grad[1], 3.0, atol=1e-10)

    def test_returns_d_arrays(self):
        """Returns one array per dimension, each same shape as theta."""
        theta = np.zeros((5, 7, 3))
        grad = compute_gradient(theta, [1.0, 1.0, 1.0])
        assert len(grad) == 3
        for g in grad:
            assert g.shape == (5, 7, 3)
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_pde_hedging.py::TestComputeGradient -v`
Expected: ImportError — `compute_gradient` does not exist yet.

**Step 3: Implement `compute_gradient` in `src/pde.py`**

Append after `H_execution_cost`:

```python
def compute_gradient(theta: np.ndarray, dy_list: List[float]) -> List[np.ndarray]:
    """Compute gradient of theta on a uniform grid via finite differences.

    Uses second-order central differences in the interior and first-order
    one-sided differences at the boundaries.

    Parameters
    ----------
    theta : array of any dimension d
    dy_list : list of d floats, grid spacing per axis

    Returns
    -------
    List of d arrays, each same shape as theta. Entry k is d(theta)/d(y_k).
    """
    d = theta.ndim
    grad = []
    for k in range(d):
        n = theta.shape[k]
        g = np.empty_like(theta)

        # Central differences for interior
        sl_plus = [slice(None)] * d
        sl_minus = [slice(None)] * d
        sl_center = [slice(None)] * d
        sl_plus[k] = slice(2, n)
        sl_minus[k] = slice(0, n - 2)
        sl_center[k] = slice(1, n - 1)
        g[tuple(sl_center)] = (theta[tuple(sl_plus)] - theta[tuple(sl_minus)]) / (2.0 * dy_list[k])

        # Forward difference at left boundary
        sl_0 = [slice(None)] * d
        sl_1 = [slice(None)] * d
        sl_0[k] = 0
        sl_1[k] = 1
        g[tuple(sl_0)] = (theta[tuple(sl_1)] - theta[tuple(sl_0)]) / dy_list[k]

        # Backward difference at right boundary
        sl_last = [slice(None)] * d
        sl_prev = [slice(None)] * d
        sl_last[k] = n - 1
        sl_prev[k] = n - 2
        g[tuple(sl_last)] = (theta[tuple(sl_last)] - theta[tuple(sl_prev)]) / dy_list[k]

        grad.append(g)
    return grad
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_pde_hedging.py::TestComputeGradient -v`
Expected: All 4 tests PASS.

**Step 5: Commit**

```bash
git add src/pde.py tests/test_pde_hedging.py
git commit -m "add compute_gradient: central FD gradient on uniform grids"
```

---

### Task 3: `HedgingSpec` + `build_hedging_spec`

**Files:**
- Modify: `src/pde.py` (append after `compute_gradient`)
- Modify: `tests/test_pde_hedging.py` (add new test class)

**Step 1: Write tests**

Append to `tests/test_pde_hedging.py`:

```python
from src.model import (
    ModelParams, PairParams, TierParams, canon_pair,
    build_paper_example_params, restrict_currencies,
)
from src.pde import HedgingSpec, build_hedging_spec


class TestBuildHedgingSpec:
    """Tests for HedgingSpec construction."""

    def _make_2ccy_mp(self):
        """Minimal 2-currency model: USD + EUR."""
        mp = build_paper_example_params()
        return restrict_currencies(mp, ["USD", "EUR"])

    def test_returns_hedging_spec(self):
        mp = self._make_2ccy_mp()
        y_grids = [np.arange(-50, 51, 1.0), np.arange(-50, 51, 1.0)]
        spec = build_hedging_spec(y_grids, mp)
        assert isinstance(spec, HedgingSpec)
        assert spec.d == 2

    def test_one_pair_one_contribution(self):
        """2-currency model has exactly one D2D pair (0,1)."""
        mp = self._make_2ccy_mp()
        y_grids = [np.arange(-50, 51, 1.0), np.arange(-50, 51, 1.0)]
        spec = build_hedging_spec(y_grids, mp)
        assert len(spec.contributions) == 1

    def test_contribution_values(self):
        """Check that psi, eta, k values are extracted correctly."""
        mp = self._make_2ccy_mp()
        y_grids = [np.arange(-50, 51, 1.0), np.arange(-50, 51, 1.0)]
        spec = build_hedging_spec(y_grids, mp)
        i, j, psi, eta, k_i, k_j = spec.contributions[0]
        key = canon_pair("USD", "EUR")
        pp = mp.pairs[key]
        assert psi == pp.psi
        assert eta == pp.eta
        assert k_i == mp.k.get("USD", 0.0)
        assert k_j == mp.k.get("EUR", 0.0)
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_pde_hedging.py::TestBuildHedgingSpec -v`
Expected: ImportError — `HedgingSpec`, `build_hedging_spec` do not exist yet.

**Step 3: Implement in `src/pde.py`**

Append after `compute_gradient`:

```python
@dataclass(frozen=True)
class HedgingSpec:
    """Precomputed specification for the hedging Hamiltonian."""
    contributions: tuple  # tuple of (i, j, psi, eta, k_i, k_j)
    d: int                # number of currencies


def build_hedging_spec(
    y_grids: List[np.ndarray],
    mp: ModelParams,
) -> HedgingSpec:
    """Precompute hedging pair info for the hedging Hamiltonian.

    Loops over canonical pairs (i, j) with i < j that exist in mp.pairs.
    """
    ccy = mp.currencies
    d = len(ccy)
    contributions = []

    for i in range(d):
        for j in range(i + 1, d):
            key = canon_pair(ccy[i], ccy[j])
            if key not in mp.pairs:
                continue
            pp = mp.pairs[key]
            k_i = mp.k.get(ccy[i], 0.0)
            k_j = mp.k.get(ccy[j], 0.0)
            contributions.append((i, j, pp.psi, pp.eta, k_i, k_j))

    return HedgingSpec(contributions=tuple(contributions), d=d)
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_pde_hedging.py::TestBuildHedgingSpec -v`
Expected: All 3 tests PASS.

**Step 5: Commit**

```bash
git add src/pde.py tests/test_pde_hedging.py
git commit -m "add HedgingSpec and build_hedging_spec for D2D pair precomputation"
```

---

### Task 4: `hedging_hamiltonian`

**Files:**
- Modify: `src/pde.py` (append after `build_hedging_spec`)
- Modify: `tests/test_pde_hedging.py` (add new test class)

**Step 1: Write tests**

Append to `tests/test_pde_hedging.py`:

```python
from src.pde import hedging_hamiltonian


class TestHedgingHamiltonian:
    """Tests for the full hedging Hamiltonian evaluation."""

    def _make_2ccy_mp(self):
        mp = build_paper_example_params()
        return restrict_currencies(mp, ["USD", "EUR"])

    def test_zero_gradient_gives_zero(self):
        """When grad_theta = 0 everywhere and k = 0, p = 0 => H = 0."""
        mp = self._make_2ccy_mp()
        # Override k to zero for clean test
        mp.k = {c: 0.0 for c in mp.currencies}
        n = 21
        y_grids = [np.arange(-10, 11, 1.0), np.arange(-10, 11, 1.0)]
        spec = build_hedging_spec(y_grids, mp)
        grad_theta = [np.zeros((n, n)), np.zeros((n, n))]
        result = hedging_hamiltonian(grad_theta, y_grids, spec)
        np.testing.assert_allclose(result, 0.0)

    def test_uniform_gradient_no_market_impact(self):
        """With k=0, p = grad_theta[i] - grad_theta[j]. Constant gradient => constant p."""
        mp = self._make_2ccy_mp()
        mp.k = {c: 0.0 for c in mp.currencies}
        n = 21
        y_grids = [np.arange(-10, 11, 1.0), np.arange(-10, 11, 1.0)]
        spec = build_hedging_spec(y_grids, mp)
        # grad_theta[0] = 1.0, grad_theta[1] = -1.0 => p = 1 - (-1) = 2
        grad_theta = [np.ones((n, n)), -np.ones((n, n))]
        result = hedging_hamiltonian(grad_theta, y_grids, spec)
        i, j, psi, eta, _, _ = spec.contributions[0]
        expected = H_execution_cost(2.0, psi, eta)
        np.testing.assert_allclose(result, expected)

    def test_output_shape(self):
        """Output has same shape as grad_theta arrays."""
        mp = self._make_2ccy_mp()
        n = 21
        y_grids = [np.arange(-10, 11, 1.0), np.arange(-10, 11, 1.0)]
        spec = build_hedging_spec(y_grids, mp)
        grad_theta = [np.zeros((n, n)), np.zeros((n, n))]
        result = hedging_hamiltonian(grad_theta, y_grids, spec)
        assert result.shape == (n, n)

    def test_market_impact_increases_p(self):
        """With k > 0 and positive inventory, |p| should be larger than without."""
        mp = self._make_2ccy_mp()
        n = 21
        y_grids = [np.arange(-10, 11, 1.0), np.arange(-10, 11, 1.0)]

        # Without market impact
        mp_no_k = self._make_2ccy_mp()
        mp_no_k.k = {c: 0.0 for c in mp_no_k.currencies}
        spec_no_k = build_hedging_spec(y_grids, mp_no_k)

        # With market impact (paper values)
        spec_k = build_hedging_spec(y_grids, mp)

        # Use a gradient that puts us outside the dead zone
        grad_theta = [np.full((n, n), 0.5), np.full((n, n), -0.5)]
        result_no_k = hedging_hamiltonian(grad_theta, y_grids, spec_no_k)
        result_k = hedging_hamiltonian(grad_theta, y_grids, spec_k)

        # At y = (0, 0), market impact terms vanish, so results should match there
        mid = n // 2
        assert result_no_k[mid, mid] == pytest.approx(result_k[mid, mid])

        # At y = (10, 0), k*y terms contribute, so H should differ
        # (exact direction depends on sign, but they should not be equal everywhere)
        assert not np.allclose(result_no_k, result_k)
```

**Step 2: Run tests to verify they fail**

Run: `.venv/bin/python -m pytest tests/test_pde_hedging.py::TestHedgingHamiltonian -v`
Expected: ImportError — `hedging_hamiltonian` does not exist yet.

**Step 3: Implement `hedging_hamiltonian` in `src/pde.py`**

Append after `build_hedging_spec`:

```python
def hedging_hamiltonian(
    grad_theta: List[np.ndarray],
    y_grids: List[np.ndarray],
    spec: HedgingSpec,
) -> np.ndarray:
    """Evaluate the hedging Hamiltonian at every grid point.

    H_hedge(y) = sum_{i<j} H^{i,j}(p^{i,j}(y))

    where p^{i,j} = grad_theta[i] - grad_theta[j]
                   + k_i * y_i * (1 + grad_theta[i])
                   - k_j * y_j * (1 + grad_theta[j]).

    Parameters
    ----------
    grad_theta : list of d arrays, each same shape as theta.
        Entry k is d(theta)/d(y_k), e.g. from compute_gradient.
    y_grids : list of d 1D arrays (grid axes).
    spec : HedgingSpec from build_hedging_spec.
    """
    shape = grad_theta[0].shape
    result = np.zeros(shape)

    # Precompute broadcasted y arrays: y_grids[k] along axis k
    d = spec.d
    y_broadcast = []
    for k in range(d):
        slices = [np.newaxis] * d
        slices[k] = slice(None)
        y_broadcast.append(y_grids[k][tuple(slices)])

    for (i, j, psi, eta, k_i, k_j) in spec.contributions:
        p = (grad_theta[i] - grad_theta[j]
             + k_i * y_broadcast[i] * (1.0 + grad_theta[i])
             - k_j * y_broadcast[j] * (1.0 + grad_theta[j]))
        result += H_execution_cost(p, psi, eta)

    return result
```

**Step 4: Run tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_pde_hedging.py::TestHedgingHamiltonian -v`
Expected: All 4 tests PASS.

**Step 5: Commit**

```bash
git add src/pde.py tests/test_pde_hedging.py
git commit -m "add hedging_hamiltonian: full hedging term evaluation on grid"
```

---

### Task 5: Update exports and run full test suite

**Files:**
- Modify: `src/__init__.py`

**Step 1: Update `src/__init__.py`**

Replace the existing `from .pde import` block with:

```python
from .pde import (
    validate_pde_grid,
    QuotingSpec,
    build_quoting_spec,
    quoting_hamiltonian_integral,
    H_execution_cost,
    compute_gradient,
    HedgingSpec,
    build_hedging_spec,
    hedging_hamiltonian,
)
```

**Step 2: Run full test suite**

Run: `.venv/bin/python -m pytest tests/ -v`
Expected: All tests PASS.

**Step 3: Commit**

```bash
git add src/__init__.py
git commit -m "export hedging Hamiltonian functions from src package"
```
