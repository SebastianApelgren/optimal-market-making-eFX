from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

import numpy as np

from .model import ModelParams, canon_pair
from .hamiltonian import H_logistic


def validate_pde_grid(y_grids: List[np.ndarray], mp: ModelParams) -> List[float]:
    """Check that each axis is uniformly spaced and all trade sizes are exact
    multiples of each relevant grid spacing.

    Returns list of dy per axis.  Raises ValueError on failure.
    """
    d = len(mp.currencies)
    if len(y_grids) != d:
        raise ValueError(
            f"Expected {d} grid axes, got {len(y_grids)}"
        )

    dy_list: List[float] = []
    for k, g in enumerate(y_grids):
        if len(g) < 2:
            raise ValueError(f"Axis {k} must have at least 2 points")
        diffs = np.diff(g)
        dy = diffs[0]
        if not np.allclose(diffs, dy, atol=1e-12):
            raise ValueError(
                f"Axis {k} ({mp.currencies[k]}) is not uniformly spaced"
            )
        dy_list.append(float(dy))

    # Check that all trade sizes are exact multiples of the relevant grid spacings.
    ccy = mp.currencies
    for i in range(d):
        for j in range(d):
            if i == j:
                continue
            key = canon_pair(ccy[i], ccy[j])
            if key not in mp.pairs:
                continue
            pp = mp.pairs[key]
            for z in pp.sizes_musd:
                for axis, sign in [(i, +1), (j, -1)]:
                    dy = dy_list[axis]
                    ratio = z / dy
                    if abs(ratio - round(ratio)) > 1e-9:
                        raise ValueError(
                            f"Trade size {z} M$ is not an exact multiple of "
                            f"grid spacing dy={dy} on axis {axis} ({ccy[axis]})"
                        )

    return dy_list


def _compute_shift_slices(
    grid_shape: Tuple[int, ...],
    shift_indices: Tuple[int, ...],
) -> Tuple[tuple, tuple]:
    """Compute array slices for a shifted lookup.

    For each axis k with grid size n and shift s:
      s >= 0: src = [s, n), dst = [0, n-s)
      s < 0:  src = [0, n+s), dst = [-s, n)

    Returns (src_slices, dst_slices) — both select sub-arrays of the same shape.
    """
    src_slices = []
    dst_slices = []
    for k, s in enumerate(shift_indices):
        n = grid_shape[k]
        if s >= 0:
            src_slices.append(slice(s, n))
            dst_slices.append(slice(0, n - s))
        else:
            src_slices.append(slice(0, n + s))
            dst_slices.append(slice(-s, n))
    return tuple(src_slices), tuple(dst_slices)


@dataclass(frozen=True)
class QuotingSpec:
    """Precomputed specification for the quoting Hamiltonian integral."""
    contributions: tuple  # tuple of (z, lam, alpha, beta, src_slices, dst_slices)
    d: int                # number of currencies


def build_quoting_spec(
    y_grids: List[np.ndarray],
    mp: ModelParams,
) -> QuotingSpec:
    """Precompute all (pair, tier, size) contributions with shift slices.

    Loops over directed pairs (i,j) with i != j, canonical pair lookup, tiers,
    and sizes — the same structure as build_M_tildeM_P in riccati.py.
    """
    dy_list = validate_pde_grid(y_grids, mp)
    grid_shape = tuple(len(g) for g in y_grids)
    ccy = mp.currencies
    d = len(ccy)

    contributions = []

    for i in range(d):
        for j in range(d):
            if i == j:
                continue
            key = canon_pair(ccy[i], ccy[j])
            if key not in mp.pairs:
                continue

            pp = mp.pairs[key]
            z_arr = pp.sizes_musd
            lam_arr = pp.lambdas_per_day

            for size_idx in range(len(z_arr)):
                z = float(z_arr[size_idx])
                lam = float(lam_arr[size_idx])

                # Direction vector: d_vec[i] = +1, d_vec[j] = -1.
                # Shift on grid: y + z * d_vec => axis i shifts by +z, axis j by -z.
                shift_indices = [0] * d
                shift_indices[i] = round(z / dy_list[i])
                shift_indices[j] = round(-z / dy_list[j])
                shift_indices = tuple(shift_indices)

                src_sl, dst_sl = _compute_shift_slices(grid_shape, shift_indices)

                for tier in pp.tiers:
                    contributions.append(
                        (z, lam, tier.alpha, tier.beta, src_sl, dst_sl)
                    )

    return QuotingSpec(contributions=tuple(contributions), d=d)


def quoting_hamiltonian_integral(
    theta: np.ndarray,
    spec: QuotingSpec,
) -> np.ndarray:
    """Evaluate the quoting Hamiltonian integral Q(y) at every grid point.

    Q(y) = sum_{n,i,j} sum_k  z_k * lambda_k * H(p_k(y), alpha_n, beta_n)

    where p_k(y) = (theta(y) - theta(y + z_k*(e_i - e_j))) / z_k.

    Points near grid edges where the shifted lookup falls outside the grid
    receive H = 0 for those out-of-bounds contributions.
    """
    Q = np.zeros_like(theta)
    slice_cache = {}

    def _slice_key(sl_tuple):
        return tuple((s.start, s.stop) for s in sl_tuple)

    for (z, lam, alpha, beta, src_sl, dst_sl) in spec.contributions:
        cache_key = (_slice_key(src_sl), _slice_key(dst_sl))
        if cache_key not in slice_cache:
            slice_cache[cache_key] = (theta[src_sl], theta[dst_sl])
        theta_shifted, theta_here = slice_cache[cache_key]

        p = (theta_here - theta_shifted) / z
        H_vals, _, _ = H_logistic(p, alpha, beta)
        Q[dst_sl] += z * lam * H_vals

    return Q


def H_execution_cost(p, psi: float, eta: float):
    """Hedging Hamiltonian H(p) = (max(|p| - psi, 0))^2 / (4*eta).

    Legendre-Fenchel transform of L(xi) = psi*|xi| + eta*xi^2.
    Returns 0 in the dead zone |p| <= psi.
    """
    p = np.asarray(p, dtype=float)
    excess = np.maximum(np.abs(p) - psi, 0.0)
    return excess * excess / (4.0 * eta)


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


def _second_deriv_diagonal(theta: np.ndarray, axis: int, dy: float) -> np.ndarray:
    """Compute d²θ/dy_i² via central differences (interior) and one-sided at boundaries."""
    d = theta.ndim
    n = theta.shape[axis]
    out = np.empty_like(theta)
    dy2 = dy * dy

    # Central: (θ_{k+1} - 2θ_k + θ_{k-1}) / Δy²
    sl_plus = [slice(None)] * d
    sl_mid = [slice(None)] * d
    sl_minus = [slice(None)] * d
    sl_plus[axis] = slice(2, n)
    sl_mid[axis] = slice(1, n - 1)
    sl_minus[axis] = slice(0, n - 2)
    out[tuple(sl_mid)] = (theta[tuple(sl_plus)] - 2.0 * theta[tuple(sl_mid)] + theta[tuple(sl_minus)]) / dy2

    # Left boundary: (θ_0 - 2θ_1 + θ_2) / Δy²
    sl_0 = [slice(None)] * d; sl_0[axis] = 0
    sl_1 = [slice(None)] * d; sl_1[axis] = 1
    sl_2 = [slice(None)] * d; sl_2[axis] = 2
    out[tuple(sl_0)] = (theta[tuple(sl_0)] - 2.0 * theta[tuple(sl_1)] + theta[tuple(sl_2)]) / dy2

    # Right boundary: (θ_{n-3} - 2θ_{n-2} + θ_{n-1}) / Δy²
    sl_n1 = [slice(None)] * d; sl_n1[axis] = n - 1
    sl_n2 = [slice(None)] * d; sl_n2[axis] = n - 2
    sl_n3 = [slice(None)] * d; sl_n3[axis] = n - 3
    out[tuple(sl_n1)] = (theta[tuple(sl_n3)] - 2.0 * theta[tuple(sl_n2)] + theta[tuple(sl_n1)]) / dy2

    return out


def _second_deriv_cross(theta: np.ndarray, axis_i: int, axis_j: int,
                        dy_i: float, dy_j: float) -> np.ndarray:
    """Compute d²θ/(dy_i dy_j) via central differences (interior), one-sided at boundaries."""
    d = theta.ndim
    n_i = theta.shape[axis_i]
    n_j = theta.shape[axis_j]
    out = np.empty_like(theta)
    dydy4 = 4.0 * dy_i * dy_j

    # Interior: central cross stencil
    # (θ_{i+1,j+1} - θ_{i+1,j-1} - θ_{i-1,j+1} + θ_{i-1,j-1}) / (4 Δy_i Δy_j)
    def _make_sl(si, sj):
        sl = [slice(None)] * d
        sl[axis_i] = si
        sl[axis_j] = sj
        return tuple(sl)

    interior_i = slice(1, n_i - 1)
    interior_j = slice(1, n_j - 1)

    out[_make_sl(interior_i, interior_j)] = (
        theta[_make_sl(slice(2, n_i), slice(2, n_j))]
        - theta[_make_sl(slice(2, n_i), slice(0, n_j - 2))]
        - theta[_make_sl(slice(0, n_i - 2), slice(2, n_j))]
        + theta[_make_sl(slice(0, n_i - 2), slice(0, n_j - 2))]
    ) / dydy4

    # Boundaries: use forward/backward differences to approximate d/dy_i and d/dy_j,
    # then combine.  We compute the full cross-derivative via
    #   d²θ/(dy_i dy_j) ≈ d/dy_i [dθ/dy_j]
    # using one-sided first derivatives at boundary rows/columns.

    # Helper: first derivative along one axis (forward or backward)
    def _d1(arr, ax, dh, side):
        """First derivative of arr along ax using forward (side=0) or backward (side=-1)."""
        s0 = [slice(None)] * d
        s1 = [slice(None)] * d
        if side == 0:  # forward
            s0[ax] = 0
            s1[ax] = 1
        else:  # backward
            s0[ax] = -1
            s1[ax] = -2
        if side == 0:
            return (arr[tuple(s1)] - arr[tuple(s0)]) / dh
        else:
            return (arr[tuple(s0)] - arr[tuple(s1)]) / dh

    # First compute dθ/dy_j everywhere (central interior, one-sided boundary)
    dtheta_dj = np.empty_like(theta)
    sl_p = [slice(None)] * d; sl_p[axis_j] = slice(2, n_j)
    sl_m = [slice(None)] * d; sl_m[axis_j] = slice(0, n_j - 2)
    sl_c = [slice(None)] * d; sl_c[axis_j] = slice(1, n_j - 1)
    dtheta_dj[tuple(sl_c)] = (theta[tuple(sl_p)] - theta[tuple(sl_m)]) / (2.0 * dy_j)
    sl_0j = [slice(None)] * d; sl_0j[axis_j] = 0
    sl_1j = [slice(None)] * d; sl_1j[axis_j] = 1
    dtheta_dj[tuple(sl_0j)] = (theta[tuple(sl_1j)] - theta[tuple(sl_0j)]) / dy_j
    sl_lj = [slice(None)] * d; sl_lj[axis_j] = n_j - 1
    sl_pj = [slice(None)] * d; sl_pj[axis_j] = n_j - 2
    dtheta_dj[tuple(sl_lj)] = (theta[tuple(sl_lj)] - theta[tuple(sl_pj)]) / dy_j

    # Now differentiate dtheta_dj along axis_i for boundary rows of axis_i
    # Left boundary of axis i (k_i = 0): forward difference
    sl_i0 = [slice(None)] * d; sl_i0[axis_i] = 0
    sl_i1 = [slice(None)] * d; sl_i1[axis_i] = 1
    out[tuple(sl_i0)] = (dtheta_dj[tuple(sl_i1)] - dtheta_dj[tuple(sl_i0)]) / dy_i

    # Right boundary of axis i (k_i = n_i-1): backward difference
    sl_il = [slice(None)] * d; sl_il[axis_i] = n_i - 1
    sl_ip = [slice(None)] * d; sl_ip[axis_i] = n_i - 2
    out[tuple(sl_il)] = (dtheta_dj[tuple(sl_il)] - dtheta_dj[tuple(sl_ip)]) / dy_i

    # Left/right boundary of axis j (but interior of axis i): already covered by
    # the dtheta_dj one-sided stencil + central diff in axis_i from the interior block.
    # However, the interior block only filled interior_i x interior_j.
    # Fill the boundary columns of axis j for interior rows of axis i.
    sl_j0_ii = [slice(None)] * d; sl_j0_ii[axis_i] = interior_i; sl_j0_ii[axis_j] = 0
    sl_j1_ii = [slice(None)] * d; sl_j1_ii[axis_i] = interior_i; sl_j1_ii[axis_j] = 1

    # d²θ/(dy_i dy_j) at j=0: central in i, forward in j
    sl_ip1_j0 = [slice(None)] * d; sl_ip1_j0[axis_i] = slice(2, n_i); sl_ip1_j0[axis_j] = 0
    sl_ip1_j1 = [slice(None)] * d; sl_ip1_j1[axis_i] = slice(2, n_i); sl_ip1_j1[axis_j] = 1
    sl_im1_j0 = [slice(None)] * d; sl_im1_j0[axis_i] = slice(0, n_i - 2); sl_im1_j0[axis_j] = 0
    sl_im1_j1 = [slice(None)] * d; sl_im1_j1[axis_i] = slice(0, n_i - 2); sl_im1_j1[axis_j] = 1

    out[tuple(sl_j0_ii)] = (
        (theta[tuple(sl_ip1_j1)] - theta[tuple(sl_ip1_j0)]) / dy_j
        - (theta[tuple(sl_im1_j1)] - theta[tuple(sl_im1_j0)]) / dy_j
    ) / (2.0 * dy_i)

    sl_jl_ii = [slice(None)] * d; sl_jl_ii[axis_i] = interior_i; sl_jl_ii[axis_j] = n_j - 1
    sl_jp_ii = [slice(None)] * d; sl_jp_ii[axis_i] = interior_i; sl_jp_ii[axis_j] = n_j - 2
    sl_ip1_jl = [slice(None)] * d; sl_ip1_jl[axis_i] = slice(2, n_i); sl_ip1_jl[axis_j] = n_j - 1
    sl_ip1_jp = [slice(None)] * d; sl_ip1_jp[axis_i] = slice(2, n_i); sl_ip1_jp[axis_j] = n_j - 2
    sl_im1_jl = [slice(None)] * d; sl_im1_jl[axis_i] = slice(0, n_i - 2); sl_im1_jl[axis_j] = n_j - 1
    sl_im1_jp = [slice(None)] * d; sl_im1_jp[axis_i] = slice(0, n_i - 2); sl_im1_jp[axis_j] = n_j - 2

    out[tuple(sl_jl_ii)] = (
        (theta[tuple(sl_ip1_jl)] - theta[tuple(sl_ip1_jp)]) / dy_j
        - (theta[tuple(sl_im1_jl)] - theta[tuple(sl_im1_jp)]) / dy_j
    ) / (2.0 * dy_i)

    return out


def diffusion_term(
    theta: np.ndarray,
    y_grids: List[np.ndarray],
    dy_list: List[float],
    Sigma: np.ndarray,
) -> np.ndarray:
    """Compute ½ Tr(D(y) Sigma D(y) D²_yy theta) at every grid point.

    Parameters
    ----------
    theta : d-dimensional array on the inventory grid.
    y_grids : list of d 1D arrays (grid axes).
    dy_list : list of d floats (grid spacing per axis).
    Sigma : (d, d) covariance matrix.

    Returns
    -------
    Array same shape as theta.
    """
    d = theta.ndim
    result = np.zeros_like(theta)

    # Precompute broadcasted y arrays
    y_bc = []
    for k in range(d):
        slices = [np.newaxis] * d
        slices[k] = slice(None)
        y_bc.append(y_grids[k][tuple(slices)])

    for i in range(d):
        for j in range(i, d):
            sij = Sigma[i, j]
            if abs(sij) < 1e-30:
                continue

            if i == j:
                d2 = _second_deriv_diagonal(theta, i, dy_list[i])
                result += 0.5 * sij * y_bc[i] * y_bc[i] * d2
            else:
                d2 = _second_deriv_cross(theta, i, j, dy_list[i], dy_list[j])
                # Factor 1: off-diagonal appears twice (i,j) and (j,i) in the
                # full double sum, but we loop i<j once, so weight is 2 * ½ = 1.
                result += sij * y_bc[i] * y_bc[j] * d2

    return result


def drift_term(
    grad_theta: List[np.ndarray],
    y_grids: List[np.ndarray],
    mu_vec: np.ndarray,
) -> np.ndarray:
    """Compute y^T mu + y^T D(mu) nabla_y theta = sum_i mu_i y_i (1 + dtheta/dy_i).

    Parameters
    ----------
    grad_theta : list of d arrays, each same shape as the grid.
        Entry k is d(theta)/d(y_k), e.g. from compute_gradient.
    y_grids : list of d 1D arrays (grid axes).
    mu_vec : length-d array of drift coefficients.

    Returns
    -------
    Array same shape as grad_theta[0].
    """
    d = len(y_grids)
    result = np.zeros_like(grad_theta[0])

    for k in range(d):
        if abs(mu_vec[k]) < 1e-30:
            continue
        slices = [np.newaxis] * d
        slices[k] = slice(None)
        y_k = y_grids[k][tuple(slices)]
        result += mu_vec[k] * y_k * (1.0 + grad_theta[k])

    return result


def running_penalty(
    y_grids: List[np.ndarray],
    Sigma: np.ndarray,
    gamma: float,
) -> np.ndarray:
    """Compute -gamma/2 * y^T Sigma y at every grid point.

    This is a known function of y (no theta dependence) and only needs
    to be computed once per grid.

    Parameters
    ----------
    y_grids : list of d 1D arrays (grid axes).
    Sigma : (d, d) covariance matrix.
    gamma : risk aversion parameter.

    Returns
    -------
    d-dimensional array on the grid.
    """
    d = len(y_grids)
    shape = tuple(len(g) for g in y_grids)
    result = np.zeros(shape)

    y_bc = []
    for k in range(d):
        slices = [np.newaxis] * d
        slices[k] = slice(None)
        y_bc.append(y_grids[k][tuple(slices)])

    for i in range(d):
        for j in range(i, d):
            sij = Sigma[i, j]
            if abs(sij) < 1e-30:
                continue
            if i == j:
                result += sij * y_bc[i] * y_bc[i]
            else:
                result += 2.0 * sij * y_bc[i] * y_bc[j]

    return -0.5 * gamma * result
