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


def terminal_condition(
    y_grids: List[np.ndarray],
    kappa: np.ndarray,
) -> np.ndarray:
    """Compute θ(T, y) = -ℓ(y) = -y^T κ y at every grid point.

    Parameters
    ----------
    y_grids : list of d 1D arrays (grid axes).
    kappa : (d, d) positive semi-definite symmetric matrix (terminal penalty).

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
            kij = kappa[i, j]
            if abs(kij) < 1e-30:
                continue
            if i == j:
                result += kij * y_bc[i] * y_bc[i]
            else:
                result += 2.0 * kij * y_bc[i] * y_bc[j]

    return -result


def _diffusion_term_cross(
    theta: np.ndarray,
    y_grids: List[np.ndarray],
    dy_list: List[float],
    Sigma: np.ndarray,
) -> np.ndarray:
    """Compute only the cross-derivative part of the diffusion term.

    Returns Σ_{i<j} Σ_{ij} y_i y_j ∂²θ/(∂y_i ∂y_j).
    For d=2 with σ_ref=0, this is zero.
    """
    d = theta.ndim
    result = np.zeros_like(theta)

    y_bc = []
    for k in range(d):
        slices = [np.newaxis] * d
        slices[k] = slice(None)
        y_bc.append(y_grids[k][tuple(slices)])

    for i in range(d):
        for j in range(i + 1, d):
            sij = Sigma[i, j]
            if abs(sij) < 1e-30:
                continue
            d2 = _second_deriv_cross(theta, i, j, dy_list[i], dy_list[j])
            result += sij * y_bc[i] * y_bc[j] * d2

    return result


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


def pde_rhs(theta: np.ndarray, spec: PDESpec) -> np.ndarray:
    """Evaluate the spatial RHS of PDE (1) at every grid point.

    Returns RHS such that dθ/dt = -RHS, i.e. the backward Euler update is
    θ^{m-1} = θ^m + Δt * RHS.
    """
    grad = compute_gradient(theta, spec.dy_list)

    rhs = spec.penalty.copy()                                             # -γ/2 y^T Σ y
    rhs += diffusion_term(theta, spec.y_grids, spec.dy_list, spec.Sigma)  # ½ Tr(...)
    rhs += drift_term(grad, spec.y_grids, spec.mu_vec)                    # y^T μ (1+∂θ)
    rhs += quoting_hamiltonian_integral(theta, spec.quoting)              # Q(y)
    rhs += hedging_hamiltonian(grad, spec.y_grids, spec.hedging)          # H_hedge(y)

    return rhs


def pde_rhs_nonlinear(theta: np.ndarray, spec: PDESpec) -> np.ndarray:
    """Evaluate the explicit (nonlinear) part of the spatial RHS.

    Includes: penalty, drift, quoting & hedging Hamiltonians, cross diffusion.
    Excludes: diagonal diffusion along each axis (handled implicitly).

    The full drift y^T μ (1 + ∂θ/∂y) is treated explicitly here. For μ=0
    (the paper's example) this is zero. For nonzero μ, moving the drift
    gradient into the implicit operator would improve stability.
    """
    grad = compute_gradient(theta, spec.dy_list)

    rhs = spec.penalty.copy()                                             # -γ/2 y^T Σ y
    rhs += _diffusion_term_cross(theta, spec.y_grids, spec.dy_list, spec.Sigma)
    rhs += drift_term(grad, spec.y_grids, spec.mu_vec)                    # y^T μ (1+∂θ)
    rhs += quoting_hamiltonian_integral(theta, spec.quoting)              # Q(y)
    rhs += hedging_hamiltonian(grad, spec.y_grids, spec.hedging)          # H_hedge(y)

    return rhs


# ---------------------------------------------------------------------------
# Thomas algorithm for tridiagonal systems
# ---------------------------------------------------------------------------

def _thomas_factor(
    lower: np.ndarray,
    diag: np.ndarray,
    upper: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Precompute forward-sweep factors for the Thomas algorithm.

    Parameters
    ----------
    lower, diag, upper : 1D arrays of length n.
        lower[0] is unused (no x[-1]).

    Returns
    -------
    c_prime : modified upper diagonal (length n).
    denom_inv : 1/denominator at each row (length n).
    """
    n = len(diag)
    c_prime = np.empty(n)
    denom_inv = np.empty(n)

    denom_inv[0] = 1.0 / diag[0]
    c_prime[0] = upper[0] * denom_inv[0]

    for i in range(1, n):
        denom = diag[i] - lower[i] * c_prime[i - 1]
        denom_inv[i] = 1.0 / denom
        c_prime[i] = upper[i] * denom_inv[i]

    return c_prime, denom_inv


def _thomas_solve_batch(
    lower: np.ndarray,
    c_prime: np.ndarray,
    denom_inv: np.ndarray,
    rhs: np.ndarray,
) -> np.ndarray:
    """Solve tridiagonal system for multiple RHS vectors.

    Parameters
    ----------
    lower : 1D, length n (lower[0] unused).
    c_prime, denom_inv : from _thomas_factor, length n.
    rhs : (batch, n) array — each row is an independent RHS.

    Returns
    -------
    x : (batch, n) solution array.
    """
    n = len(c_prime)

    # Forward sweep: compute modified RHS d'
    d_prime = np.empty_like(rhs)
    d_prime[:, 0] = rhs[:, 0] * denom_inv[0]
    for i in range(1, n):
        d_prime[:, i] = (rhs[:, i] - lower[i] * d_prime[:, i - 1]) * denom_inv[i]

    # Back substitution
    x = np.empty_like(rhs)
    x[:, -1] = d_prime[:, -1]
    for i in range(n - 2, -1, -1):
        x[:, i] = d_prime[:, i] - c_prime[i] * x[:, i + 1]

    return x


# ---------------------------------------------------------------------------
# Implicit diffusion operator
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class _ImplicitAxisOp:
    """Precomputed Thomas factors for (I - dt*L_k) along one grid axis."""
    axis: int
    lower: np.ndarray       # length n, lower[0] unused
    c_prime: np.ndarray     # from _thomas_factor
    denom_inv: np.ndarray   # from _thomas_factor


def _build_implicit_ops(
    spec: PDESpec,
    dt: float,
) -> Tuple['_ImplicitAxisOp', ...]:
    """Build tridiagonal implicit operators for each axis with nonzero diffusion.

    The operator along axis k discretises:
        L_k[θ]_j = ½ Σ_{kk} y_j² (θ_{j+1} - 2θ_j + θ_{j-1}) / Δy²

    and returns the Thomas factors for (I - dt * L_k).
    Boundary points (j=0, j=n-1) are identity rows.
    """
    ops = []
    d = len(spec.y_grids)

    for k in range(d):
        sigma_kk = spec.Sigma[k, k]
        if sigma_kk < 1e-30:
            continue

        y = spec.y_grids[k]
        dy = spec.dy_list[k]
        n = len(y)
        dy2 = dy * dy

        lower = np.zeros(n)
        diag = np.ones(n)
        upper = np.zeros(n)

        # Interior points j = 1, ..., n-2: central differences
        y_j = y[1:-1]
        diff_c = 0.5 * sigma_kk * y_j * y_j / dy2  # ½ Σ_kk y_j² / Δy²

        # (I - dt * L_k):
        #   lower[j] = -dt * diff_c   (coeff of θ_{j-1})
        #   diag[j]  = 1 + 2*dt * diff_c  (coeff of θ_j)
        #   upper[j] = -dt * diff_c   (coeff of θ_{j+1})
        lower[1:-1] = -dt * diff_c
        diag[1:-1] = 1.0 + 2.0 * dt * diff_c
        upper[1:-1] = -dt * diff_c

        c_prime, denom_inv = _thomas_factor(lower, diag, upper)
        ops.append(_ImplicitAxisOp(
            axis=k, lower=lower, c_prime=c_prime, denom_inv=denom_inv,
        ))

    return tuple(ops)


def _apply_implicit_step(
    theta: np.ndarray,
    ops: Tuple[_ImplicitAxisOp, ...],
) -> np.ndarray:
    """Solve (I - dt*L) θ_new = θ by sequential tridiagonal solves per axis."""
    result = theta
    for op in ops:
        # Move target axis to last position, flatten batch dims
        moved = np.moveaxis(result, op.axis, -1)
        shape = moved.shape
        n = shape[-1]
        rhs_2d = moved.reshape(-1, n)

        sol_2d = _thomas_solve_batch(op.lower, op.c_prime, op.denom_inv, rhs_2d)

        result = np.moveaxis(sol_2d.reshape(shape), -1, op.axis)
    return result


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
    from tqdm import trange

    for m in trange(n_steps, desc="HJB PDE solve", unit="step"):
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


def solve_hjb_semi_implicit(
    y_grids: List[np.ndarray],
    mp: ModelParams,
    n_steps: int = 1000,
    snapshot_times: List[float] | None = None,
) -> dict:
    """Solve PDE (1) backward from T to 0 via IMEX Euler.

    Diagonal diffusion is treated implicitly (tridiagonal solve per axis).
    All other terms (quoting, hedging, drift, penalty, cross diffusion) are
    treated explicitly.

    Parameters and return value are identical to solve_hjb_explicit.
    """
    d = len(mp.currencies)
    kappa = mp.kappa if mp.kappa is not None else np.zeros((d, d))
    T = mp.T_days
    dt = T / n_steps

    spec = build_pde_spec(y_grids, mp)
    theta = terminal_condition(y_grids, kappa)

    # Precompute implicit operator (Thomas factors, done once)
    implicit_ops = _build_implicit_ops(spec, dt)

    # Prepare snapshot collection
    if snapshot_times is None:
        snapshot_times = []
    snap_steps = {}
    for t_snap in snapshot_times:
        step_idx = round((T - t_snap) / dt)
        step_idx = max(0, min(n_steps, step_idx))
        snap_steps[step_idx] = t_snap
    snapshots = {}

    # Backward march in τ = T - t
    from tqdm import trange

    for m in trange(n_steps, desc="HJB PDE solve (IMEX)", unit="step"):
        # Explicit step: b = θ^m + Δt * N[θ^m]
        rhs_N = pde_rhs_nonlinear(theta, spec)
        b = theta + dt * rhs_N

        # Implicit step: solve (I - Δt L) θ^{m+1} = b
        theta = _apply_implicit_step(b, implicit_ops)

        if (m + 1) in snap_steps:
            t_snap = snap_steps[m + 1]
            snapshots[t_snap] = theta.copy()

    return {
        'theta_0': theta,
        'snapshots': snapshots,
        'spec': spec,
        'dt': dt,
    }


# ---------------------------------------------------------------------------
# Fully implicit Euler with policy iteration (Howard's algorithm)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class QuotingControl:
    """Optimal quoting control for one contribution."""
    delta_star: np.ndarray  # on overlap sub-grid
    f_star: np.ndarray      # on overlap sub-grid
    z: float
    lam: float
    src_sl: tuple
    dst_sl: tuple


@dataclass(frozen=True)
class HedgingControl:
    """Optimal hedging control for one canonical pair."""
    xi_star: np.ndarray  # on full grid
    psi: float
    eta: float
    i: int
    j: int
    k_i: float
    k_j: float


def extract_quoting_controls(
    theta: np.ndarray,
    spec: PDESpec,
) -> List[QuotingControl]:
    """Compute optimal markup and fill probability for each quoting contribution.

    For each (pair, tier, size) contribution, evaluates the momentum
    p(y) = (θ(y) - θ(y + z·d)) / z and computes δ*(y), f*(y) via Lambert W.
    """
    controls = []
    for (z, lam, alpha, beta, src_sl, dst_sl) in spec.quoting.contributions:
        p = (theta[dst_sl] - theta[src_sl]) / z
        _, delta_star, f_star = H_logistic(p, alpha, beta)
        controls.append(QuotingControl(
            delta_star=delta_star, f_star=f_star,
            z=z, lam=lam, src_sl=src_sl, dst_sl=dst_sl,
        ))
    return controls


def extract_hedging_controls(
    grad_theta: List[np.ndarray],
    y_grids: List[np.ndarray],
    spec: PDESpec,
) -> List[HedgingControl]:
    """Compute optimal hedge rate for each hedging pair.

    For each canonical pair (i, j), evaluates the momentum
    p = (1+k_i y_i)∇θ_i - (1+k_j y_j)∇θ_j + k_i y_i - k_j y_j
    and computes ξ* = sign(p)·max(|p|-ψ, 0) / (2η).
    """
    d = spec.hedging.d
    y_bc = []
    for k in range(d):
        sl = [np.newaxis] * d
        sl[k] = slice(None)
        y_bc.append(y_grids[k][tuple(sl)])

    controls = []
    for (i, j, psi, eta, k_i, k_j) in spec.hedging.contributions:
        p = (grad_theta[i] - grad_theta[j]
             + k_i * y_bc[i] * (1.0 + grad_theta[i])
             - k_j * y_bc[j] * (1.0 + grad_theta[j]))
        excess = np.maximum(np.abs(p) - psi, 0.0)
        xi_star = np.sign(p) * excess / (2.0 * eta)
        controls.append(HedgingControl(
            xi_star=xi_star, psi=psi, eta=eta,
            i=i, j=j, k_i=k_i, k_j=k_j,
        ))
    return controls


def assemble_implicit_system(
    quoting_controls: List[QuotingControl],
    hedging_controls: List[HedgingControl],
    spec: PDESpec,
    dt: float,
) -> Tuple:
    """Build sparse matrix (I - dt·A) and source vector for the linearized system.

    With fixed controls (δ*, ξ*), the HJB becomes the linear system:
        (I - dt·A) θ^{k+1} = θ^m + dt·source

    Returns (A_csr, source_flat) where A_csr is a scipy CSR matrix of shape
    (N, N) and source_flat is a 1D array of length N.
    """
    import scipy.sparse

    grid_shape = tuple(len(g) for g in spec.y_grids)
    d = len(grid_shape)
    N = int(np.prod(grid_shape))
    strides = np.array([int(np.prod(grid_shape[k+1:])) for k in range(d)])
    flat_idx = np.arange(N).reshape(grid_shape)

    rows_list: List[np.ndarray] = []
    cols_list: List[np.ndarray] = []
    vals_list: List[np.ndarray] = []
    source = np.zeros(N)

    # --- Identity diagonal ---
    all_idx = np.arange(N)
    rows_list.append(all_idx)
    cols_list.append(all_idx)
    vals_list.append(np.ones(N))

    # --- Penalty source (static, no θ dependence) — all points incl. boundary ---
    source += dt * spec.penalty.ravel()

    # Interior mask: True at points where all axes are away from the boundary.
    # Non-penalty source terms are restricted to interior points (Option A).
    interior = np.ones(grid_shape, dtype=bool)
    for k in range(d):
        sl = [slice(None)] * d
        sl[k] = 0
        interior[tuple(sl)] = False
        sl = [slice(None)] * d
        sl[k] = grid_shape[k] - 1
        interior[tuple(sl)] = False
    interior_flat = interior.ravel()

    # --- Drift source: Σ_k μ_k y_k (interior only) ---
    for k in range(d):
        if abs(spec.mu_vec[k]) < 1e-30:
            continue
        sl = [np.newaxis] * d
        sl[k] = slice(None)
        y_k = spec.y_grids[k][tuple(sl)] * np.ones(grid_shape)
        source[interior_flat] += dt * spec.mu_vec[k] * y_k.ravel()[interior_flat]

    # --- Diagonal diffusion: ½ Σ_kk y_k² ∂²θ/∂y_k² ---
    for k in range(d):
        sigma_kk = spec.Sigma[k, k]
        if sigma_kk < 1e-30:
            continue
        dy2 = spec.dy_list[k] ** 2
        n_k = grid_shape[k]

        # Interior points along axis k
        int_sl = [slice(None)] * d
        int_sl[k] = slice(1, n_k - 1)
        idx_c = flat_idx[tuple(int_sl)].ravel()

        # Coefficient: 0.5 * Σ_kk * y_k² / Δy²
        y_k_int = spec.y_grids[k][1:-1]
        bc_shape = [grid_shape[l] if l != k else n_k - 2 for l in range(d)]
        sl = [np.newaxis] * d
        sl[k] = slice(None)
        coeff = (0.5 * sigma_kk / dy2) * (
            y_k_int[tuple(sl)] ** 2 * np.ones(bc_shape)
        ).ravel()

        # (I - dt·A) entries: diag += 2·dt·coeff, off-diag = -dt·coeff
        rows_list.append(idx_c)
        cols_list.append(idx_c)
        vals_list.append(2.0 * dt * coeff)

        rows_list.append(idx_c)
        cols_list.append(idx_c + strides[k])
        vals_list.append(-dt * coeff)

        rows_list.append(idx_c)
        cols_list.append(idx_c - strides[k])
        vals_list.append(-dt * coeff)

    # --- Cross diffusion: Σ_{ij} y_i y_j ∂²θ/(∂y_i ∂y_j) for i < j ---
    for ki in range(d):
        for kj in range(ki + 1, d):
            sij = spec.Sigma[ki, kj]
            if abs(sij) < 1e-30:
                continue
            n_i, n_j = grid_shape[ki], grid_shape[kj]
            int_sl = [slice(None)] * d
            int_sl[ki] = slice(1, n_i - 1)
            int_sl[kj] = slice(1, n_j - 1)
            idx_c = flat_idx[tuple(int_sl)].ravel()

            bc_shape = [
                grid_shape[l] if l != ki and l != kj
                else (n_i - 2 if l == ki else n_j - 2)
                for l in range(d)
            ]
            sl_i = [np.newaxis] * d
            sl_i[ki] = slice(None)
            sl_j = [np.newaxis] * d
            sl_j[kj] = slice(None)
            coeff = (sij / (4.0 * spec.dy_list[ki] * spec.dy_list[kj])) * (
                spec.y_grids[ki][1:-1][tuple(sl_i)]
                * spec.y_grids[kj][1:-1][tuple(sl_j)]
                * np.ones(bc_shape)
            ).ravel()

            si, sj = strides[ki], strides[kj]
            rows_list.append(idx_c)
            cols_list.append(idx_c + si + sj)
            vals_list.append(-dt * coeff)
            rows_list.append(idx_c)
            cols_list.append(idx_c + si - sj)
            vals_list.append(dt * coeff)
            rows_list.append(idx_c)
            cols_list.append(idx_c - si + sj)
            vals_list.append(dt * coeff)
            rows_list.append(idx_c)
            cols_list.append(idx_c - si - sj)
            vals_list.append(-dt * coeff)

    # --- Drift gradient (linear: μ_k y_k ∂θ/∂y_k) ---
    for k in range(d):
        mu_k = spec.mu_vec[k]
        if abs(mu_k) < 1e-30:
            continue
        n_k = grid_shape[k]
        int_sl = [slice(None)] * d
        int_sl[k] = slice(1, n_k - 1)
        idx_c = flat_idx[tuple(int_sl)].ravel()

        y_k_int = spec.y_grids[k][1:-1]
        bc_shape = [grid_shape[l] if l != k else n_k - 2 for l in range(d)]
        sl = [np.newaxis] * d
        sl[k] = slice(None)
        coeff = (mu_k / (2.0 * spec.dy_list[k])) * (
            y_k_int[tuple(sl)] * np.ones(bc_shape)
        ).ravel()

        rows_list.append(idx_c)
        cols_list.append(idx_c + strides[k])
        vals_list.append(-dt * coeff)
        rows_list.append(idx_c)
        cols_list.append(idx_c - strides[k])
        vals_list.append(dt * coeff)

    # --- Quoting (fixed controls) ---
    for qc in quoting_controls:
        idx_here = flat_idx[qc.dst_sl].ravel()
        idx_shifted = flat_idx[qc.src_sl].ravel()
        lam_f = qc.lam * qc.f_star.ravel()

        # (I - dt·A): diagonal += dt·λ·f*, off-diag = -dt·λ·f*
        rows_list.append(idx_here)
        cols_list.append(idx_here)
        vals_list.append(dt * lam_f)

        rows_list.append(idx_here)
        cols_list.append(idx_shifted)
        vals_list.append(-dt * lam_f)

        # Source: dt · z · λ · f* · δ* (interior only)
        contrib = dt * qc.z * qc.lam * qc.f_star.ravel() * qc.delta_star.ravel()
        mask = interior_flat[idx_here]
        np.add.at(source, idx_here[mask], contrib[mask])

    # --- Hedging (fixed controls) ---
    for hc in hedging_controls:
        # Source (interior only): ξ*(k_i y_i - k_j y_j) - ψ|ξ*| - η(ξ*)²
        xi_flat = hc.xi_star.ravel()
        sl_i = [np.newaxis] * d
        sl_i[hc.i] = slice(None)
        sl_j = [np.newaxis] * d
        sl_j[hc.j] = slice(None)
        y_i_full = (spec.y_grids[hc.i][tuple(sl_i)] * np.ones(grid_shape)).ravel()
        y_j_full = (spec.y_grids[hc.j][tuple(sl_j)] * np.ones(grid_shape)).ravel()

        hedge_src = dt * (
            xi_flat * (hc.k_i * y_i_full - hc.k_j * y_j_full)
            - hc.psi * np.abs(xi_flat)
            - hc.eta * xi_flat * xi_flat
        )
        source[interior_flat] += hedge_src[interior_flat]

        # Matrix entries at interior points (central differences for gradient)
        n_i, n_j = grid_shape[hc.i], grid_shape[hc.j]
        int_sl = [slice(None)] * d
        int_sl[hc.i] = slice(1, n_i - 1)
        int_sl[hc.j] = slice(1, n_j - 1)
        int_sl_t = tuple(int_sl)

        idx_c = flat_idx[int_sl_t].ravel()
        xi_int = hc.xi_star[int_sl_t].ravel()

        bc_shape = [
            grid_shape[l] if l != hc.i and l != hc.j
            else (n_i - 2 if l == hc.i else n_j - 2)
            for l in range(d)
        ]
        sl_i2 = [np.newaxis] * d
        sl_i2[hc.i] = slice(None)
        sl_j2 = [np.newaxis] * d
        sl_j2[hc.j] = slice(None)
        y_i_int = (spec.y_grids[hc.i][1:-1][tuple(sl_i2)] * np.ones(bc_shape)).ravel()
        y_j_int = (spec.y_grids[hc.j][1:-1][tuple(sl_j2)] * np.ones(bc_shape)).ravel()

        # ∂θ/∂y_i coupling (positive sign in p) — central differences
        coeff_i = (xi_int * (1.0 + hc.k_i * y_i_int)
                   / (2.0 * spec.dy_list[hc.i]))
        rows_list.append(idx_c)
        cols_list.append(idx_c + strides[hc.i])
        vals_list.append(-dt * coeff_i)
        rows_list.append(idx_c)
        cols_list.append(idx_c - strides[hc.i])
        vals_list.append(dt * coeff_i)

        # ∂θ/∂y_j coupling (negative sign in p) — central differences
        coeff_j = (-xi_int * (1.0 + hc.k_j * y_j_int)
                   / (2.0 * spec.dy_list[hc.j]))
        rows_list.append(idx_c)
        cols_list.append(idx_c + strides[hc.j])
        vals_list.append(-dt * coeff_j)
        rows_list.append(idx_c)
        cols_list.append(idx_c - strides[hc.j])
        vals_list.append(dt * coeff_j)

    # --- Assemble sparse matrix ---
    all_rows = np.concatenate(rows_list)
    all_cols = np.concatenate(cols_list)
    all_vals = np.concatenate(vals_list)

    A = scipy.sparse.coo_matrix(
        (all_vals, (all_rows, all_cols)), shape=(N, N),
    ).tocsr()

    return A, source


def solve_hjb_implicit(
    y_grids: List[np.ndarray],
    mp: ModelParams,
    n_steps: int = 1000,
    snapshot_times: List[float] | None = None,
    pi_tol: float = 1e-10,
    pi_max_iter: int = 20,
) -> dict:
    """Solve PDE (1) backward from T to 0 via fully implicit Euler
    with policy iteration (Howard's algorithm).

    At each time step, the nonlinear implicit equation is solved by
    iterating: (1) fix controls from current θ, (2) solve the resulting
    linear system for the new θ, until convergence.

    Parameters
    ----------
    y_grids : list of d 1D arrays (grid axes).
    mp : ModelParams with T_days, gamma, kappa, etc.
    n_steps : number of time steps.
    snapshot_times : optional list of times (in days) at which to save θ.
    pi_tol : policy iteration convergence tolerance (relative l-inf).
    pi_max_iter : maximum policy iterations per time step.

    Returns
    -------
    dict with keys:
        'theta_0'   : θ(0, y), d-dimensional array on the grid.
        'snapshots' : dict mapping time -> θ(t, y) for requested snapshot_times.
        'spec'      : the PDESpec used.
        'dt'        : time step size.
        'pi_iters'  : list of policy iteration counts per time step.
    """
    import scipy.sparse.linalg

    d = len(mp.currencies)
    kappa = mp.kappa if mp.kappa is not None else np.zeros((d, d))
    T = mp.T_days
    dt = T / n_steps

    spec = build_pde_spec(y_grids, mp)
    theta = terminal_condition(y_grids, kappa)
    grid_shape = theta.shape

    # Prepare snapshot collection
    if snapshot_times is None:
        snapshot_times = []
    snap_steps = {}
    for t_snap in snapshot_times:
        step_idx = round((T - t_snap) / dt)
        step_idx = max(0, min(n_steps, step_idx))
        snap_steps[step_idx] = t_snap
    snapshots = {}

    pi_iters = []

    from tqdm import trange

    for m in trange(n_steps, desc="HJB PDE solve (implicit)", unit="step"):
        theta_old = theta.copy()

        for k_pi in range(pi_max_iter):
            # Step A: extract optimal controls from current θ
            grad = compute_gradient(theta, spec.dy_list)
            q_ctrl = extract_quoting_controls(theta, spec)
            h_ctrl = extract_hedging_controls(grad, spec.y_grids, spec)

            # Step B: assemble and solve linearized system
            A, source = assemble_implicit_system(q_ctrl, h_ctrl, spec, dt)
            b = theta_old.ravel() + source
            theta_new = scipy.sparse.linalg.spsolve(A, b).reshape(grid_shape)

            # Check convergence
            diff = np.max(np.abs(theta_new - theta))
            scale = max(1.0, np.max(np.abs(theta_new)))
            theta = theta_new

            if diff / scale < pi_tol:
                break

        pi_iters.append(k_pi + 1)

        if (m + 1) in snap_steps:
            t_snap = snap_steps[m + 1]
            snapshots[t_snap] = theta.copy()

    return {
        'theta_0': theta,
        'snapshots': snapshots,
        'spec': spec,
        'dt': dt,
        'pi_iters': pi_iters,
    }
