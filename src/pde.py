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
