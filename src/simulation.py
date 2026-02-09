from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

from .model import ModelParams, canon_pair, DAY_SECONDS
from .hamiltonian import logistic_f
from .policy import MMResult


def client_trade_intensity_per_sec(mp: ModelParams,
                                  tier_idx: int,
                                  ccy_pay: str, ccy_sell: str,
                                  z_musd: float,
                                  delta: float) -> float:
    """Return intensity (1/sec) for a given tier, direction, and size."""
    key = canon_pair(ccy_pay, ccy_sell)
    pp = mp.pairs[key]
    tier = pp.tiers[tier_idx]

    k = int(np.where(np.isclose(pp.sizes_musd, z_musd))[0][0])
    lam_day = pp.lambdas_per_day[k]

    f = logistic_f(delta, tier.alpha, tier.beta)
    return (lam_day * f) / DAY_SECONDS


def simulate_inventory_path_tau_leap(
    res: MMResult,
    T_sec: float,
    dt_sec: float,
    seed: int = 123,
    include_hedging: bool = True,
    tiers_to_use: Optional[List[int]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Simulate inventories Y_t (M$) under the approximate policy via tau-leaping.

    Returns (times, Y_path) where Y_path has shape (n_steps+1, d).
    """
    mp = res.mp
    rng = np.random.default_rng(seed)

    d = len(mp.currencies)
    n_steps = int(np.ceil(T_sec / dt_sec))
    times = np.linspace(0.0, n_steps * dt_sec, n_steps + 1)

    Y = np.zeros((n_steps + 1, d), dtype=float)

    if tiers_to_use is None:
        tiers_to_use = list(range(len(next(iter(mp.pairs.values())).tiers)))

    unordered_pairs = list(mp.pairs.keys())

    for t in range(n_steps):
        y = Y[t].copy()

        for (a, b) in unordered_pairs:
            pp = mp.pairs[(a, b)]
            for z in pp.sizes_musd:
                for tier_idx in tiers_to_use:
                    delta_ab = res.markup(tier_idx, ccy_pay=a, ccy_sell=b, z_musd=z, y=y)
                    lam_ab = client_trade_intensity_per_sec(mp, tier_idx, a, b, z, delta_ab)

                    delta_ba = res.markup(tier_idx, ccy_pay=b, ccy_sell=a, z_musd=z, y=y)
                    lam_ba = client_trade_intensity_per_sec(mp, tier_idx, b, a, z, delta_ba)

                    n_ab = rng.poisson(lam_ab * dt_sec)
                    n_ba = rng.poisson(lam_ba * dt_sec)

                    if n_ab > 0:
                        ia = mp.currencies.index(a)
                        ib = mp.currencies.index(b)
                        y[ia] += n_ab * z
                        y[ib] -= n_ab * z

                    if n_ba > 0:
                        ia = mp.currencies.index(b)
                        ib = mp.currencies.index(a)
                        y[ia] += n_ba * z
                        y[ib] -= n_ba * z

        if include_hedging:
            for (a, b) in unordered_pairs:
                xi_day = res.hedge_rate(a, b, y=y)
                xi_sec = xi_day / DAY_SECONDS
                ia = mp.currencies.index(a)
                ib = mp.currencies.index(b)
                y[ia] += xi_sec * dt_sec
                y[ib] -= xi_sec * dt_sec

        Y[t + 1] = y

    return times, Y


def autocorr(x: np.ndarray, max_lag: int) -> np.ndarray:
    """Simple autocorrelation function for 1D series."""
    x = x - np.mean(x)
    denom = np.dot(x, x)
    if denom <= 0:
        return np.zeros(max_lag + 1)
    ac = np.empty(max_lag + 1)
    ac[0] = 1.0
    for k in range(1, max_lag + 1):
        ac[k] = np.dot(x[:-k], x[k:]) / denom
    return ac
