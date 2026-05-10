"""3-currency sensitivity analysis: parameter builder and QoI evaluation.

Moved from the sensitivity_analysis.ipynb notebook so that the evaluate_qois
function is importable by multiprocessing workers (macOS spawn requires this).
"""
from __future__ import annotations

import numpy as np

from .model import ModelParams, PairParams, TierParams, BP, build_paper_example_params, restrict_currencies
from .hamiltonian import logistic_f
from .policy import run_multicurrency_mm


# Build base 3-currency params once at import time
_BASE_3CCY = restrict_currencies(build_paper_example_params(), ["USD", "EUR", "GBP"])
_BASE_EU = _BASE_3CCY.pairs[("EUR", "USD")]
_BASE_GU = _BASE_3CCY.pairs[("GBP", "USD")]
_BASE_EG = _BASE_3CCY.pairs[("EUR", "GBP")]


def build_modified_params(params):
    """Build 3-currency (USD, EUR, GBP) ModelParams from 20-element parameter vector.

    Parameters in human-readable units (sigma in bps, beta in 1/bps, eta in bps).
    Converted to internal units (decimal) inside this function.
    """
    (sigma_eur, sigma_gbp, rho, gamma, lam_scale,
     a1_eu, a2_eu, b1_eu, b2_eu, eta_eu,
     a1_gu, a2_gu, b1_gu, b2_gu, eta_gu,
     a1_eg, a2_eg, b1_eg, b2_eg, eta_eg) = params

    # EUR/USD pair
    tiers_eu = [
        TierParams(alpha=a1_eu, beta=b1_eu * 1e4),
        TierParams(alpha=a2_eu, beta=b2_eu * 1e4),
    ]
    pp_eu = PairParams(
        pair=("EUR", "USD"),
        sizes_musd=_BASE_EU.sizes_musd,
        lambdas_per_day=_BASE_EU.lambdas_per_day * lam_scale,
        tiers=tiers_eu,
        psi=_BASE_EU.psi,
        eta=eta_eu * BP,
    )

    # GBP/USD pair
    tiers_gu = [
        TierParams(alpha=a1_gu, beta=b1_gu * 1e4),
        TierParams(alpha=a2_gu, beta=b2_gu * 1e4),
    ]
    pp_gu = PairParams(
        pair=("GBP", "USD"),
        sizes_musd=_BASE_GU.sizes_musd,
        lambdas_per_day=_BASE_GU.lambdas_per_day * lam_scale,
        tiers=tiers_gu,
        psi=_BASE_GU.psi,
        eta=eta_gu * BP,
    )

    # EUR/GBP pair
    tiers_eg = [
        TierParams(alpha=a1_eg, beta=b1_eg * 1e4),
        TierParams(alpha=a2_eg, beta=b2_eg * 1e4),
    ]
    pp_eg = PairParams(
        pair=("EUR", "GBP"),
        sizes_musd=_BASE_EG.sizes_musd,
        lambdas_per_day=_BASE_EG.lambdas_per_day * lam_scale,
        tiers=tiers_eg,
        psi=_BASE_EG.psi,
        eta=eta_eg * BP,
    )

    return ModelParams(
        currencies=["USD", "EUR", "GBP"],
        ref_ccy="USD",
        sigma={"USD": 0.0, "EUR": sigma_eur * BP, "GBP": sigma_gbp * BP},
        corr={("EUR", "GBP"): rho},
        k=_BASE_3CCY.k,
        mu={"USD": 0.0, "EUR": 0.0, "GBP": 0.0},
        gamma=gamma,
        kappa=np.zeros((3, 3)),
        T_days=_BASE_3CCY.T_days,
        pairs={("EUR", "USD"): pp_eu, ("GBP", "USD"): pp_gu, ("EUR", "GBP"): pp_eg},
    )


def _revenue_rate(mp, res, y):
    """Instantaneous expected revenue rate ($/day) at inventory y.

    R = sum over all pairs, directions, tiers, sizes of: lambda_j * f(delta*) * delta* * z_j
    Converted from M$/day to $/day (* 1e6).
    """
    revenue = 0.0
    for (a, b), pp in mp.pairs.items():
        for ccy_pay, ccy_sell in [(a, b), (b, a)]:
            for t_idx, tier in enumerate(pp.tiers):
                for z, lam in zip(pp.sizes_musd, pp.lambdas_per_day):
                    delta = res.markup(t_idx, ccy_pay, ccy_sell, z, y)
                    f = logistic_f(delta, tier.alpha, tier.beta)
                    revenue += lam * f * delta * z
    return revenue * 1e6


def evaluate_qois(params, y_magnitude=10.0):
    """Solve ODE and return 6 QoIs for the 3-currency model.

    Set A (quoting):  [tier spread diff (bps), own-inv skew (bps), cross-inv skew (bps)]
    Set B (hedging):  [own-pair hedge (M$/day), cross-hedge momentum (decimal), net revenue ($/day)]

    `y_magnitude` sets the inventory (in M$) used for the skew, hedge, and revenue QoIs.
    """
    mp = build_modified_params(params)
    res = run_multicurrency_mm(mp, n_steps=500)

    y_flat = np.zeros(3)
    y_long_eur = np.array([0.0, y_magnitude, 0.0])

    # --- Set A: quoting policy ---
    delta_t1 = res.markup(0, "EUR", "USD", 1.0, y_flat)
    delta_t2 = res.markup(1, "EUR", "USD", 1.0, y_flat)
    tier_diff = (delta_t1 - delta_t2) / BP

    delta_t1_long_eur = res.markup(0, "EUR", "USD", 1.0, y_long_eur)
    own_skew = (delta_t1_long_eur - delta_t1) / BP

    y_long_gbp = np.array([0.0, 0.0, y_magnitude])
    delta_t1_long_gbp = res.markup(0, "EUR", "USD", 1.0, y_long_gbp)
    cross_skew = (delta_t1_long_gbp - delta_t1) / BP

    # --- Set B: hedging and economics ---
    xi_own = res.hedge_rate("EUR", "USD", y_long_eur)

    ccy = mp.currencies
    i_eur, i_gbp = ccy.index("EUR"), ccy.index("GBP")
    d_eg = np.zeros(3)
    d_eg[i_eur] = 1.0
    d_eg[i_gbp] = -1.0
    AyB = 2.0 * res.A0 @ y_long_eur + res.B0
    p_cross = -(AyB @ d_eg)

    R_long = _revenue_rate(mp, res, y_long_eur)
    L_total = 0.0
    for (a, b), pp in mp.pairs.items():
        xi = res.hedge_rate(a, b, y_long_eur)
        L_total += (pp.psi * abs(xi) + pp.eta * xi**2) * 1e6
    net_revenue = R_long - L_total

    return np.array([tier_diff, own_skew, cross_skew, xi_own, p_cross, net_revenue])
