"""PDE vs ODE comparison: worker functions for parallel computation."""
from __future__ import annotations
import os
import numpy as np
from .model import ModelParams, PairParams, TierParams, BP, build_paper_example_params, restrict_currencies
from .hamiltonian import optimal_delta_logistic, logistic_f
from .policy import Hprime_execution_cost, run_multicurrency_mm
from .pde import solve_hjb_implicit

# Reuse the base 2-currency params (built once at import time)
_BASE_2CCY = restrict_currencies(build_paper_example_params(), ["USD", "EUR"])
_BASE_PP = _BASE_2CCY.pairs[("EUR", "USD")]


# ---------------------------------------------------------------------------
# Parameter definitions
# ---------------------------------------------------------------------------

PARAM_LABELS = [
    "sigma_EUR", "gamma", "lambda_scale",
    "alpha_1", "alpha_2", "beta_1", "beta_2", "eta"
]

NOMINAL = np.array([80.0, 20.0, 1.0, -1.9, -0.3, 11.0, 3.5, 1e-5])

RANGES = np.array([
    [56.0, 104.0],       # sigma_EUR (bps), +-30%
    [10.0, 40.0],        # gamma, +-50%
    [0.5, 1.5],          # lambda_scale, +-50%
    [-2.5, -1.3],        # alpha_1, +-30%
    [-0.9, 0.3],         # alpha_2, +-100%
    [7.7, 14.3],         # beta_1 (1/bps), +-30%
    [2.45, 4.55],        # beta_2 (1/bps), +-30%
    [0.5e-5, 2.0e-5],    # eta (bps), +-50%
])

N_PARAMS = len(PARAM_LABELS)

QOI_NAMES = [
    "tier_spread_diff",    # delta*(tier1) - delta*(tier2) at y=0 [bps]
    "own_inv_skew",        # delta*(y=10) - delta*(y=0), tier 1 [bps]
    "hedge_rate",          # xi*(y=10) [M$/day]
    "net_revenue",         # R(y=10) - L(xi*) [$/day]
    "neutral_spread",      # delta*(y=0), tier 1 [bps]  (diagnostic)
]

QOI_LABELS = [
    r"$\delta^*_1(0) - \delta^*_2(0)$ [bps]",
    r"$\delta^*(10) - \delta^*(0)$ [bps]",
    r"$\xi^*(10)$ [M\$/day]",
    r"$R(10) - L(10)$ [\$/day]",
    r"$\delta^*(0)$ [bps]",
]

N_QOIS = len(QOI_NAMES)


# ---------------------------------------------------------------------------
# Parameter builder
# ---------------------------------------------------------------------------

def build_modified_params_2ccy(params):
    """Build 2-currency (USD, EUR) ModelParams from 8-element parameter vector.

    params: [sigma_EUR(bps), gamma, lambda_scale, alpha_1, alpha_2,
             beta_1(1/bps), beta_2(1/bps), eta(bps)]
    """
    sigma_eur, gamma, lam_scale, a1, a2, b1, b2, eta = params

    tiers = [
        TierParams(alpha=a1, beta=b1 * 1e4),
        TierParams(alpha=a2, beta=b2 * 1e4),
    ]

    pp = PairParams(
        pair=("EUR", "USD"),
        sizes_musd=_BASE_PP.sizes_musd,
        lambdas_per_day=_BASE_PP.lambdas_per_day * lam_scale,
        tiers=tiers,
        psi=_BASE_PP.psi,
        eta=eta * BP,
    )

    return ModelParams(
        currencies=["USD", "EUR"],
        ref_ccy="USD",
        sigma={"USD": 0.0, "EUR": sigma_eur * BP},
        corr={},
        k=_BASE_2CCY.k,
        mu={"USD": 0.0, "EUR": 0.0},
        gamma=gamma,
        kappa=np.zeros((2, 2)),
        T_days=_BASE_2CCY.T_days,
        pairs={("EUR", "USD"): pp},
    )


# ---------------------------------------------------------------------------
# ODE QoI evaluation
# ---------------------------------------------------------------------------

def evaluate_ode_qois(params):
    """Solve ODE and return 5 QoIs.

    Returns array of:
      [0] tier spread diff (bps): delta*(tier1, y=0) - delta*(tier2, y=0), z=1 M$
      [1] own-inventory skew (bps): delta*(y=10) - delta*(y=0), tier 1, z=1
      [2] hedge rate (M$/day): xi*(EUR/USD) at y=[0, 10]
      [3] net revenue ($/day): R(y=[0,10]) - L(xi*) at y=[0,10]
      [4] neutral spread (bps): delta*(tier1, y=0), z=1 M$
    """
    mp = build_modified_params_2ccy(params)
    res = run_multicurrency_mm(mp, n_steps=500)

    y_flat = np.zeros(2)
    y_long = np.array([0.0, 10.0])  # 10 M$ long EUR

    # QoI 4 (neutral spread) and QoI 0 (tier diff)
    delta_t1_flat = res.markup(0, "EUR", "USD", 1.0, y_flat)
    delta_t2_flat = res.markup(1, "EUR", "USD", 1.0, y_flat)
    tier_diff = (delta_t1_flat - delta_t2_flat) / BP
    neutral_spread = delta_t1_flat / BP

    # QoI 1 (own-inventory skew)
    delta_t1_long = res.markup(0, "EUR", "USD", 1.0, y_long)
    own_skew = (delta_t1_long - delta_t1_flat) / BP

    # QoI 2 (hedge rate)
    xi = res.hedge_rate("EUR", "USD", y_long)

    # QoI 3 (net revenue)
    net_rev = _revenue_rate_2ccy(mp, res, y_long) - _hedging_cost_2ccy(mp, res, y_long)

    return np.array([tier_diff, own_skew, xi, net_rev, neutral_spread])


# ---------------------------------------------------------------------------
# Revenue / hedging cost helpers
# ---------------------------------------------------------------------------

def _revenue_rate_2ccy(mp, res, y):
    """Instantaneous expected revenue rate ($/day) at inventory y.

    R = sum over directions, tiers, sizes of: lambda * f(delta*) * delta* * z
    Converted from M$/day to $/day (* 1e6).
    """
    revenue = 0.0
    pp = mp.pairs[("EUR", "USD")]
    for ccy_pay, ccy_sell in [("EUR", "USD"), ("USD", "EUR")]:
        for t_idx, tier in enumerate(pp.tiers):
            for z, lam in zip(pp.sizes_musd, pp.lambdas_per_day):
                delta = res.markup(t_idx, ccy_pay, ccy_sell, z, y)
                f = logistic_f(delta, tier.alpha, tier.beta)
                revenue += lam * f * delta * z
    return revenue * 1e6


def _hedging_cost_2ccy(mp, res, y):
    """Instantaneous hedging cost rate ($/day) at inventory y."""
    pp = mp.pairs[("EUR", "USD")]
    xi = res.hedge_rate("EUR", "USD", y)
    return (pp.psi * abs(xi) + pp.eta * xi**2) * 1e6


# ---------------------------------------------------------------------------
# PDE QoI extraction
# ---------------------------------------------------------------------------

def pde_extract_qois(theta_0, y_grids, mp):
    """Extract 5 QoIs from PDE solution theta(0, y).

    Same QoIs as evaluate_ode_qois but computed from the PDE grid.
    """
    idx0 = len(y_grids[0]) // 2  # index of y=0
    dy = y_grids[0][1] - y_grids[0][0]
    y_eur_idx = int(round(10.0 / dy))  # index offset for 10 M$ EUR

    pp = mp.pairs[("EUR", "USD")]
    tier1 = pp.tiers[0]
    tier2 = pp.tiers[1]
    z = 1.0
    shift = int(round(z / dy))

    # Quoting momentum p at y=0 for tier evaluation
    # EUR pays (i=1), USD sells (j=0) -> dvec = [-1, +1]
    # Shifted point: [y_USD - z, y_EUR + z]
    p_flat = (theta_0[idx0, idx0] - theta_0[idx0 - shift, idx0 + shift]) / z

    # QoI 4: neutral spread (tier 1)
    delta_t1_flat = optimal_delta_logistic(p_flat, tier1.alpha, tier1.beta)
    neutral_spread = delta_t1_flat / BP

    # QoI 0: tier spread differential
    delta_t2_flat = optimal_delta_logistic(p_flat, tier2.alpha, tier2.beta)
    tier_diff = (delta_t1_flat - delta_t2_flat) / BP

    # QoI 1: own-inventory skew
    p_long = (theta_0[idx0, idx0 + y_eur_idx] -
              theta_0[idx0 - shift, idx0 + y_eur_idx + shift]) / z
    delta_t1_long = optimal_delta_logistic(p_long, tier1.alpha, tier1.beta)
    own_skew = (delta_t1_long - delta_t1_flat) / BP

    # QoI 2: hedge rate at y=[0, 10]
    # Gradient via central differences
    grad_usd = (theta_0[idx0 + 1, idx0 + y_eur_idx] -
                theta_0[idx0 - 1, idx0 + y_eur_idx]) / (2 * dy)
    grad_eur = (theta_0[idx0, idx0 + y_eur_idx + 1] -
                theta_0[idx0, idx0 + y_eur_idx - 1]) / (2 * dy)
    # Hedge EUR buy, USD sell: base = grad_EUR - grad_USD
    base = grad_eur - grad_usd
    # Market impact (small but included)
    k_eur = mp.k.get("EUR", 0.0)
    k_usd = mp.k.get("USD", 0.0)
    impact = k_eur * 10.0 * (1.0 + grad_eur) - k_usd * 0.0 * (1.0 + grad_usd)
    xi = Hprime_execution_cost(base + impact, pp.psi, pp.eta)

    # QoI 3: net revenue at y=[0, 10]
    # Revenue: sum over directions, tiers, sizes of lambda * f(delta*) * delta* * z
    revenue = 0.0
    for ccy_pay, ccy_sell, sign_usd, sign_eur in [
        ("EUR", "USD", -1, +1),  # EUR pays: dvec = [-1, +1]
        ("USD", "EUR", +1, -1),  # USD pays: dvec = [+1, -1]
    ]:
        for t_idx, tier in enumerate(pp.tiers):
            for zk, lam in zip(pp.sizes_musd, pp.lambdas_per_day):
                sk = int(round(zk / dy))
                # Check bounds
                usd_idx = idx0 + sign_usd * sk
                eur_idx = idx0 + y_eur_idx + sign_eur * sk
                if (0 <= usd_idx < theta_0.shape[0] and
                    0 <= eur_idx < theta_0.shape[1]):
                    p_k = (theta_0[idx0, idx0 + y_eur_idx] -
                           theta_0[usd_idx, eur_idx]) / zk
                    delta_k = optimal_delta_logistic(p_k, tier.alpha, tier.beta)
                    f_k = logistic_f(delta_k, tier.alpha, tier.beta)
                    revenue += lam * f_k * delta_k * zk
    revenue_dollars = revenue * 1e6
    hedging_cost = (pp.psi * abs(xi) + pp.eta * xi**2) * 1e6
    net_rev = revenue_dollars - hedging_cost

    return np.array([tier_diff, own_skew, xi, net_rev, neutral_spread])


# ---------------------------------------------------------------------------
# PDE grid constants and worker
# ---------------------------------------------------------------------------

Y_MAX = 150
DY = 1.0
PDE_N_STEPS = 20
PDE_PI_TOL = 1e-4
PDE_PI_MAX_ITER = 50


def _init_worker():
    """Initializer for pool workers: limit BLAS threads and suppress tqdm."""
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["TQDM_DISABLE"] = "1"


def run_single_pde(config):
    """Worker: solve one PDE configuration and return results.

    config: (name, params_array)
    Returns: dict with name, params, ode_qois, pde_qois, theta_0
    """
    name, params = config

    # ODE (fast)
    ode_qois = evaluate_ode_qois(params)

    # PDE (slow)
    mp = build_modified_params_2ccy(params)
    y_grids = [np.arange(-Y_MAX, Y_MAX + DY, DY) for _ in range(2)]

    result = solve_hjb_implicit(
        y_grids, mp,
        n_steps=PDE_N_STEPS,
        pi_tol=PDE_PI_TOL,
        pi_max_iter=PDE_PI_MAX_ITER,
    )
    theta_0 = result['theta_0']
    pde_qois = pde_extract_qois(theta_0, y_grids, mp)

    return {
        'name': name,
        'params': params,
        'ode_qois': ode_qois,
        'pde_qois': pde_qois,
        'theta_0': theta_0,
    }


# ---------------------------------------------------------------------------
# Config generators
# ---------------------------------------------------------------------------

def generate_oat_configs():
    """Generate OAT parameter configurations: nominal + 2 edges per parameter."""
    configs = [("nominal", NOMINAL.copy())]

    for i, label in enumerate(PARAM_LABELS):
        p_lo = NOMINAL.copy()
        p_lo[i] = RANGES[i, 0]
        configs.append((f"{label}_low", p_lo))

        p_hi = NOMINAL.copy()
        p_hi[i] = RANGES[i, 1]
        configs.append((f"{label}_high", p_hi))

    return configs


def generate_lhs_configs(n_samples=20, seed=42):
    """Generate Latin Hypercube sample configurations."""
    from scipy.stats.qmc import LatinHypercube

    sampler = LatinHypercube(d=N_PARAMS, seed=seed)
    unit_samples = sampler.random(n=n_samples)

    # Scale from [0,1] to parameter ranges
    lo = RANGES[:, 0]
    hi = RANGES[:, 1]
    samples = lo + unit_samples * (hi - lo)

    configs = []
    for i, s in enumerate(samples):
        configs.append((f"lhs_{i:02d}", s))

    return configs
