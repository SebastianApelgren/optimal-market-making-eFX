from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import ModelParams, canon_pair
from .hamiltonian import optimal_delta_logistic
from .riccati import build_Sigma, build_M_tildeM_P, solve_AB_euler


def quote_p_scalar(y: np.ndarray, i: int, j: int, z: float,
                   A: np.ndarray, B: np.ndarray) -> float:
    """Compute the 'p' input used by delta_bar in the quote formula (p.6)."""
    d = y.shape[0]
    dvec = np.zeros(d)
    dvec[i] = 1.0
    dvec[j] = -1.0
    return ((2.0 * y + z * dvec) @ A + B) @ dvec


def optimal_client_markup(mp: ModelParams, tier_idx: int,
                          ccy_pay: str, ccy_sell: str,
                          z_musd: float,
                          y: np.ndarray,
                          A: np.ndarray, B: np.ndarray) -> float:
    """Compute optimal client markup delta in decimal."""
    ccy = mp.currencies
    i = ccy.index(ccy_pay)
    j = ccy.index(ccy_sell)

    key = canon_pair(ccy_pay, ccy_sell)
    if key not in mp.pairs:
        return 0.0

    pp = mp.pairs[key]
    tier = pp.tiers[tier_idx]

    p = quote_p_scalar(y, i, j, z_musd, A, B)
    return optimal_delta_logistic(p, tier.alpha, tier.beta)


def Hprime_execution_cost(p: float, psi: float, eta: float) -> float:
    """H'(p) for L(xi) = psi|xi| + eta xi^2."""
    if p > psi:
        return (p - psi) / (2.0 * eta)
    if p < -psi:
        return (p + psi) / (2.0 * eta)
    return 0.0


def optimal_hedge_rate(mp: ModelParams,
                       ccy_buy: str, ccy_sell: str,
                       y: np.ndarray,
                       A: np.ndarray, B: np.ndarray) -> float:
    """Compute optimal hedging rate xi (M$/day). Positive = buy ccy_buy, sell ccy_sell."""
    ccy = mp.currencies
    i = ccy.index(ccy_buy)
    j = ccy.index(ccy_sell)

    key = canon_pair(ccy_buy, ccy_sell)
    if key not in mp.pairs:
        return 0.0

    pp = mp.pairs[key]

    d = len(ccy)
    dvec = np.zeros(d); dvec[i] = 1.0; dvec[j] = -1.0

    AyB = 2.0 * A @ y + B
    base = -(AyB @ dvec)

    factor_i = 1.0 - AyB[i]
    factor_j = 1.0 - AyB[j]

    k_i = mp.k.get(ccy_buy, 0.0)
    k_j = mp.k.get(ccy_sell, 0.0)

    impact = k_i * y[i] * factor_i - k_j * y[j] * factor_j

    p_arg = base + impact
    return Hprime_execution_cost(p_arg, pp.psi, pp.eta)


@dataclass
class MMResult:
    mp: ModelParams
    Sigma: np.ndarray
    M: np.ndarray
    Mtilde: np.ndarray
    P: np.ndarray
    A0: np.ndarray
    B0: np.ndarray

    def markup(self, tier_idx: int, ccy_pay: str, ccy_sell: str, z_musd: float, y: np.ndarray) -> float:
        return optimal_client_markup(self.mp, tier_idx, ccy_pay, ccy_sell, z_musd, y, self.A0, self.B0)

    def hedge_rate(self, ccy_buy: str, ccy_sell: str, y: np.ndarray) -> float:
        return optimal_hedge_rate(self.mp, ccy_buy, ccy_sell, y, self.A0, self.B0)


def run_multicurrency_mm(mp: ModelParams,
                         eps_p: float = 1e-8,
                         n_steps: int = 2000) -> MMResult:
    Sigma = build_Sigma(mp)
    M, Mtilde, P = build_M_tildeM_P(mp, eps_p=eps_p)
    A0, B0 = solve_AB_euler(mp, M, Mtilde, P, Sigma, n_steps=n_steps)
    return MMResult(mp=mp, Sigma=Sigma, M=M, Mtilde=Mtilde, P=P, A0=A0, B0=B0)
