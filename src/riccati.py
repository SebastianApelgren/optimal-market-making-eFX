from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np

from .model import ModelParams, canon_pair
from .hamiltonian import quadratic_coeffs_H_logistic


def build_Sigma(mp: ModelParams) -> np.ndarray:
    """Build covariance matrix from volatilities and correlations (decimal units)."""
    ccy = mp.currencies
    d = len(ccy)

    sig = np.array([mp.sigma.get(c, 0.0) for c in ccy], dtype=float)
    Sigma = np.zeros((d, d), dtype=float)

    for i in range(d):
        for j in range(d):
            if i == j:
                Sigma[i, i] = sig[i] * sig[i]
            else:
                ci, cj = ccy[i], ccy[j]
                if ci == mp.ref_ccy or cj == mp.ref_ccy:
                    rho = 0.0
                else:
                    rho = mp.corr.get(canon_pair(ci, cj), 0.0)
                Sigma[i, j] = rho * sig[i] * sig[j]
    return Sigma


def build_M_tildeM_P(mp: ModelParams,
                     eps_p: float = 1e-8) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build (M, Mtilde, P) as dxd matrices."""
    ccy = mp.currencies
    d = len(ccy)

    M = np.zeros((d, d), dtype=float)
    Mtilde = np.zeros((d, d), dtype=float)
    P = np.zeros((d, d), dtype=float)

    coeff_cache: Dict[Tuple[str, str], List[Tuple[float, float, float]]] = {}

    for i in range(d):
        for j in range(d):
            if i == j:
                continue
            key = canon_pair(ccy[i], ccy[j])
            if key not in mp.pairs:
                continue

            pp = mp.pairs[key]

            if key not in coeff_cache:
                coeffs_per_tier = []
                for tier in pp.tiers:
                    coeffs_per_tier.append(quadratic_coeffs_H_logistic(tier.alpha, tier.beta, eps_p=eps_p))
                coeff_cache[key] = coeffs_per_tier

            z = pp.sizes_musd
            lam = pp.lambdas_per_day

            for (a0, a1, a2) in coeff_cache[key]:
                M[i, j] += np.sum(a2 * z * lam)
                Mtilde[i, j] += np.sum(a1 * z * lam)
                P[i, j] += np.sum(a2 * (z ** 2) * lam)

    return M, Mtilde, P


def _V_of_A(A: np.ndarray, P: np.ndarray) -> np.ndarray:
    """V(A) = D(A)P + P D(A) - 2(P . A)."""
    DA = np.diag(np.diag(A))
    return DA @ P + P @ DA - 2.0 * (P * A)


def solve_AB_euler(mp: ModelParams,
                   M: np.ndarray,
                   Mtilde: np.ndarray,
                   P: np.ndarray,
                   Sigma: np.ndarray,
                   n_steps: int = 2000) -> Tuple[np.ndarray, np.ndarray]:
    """Solve Eq. (5) backward with explicit Euler. Returns (A0, B0)."""
    ccy = mp.currencies
    d = len(ccy)
    U = np.ones(d)

    symM = M + M.T
    M_big = np.diag(symM @ U) - symM
    V = (Mtilde - Mtilde.T) @ U

    T = mp.T_days
    dt = T / n_steps

    kappa = mp.kappa if mp.kappa is not None else np.zeros((d, d))
    A = kappa.copy()
    B = np.zeros(d)

    mu_vec = np.array([mp.mu.get(c, 0.0) for c in ccy], dtype=float)
    Dmu = np.diag(mu_vec)

    for _ in range(n_steps):
        Vtilde = (_V_of_A(A, P) - _V_of_A(A, P).T) @ U

        A_dot = 2.0 * A @ M_big @ A - (Sigma * A) - 2.0 * Dmu @ A - 0.5 * mp.gamma * Sigma
        B_dot = mu_vec - Dmu @ B + 2.0 * A @ V + 2.0 * A @ Vtilde + 2.0 * A @ M_big @ B

        A = A - dt * A_dot
        B = B - dt * B_dot

        A = 0.5 * (A + A.T)

    return A, B
