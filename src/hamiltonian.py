from __future__ import annotations

from typing import Tuple

import numpy as np


def logistic_f(delta: float, alpha: float, beta: float) -> float:
    return 1.0 / (1.0 + np.exp(alpha + beta * delta))


def _g_prime(delta: float, p: float, alpha: float, beta: float) -> float:
    """Derivative of g(delta)=f(delta)*(delta-p) wrt delta."""
    f = logistic_f(delta, alpha, beta)
    fp = -beta * f * (1.0 - f)
    return fp * (delta - p) + f


def optimal_delta_logistic(p: float, alpha: float, beta: float,
                           tol: float = 1e-14,
                           max_iter: int = 200) -> float:
    """Compute delta_bar(p) = argmax_delta f(delta)*(delta-p) for logistic f."""
    lower = float(p)
    step = 1.0 / max(abs(beta), 1.0)
    upper = lower + step

    gp_upper = _g_prime(upper, p, alpha, beta)
    it = 0
    while gp_upper > 0.0 and it < 100:
        step *= 2.0
        upper = lower + step
        gp_upper = _g_prime(upper, p, alpha, beta)
        it += 1

    for _ in range(max_iter):
        mid = 0.5 * (lower + upper)
        gp_mid = _g_prime(mid, p, alpha, beta)
        if gp_mid > 0.0:
            lower = mid
        else:
            upper = mid
        if (upper - lower) < tol:
            break

    return 0.5 * (lower + upper)


def H_logistic(p: float, alpha: float, beta: float) -> Tuple[float, float, float]:
    """Return (H(p), delta_bar(p), f(delta_bar(p)))."""
    delta_star = optimal_delta_logistic(p, alpha, beta)
    f_star = logistic_f(delta_star, alpha, beta)
    H = f_star * (delta_star - p)
    return H, delta_star, f_star


def quadratic_coeffs_H_logistic(alpha: float, beta: float,
                                eps_p: float = 1e-8) -> Tuple[float, float, float]:
    """Compute (alpha0, alpha1, alpha2) for H(p) around p=0."""
    H0, delta0, f0 = H_logistic(0.0, alpha, beta)
    alpha0 = H0
    alpha1 = -f0

    _, _, f_plus = H_logistic(+eps_p, alpha, beta)
    _, _, f_minus = H_logistic(-eps_p, alpha, beta)
    Hp_plus = -f_plus
    Hp_minus = -f_minus
    alpha2 = (Hp_plus - Hp_minus) / (2.0 * eps_p)

    return alpha0, alpha1, alpha2
