from __future__ import annotations

from typing import Tuple

import numpy as np


def logistic_f(delta: float, alpha: float, beta: float) -> float:
    return 1.0 / (1.0 + np.exp(alpha + beta * delta))


# ---------------------------------------------------------------------------
# Lambert-W via Newton iteration (replaces bisection)
# ---------------------------------------------------------------------------

def _lambert_w0_newton(x: float, tol: float = 1e-12, max_iter: int = 8) -> float:
    """Compute W_0(x) for x >= 0 via Newton iteration on w*e^w = x."""
    if x < 1.0:
        w = x
    else:
        w = np.log(1.0 + x)
    for _ in range(max_iter):
        ew = np.exp(w)
        dw = (w * ew - x) / (ew * (w + 1.0))
        w = w - dw
        if abs(dw) < tol:
            break
    return w


def optimal_delta_logistic(p: float, alpha: float, beta: float) -> float:
    """Compute delta_bar(p) = argmax_delta f(delta)*(delta-p) via Lambert W."""
    x = np.exp(-(1.0 + alpha + beta * p))
    w = _lambert_w0_newton(x)
    return p + (w + 1.0) / beta


def H_logistic(p: float, alpha: float, beta: float) -> Tuple[float, float, float]:
    """Return (H(p), delta_bar(p), f(delta_bar(p))) via Lambert W."""
    x = np.exp(-(1.0 + alpha + beta * p))
    w = _lambert_w0_newton(x)
    delta_star = p + (w + 1.0) / beta
    f_star = w / (w + 1.0)
    H = w / beta
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
