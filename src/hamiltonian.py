from __future__ import annotations

from typing import Tuple

import numpy as np


def logistic_f(delta: float, alpha: float, beta: float) -> float:
    return 1.0 / (1.0 + np.exp(alpha + beta * delta))


# ---------------------------------------------------------------------------
# Lambert-W via Newton iteration (replaces bisection)
# ---------------------------------------------------------------------------

def _lambert_w0_newton(x, tol=1e-12, max_iter=8):
    """Compute W_0(x) for x >= 0 via Newton iteration on w*e^w = x.

    Works on both scalars and numpy arrays.
    For very large x (> e^500), uses the asymptotic approximation W(x) ≈ ln(x).
    """
    x = np.asarray(x, dtype=float)
    BIG = 1e200
    big = x > BIG
    w = np.where(x < 1.0, x, np.log(1.0 + x))
    if np.any(big):
        w = np.where(big, np.where(np.isinf(x), x, np.log(x)), w)
    for _ in range(max_iter):
        if np.all(big):
            break
        ew = np.exp(w)
        denom = ew * (w + 1.0)
        safe = (~big) & (denom != 0.0)
        dw = np.where(safe, (w * ew - x) / denom, 0.0)
        w = w - dw
        if np.all(np.abs(dw) < tol):
            break
    return w


def optimal_delta_logistic(p, alpha, beta):
    """Compute delta_bar(p) = argmax_delta f(delta)*(delta-p) via Lambert W.

    Overflow-safe: for very negative p, arg = -(1+α+βp) is large and
    exp(arg) overflows.  In that regime W(x) ≈ arg, so
    δ* = p + (W+1)/β ≈ p + (arg+1)/β → -α/β.
    """
    arg = -(1.0 + alpha + beta * np.asarray(p, dtype=float))
    safe_arg = np.minimum(arg, 700.0)
    x = np.exp(safe_arg)
    w = _lambert_w0_newton(x)
    overflow = arg > 700.0
    w = np.where(overflow, arg, w)
    return p + (w + 1.0) / beta


def H_logistic(p, alpha, beta):
    """Return (H(p), delta_bar(p), f(delta_bar(p))) via Lambert W.

    Overflow-safe: for very negative p, arg = -(1+α+βp) is large and
    exp(arg) overflows.  In that regime W(x) ≈ arg, so:
      H = W/β ≈ arg/β
      δ* = p + (W+1)/β ≈ p + (arg+1)/β
      f* = W/(W+1) → 1  (fill probability saturates)
    """
    p = np.asarray(p, dtype=float)
    arg = -(1.0 + alpha + beta * p)
    safe_arg = np.minimum(arg, 700.0)
    x = np.exp(safe_arg)
    w = _lambert_w0_newton(x)
    overflow = arg > 700.0
    w = np.where(overflow, arg, w)
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
