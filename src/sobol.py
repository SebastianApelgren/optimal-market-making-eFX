"""Saltelli sampling and Sobol index computation.

Implements the Saltelli (2010) Monte Carlo estimator for first-order and
total-effect Sobol sensitivity indices. The three functions are generic
and do not depend on any model-specific code.
"""

from __future__ import annotations

import numpy as np
from tqdm import tqdm


def saltelli_sample(N, ranges, seed=42):
    """Generate two independent uniform sample matrices A, B scaled to parameter ranges."""
    rng = np.random.default_rng(seed)
    k = ranges.shape[0]
    lo, hi = ranges[:, 0], ranges[:, 1]
    A = lo + rng.random((N, k)) * (hi - lo)
    B = lo + rng.random((N, k)) * (hi - lo)
    return A, B


def evaluate_saltelli_samples(model_func, A, B):
    """Evaluate model at all Saltelli sample points.

    Returns f_A (N, n_qoi), f_B (N, n_qoi), f_C (k, N, n_qoi).
    """
    N, k = A.shape
    total = N * (k + 2)

    # Stack all parameter sets: A, B, C_0, ..., C_{k-1}
    all_params = np.empty((total, k))
    all_params[:N] = A
    all_params[N:2*N] = B
    for i in range(k):
        C_i = A.copy()
        C_i[:, i] = B[:, i]
        all_params[(2 + i) * N:(3 + i) * N] = C_i

    # Evaluate all samples
    n_qoi = len(model_func(A[0]))
    results = np.empty((total, n_qoi))

    for j in tqdm(range(total), desc="Saltelli evaluations", mininterval=2.0):
        results[j] = model_func(all_params[j])

    # Unpack
    f_A = results[:N]
    f_B = results[N:2*N]
    f_C = np.empty((k, N, n_qoi))
    for i in range(k):
        f_C[i] = results[(2 + i) * N:(3 + i) * N]

    return f_A, f_B, f_C


def compute_sobol_indices(f_A, f_B, f_C):
    """Compute first-order and total-effect Sobol indices.

    First-order (Saltelli 2010):  S_i  = mean(f_B * (f_Ci - f_A)) / Var(Y)
    Total-effect (Jansen 1999):  S_Ti = mean((f_A - f_Ci)^2) / (2 * Var(Y))
    """
    N, n_qoi = f_A.shape
    k = f_C.shape[0]

    f_all = np.concatenate([f_A, f_B], axis=0)
    var_Y = np.var(f_all, axis=0)

    S = np.empty((k, n_qoi))
    S_T = np.empty((k, n_qoi))

    for i in range(k):
        S[i] = np.mean(f_B * (f_C[i] - f_A), axis=0) / var_Y
        S_T[i] = np.mean((f_A - f_C[i])**2, axis=0) / (2.0 * var_Y)

    return S, S_T
