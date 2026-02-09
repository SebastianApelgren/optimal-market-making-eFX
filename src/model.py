from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional

import numpy as np

BP = 1e-4
DAY_SECONDS = 86400.0


@dataclass(frozen=True)
class TierParams:
    alpha: float
    beta: float


@dataclass(frozen=True)
class PairParams:
    pair: Tuple[str, str]
    sizes_musd: np.ndarray
    lambdas_per_day: np.ndarray
    tiers: List[TierParams]
    psi: float
    eta: float


@dataclass
class ModelParams:
    currencies: List[str]
    ref_ccy: str = "USD"
    sigma: Dict[str, float] = field(default_factory=dict)
    corr: Dict[Tuple[str, str], float] = field(default_factory=dict)
    k: Dict[str, float] = field(default_factory=dict)
    mu: Dict[str, float] = field(default_factory=dict)
    gamma: float = 20.0
    kappa: Optional[np.ndarray] = None
    T_days: float = 0.05
    pairs: Dict[Tuple[str, str], PairParams] = field(default_factory=dict)


def canon_pair(a: str, b: str) -> Tuple[str, str]:
    return tuple(sorted((a, b)))


def build_paper_example_params() -> ModelParams:
    """Reproduce the paper's 5-currency numerical example (Table 1, p.7)."""
    currencies = ["USD", "EUR", "JPY", "GBP", "CHF"]

    sigma = {
        "USD": 0.0,
        "EUR": 80.0 * BP,
        "GBP": 70.0 * BP,
        "CHF": 60.0 * BP,
        "JPY": 60.0 * BP,
    }

    corr_raw = {
        ("EUR", "GBP"): 0.6,
        ("EUR", "CHF"): 0.5,
        ("EUR", "JPY"): 0.3,
        ("GBP", "CHF"): 0.3,
        ("GBP", "JPY"): 0.2,
        ("CHF", "JPY"): 0.4,
    }
    corr = {tuple(sorted(k)): v for k, v in corr_raw.items()}

    k = {
        "USD": 0.0,
        "EUR": 5e-3 * BP,
        "GBP": 7e-3 * BP,
        "CHF": 8e-3 * BP,
        "JPY": 6e-3 * BP,
    }

    mu = {c: 0.0 for c in currencies}

    sizes = np.array([1.0, 5.0, 10.0, 20.0, 50.0])

    table = {
        ("EUR", "USD"): dict(lambdas=[900, 540, 234, 90, 36], alpha=[-1.9, -0.3], beta=[11, 3.5], psi_bps=0.1, eta_bps=1e-5),
        ("GBP", "USD"): dict(lambdas=[600, 200, 150, 40, 10], alpha=[-1.4,  0.0], beta=[5.5, 2.0], psi_bps=0.15, eta_bps=1.5e-5),
        ("CHF", "USD"): dict(lambdas=[420, 140, 105, 28,  7], alpha=[-1.2,  0.0], beta=[4.5, 1.9], psi_bps=0.25, eta_bps=2.5e-5),
        ("JPY", "USD"): dict(lambdas=[825, 375, 180,105, 15], alpha=[-1.6, -0.1], beta=[9.0, 3.0], psi_bps=0.1, eta_bps=1.5e-5),
        ("EUR", "GBP"): dict(lambdas=[400, 50, 25, 20, 5], alpha=[-0.5, 0.5], beta=[3.5, 2.5], psi_bps=0.25, eta_bps=3e-5),
        ("EUR", "CHF"): dict(lambdas=[400, 50, 25, 20, 5], alpha=[-0.5, 0.5], beta=[3.5, 2.5], psi_bps=0.25, eta_bps=3e-5),
        ("EUR", "JPY"): dict(lambdas=[400, 50, 25, 20, 5], alpha=[-0.5, 0.5], beta=[3.5, 2.5], psi_bps=0.25, eta_bps=3e-5),
        ("GBP", "CHF"): dict(lambdas=[160, 20, 10,  8, 2], alpha=[-0.5, 0.5], beta=[3.5, 2.5], psi_bps=0.4, eta_bps=5e-5),
        ("GBP", "JPY"): dict(lambdas=[160, 20, 10,  8, 2], alpha=[-0.5, 0.5], beta=[3.5, 2.5], psi_bps=0.4, eta_bps=5e-5),
        ("CHF", "JPY"): dict(lambdas=[ 80, 10,  5,  4, 1], alpha=[-0.5, 0.5], beta=[3.5, 2.5], psi_bps=0.4, eta_bps=5e-5),
    }

    pairs: Dict[Tuple[str, str], PairParams] = {}
    for (a, b), p in table.items():
        tiers = [
            TierParams(alpha=p["alpha"][0], beta=p["beta"][0] * 1e4),
            TierParams(alpha=p["alpha"][1], beta=p["beta"][1] * 1e4),
        ]
        key = canon_pair(a, b)
        pairs[key] = PairParams(
            pair=key,
            sizes_musd=sizes,
            lambdas_per_day=np.array(p["lambdas"], dtype=float),
            tiers=tiers,
            psi=p["psi_bps"] * BP,
            eta=p["eta_bps"] * BP,
        )

    d = len(currencies)
    kappa = np.zeros((d, d))

    return ModelParams(
        currencies=currencies,
        ref_ccy="USD",
        sigma=sigma,
        corr=corr,
        k=k,
        mu=mu,
        gamma=20.0,
        kappa=kappa,
        T_days=0.05,
        pairs=pairs,
    )


def restrict_currencies(mp: ModelParams, keep: List[str]) -> ModelParams:
    """Return a new ModelParams restricted to a subset of currencies."""
    keep_set = set(keep)
    if mp.ref_ccy not in keep_set:
        raise ValueError(f"ref_ccy={mp.ref_ccy} must be included in keep={keep}")
    new_currencies = [c for c in mp.currencies if c in keep_set]

    new_pairs = {
        k: v for k, v in mp.pairs.items()
        if (k[0] in keep_set and k[1] in keep_set)
    }

    def _filt(dct):
        return {k: v for k, v in dct.items() if k in keep_set}

    new_corr = {
        k: v for k, v in mp.corr.items()
        if (k[0] in keep_set and k[1] in keep_set)
    }

    d = len(new_currencies)
    kappa = mp.kappa
    if kappa is None:
        kappa = np.zeros((d, d))
    else:
        idx_old = {c: i for i, c in enumerate(mp.currencies)}
        idx_new = [idx_old[c] for c in new_currencies]
        kappa = kappa[np.ix_(idx_new, idx_new)]

    return ModelParams(
        currencies=new_currencies,
        ref_ccy=mp.ref_ccy,
        sigma=_filt(mp.sigma),
        corr=new_corr,
        k=_filt(mp.k),
        mu=_filt(mp.mu),
        gamma=mp.gamma,
        kappa=kappa,
        T_days=mp.T_days,
        pairs=new_pairs,
    )
