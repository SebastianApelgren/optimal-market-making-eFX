from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt

from .model import BP, DAY_SECONDS
from .policy import MMResult


def _y_vector(mp, overrides_musd: Dict[str, float]) -> np.ndarray:
    y = np.zeros(len(mp.currencies))
    for c, v in overrides_musd.items():
        y[mp.currencies.index(c)] = v
    return y


def plot_top_of_book_quotes_vs_inventory(res: MMResult,
                                         tier_idx: int,
                                         z_musd: float,
                                         inventory_ccy: str,
                                         inventory_grid_musd: np.ndarray,
                                         pairs_to_plot: List[Tuple[str, str]]) -> None:
    """Plot bid/ask quotes vs inventory for selected currency pairs."""
    mp = res.mp

    plt.figure(figsize=(9, 5))

    for X, Y in pairs_to_plot:
        bid_bps = []
        ask_bps = []
        for inv in inventory_grid_musd:
            y = _y_vector(mp, {inventory_ccy: inv})

            delta_bid = res.markup(tier_idx, ccy_pay=X, ccy_sell=Y, z_musd=z_musd, y=y)
            delta_ask = -res.markup(tier_idx, ccy_pay=Y, ccy_sell=X, z_musd=z_musd, y=y)

            bid_bps.append(delta_bid / BP)
            ask_bps.append(delta_ask / BP)

        bid_bps = np.array(bid_bps)
        ask_bps = np.array(ask_bps)

        plt.plot(inventory_grid_musd, bid_bps, label=f"{X}{Y} bid")
        plt.plot(inventory_grid_musd, ask_bps, linestyle="--", label=f"{X}{Y} ask")

    plt.axhline(0.0, linewidth=0.8)
    plt.xlabel(f"{inventory_ccy} Inventory (M$)")
    plt.ylabel("Bid and Ask Quotes (bps)")
    plt.title("Approximate Optimal Top-of-Book Quotes")
    plt.legend(ncol=2)
    plt.grid(True, alpha=0.3)
    plt.show()


def plot_hedge_rates_vs_inventory(res: MMResult,
                                 inventory_ccy: str,
                                 inventory_grid_musd: np.ndarray,
                                 ordered_pairs_to_plot: List[Tuple[str, str]]) -> None:
    """Plot hedging rates (M$/s) vs inventory for selected currency directions."""
    mp = res.mp

    plt.figure(figsize=(9, 5))

    for buy_ccy, sell_ccy in ordered_pairs_to_plot:
        rates = []
        for inv in inventory_grid_musd:
            y = _y_vector(mp, {inventory_ccy: inv})
            xi_per_day = res.hedge_rate(buy_ccy, sell_ccy, y=y)
            rates.append(xi_per_day / DAY_SECONDS)

        plt.plot(inventory_grid_musd, np.array(rates), label=f"{buy_ccy}{sell_ccy}")

    plt.axhline(0.0, linewidth=0.8)
    plt.xlabel(f"{inventory_ccy} Inventory (M$)")
    plt.ylabel("Execution Rate (M$ / s)")
    plt.title("Approximate Optimal Externalization (Hedging) Rates")
    plt.legend(ncol=2)
    plt.grid(True, alpha=0.3)
    plt.show()
