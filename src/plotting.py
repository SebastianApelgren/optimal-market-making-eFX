from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import matplotlib.pyplot as plt

from .model import BP, DAY_SECONDS, canon_pair
from .policy import MMResult


# ---------------------------------------------------------------------------
# Style configuration (matches the publication-quality settings from the
# original ode_chapter_figures.py script)
# ---------------------------------------------------------------------------
_RC_PARAMS = {
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "axes.grid": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "lines.linewidth": 1.4,
    "text.usetex": False,
    "mathtext.fontset": "dejavusans",
}

# Colorblind-friendly palette (Wong 2011 / IBM Design Library style).
COLORS = [
    "#0072B2",  # blue
    "#D55E00",  # vermilion
    "#009E73",  # bluish green
    "#E69F00",  # orange
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
    "#882255",  # wine
    "#44AA99",  # teal
]


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
                                         pairs_to_plot: List[Tuple[str, str]],
                                         save_path: str | None = None) -> None:
    """Plot bid/ask quotes vs inventory for selected currency pairs."""
    mp = res.mp

    with plt.rc_context(_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(6, 4))

        for idx, (X, Y) in enumerate(pairs_to_plot):
            bid_bps = np.empty(len(inventory_grid_musd))
            ask_bps = np.empty(len(inventory_grid_musd))

            for k, inv in enumerate(inventory_grid_musd):
                y = _y_vector(mp, {inventory_ccy: inv})

                # Quote relative to mid: bid = mid - δ_bid, ask = mid + δ_ask
                delta_bid = res.markup(tier_idx, ccy_pay=X, ccy_sell=Y, z_musd=z_musd, y=y)
                delta_ask = res.markup(tier_idx, ccy_pay=Y, ccy_sell=X, z_musd=z_musd, y=y)

                bid_bps[k] = -delta_bid / BP
                ask_bps[k] = delta_ask / BP

            color = COLORS[idx % len(COLORS)]
            ax.plot(inventory_grid_musd, ask_bps, color=color, linestyle="-",
                    label=f"{X}{Y} ask")
            ax.plot(inventory_grid_musd, bid_bps, color=color, linestyle="--",
                    label=f"{X}{Y} bid")

        ax.axhline(0.0, color="grey", linewidth=0.6, zorder=0)
        ax.set_xlabel(f"{inventory_ccy} inventory (M\\$)")
        ax.set_ylabel("Quote relative to mid (bps)")
        ax.legend(ncol=2, frameon=True, fancybox=False, edgecolor="0.8",
                  loc="best")

        if save_path is not None:
            fig.savefig(save_path)
        plt.show()


def plot_hedge_rates_vs_inventory(res: MMResult,
                                 inventory_ccy: str,
                                 inventory_grid_musd: np.ndarray,
                                 ordered_pairs_to_plot: List[Tuple[str, str]],
                                 save_path: str | None = None) -> None:
    """Plot hedging rates (M$/s) vs inventory for selected currency directions."""
    mp = res.mp

    with plt.rc_context(_RC_PARAMS):
        fig, ax = plt.subplots(figsize=(6, 4))

        for idx, (buy_ccy, sell_ccy) in enumerate(ordered_pairs_to_plot):
            rates = np.empty(len(inventory_grid_musd))
            for k, inv in enumerate(inventory_grid_musd):
                y = _y_vector(mp, {inventory_ccy: inv})
                xi_per_day = res.hedge_rate(buy_ccy, sell_ccy, y=y)
                rates[k] = xi_per_day / DAY_SECONDS

            color = COLORS[idx % len(COLORS)]
            ax.plot(inventory_grid_musd, rates, color=color,
                    label=f"{buy_ccy}{sell_ccy}")

        ax.axhline(0.0, color="grey", linewidth=0.6, zorder=0)
        ax.set_xlabel(f"{inventory_ccy} inventory (M\\$)")
        ax.set_ylabel("Execution rate (M\\$/s)")
        ax.legend(ncol=2, frameon=True, fancybox=False, edgecolor="0.8",
                  fontsize=7.5, loc="lower left")

        # --- Correlation matrix inset ---
        non_ref = [c for c in mp.currencies if c != mp.ref_ccy]
        d = len(non_ref)
        corr_mat = np.eye(d)
        for i in range(d):
            for j in range(i + 1, d):
                rho = mp.corr.get(canon_pair(non_ref[i], non_ref[j]), 0.0)
                corr_mat[i, j] = rho
                corr_mat[j, i] = rho

        ax_inset = fig.add_axes([0.64, 0.58, 0.28, 0.34])
        ax_inset.imshow(corr_mat, cmap="RdBu_r", vmin=-1, vmax=1,
                        aspect="equal")
        ax_inset.set_xticks(range(d))
        ax_inset.set_xticklabels(non_ref, fontsize=7)
        ax_inset.set_yticks(range(d))
        ax_inset.set_yticklabels(non_ref, fontsize=7)
        ax_inset.tick_params(length=0)
        for i in range(d):
            for j in range(d):
                val = corr_mat[i, j]
                text_color = "white" if abs(val) > 0.6 else "black"
                ax_inset.text(j, i, f"{val:.1f}", ha="center", va="center",
                              fontsize=7, color=text_color)
        ax_inset.set_title("Correlation", fontsize=8, pad=3)

        if save_path is not None:
            fig.savefig(save_path)
        plt.show()
