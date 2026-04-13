"""
Generate publication-quality figures for the ODE chapter of the thesis report.

Reproduces Figures 1, 2, and 4 from Barzykin, Bergault & Gueant (2023),
arXiv:2207.04100v4, using the Table 1 parameters.

Usage:
    .venv/Scripts/python report/figures/ode_chapter_figures.py
"""
from __future__ import annotations

import sys
import os

# Ensure the project root is on the path so `from src import ...` works.
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, PROJECT_ROOT)

import logging
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for PDF output
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt

from src import (
    BP, DAY_SECONDS,
    build_paper_example_params, restrict_currencies,
    run_multicurrency_mm,
    simulate_inventory_path_tau_leap,
    build_Sigma,
    canon_pair,
)

# ---------------------------------------------------------------------------
# Output directory
# ---------------------------------------------------------------------------
OUTDIR = os.path.join(PROJECT_ROOT, "report", "figures")
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Style configuration
# ---------------------------------------------------------------------------
plt.rcParams.update({
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
})

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


def _y_vector(mp, overrides_musd):
    """Build an inventory vector with specified overrides."""
    y = np.zeros(len(mp.currencies))
    for c, v in overrides_musd.items():
        y[mp.currencies.index(c)] = v
    return y


# ===========================================================================
# Figure 1: Top-of-book quotes vs GBP inventory
# ===========================================================================
def figure_quotes_vs_inventory(res, outpath):
    """Bid/ask quotes for GBPUSD, EURUSD, EURGBP vs GBP inventory."""
    mp = res.mp
    inv_grid = np.linspace(-100, 100, 401)

    pairs_to_plot = [("GBP", "USD"), ("EUR", "USD"), ("EUR", "GBP")]
    pair_labels = ["GBPUSD", "EURUSD", "EURGBP"]

    fig, ax = plt.subplots(figsize=(6, 4))

    for idx, (X, Y) in enumerate(pairs_to_plot):
        bid_bps = np.empty_like(inv_grid)
        ask_bps = np.empty_like(inv_grid)

        for k, inv in enumerate(inv_grid):
            y = _y_vector(mp, {"GBP": inv})
            delta_bid = res.markup(0, ccy_pay=X, ccy_sell=Y, z_musd=1.0, y=y)
            delta_ask = -res.markup(0, ccy_pay=Y, ccy_sell=X, z_musd=1.0, y=y)
            bid_bps[k] = delta_bid / BP
            ask_bps[k] = delta_ask / BP

        color = COLORS[idx]
        ax.plot(inv_grid, bid_bps, color=color, linestyle="-",
                label=f"{pair_labels[idx]} bid")
        ax.plot(inv_grid, ask_bps, color=color, linestyle="--",
                label=f"{pair_labels[idx]} ask")

    ax.axhline(0.0, color="grey", linewidth=0.6, zorder=0)
    ax.set_xlabel("GBP inventory (M\\$)")
    ax.set_ylabel("Quote (bps)")
    ax.legend(ncol=2, frameon=True, fancybox=False, edgecolor="0.8",
              loc="best")

    fig.savefig(outpath, format="pdf")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ===========================================================================
# Figure 2: Hedge rates vs EUR inventory (with correlation inset)
# ===========================================================================
def figure_hedge_rates_vs_inventory(res, outpath):
    """Optimal hedging rates for all 10 currency pairs vs EUR inventory."""
    mp = res.mp
    inv_grid = np.linspace(-100, 100, 401)

    ordered_pairs = [
        ("EUR", "USD"), ("GBP", "USD"), ("CHF", "USD"), ("JPY", "USD"),
        ("EUR", "GBP"), ("EUR", "CHF"), ("EUR", "JPY"),
        ("GBP", "CHF"), ("GBP", "JPY"),
        ("CHF", "JPY"),
    ]

    fig, ax = plt.subplots(figsize=(6, 4))

    for idx, (buy_ccy, sell_ccy) in enumerate(ordered_pairs):
        rates = np.empty_like(inv_grid)
        for k, inv in enumerate(inv_grid):
            y = _y_vector(mp, {"EUR": inv})
            xi_day = res.hedge_rate(buy_ccy, sell_ccy, y=y)
            rates[k] = xi_day / DAY_SECONDS

        color = COLORS[idx % len(COLORS)]
        ax.plot(inv_grid, rates, color=color, label=f"{buy_ccy}{sell_ccy}")

    ax.axhline(0.0, color="grey", linewidth=0.6, zorder=0)
    ax.set_xlabel("EUR inventory (M\\$)")
    ax.set_ylabel("Execution rate (M\\$/s)")
    ax.legend(ncol=2, frameon=True, fancybox=False, edgecolor="0.8",
              fontsize=7.5, loc="lower left")

    # --- Correlation matrix inset ---
    # Build the correlation matrix for the non-reference currencies.
    ccy = mp.currencies
    d = len(ccy)
    non_ref = [c for c in ccy if c != mp.ref_ccy]
    d_nr = len(non_ref)
    corr_matrix = np.eye(d_nr)
    for i in range(d_nr):
        for j in range(i + 1, d_nr):
            rho = mp.corr.get(canon_pair(non_ref[i], non_ref[j]), 0.0)
            corr_matrix[i, j] = rho
            corr_matrix[j, i] = rho

    # Position inset in upper-right area (no overlap with lower-left legend)
    ax_inset = fig.add_axes([0.64, 0.58, 0.28, 0.34])
    im = ax_inset.imshow(corr_matrix, cmap="RdBu_r", vmin=-1, vmax=1,
                         aspect="equal")
    ax_inset.set_xticks(range(d_nr))
    ax_inset.set_xticklabels(non_ref, fontsize=7)
    ax_inset.set_yticks(range(d_nr))
    ax_inset.set_yticklabels(non_ref, fontsize=7)
    ax_inset.tick_params(length=0)

    # Annotate correlation values
    for i in range(d_nr):
        for j in range(d_nr):
            val = corr_matrix[i, j]
            text_color = "white" if abs(val) > 0.6 else "black"
            ax_inset.text(j, i, f"{val:.1f}", ha="center", va="center",
                          fontsize=7, color=text_color)

    ax_inset.set_title("Correlation", fontsize=8, pad=3)

    fig.savefig(outpath, format="pdf")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ===========================================================================
# Figure 3: Inventory distribution from MC simulation (3-currency case)
# ===========================================================================
def figure_inventory_distribution(outpath, T_sec=500_000, dt_sec=5.0, seed=42):
    """2D histogram of EUR/GBP inventories with risk contours."""
    # 3-currency model
    mp_full = build_paper_example_params()
    mp3 = restrict_currencies(mp_full, ["USD", "EUR", "GBP"])
    res3 = run_multicurrency_mm(mp3)

    print(f"  Running MC simulation: T={T_sec:.0f}s, dt={dt_sec}s ...")
    times, Ypath = simulate_inventory_path_tau_leap(
        res3, T_sec=T_sec, dt_sec=dt_sec, seed=seed, include_hedging=True,
    )

    eur_idx = mp3.currencies.index("EUR")
    gbp_idx = mp3.currencies.index("GBP")

    eur_inv = Ypath[:, eur_idx]
    gbp_inv = Ypath[:, gbp_idx]

    # Discard burn-in: first 10% of samples
    burn = len(eur_inv) // 10
    eur_inv = eur_inv[burn:]
    gbp_inv = gbp_inv[burn:]

    # Build the Sigma matrix for risk contours
    Sigma = build_Sigma(mp3)
    gamma = mp3.gamma

    fig, ax = plt.subplots(figsize=(6, 5.5))

    # 2D histogram
    # Determine range from data for the histogram
    lim = max(np.percentile(np.abs(eur_inv), 99.5),
              np.percentile(np.abs(gbp_inv), 99.5))
    lim = max(lim, 50)  # at least +/-50

    h, xedges, yedges, im = ax.hist2d(
        eur_inv, gbp_inv,
        bins=120,
        range=[[-lim, lim], [-lim, lim]],
        cmap="Blues",
        density=True,
        rasterized=True,
    )

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label("Density", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    # Risk contours: gamma/2 * y^T Sigma y = const
    # We plot in the (EUR, GBP) subspace.
    # The quadratic form is y^T (gamma/2 * Sigma) y.
    # For the non-reference currencies only: Sigma_nr is the
    # (EUR, GBP) block of the Sigma matrix.
    nr_idx = [eur_idx, gbp_idx]
    Sigma_nr = Sigma[np.ix_(nr_idx, nr_idx)]
    Q = (gamma / 2.0) * Sigma_nr

    # Contour grid
    yy = np.linspace(-lim, lim, 300)
    YY_eur, YY_gbp = np.meshgrid(yy, yy)
    risk_vals = (Q[0, 0] * YY_eur**2 + 2.0 * Q[0, 1] * YY_eur * YY_gbp
                 + Q[1, 1] * YY_gbp**2)

    # Choose contour levels that span the distribution.
    # At y=20 M$ the quadratic form is O(0.1-0.3), so use levels
    # that cover the range where most of the density lies.
    contour_levels = np.array([0.01, 0.05, 0.1, 0.2, 0.5])
    cs = ax.contour(YY_eur, YY_gbp, risk_vals, levels=contour_levels,
                    colors="0.3", linewidths=0.8, linestyles="--")
    ax.clabel(cs, fmt="%.2f", fontsize=7, inline=True)

    ax.set_xlabel("EUR inventory (M\\$)")
    ax.set_ylabel("GBP inventory (M\\$)")
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")

    # Add a small annotation for the contours
    ax.annotate(
        r"$\gamma \mathbf{y}^T \Sigma \mathbf{y} / 2$",
        xy=(0.03, 0.97), xycoords="axes fraction",
        fontsize=9, ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9),
    )

    fig.savefig(outpath, format="pdf")
    plt.close(fig)
    print(f"  Saved {outpath}")


# ===========================================================================
# Main
# ===========================================================================
def main():
    print("Building 5-currency model and solving Riccati ODE ...")
    mp = build_paper_example_params()
    res = run_multicurrency_mm(mp)
    print(f"  Currencies: {res.mp.currencies}")
    print(f"  Max |A0| = {np.max(np.abs(res.A0)):.6e}")

    print("\nFigure 1: Quotes vs GBP inventory ...")
    figure_quotes_vs_inventory(
        res, os.path.join(OUTDIR, "fig_ode_quotes_vs_inventory.pdf"))

    print("\nFigure 2: Hedge rates vs EUR inventory ...")
    figure_hedge_rates_vs_inventory(
        res, os.path.join(OUTDIR, "fig_ode_hedge_rates_vs_inventory.pdf"))

    print("\nFigure 3: Inventory distribution (MC simulation) ...")
    figure_inventory_distribution(
        os.path.join(OUTDIR, "fig_ode_inventory_distribution.pdf"),
        T_sec=500_000, dt_sec=5.0, seed=42)

    print("\nAll figures generated successfully.")


if __name__ == "__main__":
    main()
