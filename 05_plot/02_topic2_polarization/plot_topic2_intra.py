"""Topic 2 (Polarization & Neutral Tail: A20, C04, B09, C22) Intramolecular Analysis Plots.
Features:
1. Intra contacts, Rg, and long-range fraction (|i-j| >= 12).
2. 1D Residue contact probability profile & Loop size distribution.
3. Intrachain topological network metrics (Largest component & Global efficiency).
Exports strictly PNG figures (300 DPI).
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PLOT_ROOT = HERE.parent
ROOT = PLOT_ROOT.parent
sys.path.insert(0, str(PLOT_ROOT))

from common.style import apply_nature_style, save_figure, COLORS_TOPIC2
from common.data_loader import (
    load_summary,
    filter_topic,
    load_residue_profile,
    load_loop_size_distribution,
)

VARIANTS = ["A20", "C04", "B09", "C22"]
SIM_MAP = {"A20": "A20_T260", "C04": "C04_T260", "B09": "B09_T260", "C22": "C22_T260"}


def plot_intra_rg_and_contacts(object_name="dense_chains", output_dir=None):
    """Plot intra contact count, Rg, and long-range fraction with ample headroom for Topic 2."""
    apply_nature_style()
    df = filter_topic(load_summary(object_name), topic=2)
    
    fig, axes = plt.subplots(1, 3, figsize=(8.6, 2.6))
    
    # Panel A: Intra Contact Count
    ax = axes[0]
    vals = df["intra.contact_count (pairs/chain)"].to_numpy()
    sems = df["intra.contact_count (pairs/chain) SEM"].to_numpy()
    bars = ax.bar(VARIANTS, vals, yerr=sems, capsize=3, color=[COLORS_TOPIC2[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Intra contacts (pairs/chain)")
    ax.set_title("A. Intrachain contacts", loc="left", fontweight="bold")
    ax.set_ylim(0, max(vals + sems) * 1.32)
    for bar, val, sem in zip(bars, vals, sems):
        ax.text(bar.get_x() + bar.get_width()/2, val + sem + max(vals)*0.03, f"{val:.1f}", ha="center", va="bottom", fontsize=7.5)

    # Panel B: Rg (nm)
    ax = axes[1]
    vals = df["intra.rg_nm (nm)"].to_numpy()
    sems = df["intra.rg_nm (nm) SEM"].to_numpy()
    bars = ax.bar(VARIANTS, vals, yerr=sems, capsize=3, color=[COLORS_TOPIC2[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax.set_ylabel(r"$R_g$ (nm)")
    ax.set_title(r"B. Radius of gyration $R_g$", loc="left", fontweight="bold")
    ax.set_ylim(0, max(vals + sems) * 1.32)
    for bar, val, sem in zip(bars, vals, sems):
        ax.text(bar.get_x() + bar.get_width()/2, val + sem + max(vals)*0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=7.5)

    # Panel C: Long-range fraction (|i-j| >= 12)
    ax = axes[2]
    vals = df["intra.long_range_fraction (dimensionless)"].to_numpy()
    bars = ax.bar(VARIANTS, vals, color=[COLORS_TOPIC2[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Long-range contact fraction")
    ax.set_title(r"C. Long-range fraction ($|i-j|\geq 12$)", loc="left", fontweight="bold")
    ax.set_ylim(0, max(vals) * 1.35)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + max(vals)*0.03, f"{val:.3f}", ha="center", va="bottom", fontsize=7.5)

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic2_intra_metrics", output_dir)
    return fig


def plot_residue_and_loop_profile(object_name="dense_chains", output_dir=None):
    """Plot 1D residue contact probability profile & loop size distribution for Topic 2."""
    apply_nature_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 2.7))
    
    # Panel 1: 1D residue profile
    for var in VARIANTS:
        df_res = load_residue_profile(SIM_MAP[var], object_name)
        ax1.plot(df_res["residue_index"], df_res["fraction"], label=var, color=COLORS_TOPIC2[var], linewidth=1.1)
        
    ax1.set_xlabel("Residue index (1–129)")
    ax1.set_ylabel("Mean contact fraction")
    ax1.set_title("A. Residue contact profile", loc="left", fontweight="bold")
    ax1.set_xlim(1, 129)
    ax1.set_ylim(0, 0.45)
    # Add annotation for C22 neutral tail (residues 1-60)
    ax1.axvspan(1, 60, color="#F4A261", alpha=0.10, label="C22 Neutral Tail (1–60)")
    ax1.legend(loc="upper right", frameon=False, fontsize=6.5)

    # Panel 2: Loop size distribution (|i-j|)
    for var in VARIANTS:
        loop_dist = load_loop_size_distribution(SIM_MAP[var], object_name)
        seq_sep = np.arange(len(loop_dist))
        ax2.plot(seq_sep[2:], loop_dist[2:], label=var, color=COLORS_TOPIC2[var], linewidth=1.1)

    ax2.set_xlabel("Sequence separation $|i - j|$")
    ax2.set_ylabel("Contact probability")
    ax2.set_title("B. Loop size distribution", loc="left", fontweight="bold")
    ax2.set_xlim(2, 128)
    ax2.set_yscale("log")
    ax2.set_ylim(1e-4, 0.5)
    ax2.legend(loc="upper right", frameon=False, fontsize=7)

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic2_residue_and_loop_profile", output_dir)
    return fig


def plot_intra_network_metrics(object_name="dense_chains", output_dir=None):
    """Plot intrachain largest component fraction & global efficiency for Topic 2."""
    apply_nature_style()
    df = filter_topic(load_summary(object_name), topic=2)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 2.5))
    
    # Panel 1: Largest component fraction
    vals = df["intra.largest_component_fraction (dimensionless)"].to_numpy()
    sems = df["intra.largest_component_fraction (dimensionless) SEM"].to_numpy()
    bars = ax1.bar(VARIANTS, vals, yerr=sems, capsize=3, color=[COLORS_TOPIC2[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax1.set_ylabel("Largest component fraction")
    ax1.set_title("A. Intrachain largest component", loc="left", fontweight="bold")
    ax1.set_ylim(0, 1.25)
    for bar, val, sem in zip(bars, vals, sems):
        ax1.text(bar.get_x() + bar.get_width()/2, val + sem + 0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=7.5)

    # Panel 2: Global efficiency
    vals = df["intra.global_efficiency (dimensionless)"].to_numpy()
    sems = df["intra.global_efficiency (dimensionless) SEM"].to_numpy()
    bars = ax2.bar(VARIANTS, vals, yerr=sems, capsize=3, color=[COLORS_TOPIC2[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax2.set_ylabel("Global efficiency")
    ax2.set_title("B. Intrachain network efficiency", loc="left", fontweight="bold")
    ax2.set_ylim(0, max(vals + sems) * 1.32)
    for bar, val, sem in zip(bars, vals, sems):
        ax2.text(bar.get_x() + bar.get_width()/2, val + sem + max(vals)*0.03, f"{val:.3f}", ha="center", va="bottom", fontsize=7.5)

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic2_intra_network_metrics", output_dir)
    return fig


def main(object_name="dense_chains"):
    out = HERE / "figures" / object_name
    print(f"Generating Topic 2 Intra figures ({object_name})...")
    plot_intra_rg_and_contacts(object_name, out)
    plot_residue_and_loop_profile(object_name, out)
    plot_intra_network_metrics(object_name, out)


if __name__ == "__main__":
    for obj in ["dense_chains", "central_chain"]:
        main(obj)
