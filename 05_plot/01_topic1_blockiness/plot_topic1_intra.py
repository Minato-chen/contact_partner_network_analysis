"""Topic 1 (Blockiness: A18, B01, B10) Intramolecular Analysis Plots.
Features:
1. Intra contacts, Rg, and long-range fraction (|i-j| >= 12).
2. 1D Residue contact probability profile & Loop size distribution (full 0-1.0 range, block SD shading, top legend).
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

from common.style import apply_nature_style, save_figure, COLORS_TOPIC1
from common.data_loader import (
    load_summary,
    filter_topic,
    load_residue_profile,
    load_loop_size_distribution,
)

VARIANTS = ["A18", "B01", "B10"]
SIM_MAP = {"A18": "A18_T260", "B01": "B01_T260", "B10": "B10_T260"}


def plot_intra_rg_and_contacts(object_name="dense_chains", output_dir=None):
    """Plot intra contact count, Rg, and long-range fraction with ample headroom."""
    apply_nature_style()
    df = filter_topic(load_summary(object_name), topic=1)
    
    fig, axes = plt.subplots(1, 3, figsize=(8.0, 2.6))
    
    # Panel A: Intra Contact Count
    ax = axes[0]
    vals = df["intra.contact_count (pairs/chain)"].to_numpy()
    sems = df["intra.contact_count (pairs/chain) SEM"].to_numpy()
    bars = ax.bar(VARIANTS, vals, yerr=sems, capsize=3, color=[COLORS_TOPIC1[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Intra contacts (pairs/chain)")
    ax.set_title("A. Intrachain contacts", loc="left", fontweight="bold")
    ax.set_ylim(0, max(vals + sems) * 1.35)
    for bar, val, sem in zip(bars, vals, sems):
        ax.text(bar.get_x() + bar.get_width()/2, val + sem + max(vals)*0.03, f"{val:.1f}", ha="center", va="bottom", fontsize=7.5)

    # Panel B: Rg (nm)
    ax = axes[1]
    vals = df["intra.rg_nm (nm)"].to_numpy()
    sems = df["intra.rg_nm (nm) SEM"].to_numpy()
    bars = ax.bar(VARIANTS, vals, yerr=sems, capsize=3, color=[COLORS_TOPIC1[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax.set_ylabel(r"$R_g$ (nm)")
    ax.set_title(r"B. Radius of gyration $R_g$", loc="left", fontweight="bold")
    ax.set_ylim(0, max(vals + sems) * 1.35)
    for bar, val, sem in zip(bars, vals, sems):
        ax.text(bar.get_x() + bar.get_width()/2, val + sem + max(vals)*0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=7.5)

    # Panel C: Long-range fraction (|i-j| >= 12)
    ax = axes[2]
    vals = df["intra.long_range_fraction (dimensionless)"].to_numpy()
    bars = ax.bar(VARIANTS, vals, color=[COLORS_TOPIC1[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax.set_ylabel("Long-range contact fraction")
    ax.set_title(r"C. Long-range fraction ($|i-j|\geq 12$)", loc="left", fontweight="bold")
    ax.set_ylim(0, max(vals) * 1.35)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + max(vals)*0.03, f"{val:.3f}", ha="center", va="bottom", fontsize=7.5)

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic1_intra_metrics", output_dir)
    return fig


def plot_residue_and_loop_profile(object_name="dense_chains", output_dir=None):
    """Plot 1D residue contact probability profile (full 0-1.0 range, block SD shading) & loop size distribution."""
    apply_nature_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 2.9))
    
    # Panel 1: 1D residue profile (0-1.0 probability range with block SD)
    for var in VARIANTS:
        df_res = load_residue_profile(SIM_MAP[var], object_name)
        x = (df_res["residue_index"] + 1).to_numpy(dtype=float)
        y = df_res["mean"].to_numpy(dtype=float)
        sd = df_res["std"].to_numpy(dtype=float)
        
        ax1.plot(x, y, label=var, color=COLORS_TOPIC1[var], linewidth=1.2)
        ax1.fill_between(x, np.maximum(0.0, y - sd), np.minimum(1.0, y + sd), color=COLORS_TOPIC1[var], alpha=0.18)
        
    ax1.set_xlabel("Residue index (1–129)")
    ax1.set_ylabel("Contact probability (mean ± block SD)")
    ax1.set_title("A. Residue contact profile", loc="left", fontweight="bold")
    ax1.set_xlim(1, 129)
    ax1.set_ylim(0, 1.0)  # Correct full probability range

    # Panel 2: Loop size distribution (|i-j|)
    for var in VARIANTS:
        loop_dist = load_loop_size_distribution(SIM_MAP[var], object_name)
        seq_sep = np.arange(len(loop_dist))
        ax2.plot(seq_sep[2:], loop_dist[2:], label=var, color=COLORS_TOPIC1[var], linewidth=1.2)

    ax2.set_xlabel("Sequence separation $|i - j|$")
    ax2.set_ylabel("Contact probability")
    ax2.set_title("B. Loop size distribution", loc="left", fontweight="bold")
    ax2.set_xlim(2, 128)
    ax2.set_yscale("log")
    ax2.set_ylim(1e-4, 0.5)

    # Place single shared legend at the top outside plots to avoid any curve collision
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=len(VARIANTS), frameon=False, fontsize=8.0)
    plt.tight_layout(rect=[0, 0, 1, 0.90])
    
    if output_dir:
        save_figure(fig, "topic1_residue_and_loop_profile", output_dir)
    return fig


def plot_intra_network_metrics(object_name="dense_chains", output_dir=None):
    """Plot intrachain largest component fraction & global efficiency."""
    apply_nature_style()
    df = filter_topic(load_summary(object_name), topic=1)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.0, 2.6))
    
    # Panel 1: Largest component fraction
    vals = df["intra.largest_component_fraction (dimensionless)"].to_numpy()
    sems = df["intra.largest_component_fraction (dimensionless) SEM"].to_numpy()
    bars = ax1.bar(VARIANTS, vals, yerr=sems, capsize=3, color=[COLORS_TOPIC1[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax1.set_ylabel("Largest component fraction")
    ax1.set_title("A. Intrachain largest component", loc="left", fontweight="bold")
    ax1.set_ylim(0, 1.25)
    for bar, val, sem in zip(bars, vals, sems):
        ax1.text(bar.get_x() + bar.get_width()/2, val + sem + 0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=7.5)

    # Panel 2: Global efficiency
    vals = df["intra.global_efficiency (dimensionless)"].to_numpy()
    sems = df["intra.global_efficiency (dimensionless) SEM"].to_numpy()
    bars = ax2.bar(VARIANTS, vals, yerr=sems, capsize=3, color=[COLORS_TOPIC1[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax2.set_ylabel("Global efficiency")
    ax2.set_title("B. Intrachain network efficiency", loc="left", fontweight="bold")
    ax2.set_ylim(0, max(vals + sems) * 1.35)
    for bar, val, sem in zip(bars, vals, sems):
        ax2.text(bar.get_x() + bar.get_width()/2, val + sem + max(vals)*0.03, f"{val:.3f}", ha="center", va="bottom", fontsize=7.5)

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic1_intra_network_metrics", output_dir)
    return fig


def main(object_name="dense_chains"):
    out = HERE / "figures" / object_name
    print(f"Generating Topic 1 Intra figures ({object_name})...")
    plot_intra_rg_and_contacts(object_name, out)
    plot_residue_and_loop_profile(object_name, out)
    plot_intra_network_metrics(object_name, out)


if __name__ == "__main__":
    for obj in ["dense_chains", "central_chain"]:
        main(obj)
