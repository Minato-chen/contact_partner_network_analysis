"""Cross-Topic Overview: Full 7-Variant Synthesis Plots.
Synthesizes Topic 1 (A18, B01, B10) and Topic 2 (A20, C04, B09, C22).

Figures:
1. Phase Behavior Overview: Slab density profiles & Csat comparison bar chart (with clean top legend).
2. Multivalent Trade-off Map: Partners vs Single-Partner Max Hetero Contacts vs Csat (with Anti-Collision).
3. Comprehensive Descriptor Comparison: Multi-panel bar plots for structural, contact, partner, and network metrics.

Strictly exports only PNG figures (300 DPI) as requested.
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

HERE = Path(__file__).resolve().parent
PLOT_ROOT = HERE.parent
ROOT = PLOT_ROOT.parent
sys.path.insert(0, str(PLOT_ROOT))

from common.style import apply_nature_style, save_figure, COLORS_ALL, VARIANT_INFO
from common.data_loader import load_summary, load_density_profile

ALL_VARIANTS = ["A18", "B01", "B10", "A20", "C04", "B09", "C22"]
SIM_MAP = {
    "A18": "A18_T260",
    "B01": "B01_T260",
    "B10": "B10_T260",
    "A20": "A20_T260",
    "C04": "C04_T260",
    "B09": "B09_T260",
    "C22": "C22_T260",
}


def plot_phase_behavior_overview(object_name="dense_chains", output_dir=None):
    """Plot slab density profile and dilute phase saturation concentration (Csat) with ample headroom and clean top legend."""
    apply_nature_style()
    df = load_summary(object_name).set_index("variant")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.4, 3.2), gridspec_kw={"width_ratios": [1.15, 1.0]})
    
    # Panel 1: Density profiles
    for var in ALL_VARIANTS:
        df_dens = load_density_profile(SIM_MAP[var])
        z = df_dens["z_nm"]
        rho = df_dens["mean_chain_mM"].rolling(window=5, center=True).mean().fillna(df_dens["mean_chain_mM"])
        ax1.plot(z, rho, label=var, color=COLORS_ALL[var], linewidth=1.1, alpha=0.9)
        
    ax1.set_xlabel("z (nm)")
    ax1.set_ylabel("Concentration (mM)")
    ax1.set_xlim(-60, 60)
    ax1.set_ylim(0, 16)
    ax1.set_title("A. Slab Density Profiles", loc="left", fontweight="bold")
    
    # Panel 2: Csat bar chart
    csat_vals = np.array([df.loc[v, "csat (mM)"] for v in ALL_VARIANTS])
    csat_errs = np.array([df.loc[v, "csat (mM) SEM"] for v in ALL_VARIANTS])
    x = np.arange(len(ALL_VARIANTS))
    
    bars = ax2.bar(x, csat_vals, yerr=csat_errs, capsize=3, color=[COLORS_ALL[v] for v in ALL_VARIANTS], edgecolor="black", linewidth=0.6, width=0.6)
    for bar, val, sem in zip(bars, csat_vals, csat_errs):
        ax2.text(bar.get_x() + bar.get_width()/2, val + sem + max(csat_vals)*0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=7.0)
        
    ax2.set_xticks(x)
    ax2.set_xticklabels(ALL_VARIANTS)
    ax2.set_ylabel(r"$C_{\mathrm{sat}}$ (mM)")
    ax2.set_title(r"B. Dilute Saturation Conc. ($C_{\mathrm{sat}}$)", loc="left", fontweight="bold")
    ax2.set_ylim(0, max(csat_vals + csat_errs) * 1.35)
    
    # Divider line between Topic 1 and Topic 2
    ax2.axvline(2.5, color="gray", linestyle="--", linewidth=0.7, alpha=0.6)
    y_top = max(csat_vals + csat_errs) * 1.25
    ax2.text(1.0, y_top, "Topic 1 (Blockiness)", ha="center", fontsize=6.8, color="#555555", fontweight="bold")
    ax2.text(4.5, y_top, "Topic 2 (Polarization)", ha="center", fontsize=6.8, color="#555555", fontweight="bold")
    
    # Place unified shared legend at the top outside subplots
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=len(ALL_VARIANTS), frameon=False, fontsize=7.5)
    plt.tight_layout(rect=[0, 0, 1, 0.88])
    if output_dir:
        save_figure(fig, "cross_topic_phase_behavior", output_dir)
    return fig


def plot_multivalent_tradeoff_map(object_name="dense_chains", output_dir=None):
    """Plot Partners vs Max Attract Contacts trade-off map colored by Csat, with anti-collision labels."""
    apply_nature_style()
    df = load_summary(object_name).set_index("simulation_id")
    
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    
    # Precise collision-free label offsets (dx, dy in points)
    offsets = {
        "A18": (10, -6),    # Right-down
        "B01": (10, 6),     # Right-up
        "B10": (10, 6),     # Right-up
        "A20": (-12, -18),  # Left-down
        "C04": (10, -14),   # Right-down
        "B09": (-28, 8),    # Left-up
        "C22": (10, 4),     # Right-up
    }
    
    for var in ALL_VARIANTS:
        sim = SIM_MAP[var]
        row = df.loc[sim]
        x = row["inter.max_attract_contacts_per_partner [threshold=1] (pairs/partner)"]
        y = row["inter.partner_count [threshold=1] (chains)"]
        csat = row["csat (mM)"]
        
        # Plot bubble point
        ax.scatter(x, y, s=130, color=COLORS_ALL[var], edgecolor="black", linewidth=0.8, zorder=5)
        
        dx, dy = offsets[var]
        ax.annotate(
            f"{var}\n({csat:.2f} mM)",
            (x, y),
            textcoords="offset points",
            xytext=(dx, dy),
            ha="right" if dx < 0 else "left",
            fontsize=7.5,
            fontweight="bold",
            color=COLORS_ALL[var],
            zorder=6,
        )

    # Highlight Fluid Multivalent LLPS region (Top-Left)
    ax.axvspan(0, 8.5, ymin=0.48, ymax=0.98, color="#E8F4F8", alpha=0.55, zorder=1)
    ax.text(4.2, 12.2, "Fluid Multivalent\nPercolation (Low Csat)", ha="center", fontsize=7.2, color="#1D70B8", fontweight="bold")

    # Highlight Electrostatic Zipper Trap region (Bottom-Right)
    ax.axvspan(15, 23.5, ymin=0.22, ymax=0.68, color="#FBEAEB", alpha=0.55, zorder=1)
    ax.text(19.5, 6.0, "Strong Dimer / Oligomer\nTrap (Elevated Csat)", ha="center", fontsize=7.2, color="#D95F02", fontweight="bold")

    # Highlight Neutral Tail Shielding (Bottom-Left)
    ax.text(5.5, 2.3, "Steric Tail Shielding & Self-Neutralization (C22)", ha="left", fontsize=7.0, color="#777777", style="italic")

    ax.set_xlabel("Single-partner max hetero contacts (pairs / partner)")
    ax.set_ylabel("Effective intermolecular partners (partners / chain)")
    ax.set_title("Multivalent Phase Separation Landscape\n(Low Csat = Fluid Network, High Csat = Trapped/Shielded)", fontsize=8.5, fontweight="bold", pad=8)
    
    ax.set_xlim(0, 24)
    ax.set_ylim(1.2, 13.0)
    plt.tight_layout()
    if output_dir:
        save_figure(fig, "cross_topic_multivalent_tradeoff_map", output_dir)
    return fig


def plot_descriptors_matrix(object_name="dense_chains", output_dir=None):
    """Plot multi-panel comparison of all key metrics across 7 variants without clipping."""
    apply_nature_style()
    df = load_summary(object_name).set_index("simulation_id")
    
    fig, axes = plt.subplots(2, 3, figsize=(9.2, 5.2), sharex=True)
    x = np.arange(len(ALL_VARIANTS))
    colors = [COLORS_ALL[v] for v in ALL_VARIANTS]
    
    has_network = "network.largest_component_fraction [threshold=1] (dimensionless)" in df.columns
    f_metric = (
        ("network.largest_component_fraction [threshold=1] (dimensionless)", "network.largest_component_fraction [threshold=1] (dimensionless) SEM", "F. Network max comp. ratio", 1, 2)
        if has_network else
        ("joint.inter_fraction (dimensionless)", "joint.inter_fraction (dimensionless) SEM", "F. Inter contact fraction", 1, 2)
    )

    metrics = [
        ("intra.rg_nm (nm)", "intra.rg_nm (nm) SEM", r"A. Chain $R_g$ (nm)", 0, 0),
        ("inter.total_contacts (pairs/chain)", "inter.total_contacts (pairs/chain) SEM", "B. Inter contacts / chain", 0, 1),
        ("inter.attract_contacts (pairs/chain)", "inter.attract_contacts (pairs/chain) SEM", "C. Inter hetero (+/-) contacts", 0, 2),
        ("inter.partner_count [threshold=1] (chains)", "inter.partner_count [threshold=1] (chains) SEM", "D. Intermolecular partners", 1, 0),
        ("inter.max_attract_contacts_per_partner [threshold=1] (pairs/partner)", "inter.max_attract_contacts_per_partner [threshold=1] (pairs/partner) SEM", "E. Max hetero / partner", 1, 1),
        f_metric,
    ]
    
    for key, sem_key, title, r, c in metrics:
        ax = axes[r, c]
        vals = np.array([df.loc[SIM_MAP[v], key] for v in ALL_VARIANTS])
        errs = np.array([df.loc[SIM_MAP[v], sem_key] for v in ALL_VARIANTS])
        ax.bar(x, vals, yerr=errs, capsize=2.5, color=colors, edgecolor="black", linewidth=0.5, width=0.6)
        ax.set_title(title, loc="left", fontweight="bold", fontsize=7.8)
        ax.set_ylim(0, max(vals + errs) * 1.30)
        # Divider line between topics
        ax.axvline(2.5, color="gray", linestyle=":", linewidth=0.6, alpha=0.5)
        
    for c in range(3):
        axes[1, c].set_xticks(x)
        axes[1, c].set_xticklabels(ALL_VARIANTS, fontsize=7.5)
        
    plt.tight_layout()
    if output_dir:
        save_figure(fig, "cross_topic_descriptors_matrix", output_dir)
    return fig


def main(object_name="dense_chains"):
    out = HERE / "figures" / object_name
    print(f"Generating Cross-Topic Overview figures ({object_name})...")
    plot_phase_behavior_overview(object_name, out)
    plot_multivalent_tradeoff_map(object_name, out)
    plot_descriptors_matrix(object_name, out)


if __name__ == "__main__":
    for obj in ["dense_chains", "central_chain"]:
        main(obj)
