"""Topic 1 (Blockiness: A18, B01, B10) Intermolecular Analysis Plots.
Features:
1. Microscopic charge decomposition (Attract, Repulse, Neutral).
2. Multivalent coordination vs single-partner peak adhesion (Partner count vs Max attract contacts).
3. Partner count probability mass function (PMF) distribution curves with block SD shading.
4. Dense network percolation and global efficiency across thresholds.
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

from common.style import apply_nature_style, save_figure, COLORS_TOPIC1, CHARGE_COLORS
from common.data_loader import load_summary, filter_topic, load_partner_distribution

VARIANTS = ["A18", "B01", "B10"]
SIM_MAP = {"A18": "A18_T260", "B01": "B01_T260", "B10": "B10_T260"}


def plot_inter_charge_decomposition(object_name="dense_chains", output_dir=None):
    """Plot total intermolecular contacts decomposed into attract/repulse/neutral, with attract fraction."""
    apply_nature_style()
    df = filter_topic(load_summary(object_name), topic=1)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.4, 1.0]})
    
    x = np.arange(len(VARIANTS))
    width = 0.55
    
    attract = df["inter.attract_contacts (pairs/chain)"].to_numpy()
    repulse = df["inter.repulse_contacts (pairs/chain)"].to_numpy()
    neutral = df["inter.neutral_contacts (pairs/chain)"].to_numpy()
    totals = df["inter.total_contacts (pairs/chain)"].to_numpy()
    fracs = df["inter.attract_fraction (dimensionless)"].to_numpy() * 100
    frac_sems = df["inter.attract_fraction (dimensionless) SEM"].to_numpy() * 100
    
    # Panel 1: Stacked Bar Chart with dedicated electrostatic colors
    ax1.bar(x, attract, width, label="Opposite-charge (+/-)", color=CHARGE_COLORS["attract"], edgecolor="black", linewidth=0.5)
    ax1.bar(x, repulse, width, bottom=attract, label="Same-charge (+/+, -/-)", color=CHARGE_COLORS["repulse"], edgecolor="black", linewidth=0.5)
    ax1.bar(x, neutral, width, bottom=attract + repulse, label="Neutral", color=CHARGE_COLORS["neutral"], edgecolor="black", linewidth=0.5)
    
    for i in range(len(VARIANTS)):
        ax1.text(x[i], totals[i] + max(totals)*0.03, f"{totals[i]:.1f}", ha="center", va="bottom", fontsize=7.5, fontweight="bold")
        
    ax1.set_xticks(x)
    ax1.set_xticklabels(VARIANTS)
    ax1.set_ylabel("Inter contacts (pairs/chain)")
    ax1.set_title("A. Contact charge decomposition", loc="left", fontweight="bold", pad=20)
    ax1.set_ylim(0, max(totals) * 1.35)
    ax1.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=3, fontsize=6.8, frameon=False)
    
    # Panel 2: Attractive Fraction
    ax2.bar(x, fracs, width, yerr=frac_sems, capsize=3, color=[COLORS_TOPIC1[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    for i in range(len(VARIANTS)):
        ax2.text(x[i], fracs[i] + frac_sems[i] + max(fracs)*0.03, f"{fracs[i]:.1f}%", ha="center", va="bottom", fontsize=7.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(VARIANTS)
    ax2.set_ylabel("Opposite-charge fraction (%)")
    ax2.set_title("B. Attractive contact fraction", loc="left", fontweight="bold")
    ax2.set_ylim(0, max(fracs + frac_sems) * 1.35)
    
    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic1_inter_charge_decomposition", output_dir)
    return fig


def plot_partners_and_max_adhesion(object_name="dense_chains", output_dir=None):
    """Plot multivalent partner count vs single-partner max attractive contacts (adhesion tradeoff)."""
    apply_nature_style()
    df = filter_topic(load_summary(object_name), topic=1)
    
    fig, axes = plt.subplots(1, 3, figsize=(8.2, 2.7))
    
    x = np.arange(len(VARIANTS))
    width = 0.55
    
    # Panel 1: Partner Count (threshold = 1)
    ax = axes[0]
    p_cnt = df["inter.partner_count [threshold=1] (chains)"].to_numpy()
    p_sem = df["inter.partner_count [threshold=1] (chains) SEM"].to_numpy()
    bars = ax.bar(x, p_cnt, width, yerr=p_sem, capsize=3, color=[COLORS_TOPIC1[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    for bar, val, sem in zip(bars, p_cnt, p_sem):
        ax.text(bar.get_x() + bar.get_width()/2, val + sem + max(p_cnt)*0.03, f"{val:.1f}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels(VARIANTS)
    ax.set_ylabel("Partner chains (count)")
    ax.set_title("A. Multivalent partners", loc="left", fontweight="bold")
    ax.set_ylim(0, max(p_cnt + p_sem) * 1.35)

    # Panel 2: Single-partner Max Attractive Contacts
    ax = axes[1]
    max_att = df["inter.max_attract_contacts_per_partner [threshold=1] (pairs/partner)"].to_numpy()
    max_sem = df["inter.max_attract_contacts_per_partner [threshold=1] (pairs/partner) SEM"].to_numpy()
    bars = ax.bar(x, max_att, width, yerr=max_sem, capsize=3, color=[COLORS_TOPIC1[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    for bar, val, sem in zip(bars, max_att, max_sem):
        ax.text(bar.get_x() + bar.get_width()/2, val + sem + max(max_att)*0.03, f"{val:.1f}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels(VARIANTS)
    ax.set_ylabel("Max attract contacts / partner")
    ax.set_title("B. Peak two-body adhesion", loc="left", fontweight="bold")
    ax.set_ylim(0, max(max_att + max_sem) * 1.35)

    # Panel 3: Trade-off Scatter (Partners vs Peak Adhesion)
    ax = axes[2]
    csat_vals = [0.949, 0.145, 1.305]  # A18, B01, B10
    for i, v in enumerate(VARIANTS):
        ax.errorbar(
            p_cnt[i],
            max_att[i],
            xerr=p_sem[i],
            yerr=max_sem[i],
            fmt="o",
            color=COLORS_TOPIC1[v],
            ecolor=COLORS_TOPIC1[v],
            markersize=7,
            capsize=3,
        )
        # Avoid collision with points
        dx = 7 if v != "B10" else -28
        dy = 2 if v != "B10" else -14
        ax.annotate(
            f"{v}\n({csat_vals[i]:.2f} mM)",
            (p_cnt[i], max_att[i]),
            textcoords="offset points",
            xytext=(dx, dy),
            fontsize=7.5,
            fontweight="bold",
            color=COLORS_TOPIC1[v],
        )
    ax.set_xlabel("Partner count (chains)")
    ax.set_ylabel("Max attract contacts / partner")
    ax.set_title("C. Multivalency vs Adhesion", loc="left", fontweight="bold")
    ax.set_xlim(min(p_cnt) - 1.5, max(p_cnt) + 1.8)
    ax.set_ylim(0, max(max_att + max_sem) * 1.35)

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic1_inter_partners_and_max_adhesion", output_dir)
    return fig


def plot_partner_distribution_pmf(object_name="dense_chains", output_dir=None):
    """Plot full partner count probability mass function (PMF) distribution curves with block SD shading and clean top legend."""
    apply_nature_style()
    thresholds = [1, 5, 10]
    fig, axes = plt.subplots(1, len(thresholds), figsize=(8.8, 2.9), sharey=True)
    
    for idx, thresh in enumerate(thresholds):
        ax = axes[idx]
        for var in VARIANTS:
            df_dist = load_partner_distribution(SIM_MAP[var], object_name)
            sub = df_dist[df_dist["threshold"] == thresh].sort_values("partner_count")
            if sub.empty:
                continue
            x = sub["partner_count"].to_numpy()
            y = sub["mean"].to_numpy()
            sd = sub["fraction_std_blocks"].to_numpy()
            
            ax.plot(x, y, label=var, color=COLORS_TOPIC1[var], linewidth=1.2)
            ax.fill_between(x, np.maximum(0, y - sd), np.minimum(1, y + sd), color=COLORS_TOPIC1[var], alpha=0.18)
            
        ax.set_xlim(0, 24)
        ax.set_xticks(np.arange(0, 25, 5))
        ax.set_xlabel("Partners per chain")
        ax.set_title(f"Threshold ≥ {thresh} contacts", fontsize=8.0, fontweight="bold")
        if idx == 0:
            ax.set_ylabel("Frame probability (PMF)")
            
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=len(VARIANTS), frameon=False, fontsize=8.0)
    plt.tight_layout(rect=[0, 0, 1, 0.90])
    
    if output_dir:
        save_figure(fig, "topic1_inter_partner_distribution_pmf", output_dir)
    return fig


def plot_network_percolation(object_name="dense_chains", output_dir=None):
    """Plot network global efficiency and largest component across thresholds with clean top legend."""
    if object_name != "dense_chains":
        return None
    apply_nature_style()
    df = filter_topic(load_summary(object_name), topic=1)
    
    thresholds = [1, 5, 10, 20, 30, 35]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.1))
    
    for i, var in enumerate(VARIANTS):
        row = df[df["variant"] == var].iloc[0]
        lc_vals = [row[f"network.largest_component_fraction [threshold={t}] (dimensionless)"] for t in thresholds]
        lc_sems = [row[f"network.largest_component_fraction [threshold={t}] (dimensionless) SEM"] for t in thresholds]
        ax1.errorbar(thresholds, lc_vals, yerr=lc_sems, marker="o", label=var, color=COLORS_TOPIC1[var], capsize=2.5, linewidth=1.1)
        
        ge_vals = [row[f"network.global_efficiency [threshold={t}] (dimensionless)"] for t in thresholds]
        ge_sems = [row[f"network.global_efficiency [threshold={t}] (dimensionless) SEM"] for t in thresholds]
        ax2.errorbar(thresholds, ge_vals, yerr=ge_sems, marker="s", label=var, color=COLORS_TOPIC1[var], capsize=2.5, linewidth=1.1)

    ax1.set_xlim(-1, 37)
    ax1.set_xticks(thresholds)
    ax1.set_xlabel("Threshold (contacts / pair)")
    ax1.set_ylabel("Largest component fraction")
    ax1.set_title("A. Network percolation", loc="left", fontweight="bold")
    ax1.set_ylim(-0.05, 1.15)

    ax2.set_xlim(-1, 37)
    ax2.set_xticks(thresholds)
    ax2.set_xlabel("Threshold (contacts / pair)")
    ax2.set_ylabel("Global efficiency")
    ax2.set_title("B. Network topological efficiency", loc="left", fontweight="bold")
    ax2.set_ylim(-0.02, 0.65)  # Extended headroom to prevent truncation (max GE ~0.55)

    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=len(VARIANTS), frameon=False, fontsize=8.0)
    plt.tight_layout(rect=[0, 0, 1, 0.88])
    
    if output_dir:
        save_figure(fig, "topic1_network_percolation", output_dir)
    return fig


def main(object_name="dense_chains"):
    out = HERE / "figures" / object_name
    print(f"Generating Topic 1 Inter figures ({object_name})...")
    plot_inter_charge_decomposition(object_name, out)
    plot_partners_and_max_adhesion(object_name, out)
    plot_partner_distribution_pmf(object_name, out)
    plot_network_percolation(object_name, out)


if __name__ == "__main__":
    for obj in ["dense_chains", "central_chain"]:
        main(obj)
