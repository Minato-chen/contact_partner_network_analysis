"""Topic 2 (Polarization & Neutral Tail: A20, C04, B09, C22) Joint Analysis Plots.
Features:
1. Intra vs Inter contacts per residue (trade-off balance).
2. Intermolecular contact fraction (inter_fraction).
3. Multivalent state landscape: 2D joint density contour of Rg vs Inter Contacts.
Exports strictly PNG figures (300 DPI).
"""
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

HERE = Path(__file__).resolve().parent
PLOT_ROOT = HERE.parent
ROOT = PLOT_ROOT.parent
sys.path.insert(0, str(PLOT_ROOT))

from common.style import apply_nature_style, save_figure, COLORS_TOPIC2
from common.data_loader import load_summary, filter_topic, load_paired_chain_frames

VARIANTS = ["A20", "C04", "B09", "C22"]
SIM_MAP = {"A20": "A20_T260", "C04": "C04_T260", "B09": "B09_T260", "C22": "C22_T260"}


def plot_joint_contacts_balance(object_name="dense_chains", output_dir=None):
    """Plot intra vs inter contacts per residue and intermolecular fraction with ample headroom."""
    apply_nature_style()
    df = filter_topic(load_summary(object_name), topic=2)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.6))
    
    x = np.arange(len(VARIANTS))
    width = 0.35
    
    # Panel 1: Intra vs Inter per residue
    intra_res = df["joint.intra_contacts_per_residue (contacts/residue)"].to_numpy()
    intra_sem = df["joint.intra_contacts_per_residue (contacts/residue) SEM"].to_numpy()
    inter_res = df["joint.inter_contacts_per_residue (contacts/residue)"].to_numpy()
    inter_sem = df["joint.inter_contacts_per_residue (contacts/residue) SEM"].to_numpy()
    
    ax1.bar(x - width/2, intra_res, width, yerr=intra_sem, capsize=3, label="Intra / res", color="#A0AAB2", edgecolor="black", linewidth=0.6)
    ax1.bar(x + width/2, inter_res, width, yerr=inter_sem, capsize=3, label="Inter / res", color=[COLORS_TOPIC2[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    ax1.set_xticks(x)
    ax1.set_xticklabels(VARIANTS)
    ax1.set_ylabel("Contacts / residue")
    ax1.set_title("A. Intra vs Inter per residue", loc="left", fontweight="bold")
    ax1.legend(frameon=False, fontsize=7)
    max_c = max(max(intra_res + intra_sem), max(inter_res + inter_sem))
    ax1.set_ylim(0, max_c * 1.30)

    # Panel 2: Inter fraction
    inter_frac = df["joint.inter_fraction (dimensionless)"].to_numpy() * 100
    inter_frac_sem = df["joint.inter_fraction (dimensionless) SEM"].to_numpy() * 100
    bars = ax2.bar(x, inter_frac, width=0.55, yerr=inter_frac_sem, capsize=3, color=[COLORS_TOPIC2[v] for v in VARIANTS], edgecolor="black", linewidth=0.6)
    for bar, val, sem in zip(bars, inter_frac, inter_frac_sem):
        ax2.text(bar.get_x() + bar.get_width()/2, val + sem + max(inter_frac)*0.03, f"{val:.1f}%", ha="center", va="bottom", fontsize=7.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(VARIANTS)
    ax2.set_ylabel("Inter contact fraction (%)")
    ax2.set_title("B. Outward contact fraction", loc="left", fontweight="bold")
    ax2.set_ylim(0, max(inter_frac + inter_frac_sem) * 1.30)
    
    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic2_joint_contacts_balance", output_dir)
    return fig


def plot_multivalent_landscape_contours(object_name="dense_chains", output_dir=None):
    """Plot 2D KDE multivalent landscapes (Rg vs Total Inter Contacts) for Topic 2 without clipping."""
    apply_nature_style()
    fig, axes = plt.subplots(1, 4, figsize=(10.5, 2.7), sharex=True, sharey=True)
    
    for i, var in enumerate(VARIANTS):
        ax = axes[i]
        df_paired = load_paired_chain_frames(SIM_MAP[var], object_name)
        df_sub = df_paired[df_paired["threshold"] == 1]
        
        sns.kdeplot(
            data=df_sub,
            x="total_contacts",
            y="rg_nm",
            ax=ax,
            color=COLORS_TOPIC2[var],
            fill=True,
            alpha=0.45,
            levels=7,
            thresh=0.05,
        )
        sns.kdeplot(
            data=df_sub,
            x="total_contacts",
            y="rg_nm",
            ax=ax,
            color=COLORS_TOPIC2[var],
            linewidths=0.9,
            levels=7,
            thresh=0.05,
        )
        
        mean_x = df_sub["total_contacts"].mean()
        mean_y = df_sub["rg_nm"].mean()
        ax.plot(mean_x, mean_y, "x", color="black", markersize=7, markeredgewidth=1.5)
        
        ax.set_title(f"{var} ({mean_x:.1f}, {mean_y:.2f} nm)", fontsize=8.0, fontweight="bold", pad=8)
        ax.set_xlabel("Inter contacts (pairs/chain)")
        if i == 0:
            ax.set_ylabel(r"$R_g$ (nm)")
            
    fig.suptitle(f"Topic 2: Multivalent Landscape (Rg vs Inter Contacts, {object_name})", fontsize=9.0, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.92])
    if output_dir:
        save_figure(fig, "topic2_multivalent_landscape", output_dir)
    return fig


def main(object_name="dense_chains"):
    out = HERE / "figures" / object_name
    print(f"Generating Topic 2 Joint figures ({object_name})...")
    plot_joint_contacts_balance(object_name, out)
    plot_multivalent_landscape_contours(object_name, out)


if __name__ == "__main__":
    for obj in ["dense_chains", "central_chain"]:
        main(obj)
