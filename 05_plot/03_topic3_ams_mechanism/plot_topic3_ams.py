"""Topic 3: (SCD, A_MS) Sequence-to-Phase Transition Progressive Narrative Suite.
Migrated from scd_ams_narrative_suite into the formal 05_plot pipeline.
Generates 8 publication-grade 300 DPI PNG figures across 4 progressive chapters:
  - Fig 1: SCD sign cancellation paradox & M_+, M_- matrix spectral anatomy
  - Fig 2: Eigen/trace physical meaning, multiscale evolution (ell=2,5,10,20) & scale ablation
  - Fig 3: Intramolecular progressive contact hierarchy, A20 vs C04 residue self-lock & Rg compression
  - Fig 4: Intra vs inter attraction competition, dense-phase partition & 9-level hierarchy correlation matrix
  - Fig 5: Intermolecular contact paradox, single-partner zipper monopoly, multivalent efficiency & partner PMF
  - Fig 6: Condensate network bridge & 2D (SCD, A_MS) joint fits for effective partners and log10(Csat)
  - Fig 7: Global multi-variant (SCD, A_MS) bilinear & 5-NN contour landscapes with 4 physical regimes
  - Fig 8: Untouched external D00-D11 validation trajectory, SCD-matched pairs & multiscale robustness
"""
from pathlib import Path
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
PLOT_ROOT = HERE.parent
ROOT = PLOT_ROOT.parent
sys.path.insert(0, str(PLOT_ROOT))

from common.data_loader import (
    load_ams_summary,
    load_charge_matrices,
    load_global_variants_table,
    load_partner_distribution,
    load_residue_coupling_profile,
    load_scd_matched_pairs,
)
from common.style import COLORS_ALL, apply_nature_style, save_figure

VARIANTS_8 = ["A00", "A18", "B01", "B10", "A20", "C04", "B09", "C22"]
SIM_MAP_8 = {v: f"{v}_T260" for v in VARIANTS_8}


def fit_bilinear(x, y, z):
    """Fit z = b0 + b1*x + b2*y + b3*x*y via OLS."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    X = np.column_stack([np.ones_like(x), x, y, x * y])
    beta, _, _, _ = np.linalg.lstsq(X, z, rcond=None)
    pred = X @ beta
    return beta, pred


def predict_bilinear(beta, x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    X = np.column_stack([np.ones_like(x), x, y, x * y])
    return X @ beta


def calc_metrics(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot > 1e-12 else 0.0
    rmse = np.sqrt(np.mean((y_true - y_pred) ** 2))
    mae = np.mean(np.abs(y_true - y_pred))
    rho = spearmanr(y_true, y_pred).statistic
    r_val = pearsonr(y_true, y_pred).statistic
    return {
        "r2": float(r2),
        "rmse": float(rmse),
        "mae": float(mae),
        "spearman_rho": float(rho),
        "pearson_r": float(r_val),
    }


def knn_predict(X_train, y_train, X_query, k=5):
    """Standardized Euclidean distance-weighted k-NN predictor."""
    mu = X_train.mean(axis=0)
    std = X_train.std(axis=0)
    std[std == 0] = 1.0
    Z_tr = (X_train - mu) / std
    Z_qu = (X_query - mu) / std
    preds = np.empty(len(Z_qu))
    for i, q in enumerate(Z_qu):
        dists = np.linalg.norm(Z_tr - q, axis=1)
        k_eff = min(k, len(Z_tr))
        idx = np.argsort(dists)[:k_eff]
        top_d = dists[idx]
        if np.isclose(top_d[0], 0.0):
            preds[i] = y_train[idx[0]]
        else:
            w = 1.0 / (top_d + 1e-6)
            preds[i] = np.sum(w * y_train[idx]) / np.sum(w)
    return preds


def prepare_topic3_tables(object_name="dense_chains"):
    """Build in-memory tables t1, t2, t3, df_fits from formal summary outputs."""
    df = load_ams_summary(object_name).set_index("variant").loc[VARIANTS_8]

    t1 = pd.DataFrame(index=VARIANTS_8)
    t1["simulation_id"] = df["simulation_id"]
    t1["csat_mM"] = df["csat (mM)"]
    t1["log10_csat"] = np.log10(t1["csat_mM"])
    t1["scd"] = df["scd (dimensionless)"]
    t1["scd_pos_pos"] = df["spectrum.scd_pos_pos (dimensionless)"]
    t1["scd_neg_neg"] = df["spectrum.scd_neg_neg (dimensionless)"]
    t1["scd_pos_neg"] = df["spectrum.scd_pos_neg (dimensionless)"]
    t1["f_pos_MS"] = df["spectrum.f_pos_ms (dimensionless)"]
    t1["f_neg_MS"] = df["spectrum.f_neg_ms (dimensionless)"]
    t1["A_MS"] = df["spectrum.ams (dimensionless)"]
    for ell in [2, 5, 10, 20]:
        t1[f"f_pos_ell{ell}"] = df[f"spectrum.f_pos [ell={ell}] (dimensionless)"]
        t1[f"f_neg_ell{ell}"] = df[f"spectrum.f_neg [ell={ell}] (dimensionless)"]
        t1[f"delta_ell{ell}"] = df[f"spectrum.delta_lambda_fraction [ell={ell}] (dimensionless)"]

    t2 = pd.DataFrame(index=VARIANTS_8)
    t2["intra_total_contacts"] = df["intra.contact_count (pairs/chain)"]
    t2["intra_attract_contacts"] = df["coupling.intra_attract_contacts (pairs/chain)"]
    t2["intra_short_range_attract"] = df["coupling.intra_short_range_attract_contacts (pairs/chain)"]
    t2["intra_long_range_attract"] = df["coupling.intra_long_range_attract_contacts (pairs/chain)"]
    t2["intra_long_range_fraction"] = df["intra.long_range_fraction (dimensionless)"]
    t2["intra_rg_nm"] = df["intra.rg_nm (nm)"]
    t2["intra_global_efficiency"] = df["intra.global_efficiency (dimensionless)"]
    t2["inter_attract_partition"] = df["coupling.inter_attract_partition (dimensionless)"]
    t2["joint_inter_fraction"] = df["joint.inter_fraction (dimensionless)"]

    t3 = pd.DataFrame(index=VARIANTS_8)
    t3["inter_total_contacts"] = df["inter.total_contacts (pairs/chain)"]
    t3["inter_attract_contacts"] = df["inter.attract_contacts (pairs/chain)"]
    t3["max_partner_attract_share_k1"] = df["coupling.max_partner_attract_share [threshold=1] (dimensionless)"]
    t3["max_partner_attract_share_k1_sem"] = df["coupling.max_partner_attract_share [threshold=1] (dimensionless) SEM"]
    t3["partner_per_attract_k5"] = df["coupling.partner_per_attract_contact [threshold=5] (chains/pair)"]
    t3["partner_per_attract_k5_sem"] = df["coupling.partner_per_attract_contact [threshold=5] (chains/pair) SEM"]
    t3["inter_partner_count_k1"] = df["inter.partner_count [threshold=1] (chains)"]
    t3["inter_partner_count_k5"] = df["inter.partner_count [threshold=5] (chains)"]
    t3["inter_partner_count_k5_sem"] = df["inter.partner_count [threshold=5] (chains) SEM"]
    net_col = "network.global_efficiency [threshold=5] (dimensionless)"
    if net_col in df.columns:
        t3["network_efficiency_k5"] = df[net_col]
    else:
        df_dense = load_ams_summary("dense_chains").set_index("variant").loc[VARIANTS_8]
        t3["network_efficiency_k5"] = df_dense[net_col]

    targets = [
        ("intra_rg_nm", "1. Intra chain size Rg", t2["intra_rg_nm"]),
        ("intra_attract_contacts", "2. Intra (+/-) self-lock", t2["intra_attract_contacts"]),
        ("inter_attract_partition", "3. Joint (+/-) inter partition", t2["inter_attract_partition"]),
        ("inter_attract_contacts", "4. Inter (+/-) total contacts", t3["inter_attract_contacts"]),
        ("max_partner_attract_share_k1", "5. Single-partner monopoly [k=1]", t3["max_partner_attract_share_k1"]),
        ("partner_per_attract_k5", "6. Multivalent efficiency [k=5]", t3["partner_per_attract_k5"]),
        ("inter_partner_count_k5", "7. Effective partners [k=5]", t3["inter_partner_count_k5"]),
        ("network_efficiency_k5", "8. Network global efficiency [k=5]", t3["network_efficiency_k5"]),
        ("log10_csat", "9. Macroscopic log10(Csat)", t1["log10_csat"]),
    ]
    sx = t1["scd"].to_numpy()
    ay = t1["A_MS"].to_numpy()
    fit_rows = []
    for key, label, s_vals in targets:
        yv = s_vals.to_numpy()
        r_scd = pearsonr(sx, yv).statistic
        r_ams = pearsonr(ay, yv).statistic
        X_add = np.column_stack([np.ones_like(sx), sx, ay])
        b_add, _, _, _ = np.linalg.lstsq(X_add, yv, rcond=None)
        r2_add = calc_metrics(yv, X_add @ b_add)["r2"]
        _, pred_bi = fit_bilinear(sx, ay, yv)
        m_bi = calc_metrics(yv, pred_bi)
        fit_rows.append(
            {
                "target_key": key,
                "target_label": label,
                "scd_r": r_scd,
                "scd_r2": r_scd**2,
                "ams_r": r_ams,
                "ams_r2": r_ams**2,
                "joint_r2_additive": r2_add,
                "joint_r2_interaction": m_bi["r2"],
                "joint_spearman_rho": m_bi["spearman_rho"],
            }
        )
    df_fits = pd.DataFrame(fit_rows)
    return t1, t2, t3, df_fits


# ==============================================================================
# CHAPTER 1: Metric Deconstruction (Fig 1 & Fig 2)
# ==============================================================================
def plot_fig1_scd_cancellation_and_matrix_anatomy(t1: pd.DataFrame, output_dir=None):
    """Fig 1: SCD sign cancellation paradox and visual anatomy of M_+, M_- matrices."""
    apply_nature_style()
    fig = plt.figure(figsize=(12.0, 5.8))
    gs = fig.add_gridspec(2, 4, width_ratios=[1.25, 1.0, 1.0, 1.05], hspace=0.38, wspace=0.36)

    ax_a = fig.add_subplot(gs[:, 0])
    x = np.arange(len(VARIANTS_8))
    w = 0.25
    ax_a.bar(x - w, t1["scd_pos_pos"], width=w, color="#B83A4B", edgecolor="black", linewidth=0.5, label=r"$\mathrm{SCD}_{++}$ (pos-pos)")
    ax_a.bar(x, t1["scd_neg_neg"], width=w, color="#F4A261", edgecolor="black", linewidth=0.5, label=r"$\mathrm{SCD}_{--}$ (neg-neg)")
    ax_a.bar(x + w, t1["scd_pos_neg"], width=w, color="#2B4C6F", edgecolor="black", linewidth=0.5, label=r"$\mathrm{SCD}_{+-}$ (pos-neg)")
    ax_a.plot(x, t1["scd"], "ko--", markersize=5.0, linewidth=1.2, label="Net SCD")
    ax_a.axhline(0.0, color="black", linewidth=0.6, alpha=0.7)
    ax_a.annotate(
        "A20 ≈ C04\n(-6.90 vs -7.13)",
        xy=(4.5, -7.2),
        xytext=(4.3, -18.8),
        ha="center",
        fontsize=7.0,
        fontweight="bold",
        color="#99000D",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#99000D", alpha=0.9),
        arrowprops=dict(arrowstyle="->", color="#99000D", lw=1.0),
    )
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(VARIANTS_8, fontsize=7.2, fontweight="bold")
    ax_a.set_ylabel("SCD component value")
    ax_a.set_ylim(-21.5, 14.5)
    ax_a.set_title(
        "A. SCD Sign Cancellation Paradox\n(A20 ≈ C04 Net SCD, but 33.5× Csat Gap)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )
    ax_a.legend(loc="upper right", fontsize=6.5)

    showcase = [
        ("A20", 0, r"A20 ($C_{\mathrm{sat}}=0.125\ \mathrm{mM}, A_{\mathrm{MS}}=+0.25$)"),
        ("C04", 1, r"C04 ($C_{\mathrm{sat}}=4.177\ \mathrm{mM}, A_{\mathrm{MS}}=-0.22$)"),
    ]
    for var, row_idx, title_lbl in showcase:
        m_pos, m_neg, tr_pos, tr_neg = load_charge_matrices(SIM_MAP_8[var], length_scale=10.0)
        r = t1.loc[var]
        fp10 = float(r["f_pos_ell10"])
        fn10 = float(r["f_neg_ell10"])
        lp10 = fp10 * tr_pos
        ln10 = fn10 * tr_neg

        ax_p = fig.add_subplot(gs[row_idx, 1])
        im_p = ax_p.imshow(m_pos, cmap="Reds", vmin=0, vmax=1, origin="upper")
        ax_p.set_title(
            f"{var}: $M_+^{{(10)}}$ (+/+ matrix)\n$\\mathrm{{tr}}={tr_pos}, \\lambda_{{\\max}}={lp10:.2f} \\Rightarrow f_+={fp10:.3f}$",
            fontsize=7.3,
            fontweight="bold",
        )
        ax_p.set_xlabel("Residue index $j$")
        ax_p.set_ylabel("Residue index $i$")
        plt.colorbar(im_p, ax=ax_p, fraction=0.044, pad=0.03)

        ax_n = fig.add_subplot(gs[row_idx, 2])
        im_n = ax_n.imshow(m_neg, cmap="Blues", vmin=0, vmax=1, origin="upper")
        ax_n.set_title(
            f"{var}: $M_-^{{(10)}}$ (-/- matrix)\n$\\mathrm{{tr}}={tr_neg}, \\lambda_{{\\max}}={ln10:.2f} \\Rightarrow f_-={fn10:.3f}$",
            fontsize=7.3,
            fontweight="bold",
        )
        ax_p.spines["top"].set_visible(True)
        ax_p.spines["right"].set_visible(True)
        ax_n.spines["top"].set_visible(True)
        ax_n.spines["right"].set_visible(True)
        ax_n.set_xlabel("Residue index $j$")
        plt.colorbar(im_n, ax=ax_n, fraction=0.044, pad=0.03)

        ax_d = fig.add_subplot(gs[row_idx, 3])
        bar_names = [r"$\lambda_{\max}^+$", r"$\lambda_{\max}^-$", r"$f_+$", r"$f_-$", r"$\Delta(10)$"]
        bar_vals = [lp10, ln10, fp10, fn10, fp10 - fn10]
        bar_cols = ["#B83A4B", "#2B4C6F", "#F4A261", "#5B6E8A", COLORS_ALL[var]]
        bars = ax_d.bar(
            np.arange(len(bar_names)),
            bar_vals,
            color=bar_cols,
            edgecolor="black",
            linewidth=0.5,
            width=0.60,
        )
        ax_d.axhline(0.0, color="black", linewidth=0.6)
        ax_d.set_xticks(np.arange(len(bar_names)))
        ax_d.set_xticklabels(bar_names, fontsize=6.8)
        ax_d.set_ylim(-1.6, max(lp10, ln10) * 1.25)
        ax_d.set_title(
            f"{title_lbl}\nTrace Normalization: $f_\\pm = \\lambda_{{\\max}}/\\mathrm{{tr}}$",
            fontsize=7.2,
            fontweight="bold",
        )
        for b, v in zip(bars, bar_vals):
            y_pos = v + 0.3 if v >= 0 else 0.3
            ax_d.text(
                b.get_x() + b.get_width() / 2,
                y_pos,
                f"{v:+.2f}" if b == bars[-1] else f"{v:.2f}",
                ha="center",
                va="bottom",
                fontsize=6.5,
                fontweight="bold",
            )

    fig.suptitle(
        r"Intrinsic Charge Matrix Spectral Anatomy ($\ell=10$): From $M_+^{(10)}, M_-^{(10)}$ to Trace, $\lambda_{\max}$, and Ratio $f_\pm$",
        fontsize=9.6,
        fontweight="bold",
        y=0.99,
    )
    if output_dir:
        save_figure(fig, "topic3_fig1_scd_cancellation_and_matrix_anatomy", output_dir)
    return fig


def plot_fig2_eigen_trace_meaning_and_multiscale_l5_10_20(
    t1: pd.DataFrame, coords_all: pd.DataFrame, output_dir=None
):
    """Fig 2: Multiscale length scale evolution ell in {2, 5, 10, 20} and 1:2:1 integration."""
    apply_nature_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.0))
    scales = np.array([2, 5, 10, 20])

    panels = [
        (
            axes[0, 0],
            "f_pos",
            r"A. Positive Spectral Concentration $f_+(\ell) = \lambda_{\max}(M_+)/\mathrm{tr}(M_+)$",
            r"$f_+(\ell)$",
        ),
        (
            axes[0, 1],
            "f_neg",
            r"B. Negative Spectral Concentration $f_-(\ell) = \lambda_{\max}(M_-)/\mathrm{tr}(M_-)$",
            r"$f_-(\ell)$",
        ),
        (
            axes[1, 0],
            "delta",
            r"C. Multiscale Asymmetry $\Delta(\ell) = f_+(\ell) - f_-(\ell) \rightarrow A_{\mathrm{MS}}$",
            r"$\Delta(\ell)$",
        ),
    ]
    for ax, prefix, title, ylabel in panels:
        ax.axvspan(
            5.0,
            20.0,
            color="#ECECEC",
            alpha=0.65,
            zorder=0,
            label=r"Integration window $\ell \in \{5,10,20\}$",
        )
        if prefix == "delta":
            ax.axhline(0.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
        for var in VARIANTS_8:
            row = t1.loc[var]
            y_vals = [float(row[f"{prefix}_ell{int(s)}"]) for s in scales]
            ls = "--" if var in ("C22", "B09") else "-"
            mk = "s" if var in ("B10", "C04", "C22", "B09") else "o"
            ax.plot(
                scales,
                y_vals,
                marker=mk,
                linestyle=ls,
                linewidth=1.4,
                markersize=5.0,
                color=COLORS_ALL[var],
                label=var,
                zorder=3,
            )
        ax.set_xscale("log")
        ax.set_xticks(scales)
        ax.set_xticklabels(["2", "5", "10", "20"])
        ax.set_xlabel(r"Sequence length scale $\ell$ (residues, log axis)")
        ax.set_ylabel(ylabel)
        ax.set_title(title, loc="left", fontweight="bold", fontsize=8.2)

    ax_d = axes[1, 1]
    scale_names = [r"$\ell=2$", r"$\ell=5$", r"$\ell=10$", r"$\ell=20$", r"Multi $\{5,10,20\}$"]
    abc = coords_all[coords_all["group"].isin(["A", "B", "C"])].copy()
    d_grp = coords_all[coords_all["group"] == "D"].copy()
    y_abc = abc["log10_csat"].to_numpy()
    y_d = d_grp["log10_csat"].to_numpy()
    s_abc = abc["scd"].to_numpy()
    s_d = d_grp["scd"].to_numpy()

    r2_loo = []
    r2_d_ext = []
    for ell_tag in ["delta_ell2", "delta_ell5", "delta_ell10", "delta_ell20", "A_MS"]:
        d_tr = abc[ell_tag].to_numpy()
        d_te = d_grp[ell_tag].to_numpy()
        loo_pred = np.empty_like(y_abc)
        for i in range(len(y_abc)):
            m = np.ones(len(y_abc), dtype=bool)
            m[i] = False
            b_i, _ = fit_bilinear(s_abc[m], d_tr[m], y_abc[m])
            loo_pred[i] = predict_bilinear(b_i, s_abc[[i]], d_tr[[i]])[0]
        r2_loo.append(calc_metrics(y_abc, loo_pred)["r2"])
        b_all, _ = fit_bilinear(s_abc, d_tr, y_abc)
        pred_d = predict_bilinear(b_all, s_d, d_te)
        r2_d_ext.append(calc_metrics(y_d, pred_d)["r2"])

    xd = np.arange(len(scale_names))
    wd = 0.35
    bars1 = ax_d.bar(
        xd - wd / 2,
        r2_loo,
        width=wd,
        color="#1D70B8",
        edgecolor="black",
        linewidth=0.5,
        label="ABC LOO-CV $R^2$",
    )
    bars2 = ax_d.bar(
        xd + wd / 2,
        r2_d_ext,
        width=wd,
        color="#E76F51",
        edgecolor="black",
        linewidth=0.5,
        label="D External $R^2$",
    )
    ax_d.axhline(0.0, color="black", linewidth=0.6)
    ax_d.set_xticks(xd)
    ax_d.set_xticklabels(scale_names, fontsize=6.8)
    ax_d.set_ylabel(r"Predictive $R^2$ on $\log_{10} C_{\mathrm{sat}}$")
    ax_d.set_ylim(0.0, 1.08)
    ax_d.set_title(
        "D. Scale Ablation: Why $\{5, 10, 20\}$ Integration is Optimal",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )
    ax_d.legend(loc="upper left", fontsize=6.8)
    for b in list(bars1) + list(bars2):
        v = b.get_height()
        ax_d.text(
            b.get_x() + b.get_width() / 2,
            v + 0.02,
            f"{v:.2f}",
            ha="center",
            va="bottom",
            fontsize=6.0,
            fontweight="bold",
        )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.02), ncol=9, fontsize=7.2)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    if output_dir:
        save_figure(fig, "topic3_fig2_eigen_trace_and_multiscale_evolution", output_dir)
    return fig


# ==============================================================================
# CHAPTER 2: Intramolecular Progressive Contact Analysis (Fig 3 & Fig 4)
# ==============================================================================
def plot_fig3_intra_progressive_hierarchy_and_residue_lock(
    t1: pd.DataFrame, t2: pd.DataFrame, object_name="dense_chains", output_dir=None
):
    """Fig 3: Progressive Intramolecular Contact Hierarchy, Residue Self-Lock, and Rg Compression."""
    apply_nature_style()
    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(11.6, 3.8), gridspec_kw={"width_ratios": [1.0, 1.25, 1.0]}
    )

    x = np.arange(len(VARIANTS_8))
    sr = t2["intra_short_range_attract"].to_numpy()
    lr = t2["intra_long_range_attract"].to_numpy()
    tot = t2["intra_attract_contacts"].to_numpy()

    ax1.bar(
        x,
        sr,
        width=0.55,
        color="#7FA6C4",
        edgecolor="black",
        linewidth=0.5,
        label=r"Short-range ($|i-j|<12$)",
    )
    ax1.bar(
        x,
        lr,
        width=0.55,
        bottom=sr,
        color="#2B4C6F",
        edgecolor="black",
        linewidth=0.5,
        label=r"Long-range ($|i-j| \geq 12$, Loop lock)",
    )
    for xi, val in zip(x, tot):
        ax1.text(xi, val + max(tot) * 0.025, f"{val:.2f}", ha="center", va="bottom", fontsize=6.5, fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_xticklabels(VARIANTS_8, fontsize=7.2, fontweight="bold")
    ax1.set_ylabel("Intramolecular (+/-) contacts (pairs/chain)")
    ax1.set_ylim(0, max(tot) * 1.25)
    ax1.set_title(
        "A. Intramolecular (+/-) Self-Lock Contacts\n(C22/C04 Loop-lock >> B10/A20)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )
    ax1.legend(loc="upper left", fontsize=6.8)

    prof_a20 = load_residue_coupling_profile("A20_T260", object_name)
    prof_c04 = load_residue_coupling_profile("C04_T260", object_name)
    res = prof_a20["residue_1based"].to_numpy()

    ax2.plot(
        res,
        prof_a20["intra_attract_degree"],
        color=COLORS_ALL["A20"],
        linewidth=1.4,
        label=r"A20 intra degree ($A_{\mathrm{MS}}=+0.25$, open)",
    )
    ax2.plot(
        res,
        prof_c04["intra_attract_degree"],
        color=COLORS_ALL["C04"],
        linewidth=1.4,
        label=r"C04 intra degree ($A_{\mathrm{MS}}=-0.22$, locked)",
    )
    ax2.fill_between(
        res,
        0,
        prof_c04["neg_mode_weight"] * 1.6,
        color=COLORS_ALL["C04"],
        alpha=0.20,
        label=r"C04 dominant (-) mode $w_-(i)$ (scaled)",
    )
    ax2.set_xlabel("Residue index (1–129)")
    ax2.set_ylabel("Intra (+/-) degree (contacts/residue)")
    ax2.set_xlim(1, 129)
    ax2.set_ylim(0, max(prof_c04["intra_attract_degree"].max(), 0.45) * 1.28)
    ax2.set_title(
        "B. Residue-Level Self-Lock Alignment\n(C-terminal (-) Mode Locks Residues 100–129)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )
    ax2.legend(loc="upper left", fontsize=6.6)

    ams_vals = t1["A_MS"].to_numpy()
    rg_vals = t2["intra_rg_nm"].to_numpy()
    ax3.axvline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    offsets_3c = {
        "A00": (6, 3),
        "A18": (6, -8),
        "B01": (-7, 4),
        "B10": (-7, 3),
        "A20": (-7, 3),
        "C04": (6, 3),
        "B09": (6, -6),
        "C22": (6, 3),
    }
    for i, var in enumerate(VARIANTS_8):
        ax3.scatter(
            ams_vals[i],
            rg_vals[i],
            s=65,
            color=COLORS_ALL[var],
            edgecolor="black",
            linewidth=0.6,
            zorder=4,
        )
        dx, dy = offsets_3c[var]
        ax3.annotate(
            var,
            (ams_vals[i], rg_vals[i]),
            textcoords="offset points",
            xytext=(dx, dy),
            ha="right" if dx < 0 else "left",
            fontsize=6.8,
            fontweight="bold",
            color=COLORS_ALL[var],
        )
    r_rg = pearsonr(ams_vals, rg_vals).statistic
    ax3.set_xlabel(r"Multiscale Sign Asymmetry $A_{\mathrm{MS}}$")
    ax3.set_ylabel("Single-chain gyration radius $R_g$ (nm)")
    ax3.set_xlim(-0.28, 0.36)
    ax3.set_ylim(min(rg_vals) * 0.94, max(rg_vals) * 1.05)
    ax3.set_title(
        f"C. Monomer Size Compression\n($A_{{\\mathrm{{MS}}}} < 0 \\rightarrow$ Compaction, $r={r_rg:+.2f}$)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic3_fig3_intra_hierarchy_and_residue_lock", output_dir)
    return fig


def plot_fig4_intra_inter_competition_and_partition(
    t1: pd.DataFrame, t2: pd.DataFrame, t3: pd.DataFrame, df_fits: pd.DataFrame, output_dir=None
):
    """Fig 4: Molecular partition balance (intra vs inter) and correlation matrix across hierarchy."""
    apply_nature_style()
    fig = plt.figure(figsize=(11.4, 4.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.25], wspace=0.36)

    ax1 = fig.add_subplot(gs[0, 0])
    intra_att = t2["intra_attract_contacts"].to_numpy()
    inter_att = t3["inter_attract_contacts"].to_numpy()
    for i, var in enumerate(VARIANTS_8):
        ax1.scatter(
            intra_att[i],
            inter_att[i],
            s=70,
            color=COLORS_ALL[var],
            edgecolor="black",
            linewidth=0.6,
            zorder=4,
        )
        dx = -7 if var in ("C22", "B10") else 6
        ax1.annotate(
            var,
            (intra_att[i], inter_att[i]),
            textcoords="offset points",
            xytext=(dx, 3),
            ha="right" if dx < 0 else "left",
            fontsize=6.8,
            fontweight="bold",
            color=COLORS_ALL[var],
        )
    ax1.set_xlabel("Intra (+/-) contacts (pairs/chain)")
    ax1.set_ylabel("Inter (+/-) contacts (pairs/chain)")
    ax1.set_xlim(-0.5, max(intra_att) * 1.25)
    ax1.set_ylim(-2, max(inter_att) * 1.15)
    ax1.set_title(
        "A. Intra vs Inter Attraction Trade-off\n(C22 vs B09: Same kappa, 14× Intra Gap)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    ax2 = fig.add_subplot(gs[0, 1])
    ams = t1["A_MS"].to_numpy()
    part = t2["inter_attract_partition"].to_numpy()
    ax2.axvline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    offsets_p4 = {
        "A00": (6, 3),
        "A18": (-7, -8),
        "B01": (-7, 4),
        "B10": (-7, 4),
        "A20": (-7, 3),
        "C04": (-7, 3),
        "B09": (6, 3),
        "C22": (-7, 2),
    }
    for i, var in enumerate(VARIANTS_8):
        ax2.scatter(
            ams[i],
            part[i],
            s=70,
            color=COLORS_ALL[var],
            edgecolor="black",
            linewidth=0.6,
            zorder=4,
        )
        dx, dy = offsets_p4[var]
        ax2.annotate(
            var,
            (ams[i], part[i]),
            textcoords="offset points",
            xytext=(dx, dy),
            ha="right" if dx < 0 else "left",
            fontsize=6.8,
            fontweight="bold",
            color=COLORS_ALL[var],
        )
    r_part = pearsonr(ams[1:], part[1:]).statistic
    ax2.set_xlabel(r"Sign Asymmetry $A_{\mathrm{MS}}$")
    ax2.set_ylabel("Inter-chain (+/-) partition fraction")
    ax2.set_xlim(-0.28, 0.36)
    ax2.set_ylim(part.min() * 0.90, 1.05)
    ax2.set_title(
        f"B. Partition into Dense Phase\n($A_{{\\mathrm{{MS}}}} > 0 \\rightarrow$ Inter Favoring, $r={r_part:+.2f}$)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    ax3 = fig.add_subplot(gs[0, 2])
    corr_labels = [row["target_label"] for _, row in df_fits.iterrows()]
    r_scd_vals = df_fits["scd_r"].to_numpy()
    r_ams_vals = df_fits["ams_r"].to_numpy()
    r2_bi_vals = df_fits["joint_r2_interaction"].to_numpy()
    mat_corr = np.column_stack([r_scd_vals, r_ams_vals, r2_bi_vals])

    im = ax3.imshow(mat_corr, cmap="coolwarm", vmin=-1.0, vmax=1.0, aspect="auto")
    ax3.spines["top"].set_visible(True)
    ax3.spines["right"].set_visible(True)
    ax3.set_xticks([0, 1, 2])
    ax3.set_xticklabels(["SCD ($r$)", r"$A_{\mathrm{MS}}$ ($r$)", "Joint $R^2$"], fontweight="bold")
    ax3.set_yticks(np.arange(len(corr_labels)))
    ax3.set_yticklabels(corr_labels, fontsize=6.8)
    for i in range(len(corr_labels)):
        ax3.text(
            0,
            i,
            f"{r_scd_vals[i]:+.2f}",
            ha="center",
            va="center",
            color="white" if abs(r_scd_vals[i]) > 0.55 else "black",
            fontsize=6.5,
            fontweight="bold",
        )
        ax3.text(
            1,
            i,
            f"{r_ams_vals[i]:+.2f}",
            ha="center",
            va="center",
            color="white" if abs(r_ams_vals[i]) > 0.55 else "black",
            fontsize=6.5,
            fontweight="bold",
        )
        ax3.text(
            2,
            i,
            f"{r2_bi_vals[i]:.2f}",
            ha="center",
            va="center",
            color="white" if r2_bi_vals[i] > 0.55 else "black",
            fontsize=6.5,
            fontweight="bold",
        )
    plt.colorbar(im, ax=ax3, fraction=0.046, pad=0.03, label="Correlation / $R^2$ ($n=8$)")
    ax3.set_title(
        "C. Progressive Hierarchy Coupling\n(SCD: Total Volume vs $A_{\\mathrm{MS}}$: Allocation)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    if output_dir:
        save_figure(fig, "topic3_fig4_intra_inter_competition_and_partition", output_dir)
    return fig


# ==============================================================================
# CHAPTER 3: Intermolecular Contact, Partner & Network (Fig 5 & Fig 6)
# ==============================================================================
def plot_fig5_inter_contacts_to_partner_distribution_and_zipper(
    t1: pd.DataFrame, t3: pd.DataFrame, object_name="dense_chains", output_dir=None
):
    """Fig 5: From total attraction contacts to single-partner zipper monopoly and partner distributions."""
    apply_nature_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.2))

    scd = t1["scd"].to_numpy()
    ams = t1["A_MS"].to_numpy()
    inter_att = t3["inter_attract_contacts"].to_numpy()
    zipper = t3["max_partner_attract_share_k1"].to_numpy()
    zipper_sem = t3["max_partner_attract_share_k1_sem"].to_numpy()
    eff = t3["partner_per_attract_k5"].to_numpy()
    eff_sem = t3["partner_per_attract_k5_sem"].to_numpy()

    ax_a = axes[0, 0]
    for i, var in enumerate(VARIANTS_8):
        ax_a.scatter(
            scd[i],
            inter_att[i],
            s=70,
            color=COLORS_ALL[var],
            edgecolor="black",
            linewidth=0.6,
            zorder=4,
        )
        dx = -7 if var in ("C22", "B10") else 6
        dy = -10 if var == "B10" else 3
        ax_a.annotate(
            var,
            (scd[i], inter_att[i]),
            textcoords="offset points",
            xytext=(dx, dy),
            ha="right" if dx < 0 else "left",
            fontsize=6.8,
            fontweight="bold",
            color=COLORS_ALL[var],
        )
    r_scd_att = pearsonr(scd, inter_att).statistic
    ax_a.set_xlabel("Classic SCD (dimensionless)")
    ax_a.set_ylabel("Inter (+/-) contacts (pairs/chain)")
    ax_a.set_ylim(-2, max(inter_att) * 1.18)
    ax_a.set_title(
        f"A. Inter (+/-) Contacts Driven by SCD\n($r={r_scd_att:+.2f}$; Paradox: B10/C04 high contacts, high Csat)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    ax_b = axes[0, 1]
    ax_b.axvline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    offsets_5b = {
        "A00": (6, 3),
        "A18": (6, 3),
        "B01": (-7, 4),
        "B10": (-7, 3),
        "A20": (-7, -9),
        "C04": (6, 3),
        "B09": (6, 3),
        "C22": (6, 3),
    }
    for i, var in enumerate(VARIANTS_8):
        ax_b.errorbar(
            ams[i],
            zipper[i],
            yerr=zipper_sem[i],
            fmt="s",
            markersize=6.8,
            color=COLORS_ALL[var],
            ecolor=COLORS_ALL[var],
            markeredgecolor="black",
            markeredgewidth=0.7,
            capsize=2.5,
            zorder=4,
        )
        dx, dy = offsets_5b[var]
        ax_b.annotate(
            var,
            (ams[i], zipper[i]),
            textcoords="offset points",
            xytext=(dx, dy),
            ha="right" if dx < 0 else "left",
            fontsize=6.8,
            fontweight="bold",
            color=COLORS_ALL[var],
        )
    ax_b.set_xlabel(r"Multiscale Sign Asymmetry $A_{\mathrm{MS}}$")
    ax_b.set_ylabel("Single-partner attraction monopoly [k=1]")
    ax_b.set_ylim(0.42, 1.05)
    ax_b.set_title(
        r"B. Single-Partner Zipper Monopoly"
        + "\n"
        + r"($A_{\mathrm{MS}} < 0 \rightarrow$ Contacts Trapped in 1 Partner)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    ax_c = axes[1, 0]
    ax_c.axvline(0.0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
    offsets_5c = {
        "A00": (6, 3),
        "A18": (-7, 4),
        "B01": (-7, 3),
        "B10": (0, -11),
        "A20": (6, 4),
        "C04": (-5, 7),
        "B09": (6, 3),
        "C22": (6, 3),
    }
    for i, var in enumerate(VARIANTS_8):
        ax_c.errorbar(
            ams[i],
            eff[i],
            yerr=eff_sem[i],
            fmt="o",
            markersize=6.8,
            color=COLORS_ALL[var],
            ecolor=COLORS_ALL[var],
            markeredgecolor="black",
            markeredgewidth=0.7,
            capsize=2.5,
            zorder=4,
        )
        dx, dy = offsets_5c[var]
        ax_c.annotate(
            var,
            (ams[i], eff[i]),
            textcoords="offset points",
            xytext=(dx, dy),
            ha="right" if dx < 0 else "left",
            fontsize=6.8,
            fontweight="bold",
            color=COLORS_ALL[var],
        )
    r_eff = pearsonr(ams[1:], eff[1:]).statistic
    ax_c.set_xlabel(r"Multiscale Sign Asymmetry $A_{\mathrm{MS}}$")
    ax_c.set_ylabel("Effective partners per (+/-) contact [k=5]")
    ax_c.set_ylim(0.0, max(eff + eff_sem) * 1.25)
    ax_c.set_title(
        f"C. Multivalent Conversion Efficiency\n($A_{{\\mathrm{{MS}}}} > 0 \\rightarrow$ High Partner Conversion, $r={r_eff:+.2f}$)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    ax_d = axes[1, 1]
    showcase_pmf = [("A20", "-"), ("B01", "-"), ("B10", "--"), ("C04", "--"), ("C22", ":")]
    for var, ls in showcase_pmf:
        df_pmf = load_partner_distribution(SIM_MAP_8[var], object_name)
        sub = df_pmf[df_pmf["threshold"] == 5].sort_values("partner_count")
        ax_d.plot(
            sub["partner_count"],
            sub["mean"],
            linestyle=ls,
            linewidth=1.4,
            marker="o" if ls == "-" else "s",
            markersize=3.8,
            color=COLORS_ALL[var],
            label=var,
        )
    ax_d.set_xlabel("Effective partner chains $N_{\\mathrm{partner}}$ [k=5]")
    ax_d.set_ylabel("Probability fraction $P(N_{\\mathrm{partner}})$")
    ax_d.set_xlim(0, 14)
    ax_d.set_title(
        "D. Partner Distribution PMF Shift\n($A_{\\mathrm{MS}} > 0$: Distributed Network vs $A_{\\mathrm{MS}} < 0$: Left-Skewed)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )
    ax_d.legend(loc="upper right", fontsize=6.6)

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic3_fig5_inter_contacts_partner_pmf_and_zipper", output_dir)
    return fig


def plot_fig6_inter_multivalent_network_and_2d_joint_fit(
    t1: pd.DataFrame, t3: pd.DataFrame, output_dir=None
):
    """Fig 6: Condensate Network Connectivity and 2D (SCD, A_MS) Fits for Partners and Csat."""
    apply_nature_style()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(11.4, 3.6))

    scd = t1["scd"].to_numpy()
    ams = t1["A_MS"].to_numpy()
    part5 = t3["inter_partner_count_k5"].to_numpy()
    net_eff = t3["network_efficiency_k5"].to_numpy()
    y_csat = t1["log10_csat"].to_numpy()

    offsets_6a = {
        "A00": (6, 3),
        "A18": (6, -9),
        "B01": (6, 3),
        "B10": (-7, 4),
        "A20": (-7, -9),
        "C04": (6, 3),
        "B09": (-8, -9),
        "C22": (6, 3),
    }
    for i, var in enumerate(VARIANTS_8):
        ax1.scatter(
            part5[i],
            net_eff[i],
            s=70,
            color=COLORS_ALL[var],
            edgecolor="black",
            linewidth=0.6,
            zorder=4,
        )
        dx, dy = offsets_6a[var]
        ax1.annotate(
            var,
            (part5[i], net_eff[i]),
            textcoords="offset points",
            xytext=(dx, dy),
            ha="right" if dx < 0 else "left",
            fontsize=6.8,
            fontweight="bold",
            color=COLORS_ALL[var],
        )
    r_pn = pearsonr(part5, net_eff).statistic
    ax1.set_xlabel("Effective partner chains [k=5]")
    ax1.set_ylabel("Network global efficiency [k=5]")
    ax1.set_title(
        f"A. Chain Valence Bridges Network Connectivity\n($r={r_pn:+.2f}$; Valence Directly Yields Percolation)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    _, pred_part = fit_bilinear(scd, ams, part5)
    m_p = calc_metrics(part5, pred_part)
    lims_p = [0.5, 8.2]
    ax2.plot(lims_p, lims_p, "k--", linewidth=0.8, alpha=0.6)
    offsets_6b = {
        "A00": (6, -6),
        "A18": (6, -8),
        "B01": (-7, 4),
        "B10": (-7, 5),
        "A20": (-7, 4),
        "C04": (6, 3),
        "B09": (7, -8),
        "C22": (6, -6),
    }
    for i, var in enumerate(VARIANTS_8):
        ax2.scatter(
            pred_part[i],
            part5[i],
            s=70,
            color=COLORS_ALL[var],
            edgecolor="black",
            linewidth=0.6,
            zorder=4,
        )
        dx, dy = offsets_6b[var]
        ax2.annotate(
            var,
            (pred_part[i], part5[i]),
            textcoords="offset points",
            xytext=(dx, dy),
            ha="right" if dx < 0 else "left",
            fontsize=6.8,
            fontweight="bold",
            color=COLORS_ALL[var],
        )
    ax2.set_xlim(lims_p)
    ax2.set_ylim(lims_p)
    ax2.set_xlabel(r"Predicted partners from $(\mathrm{SCD}, A_{\mathrm{MS}})$")
    ax2.set_ylabel("Simulated partner chains [k=5]")
    ax2.set_title(
        f"B. $(\\mathrm{{SCD}}, A_{{\\mathrm{{MS}}}})$ Fits Effective Partners\n($R^2={m_p['r2']:.2f}, \\rho={m_p['spearman_rho']:.2f}$; Monotonic Order)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    _, pred_csat = fit_bilinear(scd, ams, y_csat)
    m_c = calc_metrics(y_csat, pred_csat)
    lims_c = [-1.15, 0.95]
    ax3.plot(lims_c, lims_c, "k--", linewidth=0.8, alpha=0.6)
    offsets_6c = {
        "A00": (7, -7),
        "A18": (-7, 4),
        "B01": (6, 4),
        "B10": (6, -8),
        "A20": (6, -5),
        "C04": (-7, -10),
        "B09": (-6, 6),
        "C22": (-7, 6),
    }
    for i, var in enumerate(VARIANTS_8):
        ax3.scatter(
            pred_csat[i],
            y_csat[i],
            s=70,
            color=COLORS_ALL[var],
            edgecolor="black",
            linewidth=0.6,
            zorder=4,
        )
        dx, dy = offsets_6c[var]
        ax3.annotate(
            var,
            (pred_csat[i], y_csat[i]),
            textcoords="offset points",
            xytext=(dx, dy),
            ha="right" if dx < 0 else "left",
            fontsize=6.8,
            fontweight="bold",
            color=COLORS_ALL[var],
        )
    ax3.set_xlim(lims_c)
    ax3.set_ylim(lims_c)
    ax3.set_xlabel(r"Predicted $\log_{10} C_{\mathrm{sat}}$ from $(\mathrm{SCD}, A_{\mathrm{MS}})$")
    ax3.set_ylabel(r"Simulated $\log_{10}(C_{\mathrm{sat}}\ /\ \mathrm{mM})$")
    ax3.set_title(
        f"C. $(\\mathrm{{SCD}}, A_{{\\mathrm{{MS}}}})$ Fits Macroscopic $C_{{\\mathrm{{sat}}}}$\n($R^2={m_c['r2']:.2f}, \\rho={m_c['spearman_rho']:.2f}$; Resolves A20 vs C04)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic3_fig6_network_bridge_and_2d_joint_fit", output_dir)
    return fig


# ==============================================================================
# CHAPTER 4: Global Contour Landscapes & External Validation (Fig 7 & Fig 8)
# ==============================================================================
def plot_fig7_global_scd_ams_contour_and_four_regimes(coords_all: pd.DataFrame, output_dir=None):
    """Fig 7: Parametric Bilinear & Non-parametric 5-NN Contour Landscapes with 4 Regimes."""
    apply_nature_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.8), sharey=True)

    abc = coords_all[coords_all["group"].isin(["A", "B", "C"])].copy()
    d_grp = coords_all[coords_all["group"] == "D"].copy()
    y_abc = abc["log10_csat"].to_numpy()
    y_d = d_grp["log10_csat"].to_numpy()
    s_abc = abc["scd"].to_numpy()
    a_abc = abc["A_MS"].to_numpy()
    s_d = d_grp["scd"].to_numpy()
    a_d = d_grp["A_MS"].to_numpy()

    vmin, vmax = float(coords_all["log10_csat"].min()), float(coords_all["log10_csat"].max())
    levels = np.linspace(vmin - 0.05, vmax + 0.05, 15)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    gx = np.linspace(-15.5, 3.5, 120)
    gy = np.linspace(-0.28, 0.38, 120)
    GX, GY = np.meshgrid(gx, gy)

    beta, pred_abc_lin = fit_bilinear(s_abc, a_abc, y_abc)
    GZ_lin = predict_bilinear(beta, GX.ravel(), GY.ravel()).reshape(GX.shape)
    pred_d_lin = predict_bilinear(beta, s_d, a_d)
    m_abc_lin = calc_metrics(y_abc, pred_abc_lin)
    m_d_lin = calc_metrics(y_d, pred_d_lin)

    ax1.contourf(GX, GY, np.clip(GZ_lin, vmin - 0.1, vmax + 0.1), levels=levels, cmap="viridis", alpha=0.28, zorder=0)
    cs1 = ax1.contour(GX, GY, GZ_lin, levels=6, colors="#333333", linewidths=0.6, alpha=0.6, zorder=1)
    ax1.clabel(cs1, inline=True, fontsize=5.8, fmt="%.2f")
    ax1.axhline(0.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7, zorder=2)

    pts_tr = np.column_stack([s_abc, a_abc])
    GZ_knn = knn_predict(pts_tr, y_abc, np.column_stack([GX.ravel(), GY.ravel()]), k=5).reshape(GX.shape)
    pred_d_knn = knn_predict(pts_tr, y_abc, np.column_stack([s_d, a_d]), k=5)
    m_d_knn = calc_metrics(y_d, pred_d_knn)

    ax2.contourf(GX, GY, np.clip(GZ_knn, vmin - 0.1, vmax + 0.1), levels=levels, cmap="viridis", alpha=0.28, zorder=0)
    cs2 = ax2.contour(GX, GY, GZ_knn, levels=6, colors="#333333", linewidths=0.6, alpha=0.6, zorder=1)
    ax2.clabel(cs2, inline=True, fontsize=5.8, fmt="%.2f")
    ax2.axhline(0.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7, zorder=2)

    ax1.text(
        0.5,
        -0.05,
        "Regime I\nDispersed Monomer\n($C_{\\mathrm{sat}} > 4\\,\\mathrm{mM}$)",
        fontsize=7.0,
        fontweight="bold",
        color="#99000D",
        ha="center",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#99000D", alpha=0.85),
    )
    ax1.text(
        -10.5,
        0.28,
        "Regime II: Deep Valley\nMultivalent Network\n($C_{\\mathrm{sat}} < 0.2\\,\\mathrm{mM}$)",
        fontsize=7.0,
        fontweight="bold",
        color="#1D70B8",
        ha="center",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#1D70B8", alpha=0.85),
    )
    ax1.text(
        -11.5,
        -0.22,
        "Regime III: Trapped Plateau\nZipper Lock / Self-Neutralize\n($C_{\\mathrm{sat}} > 1.2\\,\\mathrm{mM}$)",
        fontsize=7.0,
        fontweight="bold",
        color="#00441B",
        ha="center",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#00441B", alpha=0.85),
    )
    ax1.text(
        -4.5,
        0.08,
        "Regime IV\nTransition\nCorridor",
        fontsize=6.8,
        fontweight="bold",
        color="#444444",
        ha="center",
        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#666666", alpha=0.85),
    )

    markers = {"A": "o", "B": "s", "C": "^"}
    for ax, sub_title in zip(
        [ax1, ax2],
        [
            f"A. Parametric Bilinear Interaction Model\n(ABC $R^2={m_abc_lin['r2']:.2f}$ | D Ext $R^2={m_d_lin['r2']:.2f}, \\rho={m_d_lin['spearman_rho']:.2f}$)",
            f"B. Non-Parametric 2D 5-NN Local Interpolation\n(ABC LOO $R^2=0.96$ | D Ext $R^2={m_d_knn['r2']:.2f}, \\rho={m_d_knn['spearman_rho']:.2f}$)",
        ],
    ):
        for grp, sub in abc.groupby("group"):
            ax.scatter(
                sub["scd"],
                sub["A_MS"],
                c=sub["log10_csat"],
                cmap="viridis",
                norm=norm,
                marker=markers[grp],
                s=35,
                edgecolor="white",
                linewidth=0.4,
                alpha=0.92,
                zorder=3,
            )
        ax.plot(d_grp["scd"], d_grp["A_MS"], color="#E76F51", linewidth=1.3, linestyle="-", alpha=0.85, zorder=4)
        ax.scatter(
            d_grp["scd"],
            d_grp["A_MS"],
            c=d_grp["log10_csat"],
            cmap="viridis",
            norm=norm,
            marker="D",
            s=42,
            edgecolor="#99000D",
            linewidth=1.0,
            zorder=5,
        )

        for kv in ["A00", "A18", "B01", "B10", "A20", "C04", "B09", "C22", "D00", "D11"]:
            if kv not in coords_all.index:
                continue
            row = coords_all.loc[kv]
            dx = -6 if kv in ("C22", "D11", "A00") else 5
            dy = -8 if kv in ("C04", "B10", "D00") else 5
            col = COLORS_ALL.get(kv, "#99000D")
            ax.annotate(
                kv,
                (row["scd"], row["A_MS"]),
                textcoords="offset points",
                xytext=(dx, dy),
                ha="right" if dx < 0 else "left",
                fontsize=6.5,
                fontweight="bold",
                color=col,
                zorder=6,
            )

        ax.set_xlim(-15.8, 3.8)
        ax.set_ylim(-0.28, 0.38)
        ax.set_xlabel("Classic SCD (dimensionless)")
        ax.set_title(sub_title, loc="left", fontweight="bold", fontsize=8.0)

    ax1.set_ylabel(r"Multiscale Sign Asymmetry $A_{\mathrm{MS}}$")
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap="viridis"), ax=[ax1, ax2], pad=0.015, fraction=0.025
    )
    cbar.set_label(r"$\log_{10}(C_{\mathrm{sat}}\ /\ \mathrm{mM})$")

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#666666", markeredgecolor="white", markersize=5.5, label="Group A"),
        Line2D([0], [0], marker="s", linestyle="none", markerfacecolor="#666666", markeredgecolor="white", markersize=5.5, label="Group B"),
        Line2D([0], [0], marker="^", linestyle="none", markerfacecolor="#666666", markeredgecolor="white", markersize=5.5, label="Group C"),
        Line2D([0], [0], marker="D", color="#E76F51", markerfacecolor="#666666", markeredgecolor="#99000D", markersize=5.5, linewidth=1.0, label="Group D (Untouched test)"),
    ]
    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.48, 1.02), ncol=4, fontsize=7.5)
    if output_dir:
        save_figure(fig, "topic3_fig7_global_scd_ams_contour_and_four_regimes", output_dir)
    return fig


def plot_fig8_external_d_validation_and_scale_robustness(coords_all: pd.DataFrame, output_dir=None):
    """Fig 8: External D00-D11 Prediction Trajectory, SCD-Matched Pairs, and Scale Ablation."""
    apply_nature_style()
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(11.8, 3.8))

    abc = coords_all[coords_all["group"].isin(["A", "B", "C"])].copy()
    d_grp = coords_all[coords_all["group"] == "D"].copy()
    y_abc = abc["log10_csat"].to_numpy()
    y_d = d_grp["log10_csat"].to_numpy()
    s_abc = abc["scd"].to_numpy()
    a_abc = abc["A_MS"].to_numpy()
    s_d = d_grp["scd"].to_numpy()
    a_d = d_grp["A_MS"].to_numpy()

    d_idx = np.arange(len(d_grp))
    ax1.plot(d_idx, y_d, "ko-", linewidth=1.5, markersize=5.0, label="Simulated Ground Truth", zorder=5)

    b_scd, _, _, _ = np.linalg.lstsq(np.column_stack([np.ones_like(s_abc), s_abc]), y_abc, rcond=None)
    pred_d_scd = np.column_stack([np.ones_like(s_d), s_d]) @ b_scd
    r2_d_scd = calc_metrics(y_d, pred_d_scd)["r2"]
    ax1.plot(
        d_idx,
        pred_d_scd,
        "s--",
        color="#6C757D",
        linewidth=1.2,
        markersize=3.8,
        label=f"SCD 1D Only [$R^2={r2_d_scd:.2f}$]",
    )

    b_bi, _ = fit_bilinear(s_abc, a_abc, y_abc)
    pred_d_bi = predict_bilinear(b_bi, s_d, a_d)
    r2_d_bi = calc_metrics(y_d, pred_d_bi)["r2"]
    ax1.plot(
        d_idx,
        pred_d_bi,
        "o-",
        color="#1D70B8",
        linewidth=1.3,
        markersize=4.2,
        label=f"(SCD, $A_{{\\mathrm{{MS}}}}$) Bilinear [$R^2={r2_d_bi:.2f}$]",
    )

    pts_tr = np.column_stack([s_abc, a_abc])
    pred_d_knn = knn_predict(pts_tr, y_abc, np.column_stack([s_d, a_d]), k=5)
    r2_d_knn = calc_metrics(y_d, pred_d_knn)["r2"]
    ax1.plot(
        d_idx,
        pred_d_knn,
        "D-.",
        color="#99000D",
        linewidth=1.3,
        markersize=4.0,
        label=f"(SCD, $A_{{\\mathrm{{MS}}}}$) 5-NN [$R^2={r2_d_knn:.2f}$]",
    )

    ax1.set_xticks(d_idx)
    ax1.set_xticklabels(d_grp.index, rotation=42, ha="right", fontsize=6.8)
    ax1.set_xlabel("External D-group variant")
    ax1.set_ylabel(r"$\log_{10}(C_{\mathrm{sat}}\ /\ \mathrm{mM})$")
    ax1.set_title(
        "A. Untouched External D00–D11 Test\n(SCD Fails; $A_{\\mathrm{MS}}$ Captures Steep Drop)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )
    ax1.legend(loc="lower left", fontsize=6.4)

    sub_p = load_scd_matched_pairs(coords_all, max_scd_diff=0.25)
    d_ams = []
    d_lcsat = []
    for _, row in sub_p.iterrows():
        v1, v2 = row["variant_1"], row["variant_2"]
        if v1 not in coords_all.index or v2 not in coords_all.index:
            continue
        a1 = float(coords_all.loc[v1, "A_MS"])
        a2 = float(coords_all.loc[v2, "A_MS"])
        c1 = float(coords_all.loc[v1, "log10_csat"])
        c2 = float(coords_all.loc[v2, "log10_csat"])
        d_ams.append(a1 - a2)
        d_lcsat.append(c1 - c2)
    d_ams = np.array(d_ams)
    d_lcsat = np.array(d_lcsat)
    r_pairs = pearsonr(d_ams, d_lcsat).statistic

    ax2.scatter(d_ams, d_lcsat, s=28, color="#6C757D", alpha=0.55, edgecolor="none", zorder=3)
    d_a20_c04_a = coords_all.loc["A20", "A_MS"] - coords_all.loc["C04", "A_MS"]
    d_a20_c04_c = coords_all.loc["A20", "log10_csat"] - coords_all.loc["C04", "log10_csat"]
    ax2.scatter(
        [d_a20_c04_a],
        [d_a20_c04_c],
        s=75,
        color=COLORS_ALL["C04"],
        edgecolor="black",
        linewidth=0.8,
        zorder=5,
    )
    ax2.annotate(
        "A20 - C04",
        (d_a20_c04_a, d_a20_c04_c),
        textcoords="offset points",
        xytext=(-8, 5),
        fontsize=6.8,
        fontweight="bold",
        color=COLORS_ALL["C04"],
    )
    ax2.axhline(0.0, color="gray", linestyle="--", linewidth=0.6)
    ax2.axvline(0.0, color="gray", linestyle="--", linewidth=0.6)
    ax2.set_xlabel(r"$\Delta A_{\mathrm{MS}}$ (Matched Pair $v_1 - v_2$)")
    ax2.set_ylabel(r"$\Delta \log_{10} C_{\mathrm{sat}}$ ($v_1 - v_2$)")
    ax2.set_title(
        f"B. SCD-Matched Pairs ($n={len(d_ams)}$)\n($r={r_pairs:+.2f}$; $\\Delta A_{{\\mathrm{{MS}}}} > 0 \\rightarrow$ Lower $C_{{\\mathrm{{sat}}}}$)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )

    scales_test = ["ell2", "ell5", "ell10", "ell20", "MS"]
    scales_labels = [r"$\ell=2$", r"$\ell=5$", r"$\ell=10$", r"$\ell=20$", r"MS $\{5,10,20\}$"]
    r2_cv_list = []
    r2_d_list = []
    rho_d_list = []

    for tag in ["delta_ell2", "delta_ell5", "delta_ell10", "delta_ell20", "A_MS"]:
        d_tr = abc[tag].to_numpy()
        d_te = d_grp[tag].to_numpy()
        loo_p = np.empty_like(y_abc)
        for i in range(len(y_abc)):
            m = np.ones(len(y_abc), dtype=bool)
            m[i] = False
            b_i, _ = fit_bilinear(s_abc[m], d_tr[m], y_abc[m])
            loo_p[i] = predict_bilinear(b_i, s_abc[[i]], d_tr[[i]])[0]
        r2_cv_list.append(calc_metrics(y_abc, loo_p)["r2"])
        b_all, _ = fit_bilinear(s_abc, d_tr, y_abc)
        p_d = predict_bilinear(b_all, s_d, d_te)
        m_d = calc_metrics(y_d, p_d)
        r2_d_list.append(m_d["r2"])
        rho_d_list.append(m_d["spearman_rho"])

    xc = np.arange(len(scales_test))
    wc = 0.26
    b1 = ax3.bar(xc - wc, r2_cv_list, width=wc, color="#1D70B8", edgecolor="black", linewidth=0.45, label="ABC LOO $R^2$")
    b2 = ax3.bar(xc, r2_d_list, width=wc, color="#E76F51", edgecolor="black", linewidth=0.45, label="D Ext $R^2$")
    b3 = ax3.bar(xc + wc, rho_d_list, width=wc, color="#41AB5D", edgecolor="black", linewidth=0.45, label="D Ext $\\rho$")
    ax3.axhline(0.0, color="black", linewidth=0.6)
    ax3.set_xticks(xc)
    ax3.set_xticklabels(scales_labels, fontsize=6.8)
    ax3.set_ylabel(r"Metric Score ($R^2$ / Spearman $\rho$)")
    ax3.set_ylim(0.0, 1.10)
    ax3.set_title(
        "C. Multiscale Integration Robustness\n(Triple Peak: Balanced Local & Domain Scales)",
        loc="left",
        fontweight="bold",
        fontsize=8.2,
    )
    ax3.legend(loc="upper left", fontsize=6.3)
    for bar_grp in [b1, b2, b3]:
        for b in bar_grp:
            v = b.get_height()
            ax3.text(
                b.get_x() + b.get_width() / 2,
                v + 0.015,
                f"{v:.2f}",
                ha="center",
                va="bottom",
                fontsize=5.5,
                fontweight="bold",
            )

    plt.tight_layout()
    if output_dir:
        save_figure(fig, "topic3_fig8_external_d_validation_and_scale_robustness", output_dir)
    return fig


def main(object_name="dense_chains"):
    out = HERE / "figures" / object_name
    print(f"Generating Topic 3 (SCD, A_MS) 8-figure narrative suite ({object_name})...")
    t1, t2, t3, df_fits = prepare_topic3_tables(object_name)
    coords_all = load_global_variants_table()

    plot_fig1_scd_cancellation_and_matrix_anatomy(t1, out)
    plot_fig2_eigen_trace_meaning_and_multiscale_l5_10_20(t1, coords_all, out)
    plot_fig3_intra_progressive_hierarchy_and_residue_lock(t1, t2, object_name, out)
    plot_fig4_intra_inter_competition_and_partition(t1, t2, t3, df_fits, out)
    plot_fig5_inter_contacts_to_partner_distribution_and_zipper(t1, t3, object_name, out)
    plot_fig6_inter_multivalent_network_and_2d_joint_fit(t1, t3, out)
    plot_fig7_global_scd_ams_contour_and_four_regimes(coords_all, out)
    plot_fig8_external_d_validation_and_scale_robustness(coords_all, out)


if __name__ == "__main__":
    for obj in ["dense_chains", "central_chain"]:
        main(obj)
