"""Data loader utilities for 05_plot, bridging 04_summary and 03_results."""
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PLOT_ROOT = HERE.parent
ROOT = PLOT_ROOT.parent
RESULTS_ROOT = ROOT / "03_results"
SUMMARY_ROOT = ROOT / "04_summary/output"


def load_summary(object_name="dense_chains"):
    """Load combined summary table (contacts + basic metrics)."""
    contacts_path = SUMMARY_ROOT / f"contacts_{object_name}.csv"
    basic_path = SUMMARY_ROOT / "basic.csv"
    
    if not contacts_path.is_file() or not basic_path.is_file():
        raise FileNotFoundError(f"Summary files not found in {SUMMARY_ROOT}")
        
    df_contacts = pd.read_csv(contacts_path)
    df_basic = pd.read_csv(basic_path)
    
    merged = pd.merge(
        df_contacts,
        df_basic.drop(columns=["temperature_K"], errors="ignore"),
        on=["simulation_id", "variant"],
        how="left",
    )
    return merged


def filter_topic(df, topic=1):
    """Filter DataFrame by topic (1 or 2)."""
    variants = ["A18", "B01", "B10"] if topic == 1 else ["A20", "C04", "B09", "C22"]
    filtered = df[df["variant"].isin(variants)].copy()
    filtered["variant"] = pd.Categorical(filtered["variant"], categories=variants, ordered=True)
    return filtered.sort_values("variant").reset_index(drop=True)


def load_chain_partners(simulation_id, object_name="dense_chains"):
    """Load intermolecular chain-level detail frames."""
    path = RESULTS_ROOT / simulation_id / f"02_contacts/02_inter/{object_name}/chain_partners.csv.gz"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    return pd.read_csv(path)


def load_paired_chain_frames(simulation_id, object_name="dense_chains"):
    """Load joint intra/inter matched chain frames."""
    path = RESULTS_ROOT / simulation_id / f"02_contacts/03_joint/{object_name}/paired_chain_frames.csv.gz"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    return pd.read_csv(path)


def load_contact_matrix(simulation_id, object_name="dense_chains"):
    """Load 2D intramolecular contact probability matrix (129, 129)."""
    path = RESULTS_ROOT / simulation_id / f"02_contacts/01_intra/{object_name}/contact_probability.npy"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    return np.load(path)


def load_residue_profile(simulation_id, object_name="dense_chains"):
    """Load 1D residue-level mean contact participation probabilities and block SD."""
    path = RESULTS_ROOT / simulation_id / f"02_contacts/01_intra/{object_name}/residue_profile.csv.gz"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path)
    bm = df.groupby(["block", "residue_index"])["fraction"].mean().reset_index()
    summary = bm.groupby("residue_index")["fraction"].agg(mean="mean", std="std").reset_index()
    return summary


def load_density_profile(simulation_id):
    """Load z-density profile from phase boundary calculation."""
    path = RESULTS_ROOT / simulation_id / "01_basic/phase_boundary/density_profile.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    return pd.read_csv(path)


def load_loop_size_distribution(simulation_id, object_name="dense_chains"):
    """Load 1D loop size (|i-j|) contact distribution (length 129)."""
    path = RESULTS_ROOT / simulation_id / f"02_contacts/01_intra/{object_name}/loop_size_distribution.npy"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    return np.load(path)


def load_partner_distribution(simulation_id, object_name="dense_chains"):
    """Load partner count probability distribution table."""
    path = RESULTS_ROOT / simulation_id / f"02_contacts/02_inter/{object_name}/partner_distribution.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    return pd.read_csv(path)


def load_ams_summary(object_name="dense_chains"):
    """Load merged basic + contacts + ams spectrum + ams coupling summary table."""
    base = load_summary(object_name)
    ams_spec_path = SUMMARY_ROOT / "ams.csv"
    ams_coup_path = SUMMARY_ROOT / f"ams_{object_name}.csv"
    if not ams_spec_path.is_file() or not ams_coup_path.is_file():
        raise FileNotFoundError(f"A_MS summary tables not found in {SUMMARY_ROOT}")
    df_spec = pd.read_csv(ams_spec_path).drop(columns=["temperature_K"], errors="ignore")
    df_coup = pd.read_csv(ams_coup_path).drop(columns=["temperature_K"], errors="ignore")
    merged = base.merge(df_spec, on=["simulation_id", "variant"], how="left")
    merged = merged.merge(df_coup, on=["simulation_id", "variant"], how="left")
    return merged


def load_residue_coupling_profile(simulation_id, object_name="dense_chains"):
    """Load per-residue A_MS mode weight and intramolecular degree profile."""
    path = RESULTS_ROOT / simulation_id / f"03_ams/02_coupling/{object_name}/residue_coupling_profile.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    return pd.read_csv(path)


def load_charge_matrices(simulation_id, length_scale=10.0):
    """Construct positive and negative charge kernel matrices M_+^(ell), M_-^(ell) from residue_modes.csv."""
    path = RESULTS_ROOT / simulation_id / "03_ams/01_spectrum/residue_modes.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path)
    q = df["charge"].to_numpy(dtype=int)
    pos = (q > 0).astype(float)
    neg = (q < 0).astype(float)
    positions = np.arange(len(q), dtype=float)
    kernel = np.exp(-np.abs(positions[:, None] - positions[None, :]) / float(length_scale))
    m_pos = pos[:, None] * kernel * pos[None, :]
    m_neg = neg[:, None] * kernel * neg[None, :]
    return m_pos, m_neg, int(pos.sum()), int(neg.sum())


def load_global_variants_table():
    """Load global multi-variant table (ABC + D groups).
    Automatically uses 04_summary/output/{basic,ams}.csv once server trajectories are summarized (>20 variants),
    or falls back to reference_data/global_89variants_reference.csv.
    """
    basic_path = SUMMARY_ROOT / "basic.csv"
    ams_path = SUMMARY_ROOT / "ams.csv"
    if basic_path.is_file() and ams_path.is_file():
        df_b = pd.read_csv(basic_path)
        df_a = pd.read_csv(ams_path).drop(columns=["temperature_K"], errors="ignore")
        if len(df_b) > 20 and len(df_a) > 20:
            merged = df_b.merge(df_a, on=["simulation_id", "variant"], how="inner")
            merged = merged.dropna(subset=["csat (mM)", "scd (dimensionless)", "spectrum.ams (dimensionless)"]).copy()
            merged["name"] = merged["variant"]
            merged["group"] = merged["name"].str[0]
            merged["c_dilute"] = merged["csat (mM)"]
            merged["c_dilute_err"] = merged.get("csat (mM) SEM", np.nan)
            merged["log10_csat"] = np.log10(merged["c_dilute"])
            merged["scd"] = merged["scd (dimensionless)"]
            merged["kappa"] = merged["kappa (dimensionless)"]
            merged["f_pos_MS"] = merged["spectrum.f_pos_ms (dimensionless)"]
            merged["f_neg_MS"] = merged["spectrum.f_neg_ms (dimensionless)"]
            merged["A_MS"] = merged["spectrum.ams (dimensionless)"]
            for ell in [2, 5, 10, 20]:
                merged[f"f_pos_ell{ell}"] = merged[f"spectrum.f_pos [ell={ell}] (dimensionless)"]
                merged[f"f_neg_ell{ell}"] = merged[f"spectrum.f_neg [ell={ell}] (dimensionless)"]
                merged[f"delta_ell{ell}"] = merged[f"spectrum.delta_lambda_fraction [ell={ell}] (dimensionless)"]
            return merged.set_index("name")

    ref_path = PLOT_ROOT / "03_topic3_ams_mechanism/reference_data/global_89variants_reference.csv"
    if not ref_path.is_file():
        raise FileNotFoundError(f"Missing global reference table: {ref_path}")
    return pd.read_csv(ref_path).set_index("name")


def load_scd_matched_pairs(coords_df=None, max_scd_diff=0.25):
    """Load SCD-matched variant pairs (|Delta SCD| < max_scd_diff)."""
    ref_pairs = PLOT_ROOT / "03_topic3_ams_mechanism/reference_data/matrix_similar_scd_pairs.csv"
    if ref_pairs.is_file() and (coords_df is None or len(coords_df) <= 90):
        df_p = pd.read_csv(ref_pairs)
        return df_p[df_p["absolute_scd_difference"] < max_scd_diff].copy()
    if coords_df is None:
        coords_df = load_global_variants_table()
    names = list(coords_df.index)
    scds = coords_df["scd"].to_numpy(dtype=float)
    rows = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            diff = abs(scds[i] - scds[j])
            if diff < max_scd_diff:
                rows.append({"variant_1": names[i], "variant_2": names[j], "absolute_scd_difference": diff})
    return pd.DataFrame(rows)


