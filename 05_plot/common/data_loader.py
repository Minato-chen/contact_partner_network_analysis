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
    """Load 1D residue-level mean contact participation probabilities."""
    path = RESULTS_ROOT / simulation_id / f"02_contacts/01_intra/{object_name}/residue_profile.csv.gz"
    if not path.is_file():
        raise FileNotFoundError(f"Missing {path}")
    df = pd.read_csv(path)
    means = df.groupby("residue_index")["fraction"].mean().reset_index()
    return means


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
