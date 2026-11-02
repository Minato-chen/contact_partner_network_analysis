"""Systematic exploration of Loop-Penalized Charge Descriptors in new_measure/.
Evaluates performance against all 7 variants in Topic 1 and Topic 2.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from importlib import import_module
read_sequence = import_module("02_basic_calculate._sequence").read_sequence
import yaml

# 1. Load Simulation Metadata and Sequence Data
with open(ROOT / "simulations.yaml") as f:
    cfg = yaml.safe_load(f)

# Poster Csat values (mM)
CSAT_POSTER = {
    "A18": 0.949,
    "B01": 0.145,
    "B10": 1.305,
    "A20": 0.125,
    "C04": 4.142,
    "B09": 1.212,
    "C22": 4.584,
}

variants_data = []
for s in cfg["simulations"]:
    var = s["variant"]
    src = ROOT / "01_center/output" / s["id"]
    seq, n_chains = read_sequence(src)
    q = np.array([1 if x in "KR" else (-1 if x in "DE" else 0) for x in seq], dtype=float)
    csat = CSAT_POSTER[var]
    variants_data.append({
        "id": s["id"],
        "variant": var,
        "sequence": seq,
        "charges": q,
        "csat": csat,
        "ln_csat": np.log(csat),
    })

df_vars = pd.DataFrame(variants_data)

# 2. Classic Descriptors
def compute_classic_scd(q):
    N = len(q)
    total = 0.0
    for m in range(1, N):
        for n in range(m):
            total += q[m] * q[n] * np.sqrt(m - n)
    return total / N

def compute_loop_penalized_scd(q, L0=30.0, alpha=1.5):
    """SCD with high-polymer loop closure damping."""
    N = len(q)
    total = 0.0
    for m in range(1, N):
        for n in range(m):
            k = m - n
            # Weighting kernel: suppresses large-k self-neutralizing pairs
            damping = 1.0 / (1.0 + (k / L0) ** alpha)
            w = np.sqrt(k) * damping
            total += q[m] * q[n] * w
    return total / N

def compute_exponential_scd(q, L0=35.0):
    """Exponential cutoff for effective intermolecular sticker availability."""
    N = len(q)
    total = 0.0
    for m in range(1, N):
        for n in range(m):
            k = m - n
            w = np.sqrt(k) * np.exp(-k / L0)
            total += q[m] * q[n] * w
    return total / N

def compute_bimodal_scd(q, k_min=5, k_max=45):
    """Bandpass SCD: Only counts medium-range stickers that don't self-inactivate."""
    N = len(q)
    total = 0.0
    for m in range(1, N):
        for n in range(m):
            k = m - n
            if k_min <= k <= k_max:
                total += q[m] * q[n] * np.sqrt(k)
    return total / N

print("Calculating descriptors for all variants...")
# Classic SCD
df_vars["classic_scd"] = [compute_classic_scd(q) for q in df_vars["charges"]]

# Evaluate Classic SCD
rho_classic, p_classic = spearmanr(df_vars["classic_scd"], df_vars["csat"])
r_classic, _ = pearsonr(df_vars["classic_scd"], df_vars["ln_csat"])

print(f"\n--- Classic SCD Performance ---")
print(f"Spearman rho with Csat: {rho_classic:.3f} (p={p_classic:.3f})")
print(f"Pearson r with ln(Csat): {r_classic:.3f} (R^2 = {r_classic**2:.3f})")
for _, r in df_vars.iterrows():
    print(f"  {r.variant:4s}: Csat = {r.csat:6.3f} mM, Classic SCD = {r.classic_scd:7.3f}")

# Parameter scan for Loop-Penalized SCD
best_L0, best_alpha, best_rho, best_r2 = None, None, -999, -999
results_grid = []

for L0 in np.linspace(10, 60, 26):
    for alpha in [1.0, 1.2, 1.5, 1.8, 2.0, 2.5]:
        vals = [compute_loop_penalized_scd(q, L0=L0, alpha=alpha) for q in df_vars["charges"]]
        # In SCD, more negative means stronger phase separation, so higher Csat should correlate positively with SCD!
        # If SCD is positive-correlated with Csat: as SCD increases (less negative), Csat increases.
        rho, _ = spearmanr(vals, df_vars["csat"])
        r, _ = pearsonr(vals, df_vars["ln_csat"])
        r2 = r ** 2
        results_grid.append({"L0": L0, "alpha": alpha, "rho": rho, "r": r, "r2": r2})
        if rho > best_rho or (rho == best_rho and r2 > best_r2):
            best_rho = rho
            best_r2 = r2
            best_L0 = L0
            best_alpha = alpha

df_grid = pd.DataFrame(results_grid)
print(f"\n--- Optimal Loop-Penalized SCD ---")
print(f"Optimal Parameters: L0 = {best_L0:.1f}, alpha = {best_alpha:.1f}")
print(f"Best Spearman rho: {best_rho:.3f} (Classic was {rho_classic:.3f})")
print(f"Best R^2 with ln(Csat): {best_r2:.3f} (Classic R^2 was {r_classic**2:.3f})")

# Calculate with optimal parameters
df_vars["scd_topo"] = [compute_loop_penalized_scd(q, L0=best_L0, alpha=best_alpha) for q in df_vars["charges"]]

print("\n--- Detailed Comparison Table ---")
print(df_vars[["variant", "csat", "classic_scd", "scd_topo"]])
