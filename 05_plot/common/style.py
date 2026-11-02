"""Common Nature-standard styling and configuration for 05_plot.
Strictly exports only PNG images (300 DPI) as requested.
"""
from pathlib import Path
import matplotlib.pyplot as plt

# Nature Journal Plotting Specifications
STYLE_CONFIG = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8.0,
    "axes.labelsize": 8.5,
    "axes.titlesize": 9.0,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
}

# Color palettes
COLORS_TOPIC1 = {
    "A18": "#5B6E8A",  # Baseline / Wild-type reference (slate grey)
    "B01": "#1D70B8",  # High-LLPS alternating charge (royal blue)
    "B10": "#D95F02",  # Moderate-LLPS blocky charge (vermilion orange)
}

COLORS_TOPIC2 = {
    "A20": "#2A9D8F",  # High-LLPS symmetric alternating (teal green)
    "C04": "#E76F51",  # Weak-LLPS bipolar dumbbell (coral orange)
    "B09": "#7B2CBF",  # Moderate-LLPS asymmetric cluster (purple)
    "C22": "#F4A261",  # Very-weak-LLPS neutral-tail collapse (sandy gold)
}

COLORS_ALL = {**COLORS_TOPIC1, **COLORS_TOPIC2}

# Charge decomposition colors
CHARGE_COLORS = {
    "attract": "#1D70B8",  # Attractive (+/-): Blue
    "repulse": "#D95F02",  # Repulsive (+/+, -/-): Orange
    "neutral": "#B0B0B0",  # Neutral: Gray
}

# Variant metadata with Csat
VARIANT_INFO = {
    "A18": {"label": "A18 (WT ref)", "csat": 0.949, "topic": 1},
    "B01": {"label": "B01 (Alternating)", "csat": 0.145, "topic": 1},
    "B10": {"label": "B10 (Blocky)", "csat": 1.305, "topic": 1},
    "A20": {"label": "A20 (Symmetric)", "csat": 0.125, "topic": 2},
    "C04": {"label": "C04 (Bipolar)", "csat": 4.142, "topic": 2},
    "B09": {"label": "B09 (Asymmetric)", "csat": 1.212, "topic": 2},
    "C22": {"label": "C22 (Neutral-tail)", "csat": 4.584, "topic": 2},
}


def apply_nature_style():
    """Apply Nature journal styling to global matplotlib rcParams."""
    plt.rcParams.update(STYLE_CONFIG)


def save_figure(fig, name, output_dir):
    """Save figure strictly as PNG (300 DPI). No other formats."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / f"{name}.png"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    print(f"  [saved PNG] {png_path}", flush=True)
    plt.close(fig)
