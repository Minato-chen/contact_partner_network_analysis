"""Read completed sensitivity results; collect diagnostics without rescaling data."""
from pathlib import Path
import numpy as np
import pandas as pd
from common import results as files


def collect(output, simulations, cutoffs):
    summaries, profiles, distributions, sources, diagnostics = [], [], [], {}, []
    for cutoff in cutoffs:
        folder = f"{cutoff:.1f}".replace(".", "p") + "nm"
        for simulation in simulations:
            name = simulation["id"]
            path = files.completed(Path(output) / folder / name)
            sources[f"{folder}/{name}"] = files.read_json(path / "metadata.json")
            for filename, target in [
                ("summary.csv", summaries),
                ("residue_probability.csv", profiles),
                ("partner_distribution.csv", distributions),
            ]:
                frame = pd.read_csv(path / filename)
                if "simulation_id" not in frame:
                    frame.insert(0, "simulation_id", name)
                frame.insert(0, "cutoff_nm", cutoff)
                target.append(frame)
                if filename == "residue_probability.csv":
                    for obj, group in frame.groupby("object"):
                        valid = group["mean"].dropna()
                        diagnostics.append(
                            dict(
                                simulation_id=name,
                                object=obj,
                                cutoff_nm=cutoff,
                                mean_participation=float(valid.mean()),
                                fraction_positions_ge_0p95=float(valid.ge(0.95).mean())
                                if len(valid)
                                else np.nan,
                            )
                        )
    return (
        pd.concat(summaries, ignore_index=True),
        pd.concat(profiles, ignore_index=True),
        pd.concat(distributions, ignore_index=True),
        pd.DataFrame(diagnostics),
        sources,
    )
