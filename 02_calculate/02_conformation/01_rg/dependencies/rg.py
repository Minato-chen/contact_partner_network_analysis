"""Mass-weighted radius of gyration from whole-chain coordinates (nm)."""
from pathlib import Path
import mdtraj as md
import numpy as np
import pandas as pd
import result_files as files


def radius_of_gyration(xyz, masses):
    center = (xyz * masses[..., None]).sum(axis=-2) / masses.sum(axis=-1)[..., None]
    squared = ((xyz - center[..., None, :]) ** 2).sum(axis=-1)
    return np.sqrt((squared * masses).sum(axis=-1) / masses.sum(axis=-1))


def compute(simulation, centered_root, output, parameters, io, mass_file):
    source = io.centered_source(simulation, centered_root)
    topology = md.load_topology(str(source / "topology.pdb"))
    chain_lists = files.read_json(source / "chains.json")
    table = pd.read_csv(mass_file).set_index("three").MW
    masses = np.array(
        [[table.loc[topology.atom(atom).residue.name] for atom in ids] for ids in chain_lists],
        dtype=float,
    )
    masses[:, 0] += parameters["n_terminal_mass_Da"]
    masses[:, -1] += parameters["c_terminal_mass_Da"]
    if not np.isfinite(masses).all() or (masses <= 0).any():
        raise ValueError("Residue masses must be finite and positive")
    pd.DataFrame(
        [
            dict(chain_id=c, residue_index=i, mass_Da=float(m))
            for c, row in enumerate(masses)
            for i, m in enumerate(row)
        ]
    ).to_csv(output / "bead_masses.csv", index=False)
    buffer = []
    for info, xyz, box, chains, *_ in io.centered_frames(
        simulation, centered_root, parameters["chunk_frames"]
    ):
        values = radius_of_gyration(xyz[chains].astype(float), masses)
        buffer.extend(info | dict(chain_id=c, rg_nm=float(rg)) for c, rg in enumerate(values))
        if (info["frame"] + 1) % parameters["chunk_frames"] == 0:
            io.append_csv(buffer, output / "chain_rg.csv.gz")
            buffer.clear()
    if buffer:
        io.append_csv(buffer, output / "chain_rg.csv.gz")
    return dict(
        centered=io.centered_reference(simulation, centered_root),
        mass_table=pd.read_csv(mass_file).to_dict("records"),
        coordinate_unit="nm",
    )
