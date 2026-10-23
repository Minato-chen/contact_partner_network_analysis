"""File contracts shared by object analyses; no trajectory scans here."""
from pathlib import Path
import numpy as np
import pandas as pd
import project
import result_files as files

STAGE = project.ROOT / "02_calculate"
CONTACTS = STAGE / "03_contacts"
KEYS = ["simulation_id", "source_frame", "chain_id"]


def prepared_source(simulation, centered_root, io):
    path = files.completed(CONTACTS / "00_prepare/output" / simulation["id"])
    saved = files.read_json(path / "metadata.json")["inputs"]
    phase = files.read_json(io.phase_source(simulation, centered_root) / "metadata.json")
    if (
        saved["centered"] != io.centered_reference(simulation, centered_root)
        or saved["phase"] != phase
    ):
        raise files.InputUnavailable("Shared contact inputs changed; rerun contacts/00_prepare")
    return path


def rg_source(simulation, centered_root, io):
    path = files.completed(STAGE / "02_conformation/01_rg/output" / simulation["id"])
    if files.read_json(path / "metadata.json")["inputs"]["centered"] != io.centered_reference(
        simulation, centered_root
    ):
        raise files.InputUnavailable("Rg uses another centered result; rerun conformation/01_rg")
    return path


def selected_chains(saved, index, object_name):
    if object_name == "central_chain":
        chain = int(saved["central_chain"][index])
        return np.array([chain], dtype=int) if chain >= 0 else np.array([], dtype=int)
    if object_name == "dense_chains":
        return np.flatnonzero(saved["dense_chains"][index])
    raise ValueError(f"Unknown object {object_name}")


def frame_info(simulation, saved, index, block_frames):
    frame = int(saved["frame"][index])
    return dict(
        simulation_id=simulation["id"],
        frame=frame,
        source_frame=int(saved["source_frame"][index]),
        time_ps=float(saved["time_ps"][index]),
        block=frame // block_frames,
    )


def reference(path):
    return files.read_json(Path(path) / "metadata.json")


def object_inputs(prepared, rg=None):
    result = dict(prepared=reference(prepared))
    if rg is not None:
        result["rg"] = reference(rg)
    return result


def distribution_summary(frame, keys, output, io):
    output.mkdir()
    io.summarize(frame, ["fraction"], keys, output)


def read_rg(path):
    data = pd.read_csv(path / "chain_rg.csv.gz")
    if data.duplicated(KEYS).any():
        raise ValueError("Duplicate Rg keys")
    return data.set_index(["source_frame", "chain_id"]).rg_nm
