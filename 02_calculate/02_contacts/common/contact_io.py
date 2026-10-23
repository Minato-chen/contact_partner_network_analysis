"""File contracts shared by object analyses; no trajectory scans here."""
from pathlib import Path
import numpy as np
import pandas as pd
from common import results as files

KEYS = ["simulation_id", "source_frame", "chain_id"]


def result_path(simulation, task, results_root=files.RESULTS):
    return Path(results_root) / simulation["id"] / "02_contacts" / task


def rg_source(simulation, centered_root, io, results_root=files.RESULTS):
    path = files.completed(Path(results_root) / simulation["id"] / "01_basic/rg")
    if reference(path)["inputs"]["centered"] != io.centered_reference(
        simulation, centered_root
    ):
        raise files.InputUnavailable("Rg inputs changed; rerun 01_basic/rg.py")
    return path


def chain_selection(com, box, dense):
    eligible = np.flatnonzero(dense)
    z = com[:, 2] - box[2] / 2
    z -= box[2] * np.round(z / box[2])
    central = int(eligible[np.argmin(np.abs(z[eligible]))]) if len(eligible) else -1
    return central, z


def validate_objects(objects):
    values = ["central_chain", "dense_chains"] if objects is None else list(objects)
    if (
        not values
        or len(set(values)) != len(values)
        or set(values) - {"central_chain", "dense_chains"}
    ):
        raise ValueError("Choose distinct central_chain / dense_chains objects")
    return values


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


def object_inputs(contacts, rg=None):
    result = dict(contacts=reference(contacts))
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
