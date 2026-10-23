from importlib import import_module

"""Pair internal, external and Rg results by frame and chain."""
from pathlib import Path
import sys
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = next(p for p in (HERE, *HERE.parents) if (p / "simulations.yaml").is_file())
sys.path.insert(0, str(ROOT))
from common import results as files

PARAMETERS = dict(block_frames=50)
OBJECTS = ["central_chain", "dense_chains"]
io = import_module("02_calculate.common.analysis_io")
store = import_module("02_calculate.02_contacts.common.contact_io")


def compute(
    simulation, centered_root, output, parameters, object_name, io, store, results_root
):
    intra = files.completed(
        store.result_path(simulation, "01_intra", results_root) / object_name
    )
    inter = files.completed(
        store.result_path(simulation, "02_inter", results_root) / object_name
    )
    current_centered = io.centered_reference(simulation, centered_root)
    current_phase = store.reference(io.phase_source(simulation, centered_root, results_root))
    rg = store.rg_source(simulation, centered_root, io, results_root)
    records = [store.reference(path) for path in (intra, inter)]
    for record, task in zip(records, ("01_intra", "02_inter")):
        saved = record["inputs"]
        common = files.completed(store.result_path(simulation, task, results_root) / "common")
        if saved["object"] != object_name or saved["contacts"] != store.reference(common):
            raise files.InputUnavailable(
                "Object contact results changed; rerun this object analysis"
            )
        upstream = saved["contacts"]["inputs"]
        if upstream["centered"] != current_centered or upstream["phase"] != current_phase:
            raise files.InputUnavailable("Contact upstream changed; rerun intra and inter")
    if (
        records[0]["parameters"]["contact_cutoff_nm"]
        != records[1]["parameters"]["contact_cutoff_nm"]
    ):
        raise files.InputUnavailable(
            "Intra/inter contact distances differ; use the same cutoff and rerun"
        )
    if records[0]["inputs"]["rg"] != store.reference(rg):
        raise files.InputUnavailable("Rg changed; rerun intramolecular analysis first")
    internal = pd.read_csv(intra / "chain_contacts.csv.gz")
    external = pd.read_csv(inter / "chain_partners.csv.gz")
    if (
        internal.duplicated(store.KEYS).any()
        or external.duplicated(store.KEYS + ["threshold"]).any()
    ):
        raise ValueError("Duplicate chain/frame keys in object results")
    fields = store.KEYS + ["contact_count", "rg_nm", "n_residues"]
    paired = external.merge(
        internal[fields], on=store.KEYS, how="left", validate="many_to_one", indicator=True
    )
    if not paired["_merge"].eq("both").all():
        raise ValueError("Some intermolecular observations have no intramolecular/Rg match")
    paired = paired.drop(columns="_merge")
    paired["block"] = paired.frame // parameters["block_frames"]
    paired["intra_contacts_per_residue"] = 2 * paired.contact_count / paired.n_residues
    paired["inter_contacts_per_residue"] = paired.total_contacts / paired.n_residues
    paired.to_csv(output / "paired_chain_frames.csv.gz", index=False)
    metrics = [
        "contact_count",
        "rg_nm",
        "total_contacts",
        "partner_count",
        "intra_contacts_per_residue",
        "inter_contacts_per_residue",
    ]
    # Preserve missing-object frames using the inter frame table as the frame grid.
    grid = pd.read_csv(inter / "frame_means.csv")[
        ["simulation_id", "frame", "source_frame", "time_ps", "threshold", "n_chains"]
    ]
    means = (
        paired.groupby(["simulation_id", "source_frame", "threshold"])[metrics]
        .mean()
        .reset_index()
    )
    frames = grid.merge(
        means,
        on=["simulation_id", "source_frame", "threshold"],
        how="left",
        validate="one_to_one",
    )
    frames["block"] = frames.frame // parameters["block_frames"]
    frames.to_csv(output / "frame_means.csv", index=False)
    io.summarize(frames, metrics, ["simulation_id", "threshold"], output)
    return dict(
        object=object_name,
        intramolecular=store.reference(intra),
        intermolecular=store.reference(inter),
    )


def summary_records(stage, object_name, parameters):
    stage = Path(stage)
    saved = parameters["intramolecular"]["parameters"]
    definition = {
        key: saved[key] for key in ("contact_cutoff_nm", "exclude_sequence_separation")
    }
    row = pd.read_csv(stage / "summary.csv").iloc[0]
    return [
        files.scalar_record(
            "joint." + metric,
            row[metric + "_mean_frames"],
            "contacts/residue",
            object=object_name,
            definition=definition,
            sem=row[metric + "_sem_blocks"],
        )
        for metric in ("intra_contacts_per_residue", "inter_contacts_per_residue")
    ]


def calculate(
    simulation,
    centered_root=ROOT / "01_center/output",
    results_root=files.RESULTS,
    objects=None,
):
    chosen = store.validate_objects(OBJECTS if objects is None else objects)
    parent = store.result_path(simulation, "03_joint", results_root)
    for obj in chosen:
        target = files.reset_output(parent, obj)
        io.validate_parameters(PARAMETERS | dict(chunk_frames=1))
        inputs = compute(
            simulation, centered_root, target, PARAMETERS, obj, io, store, results_root
        )
        files.finish(
            target, simulation, PARAMETERS, inputs, scalars=summary_records(target, obj, inputs)
        )


if __name__ == "__main__":
    raise SystemExit(import_module("02_calculate.02_contacts.03_joint.run").main())
