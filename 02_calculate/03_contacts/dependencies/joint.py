"""Pair intra/inter/Rg by simulation, original frame and chain identity."""
from pathlib import Path
import pandas as pd
import result_files as files


def compute(simulation, centered_root, output, parameters, object_name, io, store):
    folder = store.CONTACTS / (
        "01_central_chain" if object_name == "central_chain" else "02_dense_chains"
    )
    intra = files.completed(folder / "01_intramolecular/output" / simulation["id"])
    inter = files.completed(folder / "02_intermolecular/output" / simulation["id"])
    prepared = store.prepared_source(simulation, centered_root, io)
    rg = store.rg_source(simulation, centered_root, io)
    for path in [intra, inter]:
        saved = store.reference(path)["inputs"]
        if saved["prepared"] != store.reference(prepared) or saved["object"] != object_name:
            raise files.InputUnavailable(
                "Object/contact generation differs; rerun both object analyses"
            )
    if store.reference(intra)["inputs"]["rg"] != store.reference(rg):
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
        paired.groupby(["simulation_id", "source_frame", "threshold"])[metrics].mean().reset_index()
    )
    frames = grid.merge(
        means, on=["simulation_id", "source_frame", "threshold"], how="left", validate="one_to_one"
    )
    frames["block"] = frames.frame // parameters["block_frames"]
    frames.to_csv(output / "frame_means.csv", index=False)
    io.summarize(frames, metrics, ["simulation_id", "threshold"], output)
    return dict(
        object=object_name,
        intramolecular=store.reference(intra),
        intermolecular=store.reference(inter),
    )
