"""Read one common amino-acid sequence from the centered topology."""
import mdtraj as md
from common import results as files


def read_sequence(source):
    topology = md.load_topology(str(source / "topology.pdb"))
    chains = files.read_json(source / "chains.json")
    sequences = []
    for chain in chains:
        codes = [topology.atom(index).residue.code for index in chain]
        if not codes or any(
            code is None or code not in "ACDEFGHIKLMNPQRSTVWY" for code in codes
        ):
            raise ValueError(
                "Sequence metrics require a nonempty sequence of canonical amino acids"
            )
        sequences.append("".join(codes))
    if not sequences or len(set(sequences)) != 1:
        raise ValueError("Sequence metrics require identical sequences in all chains")
    return sequences[0], len(chains)
