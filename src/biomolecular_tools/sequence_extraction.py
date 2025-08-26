import logging
from collections.abc import Mapping

import biotite.structure.io
from biotite.sequence import ProteinSequence
from biotite.structure import AtomArray, filter_amino_acids, get_chains, get_residue_starts

logger = logging.getLogger(__name__)


def protein_sequences_ex_structure_file(filepath: str) -> Mapping[str, ProteinSequence]:
    structure = biotite.structure.io.load_structure(filepath)
    structure = structure[filter_amino_acids(structure)]
    sequences: Mapping[str, ProteinSequence] = {}
    for chain_id in get_chains(structure):
        chain = structure[structure.chain_id == chain_id]
        _sequence = get_residue_starts(chain, add_exclusive_stop=False)
        res_names = tuple(chain.res_name[_sequence].tolist())
        sequence = ProteinSequence(res_names)
        assert sequence.is_valid(), f"Invalid sequence for chain {chain_id} in {filepath}"
        sequences[chain_id.item()] = sequence

    return sequences


def extract_sequences(pdb_file: str) -> Mapping[str, str]:
    _sequences = protein_sequences_ex_structure_file(pdb_file)
    sequences = {c: "".join(s.symbols.tolist()) for c, s in _sequences.items()}
    logger.debug(f"Extracted sequences from {pdb_file}: {sequences}")
    return sequences


def main():
    import configargparse

    parser = configargparse.ArgumentParser()
    parser.add_argument("pdb_file", help="Path to the PDB file")
    args = parser.parse_args()

    sequences = extract_sequences(args.pdb_file)
    for chain_id, seq in sequences.items():
        print(f"Chain {chain_id}: {seq}")


if __name__ == "__main__":
    main()
