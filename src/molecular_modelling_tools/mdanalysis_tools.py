import MDAnalysis as mda
from MDAnalysis.guesser.default_guesser import DefaultGuesser


def enriched_universe_chain_ids_elements_(unvrs: mda.Universe) -> mda.Universe:
    chain_ids = [segid[0].upper() if segid else "X" for segid in unvrs.atoms.segids]
    unvrs.add_TopologyAttr("chainID", chain_ids)

    guesser = DefaultGuesser(unvrs)
    elements = [guesser.guess_atom_element(a) for a in unvrs.atoms.names]
    unvrs.add_TopologyAttr("element", elements)
    return unvrs
