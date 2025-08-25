import logging
import os
import os.path as osp
from collections.abc import Iterable, Mapping, Sequence
from functools import partial as prt

import MDAnalysis as mda
import mdtraj as mdt
from MDAnalysis.guesser.default_guesser import DefaultGuesser

from .mdanalysis_tools import enriched_universe_chain_ids_elements_

logger = logging.getLogger(__name__)
logging.getLogger("MDAnalysis").setLevel(logging.WARNING)


def extract_and_save_trajectory_frame_mdt(
    trajectory_filepath: str,
    topology_filepath: str,
    output_filestem: str,
    frame_index: int,
    validate_pdb_atom_ids: bool = False,
) -> None:
    """extract based on both the frame and the chain"""
    topology: mdt.Topology = mdt.load_topology(topology_filepath)
    frame: mdt.Trajectory = mdt.load_frame(trajectory_filepath, frame_index, top=topology)

    os.makedirs(osp.dirname(output_filestem), exist_ok=True)

    pdb_filepath = output_filestem + ".pdb"
    cif_filepath = output_filestem + ".cif"
    frame.save_pdb(pdb_filepath)
    frame.save_cif(cif_filepath)
    logger.info(f"Saved frame {frame_index} from {trajectory_filepath} to {pdb_filepath}")

    if validate_pdb_atom_ids:  # TODO
        raise NotImplementedError("PDB atom ID validation not implemented")


def extract_and_save_trajectory_frame_mda(
    trajectory_filepath: str,
    topology_filepath: str,
    output_filepath: str,
    frame_index: int,
    validate_pdb_atom_ids: bool = False,
    overwrite: bool = False,
) -> None:
    os.makedirs(osp.dirname(output_filepath), exist_ok=True)
    if osp.isfile(output_filepath) and not overwrite:
        logger.info(f"File {output_filepath} already exists and overwrite is False. Skipping.")
        return

    unvrs = mda.Universe(topology_filepath, trajectory_filepath)
    unvrs.trajectory[frame_index]  # set to the desired frame

    unvrs = enriched_universe_chain_ids_elements_(unvrs)

    all_atoms = unvrs.select_atoms("all")
    all_atoms.write(output_filepath)

    logger.info(f"Saved frame {frame_index} from {trajectory_filepath} to {output_filepath}")

    if validate_pdb_atom_ids:  # TODO
        raise NotImplementedError("PDB atom ID validation not implemented")


def main() -> None:
    import configargparse
    from ectools.logging import set_root_logger

    parser = configargparse.ArgumentParser(description="Extract frames from trajectory")
    parser.add_argument("trajectory", type=str, help="Path to the trajectory file")
    parser.add_argument("topology", type=str, help="Path to the topology file")
    parser.add_argument("output", type=str, help="Path to the output file stem")
    parser.add_argument("frame", type=int, help="Frame index to extract")
    args = parser.parse_args()
    set_root_logger(**getattr(args, "logging", {}))

    extract_and_save_trajectory_frame_mda(args.trajectory, args.topology, args.output, args.frame)


# python -m src.frame_from_traj
if __name__ == "__main__":
    main()
