import logging

from ectools.string import string_diff
from openmm.app import PDBFile
from pdbfixer import PDBFixer

logger = logging.getLogger(__name__)


def fixed_pdb_(in_filepath: str, out_filepath: str) -> None:
    with open(in_filepath, "r") as f:
        content_in = f.read()
    fixer = PDBFixer(filename=in_filepath)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()  # find missing atoms (this standardizes naming) but don't add them
    PDBFile.writeFile(fixer.topology, fixer.positions, out_filepath)
    with open(out_filepath, "r") as f:
        content_out = f.read()
    logger.debug(f"fixed:\n{string_diff(content_in, content_out)}")
