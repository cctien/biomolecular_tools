import os
import os.path as osp

from parmed.charmm.psf import CharmmPsfFile


def number_enriched_(psf_file: CharmmPsfFile) -> None:
    for atom in psf_file.atoms:
        atom.number = atom.idx + 1


def chain_enriched_(psf_file: CharmmPsfFile) -> None:
    for residue in psf_file.residues:
        residue.chain = residue.chain[0]
