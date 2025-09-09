from bidict import frozenbidict

_amino_acid_1to3_unambiguous = {
    "A": "ALA",  # Alanine
    "C": "CYS",  # Cysteine
    "D": "ASP",  # Aspartic acid
    "E": "GLU",  # Glutamic acid
    "F": "PHE",  # Phenylalanine
    "G": "GLY",  # Glycine
    "H": "HIS",  # Histidine
    "I": "ILE",  # Isoleucine
    "K": "LYS",  # Lysine
    "L": "LEU",  # Leucine
    "M": "MET",  # Methionine
    "N": "ASN",  # Asparagine
    "P": "PRO",  # Proline
    "Q": "GLN",  # Glutamine
    "R": "ARG",  # Arginine
    "S": "SER",  # Serine
    "T": "THR",  # Threonine
    "V": "VAL",  # Valine
    "W": "TRP",  # Tryptophan
    "Y": "TYR",  # Tyrosine
    "O": "PYL",  # Pyrrolysine (the 22nd proteinogenic amino acid)
    "U": "SEC",  # Selenocysteine (the 21st proteinogenic amino acid)
}
_amino_acid_1to3_ambiguous = {
    **_amino_acid_1to3_unambiguous,
    "B": "ASX",  # Aspartic acid or Asparagine (ambiguous code)
    "J": "XLE",  # Leucine or Isoleucine (ambiguous code)
    "X": "XAA",  # Any amino acid (ambiguous code)
    "Z": "GLX",  # Glutamic acid or Glutamine (ambiguous code)
}
amino_acid_code_1to3_unambiguous = frozenbidict(_amino_acid_1to3_unambiguous)
amino_acid_code_1to3 = frozenbidict(_amino_acid_1to3_ambiguous)


def amino_acid_code_1_to_3(one_letter_code: str) -> str | None:
    return amino_acid_code_1to3.get(one_letter_code.upper())


def amino_acid_code_3_to_1(three_letter_code: str) -> str | None:
    return amino_acid_code_1to3.inverse.get(three_letter_code.upper())


def _test():
    assert amino_acid_code_1_to_3("A") == "ALA"
    assert amino_acid_code_1_to_3("C") == "CYS"
    assert amino_acid_code_1_to_3("D") == "ASP"
    assert amino_acid_code_1_to_3("E") == "GLU"
    assert amino_acid_code_1_to_3("F") == "PHE"
    assert amino_acid_code_1_to_3("G") == "GLY"
    assert amino_acid_code_1_to_3("H") == "HIS"
    assert amino_acid_code_1_to_3("I") == "ILE"
    assert amino_acid_code_1_to_3("K") == "LYS"
    assert amino_acid_code_1_to_3("L") == "LEU"
    assert amino_acid_code_1_to_3("M") == "MET"
    assert amino_acid_code_1_to_3("N") == "ASN"
    assert amino_acid_code_1_to_3("P") == "PRO"
    assert amino_acid_code_1_to_3("Q") == "GLN"
    assert amino_acid_code_1_to_3("R") == "ARG"
    assert amino_acid_code_1_to_3("S") == "SER"
    assert amino_acid_code_1_to_3("T") == "THR"
    assert amino_acid_code_1_to_3("V") == "VAL"
    assert amino_acid_code_1_to_3("W") == "TRP"
    assert amino_acid_code_1_to_3("Y") == "TYR"
    # assert amino_acid_code_1_to_3("B") == "ASX"
    # assert amino_acid_code_1_to_3("Z") == "GLX"
    assert amino_acid_code_1_to_3("U") == "SEC"
    assert amino_acid_code_1_to_3("O") == "PYL"
    assert amino_acid_code_1_to_3("X") == "XAA"
    assert amino_acid_code_1_to_3("Z") is None

    assert amino_acid_code_3_to_1("Ala") == "A"
    assert amino_acid_code_3_to_1("Cys") == "C"
    assert amino_acid_code_3_to_1("Asp") == "D"
    assert amino_acid_code_3_to_1("Glu") == "E"
    assert amino_acid_code_3_to_1("Phe") == "F"
    assert amino_acid_code_3_to_1("Gly") == "G"
    assert amino_acid_code_3_to_1("His") == "H"
    assert amino_acid_code_3_to_1("Ile") == "I"
    assert amino_acid_code_3_to_1("Lys") == "K"
    assert amino_acid_code_3_to_1("Leu") == "L"
    assert amino_acid_code_3_to_1("Met") == "M"
    assert amino_acid_code_3_to_1("Asn") == "N"
    assert amino_acid_code_3_to_1("Pro") == "P"
    assert amino_acid_code_3_to_1("Gln") == "Q"
    assert amino_acid_code_3_to_1("Arg") == "R"
    assert amino_acid_code_3_to_1("Ser") == "S"
    assert amino_acid_code_3_to_1("Thr") == "T"
    assert amino_acid_code_3_to_1("Val") == "V"
    assert amino_acid_code_3_to_1("Trp") == "W"
    assert amino_acid_code_3_to_1("Tyr") == "Y"
    # assert amino_acid_code_3_to_1("Asx") == "B"
    # assert amino_acid_code_3_to_1("Glx") == "Z"
    assert amino_acid_code_3_to_1("Sec") == "U"
    assert amino_acid_code_3_to_1("Pyl") == "O"
    assert amino_acid_code_3_to_1("Foo") is None


# python -m src.biomolecular_tools.amino_acids
if __name__ == "__main__":
    _test()
