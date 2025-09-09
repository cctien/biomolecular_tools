from biotite.sequence import ProteinSequence

_dict_1to3 = {
    "A": "ALA",
    "B": "ASX",  # Aspartic acid or Asparagine ; ambiguous code
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "J": "XLE",  # Leucine or Isoleucine ; ambiguous code
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "O": "PYL",  # Pyrrolysine ; non-standard amino acid
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "U": "SEC",  # Selenocysteine ; non-standard amino acid
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
    "Z": "GLX",  # Glutamic acid or Glutamine ; ambiguous code
}
_dict_3to1 = {v: k for k, v in _dict_1to3.items()}


def amino_acid_code_one_to_three(one_letter_code: str) -> str | None:
    return _dict_1to3.get(one_letter_code.upper())


def amino_acid_code_three_to_one(three_letter_code: str) -> str | None:
    return _dict_3to1.get(three_letter_code.upper())


def _test():
    assert amino_acid_code_one_to_three("A") == "ALA"
    assert amino_acid_code_one_to_three("C") == "CYS"
    assert amino_acid_code_one_to_three("D") == "ASP"
    assert amino_acid_code_one_to_three("E") == "GLU"
    assert amino_acid_code_one_to_three("F") == "PHE"
    assert amino_acid_code_one_to_three("G") == "GLY"
    assert amino_acid_code_one_to_three("H") == "HIS"
    assert amino_acid_code_one_to_three("I") == "ILE"
    assert amino_acid_code_one_to_three("K") == "LYS"
    assert amino_acid_code_one_to_three("L") == "LEU"
    assert amino_acid_code_one_to_three("M") == "MET"
    assert amino_acid_code_one_to_three("N") == "ASN"
    assert amino_acid_code_one_to_three("P") == "PRO"
    assert amino_acid_code_one_to_three("Q") == "GLN"
    assert amino_acid_code_one_to_three("R") == "ARG"
    assert amino_acid_code_one_to_three("S") == "SER"
    assert amino_acid_code_one_to_three("T") == "THR"
    assert amino_acid_code_one_to_three("V") == "VAL"
    assert amino_acid_code_one_to_three("W") == "TRP"
    assert amino_acid_code_one_to_three("Y") == "TYR"
    assert amino_acid_code_one_to_three("B") == "ASX"
    assert amino_acid_code_one_to_three("Z") == "GLX"
    assert amino_acid_code_one_to_three("U") == "SEC"
    assert amino_acid_code_one_to_three("O") == "PYL"
    assert amino_acid_code_one_to_three("X") is None

    assert amino_acid_code_three_to_one("Ala") == "A"
    assert amino_acid_code_three_to_one("Cys") == "C"
    assert amino_acid_code_three_to_one("Asp") == "D"
    assert amino_acid_code_three_to_one("Glu") == "E"
    assert amino_acid_code_three_to_one("Phe") == "F"
    assert amino_acid_code_three_to_one("Gly") == "G"
    assert amino_acid_code_three_to_one("His") == "H"
    assert amino_acid_code_three_to_one("Ile") == "I"
    assert amino_acid_code_three_to_one("Lys") == "K"
    assert amino_acid_code_three_to_one("Leu") == "L"
    assert amino_acid_code_three_to_one("Met") == "M"
    assert amino_acid_code_three_to_one("Asn") == "N"
    assert amino_acid_code_three_to_one("Pro") == "P"
    assert amino_acid_code_three_to_one("Gln") == "Q"
    assert amino_acid_code_three_to_one("Arg") == "R"
    assert amino_acid_code_three_to_one("Ser") == "S"
    assert amino_acid_code_three_to_one("Thr") == "T"
    assert amino_acid_code_three_to_one("Val") == "V"
    assert amino_acid_code_three_to_one("Trp") == "W"
    assert amino_acid_code_three_to_one("Tyr") == "Y"
    assert amino_acid_code_three_to_one("Asx") == "B"
    assert amino_acid_code_three_to_one("Glx") == "Z"
    assert amino_acid_code_three_to_one("Sec") == "U"
    assert amino_acid_code_three_to_one("Pyl") == "O"
    assert amino_acid_code_three_to_one("Foo") is None


# python -m src.biomolecular_tools.amino_acids
if __name__ == "__main__":
    _test()
