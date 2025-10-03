def is_nucleic_acid(seq: str, alphabet: set = set("AaTtUuCcGg")) -> bool:
    """
    Checks whether the given sequence is a nucleic acid

    Arguments:
    - seq: a string of both UPPERCASE and lowercase letters
    - alphabet: a set of acceptable nucleotides. Default is set("AaTtUuCcGg")

    Returns bool. 
    """
    seq = set(seq)
    if any(item in set("Tt") for item in seq) and any(
        item in set("Uu") for item in seq
    ):
        return False
    else:
        return seq <= alphabet