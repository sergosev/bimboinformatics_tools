def complement(seq: str) -> str:
    """
    Give a complement version of a given string

    Arguments:
    - seq: a string of both UPPERCASE and lowercase letters

    Returns string.
    """
    
    if any(item in set("Uu") for item in set(seq)):
        letters = {
            "A": "U",
            "a": "u",
            "U": "A",
            "u": "a",
            "C": "G",
            "c": "g",
            "G": "C",
            "g": "c",
        }
    else:
        letters = {
            "A": "T",
            "a": "t",
            "T": "A",
            "t": "a",
            "C": "G",
            "c": "g",
            "G": "C",
            "g": "c",
        }

    comp = ""
    for i in seq:
        comp += letters[i]

    return comp