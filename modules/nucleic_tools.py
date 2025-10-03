ALPHABET = set("AaTtUuCcGg")
T_SET = {"T", "t"}
U_SET = {"U", "u"}

def is_nucleic_acid(seq: str, alphabet: set = ALPHABET) -> bool:
    """
    Checks whether the given sequence is a nucleic acid
     
    Arguments:
    - seq: a string of both UPPERCASE and lowercase letters
    - alphabet: a set of acceptable nucleotides. Default is set("AaTtUuCcGg")
     
    Returns False if the sequence contains both T and U.
    Returns bool. 
    """
    seq_set = set(seq)
    return seq_set <= alphabet and not (seq_set & T_SET and seq_set & U_SET)

    
def reverse(seq: str) -> str:
    """
    Gives a reversed version of a given string

    Arguments:
    - seq: a string of both UPPERCASE and lowercase letters

    Returns string.
    """
    return seq[::-1]
    
def transcribe(seq: str) -> str:
    """
    Creates a transcribed version of the given DNA sequence

    Arguments:
    - seq: a string of both UPPERCASE and lowercase letters

    Returns string.
    """ 
    trans_table = str.maketrans({"T": "U", "t": "u"})  #  making a translation table
    return seq.translate(trans_table)  #  returning transcribed sequence

def reverse_transcribe(seq: str) -> str:
    """ 
    Gives a reverse transcribed version of a given RNA sequence

    Arguments:
    - seq: a RNA string containig UPPERCASE or lowercase letters

    Return string.
    """
    trans_table = str.maketrans({"U": "T", "u": "t"})  #  making a translation table
    return seq.translate(trans_table)  #  returning transcribed sequence

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

def reverse_complement(seq: str):
    """
    Gives a reverse complement of a given string
    
    Arguments:
    - seq: a string of both UPPERCASE and lowercase letters
    
    Returns string.
    """
    return complement(seq)[::-1]