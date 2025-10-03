def reverse_transcribe(seq: str) -> str:
    """ 
    Gives a reverse transcribed version of a given RNA sequence

    Arguments:
    - seq: a RNA string containig UPPERCASE or lowercase letters

    Return string.
    """
    trans_table = str.maketrans({"U": "T", "u": "t"})  #  making a translation table
    return seq.translate(trans_table)  #  returning transcribed sequence