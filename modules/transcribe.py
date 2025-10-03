def transcribe(seq: str) -> str:
    """
    Create a transcribed version of the given sequence

    Arguments:
    - seq: a string of both UPPERCASE and lowercase letters

    Returns string.
    """ 
    trans_table = str.maketrans({"T": "U", "t": "u"})  #  making a translation table
    return seq.translate(trans_table)  #  returning transcribed sequence