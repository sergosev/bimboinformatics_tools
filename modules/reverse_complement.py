def reverse_complement(seq: str):
    """
    Gives a reverse complement of a given string
    
    Arguments:
    - seq: a string of both UPPERCASE and lowercase letters
    
    Returns string.
    """
    from complement import complement
    return complement(seq)[::-1]