from typing import Union

def gc_count(seq: str) -> float:
    """
    Counts the GC content of a given string.

    Arguments:
    - seq: a nucleic acid string

    Returns a float number.
    """
    
    seq = seq.lower()
    return (seq.count("c") + seq.count("g"))/len(seq)

def gc_filter(seq: str, gc_bounds: tuple = (0, 100)) -> bool:
    """
    Checks whether the given string has acceptable GC count.

    Arguments:
    - seq: a nucleic acid string
    - gc_bounds: a tuple with lower and upper GC count boundaries

    Return bool.
    """

    return gc_bounds[0] <= gc_count(seq)*100 <= gc_bounds[1]

def len_filter(seq: str, len_bounds: tuple = (0, 2**32)) -> bool:
    """
    Checks whether the given string is of acceptable length

    Arguments:
    - seq: a nucleic acid string
    - len_bounds: a tuple with lower and upper length boundaries

    Returns bool
    """

    return len(seq) >= len_bounds[0] and len(seq) <= len_bounds[1]

def quality_filter(seq: str, threshold: Union[int, float] = 0) -> bool:
    """
    Checks if the mean quality of a read is acceptable

    Arguments:
    - seq: a string of phed33 quality scores per each nucleotide
    - threshold: an int or float number, lower boundary for mean quality

    Returns bool
    """

    seq_qual = sum([ord(i) for i in seq])/len(seq)
    return seq_qual >= threshold
