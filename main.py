import sys
import os

# Setting the path to modules folder
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
modules_path = os.path.join(script_dir, "modules")
sys.path.append(modules_path)

from typing import Union
import nucleic_tools as nt # importing nucleic tools module
import fastq_tools as ft # importing fastq tools module


def run_dna_rna_tools(*seqs: str):
    """
    Performs certain procedures need for work with nucleic acids.

    Arguments:
    - seqs - a series of strings containing DNA or RNA sequences, separated by a coma.
    
    Last string of the series must be a procedure:
    - is_nucleic_acid: checks whether give strings are nucleic acids or not. Returns bool
    - reverse: reverts the given strings
    - transcribe: returns transcribed (DNA to RNA) versions of given strings
    - reverse_transcribe: returns reversely transcribed (RNA to DNA) versions of given strings
    - complement: returns complement vesions of given sctrings
    - reverse_complement: returns reversed complement versions of the given strings

    If a string contains both T and U (i.e. is not a nucleic acid) - results in False. 
    Otherwise returns a resulting string or bool.
    """

    command = seqs[-1]  #  saving the procedure name
    sequences = seqs[:-1]  #  saving the list of sequences

    #  creating a dictionary for procedures
    procedures = {
        "is_nucleic_acid": nt.is_nucleic_acid,
        "transcribe": nt.transcribe,
        "reverse_transcribe": nt.reverse_transcribe,
        "reverse": nt.reverse,
        "complement": nt.complement,
        "reverse_complement": nt.reverse_complement,
    }

    #  for one give sequence return a string
    #  more than one - a list of strings
    if len(sequences) == 1:
        seq = sequences[0]
        nuc_status = nt.is_nucleic_acid(seq)
        return procedures[command](seq) if nuc_status else nuc_status
    else:
        result = []
        for seq in sequences:
            nuc_status = nt.is_nucleic_acid(seq)
            if nuc_status:
                result.append(procedures[command](seq))
            else:
                result.append(nuc_status)

        return result
    
def filter_fastq(
        seqs: dict, 
        gc_bounds: tuple[Union[int, float]] = (0, 100), 
        length_bounds: tuple[int] = (0, 2**21), 
        quality_threshold: Union[int, float] = 0
) -> dict:
    """
    Filters a dictionary with fastq nucleic acid sequences.

    Arguments:
    - seqs: a dictionary of fastq nucleic acid sequences
    - gc_bounds: a tuple with GC percentage boundaries (integer or float). Default is (0, 100)
    - length_bounds: a tuple with length boundaries (only integer) Default is (0, 2**32)
    - quality_threshold: an integer or float number, lower boundary for mean quality. Default is 0.

    Returns a new dictionary containing sequences that correspond to the given filters.
    For valid results check if your sequences are nucleic acids using nucleic_tools module.
    """

    filtered_seq = {}
    for key in seqs:
        if (
            ft.gc_filter(seqs[key][0], gc_bounds=gc_bounds) and
            ft.len_filter(seqs[key][0], len_bounds=length_bounds) and
            ft.quality_filter(seqs[key][1], threshold=quality_threshold)
            ):
            filtered_seq[key] = seqs[key]
    
    return filtered_seq