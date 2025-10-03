import sys
import os


# Setting the path to modules folder
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
modules_path = os.path.join(script_dir, "modules")
sys.path.append(modules_path)


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

    Returns False if a string is not a nucleic acid. 
    Otherwise returns a resulting string or bool.
    """
    
    # importing nucleic tools module
    import nucleic_tools as nt

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