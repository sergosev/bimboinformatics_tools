import sys
import os


# Setting the paths to modules folder
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
    
    # importing all needed modules
    from is_nucleic_acid import is_nucleic_acid
    from transcribe import transcribe
    from reverse_transcribe import reverse_transcribe
    from reverse import reverse
    from complement import complement
    from reverse_complement import reverse_complement

    command = seqs[-1]  #  saving the procedure name
    sequences = seqs[:-1]  #  saving the list of sequences

    #  creating a dictionary for procedures
    procedures = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse_transcribe": reverse_transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    #  for one give sequence return a string
    #  more than one - a list of strings
    if len(sequences) == 1:
        seq = sequences[0]
        return procedures[command](seq) if is_nucleic_acid(seq) else f'{is_nucleic_acid(seq)}: not a nucleic acid'
    else:
        result = []
        for seq in sequences:
            if is_nucleic_acid(seq):
                result.append(procedures[command](seq))
            else:
                result.append(f'{is_nucleic_acid(seq)}: not a nucleic acid')

        return result