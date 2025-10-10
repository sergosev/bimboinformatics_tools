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
        input_file: str, 
        gc_bounds: tuple[Union[int, float], Union[int, float]] = (0, 100), 
        length_bounds: tuple[int] = (0, 2**21), 
        quality_threshold: Union[int, float] = 0,
        output_file: str = "output_fastq.fastq"
) -> dict:
    """
    Filters a fastq file with nucleic acid sequences.

    Arguments:
    - input_file: a string containing a path to input fastq file
    - gc_bounds: a tuple with GC percentage boundaries (integer or float). Default is (0, 100)
    - length_bounds: a tuple with length boundaries (only integer) Default is (0, 2**32)
    - quality_threshold: an integer or float number, lower boundary for mean quality. Default is 0.
    - output_file: a string containin the name of the output file

    Returns a new dictionary containing sequences that correspond to the given filters.
    Saves the result to "filtered" directory in an output fastq file.
    For valid results check if your sequences are nucleic acids using nucleic_tools module.
    """
    import os
    import sys

    #setting up directories
    work_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    if not os.path.exists(os.path.join(work_dir, 'filtered')):
        os.mkdir(os.path.join(work_dir, "filtered"))
    output_path = os.path.join(work_dir, "filtered", output_file)

    # filtering on the go
    filtered_seqs = {}
    counter = 0
    passed = 0
    with (open(input_file, mode="r") as input_fastq,
          open(output_path, mode="w") as output_fastq):
        
        for line in input_fastq:
            if line.startswith("@"):
                counter += 1
                key = line
                seq = input_fastq.readline().strip()
                next(input_fastq)
                qual_score = input_fastq.readline().strip()

                
                if (ft.gc_filter(seq=seq, gc_bounds=gc_bounds) and
                    ft.len_filter(seq=seq, len_bounds=length_bounds) and
                    ft.quality_filter(seq=qual_score,
                                      threshold=quality_threshold)):
                    passed += 1
                    filtered_seqs[key] = [seq, qual_score]
                    output_fastq.write(key)
                    output_fastq.write(seq+"\n")
                    output_fastq.write("+"+key[1:])
                    output_fastq.write(qual_score+"\n")

    print(f'Received {counter} sequences.')
    print(f'Returned {passed} sequences.')
    print(f'Filtered sequences saved to {output_path}')
    print(f'Filtered out {counter - passed} sequences.')
    return filtered_seqs