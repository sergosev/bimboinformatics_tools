from typing import Union
import modules.fastq_tools as ft 
import os
import sys
from abc import ABC, abstractmethod

# =============================== bimbo tools class refactoring ===================================
class BiologicalSequence(ABC):
    def __init__(self, seq: str = None)
        if not isinstance(seq, str):
            raise TypeError(f"Sequence must be string, got {type(seq).__name__}")
        if not seq:
            raise ValueError("Sequence cannot be empty")
        
        self.seq = seq

    def __len__(self):
        return len(self.seq)
    
    def __str__(self):
        return self.seq

    def __getitem__(self, key):
        return self.seq[key]
    
    @abstractmethod
    def check_alphabet(self, alphabet: set) -> bool:
        pass
    
class NucleicAcidSequence():
    pass

class DNASequence():
    pass

class RNASequence():
    pass

class AminoAcidSequence():
    pass

# =============================== bimbo tools original module scripts ===================================
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
    for nucl in seq:
        comp += letters[nucl]

    return comp


def reverse_complement(seq: str):
    """
    Gives a reverse complement of a given string

    Arguments:
    - seq: a string of both UPPERCASE and lowercase letters

    Returns string.
    """
    return complement(seq)[::-1]


# =============================== bimbo tools original main func ===================================
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
        nuc_status = is_nucleic_acid(seq)
        return procedures[command](seq) if nuc_status else nuc_status
    else:
        result = []
        for seq in sequences:
            nuc_status = is_nucleic_acid(seq)
            if nuc_status:
                result.append(procedures[command](seq))
            else:
                result.append(nuc_status)

        return result
    

# =============================== fastq filtrator original script ===================================
def filter_fastq(
        input_file: str, 
        gc_bounds: tuple[Union[int, float], Union[int, float]] = (0, 100), 
        length_bounds: tuple[int] = (0, 2**21), 
        quality_threshold: Union[int, float] = 0,
        save_result: bool = True,
        output_file: str = "output_fastq.fastq"
) -> dict:
    """
    Filters a fastq file with nucleic acid sequences.

    Arguments:
    - input_file: a string containing a path to input fastq file
    - gc_bounds: a tuple with GC percentage boundaries (integer or float). Default is (0, 100). Can take a single value as an upper threshold
    - length_bounds: a tuple with length boundaries (only integer) Default is (0, 2**32). length_bounds
    - quality_threshold: an integer or float number, lower boundary for mean quality. Default is 0.
    - save_result: a boolean value saying if the filtering result should be saved or not
    - output_file: a string containin the name of the output file

    Returns a new dictionary containing sequences that correspond to the given filters.
    Saves the result to "filtered" directory in an output fastq file.
    For valid results check if your sequences are nucleic acids using nucleic_tools module.
    """
    
    if isinstance(gc_bounds, int) or isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)

    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)  

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
                    if save_result:
                        output_fastq.write(key)
                        output_fastq.write(seq+"\n")
                        output_fastq.write("+"+key[1:])
                        output_fastq.write(qual_score+"\n")

    print(f'Received {counter} sequences.')
    print(f'Returned {passed} sequences.')
    print(f'Filtered sequences saved to {output_path}')
    print(f'Filtered out {counter - passed} sequences.')
    return filtered_seqs