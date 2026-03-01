from typing import Union
import os
import sys
from abc import ABC, abstractmethod

# =============================== bimbo tools class refactoring ===================================
class BiologicalSequence(ABC):
    def __init__(self, seq: str = None):
        if not isinstance(seq, str):
            raise TypeError(f"Sequence must be string, got {type(seq).__name__}")
        if not seq:
            raise ValueError("Sequence cannot be empty")
        
        self.seq = seq

        if not self.check_alphabet():
            raise ValueError(f"Invalid characters in sequence: {set(self.seq) - set('AaTtCcUuGg')}")

    def __len__(self):
        return len(self.seq)
    
    def __str__(self):
        return self.seq

    def __getitem__(self, key):
        return self.seq[key]
    
    @abstractmethod
    def check_alphabet(self) -> bool:
        pass
    
class NucleicAcidSequence(BiologicalSequence):
    def complement(self):
        """
        Returns a complement NucleicAcidSequence object.

        Returns an object of a given class (NucleicAcid, DNA or RNA Sequence)
        """
        if any(item in set("Uu") for item in set(self.seq)):
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
        for nucl in self.seq:
            comp += letters[nucl]

        return type(self)(comp)

    def reverse(self):
        """
        Returns reversed NucleicAcidSequence object.

        Returns an object of a given class (NucleicAcid, DNA or RNA Sequence)
        """
        return type(self)(self.seq[::-1])

    def reverse_complement(self):
        """
        Return reversed complement NucleicAcidSequence object.

        Returns an object of a given class (NucleicAcid, DNA or RNA Sequence)
        """
        return self.complement().reverse()

    def check_alphabet(self):
        """
        Checks whether the given sequence is a nucleic acid.

        Returns False if the sequence does not correspond to asigned class.
        Returns bool.
        """
        if type(self) == DNASequence:
            return set(self.seq) <= set("AaTtCcGg")
        elif type(self) == RNASequence:
            return set(self.seq) <= set("AaUuCcGg")
        else:
            return set(self.seq) <= set("AaTtCcUuGg")

class DNASequence(NucleicAcidSequence):
    def transcribe(self):
        """
        Returns:
            RNASequence object: transciribed RNA Sequence.
        """
        transcription = self.seq.replace('t', 'u').replace('T', 'U')
        return RNASequence(transcription)

class RNASequence(NucleicAcidSequence):
    pass

class AminoAcidSequence(BiologicalSequence):
    def check_alphabet(self):
        """
        Check whether the given sequence is a protein/peptide.

        Returns False if the sequence does not correspond to AminoAcidSequence.
        Returns bool.
        """
        return set(self.seq) <= set('ACDEFGHIKLMNPQRSTVWY')
    
    def molecular_weight(self):
        """
        Calculate the molecular weight of the amino acid sequence
        Amino acid weights are for free amino acids; water (18 Da) 
        is subtracted per peptide bond.
    
        Returns:
            float: Molecular weight in Da.
        """
        weights = {'A': 89.09, 'R': 174.2, 'N': 132.12, 'D': 133.1,
                   'C': 121.16, 'Q': 146.15, 'E': 147.13, 'G': 75.07,
                   'H': 155.16, 'I': 131.17, 'L': 131.17, 'K': 146.19,
                   'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09,
                   'T': 119.12, 'W': 204.23, 'Y': 181.19, 'V': 117.15}
        
        return sum([weights[aa] for aa in self.seq]) - 18 * (len(self.seq) - 1)
    
    def hydrophobicity_score(self):
        """
        Calculate the proportion of hydrophobic AAs in the amino acid sequence.

        Returns:
            float: proportion of hydrophobic AAs.
        """
        hydrophobic_aas = 'AVLIPFMW'
        hydrophobic_count = 0
        for aa in self.seq:
            hydrophobic_count += aa in hydrophobic_aas
        
        return hydrophobic_count / len(self.seq)

# =============================== fastq filtrator via BioPython ===================================


# =============================== fastq filtrator module scripts ===================================
def gc_count(seq: str) -> float:
    """
    Counts the GC content of a given string.

    Arguments:
    - seq: a nucleic acid string

    Returns a float number. If the sequence length is 0 returns a "Zero lengh" string.
    """
    n = len(seq)
    seq = seq.lower()
    return (seq.count("c") + seq.count("g")) / n if n > 0 else "Zero lengh"


def gc_filter(seq: str, gc_bounds: tuple = (0, 100)) -> bool:
    """
    Checks whether the given string has acceptable GC count.

    Arguments:
    - seq: a nucleic acid string
    - gc_bounds: a tuple with lower and upper GC count boundaries

    Return bool.
    """

    return gc_bounds[0] <= gc_count(seq) * 100 <= gc_bounds[1]


def len_filter(seq: str, len_bounds: tuple = (0, 2**32)) -> bool:
    """
    Checks whether the given string is of acceptable length

    Arguments:
    - seq: a nucleic acid string
    - len_bounds: a tuple with lower and upper length boundaries

    Returns bool
    """

    return len_bounds[0] <= len(seq) <= len_bounds[1]


def quality_filter(seq: str, threshold: Union[int, float] = 0) -> bool:
    """
    Checks if the mean quality of a read is acceptable

    Arguments:
    - seq: a string of phed33 quality scores per each nucleotide
    - threshold: an int or float number, lower boundary for mean quality

    Returns bool
    """

    seq_qual = sum([ord(i)-33 for i in seq]) / len(seq)
    return seq_qual >= threshold

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

                
                if (gc_filter(seq=seq, gc_bounds=gc_bounds) and
                    len_filter(seq=seq, len_bounds=length_bounds) and
                    quality_filter(seq=qual_score,
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