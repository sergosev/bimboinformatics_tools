from typing import Union
from Bio import SeqIO, SeqUtils, SeqRecord
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
            raise ValueError(
                f"Invalid characters in sequence: {set(self.seq) - set('AaTtCcUuGg')}"
            )

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
        transcription = self.seq.replace("t", "u").replace("T", "U")
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
        return set(self.seq) <= set("ACDEFGHIKLMNPQRSTVWY")

    def molecular_weight(self):
        """
        Calculate the molecular weight of the amino acid sequence
        Amino acid weights are for free amino acids; water (18 Da)
        is subtracted per peptide bond.

        Returns:
            float: Molecular weight in Da.
        """
        weights = {
            "A": 89.09,
            "R": 174.2,
            "N": 132.12,
            "D": 133.1,
            "C": 121.16,
            "Q": 146.15,
            "E": 147.13,
            "G": 75.07,
            "H": 155.16,
            "I": 131.17,
            "L": 131.17,
            "K": 146.19,
            "M": 149.21,
            "F": 165.19,
            "P": 115.13,
            "S": 105.09,
            "T": 119.12,
            "W": 204.23,
            "Y": 181.19,
            "V": 117.15,
        }

        return sum([weights[aa] for aa in self.seq]) - 18 * (len(self.seq) - 1)

    def hydrophobicity_score(self):
        """
        Calculate the proportion of hydrophobic AAs in the amino acid sequence.

        Returns:
            float: proportion of hydrophobic AAs.
        """
        hydrophobic_aas = "AVLIPFMW"
        hydrophobic_count = 0
        for aa in self.seq:
            hydrophobic_count += aa in hydrophobic_aas

        return hydrophobic_count / len(self.seq)


# =============================== fastq filtrator via BioPython ===================================
def filter_fastq(
    input_file: str,
    gc_bounds: tuple[Union[int, float], Union[int, float]] = (0, 100),
    len_bounds: tuple[int] = (0, 2**21),
    quality_threshold: Union[int, float] = 0,
    save_result: bool = True,
    output_file: str = "output_fastq.fastq",
) -> list:
    """
    Filters a fastq file with nucleic acid sequences.

    Arguments:
    - input_file: a string containing a path to input fastq file
    - gc_bounds: a tuple with GC percentage boundaries (integer or float). Default is (0, 100). Can take a single value as an upper threshold
    - length_bounds: a tuple with length boundaries (only integer) Default is (0, 2**32). Can take a single value as an upper threshold
    - quality_threshold: an integer or float number, lower boundary for mean quality. Default is 0.
    - save_result: a boolean value saying if the filtering result should be saved or not
    - output_file: a string containin the name of the output file

    Returns a list containing SeqRecord objects.
    Saves the result to "filtered" directory in an output fastq file.
    """
    gc_bounds = (0, gc_bounds) if isinstance(gc_bounds, (int, float)) else gc_bounds
    len_bounds = (0, len_bounds) if isinstance(len_bounds, (int, float)) else len_bounds

    filtered = []
    counter = 0
    filt_count = 0
    for record in SeqIO.parse(input_file, "fastq"):
        counter += 1
        GC = SeqUtils.gc_fraction(record) * 100
        length = len(record)
        qualities = record.letter_annotations["phred_quality"]
        mean_qual = sum(qualities) / length

        gc_flag = gc_bounds[0] < GC < gc_bounds[1]
        len_flag = len_bounds[0] < length < len_bounds[1]
        qual_flag = mean_qual > quality_threshold

        if gc_flag and len_flag and qual_flag:
            filtered.append(record)
            filt_count += 1

    print(f"Received {counter} sequences")
    print(f"Saved {filt_count} sequences")
    print(f"Deleted {counter-filt_count} sequences")

    if save_result:
        print(f"Saving result to {output_file}")
        with open(output_file, "w") as output:
            SeqIO.write(filtered, output, "fastq")
        return None

    return filtered