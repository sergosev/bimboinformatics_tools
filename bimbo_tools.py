from typing import Union
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