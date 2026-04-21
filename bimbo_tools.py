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
                f"Invalid characters in sequence: {set(self.seq) - self._alphabet}"
            )

    def __len__(self):
        return len(self.seq)

    def __str__(self):
        return self.seq

    def __getitem__(self, key):
        result = self.seq[key]
        if isinstance(key, slice):
            return type(self)(result)
        return result

    _alphabet = set()
    def check_alphabet(self) -> bool:
        if not self._alphabet:
            return True  # если алфавит не определён — не проверяем
        return set(self.seq) <= self._alphabet


class NucleicAcidSequence(BiologicalSequence):
    _complement_map = {}

    def complement(self):
        """
        Returns a complement NucleicAcidSequence object.

        Returns an object of a given class (NucleicAcid, DNA or RNA Sequence)
        """
        return type(self)(''.join(self._complement_map[base] for base in self.seq))

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

class DNASequence(NucleicAcidSequence):
    _complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                       'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    _alphabet = set('AaTtCcGg')

    def transcribe(self):
        """
        Returns:
            RNASequence object: transciribed RNA Sequence.
        """
        transcription = self.seq.replace("t", "u").replace("T", "U")
        return RNASequence(transcription)


class RNASequence(NucleicAcidSequence):
    _complement_map = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G',
                       'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}
    _alphabet = set('AaUuCcGg')

class AminoAcidSequence(BiologicalSequence):
    _alphabet = set('ACDEFGHIKLMNPQRSTVWY')

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