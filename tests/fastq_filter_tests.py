import pytest
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import fastq_filtrator

SEQUENCES = [
    ("seq1", "AGTCCGA", "IIIIIIII"),
    ("seq2", "ATATATTA", "!IABC?&!"),
    ("seq3", "GGCCAGGC", "!!!!!!!!")
]

def make_fastq(tmp_path, sequences):
    """
    Creates temporary .fastq file to test the fastq filtering tool.

        - sequences: list of (name, seq, quality_string) tuples

    Returns a path to the file
    """
    file = tmp_path / "test.fastq"
    lines = []
    for name, seq, qual in sequences:
        lines += [f"@{name}", seq, "+", qual]
    file.write_text("\n".join(lines))
    return str(file)

def test_make_fastq(tmp_path):
    f = make_fastq(tmp_path, SEQUENCES)
    with open(f) as file:
        print(file.read())

def test_file_reading(tmp_path):
    pass

def test_file_writing(tmp_path):
    pass

def test_GC_filtering(tmp_path):
    pass

def test_len_filtering(tmp_path):
    pass

def test_quality_filtering(tmp_path):
    pass

def test_single_gc_boundary(tmp_path):
    pass

def test_single_len_boundary(tmp_path):
    pass



