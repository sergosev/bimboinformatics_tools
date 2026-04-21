from Bio import SeqIO
import pytest
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import fastq_filtrator

SEQUENCES = [
    ("seq1", "AGTCCGAT", "IIIIIIII"),
    ("seq2", "ATATATTA", "!IABIII!"),
    ("seq3", "GGCCAGGC", "!III!!H!"),
    ("seq4", "ATATACTCATAG", "IIIIIIIIIIII")
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

@pytest.fixture
def fastq_file(tmp_path):
    """
    A pytest fixture of a temporary fastq_file for subsequent testing
    """
    return make_fastq(tmp_path, SEQUENCES)

class TestFileIO():
    def test_make_fastq(self, fastq_file):
        """
        Testing if make_fastq() works
        """
        with open(fastq_file) as file:
            print(file.read())

    def test_file_reading(self, fastq_file):
        """
        Testing for file reading
        """
        result = fastq_filtrator.filter_fastq(input_file=fastq_file)
        assert len(result) == len(SEQUENCES)

    def test_file_writing(self, fastq_file, tmp_path):
        """
        Testing for file writing
        """
        output_path = tmp_path / "test_output.fastq"
        fastq_filtrator.filter_fastq(input_file=fastq_file, output_file=output_path)
        assert len(SEQUENCES) == len(list(SeqIO.parse(output_path, "fastq")))

class TestFiltering():
    def test_GC_filtering(self, fastq_file):
        """
        Testing for proper GC percentage filtering
        """
        result = fastq_filtrator.filter_fastq(input_file=fastq_file, gc_bounds=(26, 74))

        assert len(result) == 1
        assert result[0].seq == "AGTCCGAT"

    def test_len_filtering(self, fastq_file):
        """
        Testing for proper length filtering
        """
        result = fastq_filtrator.filter_fastq(input_file=fastq_file, len_bounds=(9, 13))
        assert len(result) == 1
        assert result[0].seq == "ATATACTCATAG"

    def test_quality_filtering(self, fastq_file):
        """
        Testing for proper quality filtering
        """
        result = fastq_filtrator.filter_fastq(input_file=fastq_file, quality_threshold=30)
        assert len(result) == 2
        assert result[0].seq == "AGTCCGAT"
        assert result[1].seq == "ATATACTCATAG"

    def test_single_gc_boundary(self, fastq_file):
        """
        Testing for filtering with one GC boundary given
        """
        result = fastq_filtrator.filter_fastq(input_file=fastq_file, gc_bounds=25)

        assert len(result) == 2
        assert result[0].seq == "ATATATTA"
        assert result[1].seq == "ATATACTCATAG"

    def test_single_len_boundary(self, fastq_file):
        """
        Testing for filtering with one length boundary given
        """
        result = fastq_filtrator.filter_fastq(input_file=fastq_file, len_bounds=8)

        assert len(result) == 3
        assert result[0].seq == "AGTCCGAT"
        assert result[1].seq == "ATATATTA"
        assert result[2].seq == "GGCCAGGC"



