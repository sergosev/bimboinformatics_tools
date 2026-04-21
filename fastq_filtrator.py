from Bio import SeqIO, SeqUtils
from typing import Union
import argparse
import logging 

def filter_fastq(
    input_file: str,
    gc_bounds: tuple[Union[int, float], Union[int, float]] = (0, 100),
    len_bounds: tuple[int] = (0, 2**21),
    quality_threshold: Union[int, float] = 0,
    output_file: str = None,
) -> list:
    """
    Filters a fastq file with nucleic acid sequences.

    Arguments:
    - input_file: a string containing a path to input fastq file
    - gc_bounds: a tuple with GC percentage boundaries (integer or float). Default is (0, 100). Can take a single value as an upper threshold
    - length_bounds: a tuple with length boundaries (only integer) Default is (0, 2**21). Can take a single value as an upper threshold
    - quality_threshold: an integer or float number, lower boundary for mean quality. Default is 0.
    - output_file: a string containin the name of the output file. Default is None. If None, outputs the result to the stdout

    Returns a list containing SeqRecord objects.
    Saves the result to an output fastq file.
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

    if output_file != None:
        print(f"Saved results to {output_file}")
        with open(output_file, "w") as output:
            SeqIO.write(filtered, output, "fastq")
        return None
    else:
        print(*[record for record in filtered], sep="\n")
        return filtered

def main():
    parser = argparse.ArgumentParser(description="Fastq files filtering tool")

    parser.add_argument("-i", "--input-file", type=str, help="Path to the fastq file to filter")
    parser.add_argument("-o", "--output-file", type=str, help="Path to the output file")

    parser.add_argument("--gc-lower", type=float, help="Lower perventage boundary of GC contents")
    parser.add_argument("--gc-upper", type=float, help="Upper perventage boundary of GC contents")

    parser.add_argument("--len-lower", type=int, help="Lower sequence length boundary")
    parser.add_argument("--len-upper", type=int, help="Upper sequence length boundary")

    parser.add_argument("-q", "--qual", type=float, help="phred33 mean quality threshold")

    args = parser.parse_args()

    gc_bounds = (args.gc_lower or 0, args.gc_upper or 100)
    len_bounds = (args.len_lower or 0, args.len_upper or 2**21)
    
    filter_fastq(
        input_file=args.input_file,
        gc_bounds=gc_bounds,
        len_bounds=len_bounds,
        quality_threshold=args.qual or 0,
        output_file=args.output_file
    )

if __name__ == "__main__":
    main()