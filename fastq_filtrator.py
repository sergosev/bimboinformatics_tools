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