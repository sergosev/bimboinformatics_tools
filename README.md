# BiMbOiNfOrMaTiCs tools 👁️ 👄 👁️ 
This is my study project where I am supposed to create a python script that performs various manipulations with nucleic acid sequences and filter fastq sequences. Unfortunately for my educators the task didn't say anything specific about the naming of the repository 💀 

The tools here don't require installation, you are free to use them from your IDE. Consider the fact that all modules and scripts here were written with Python 3.12.3.
## bimbo_tools.py
This script is the entry point. It contains 2 main functions: `run_dna_rna_tools()` and `filter_fastq()`. The supplementary functions needed for these 2 to work are located in **nucleic_tools.py** and **fastq_tools.py** modules in `./modules/` directory.

### run_dna_rna_tools()
This function takes a series of strings as its argument. CAUTION: the last string in the series must be one of the procedures you'd like to perform with your strings:
- **is_nucleic_acid**: checks whether give strings are nucleic acids or not. Returns bool
- **reverse**: reverts the given strings
- **transcribe**: returns transcribed (DNA to RNA) versions of given strings
- **reverse_transcribe**: returns reversely transcribed (RNA to DNA) versions of given strings
- **complement**: returns complement vesions of given sctrings
- **reverse_complement**: returns reversed complement versions of the given strings

All procedures are case-independent so you can use both UPPERCASE and lowercase letters for nucleic acid sequences.

NOTE: if one of your sequences contains both T and U (i.e. is not a nucleic acid) the result for this string will be `False` no matter what procedure you typed in.

**Examples of use**:
```{python3}
print(run_dna_rna_tools("ATUG", "ATG", "CGaT", "transcribe")) # [False, 'AUG', 'CGaU']

print(run_dna_rna_tools("CAG", "reverse_complement")) # CTG

print(run_dna_rna_tools("GAuaGuaCCucA", "is_nucleic_acid")) # True
```

### filter_fastq()
This function takes in a fastq file. For now the function is capable of taking these parameters for filtering as its arguments:
- gc_bounds: a tuple with GC percentage boundaries (integer or float). Default is (0, 100)
- length_bounds: a tuple with length boundaries (only integer) Default is (0, $2^{32}$)
- quality_threshold: an integer or float number, lower boundary for mean quality. Default is 0.

The function returns a new, filtered dictionary, prints numbers of taken and filtered sequences, saves the filtered result to a file in a `./filtered` directory. There is a flag parameter save_ouput to mark whether you want to save results to a file.

**Example of use**
```{python3}
filter_fastq(example_dict,
			gc_bounds=(20, 80.5),
			length_bounds=(50, 5000),
			quality_threshold=99.9)
```

## bio_files_processor.py
This scripts is set to read some bioinformatics file formats. For now there are 2 functions: `select_from_gbk_to_fasta()` and `parse_blast_output()`.

### select_from_gbk_to_fasta()
Parses through gbk file and selects genes flanking the genes of interest (GoIs)

Arguments:
- input_gbk: a string containing the path to gbk file
- genes: list of GoIs as strings or a string with one GoI
- n_before: an int number of flanking genes before GoI
- n_after: an int number of flanking genes after GoI
- output_fasta: a string containing the name of output fasta file

Returns None. Saves flanking genes' names and their translations to a fasta file.
The fasta file is saved to the `./processor_output` directory.

CURRENT PROBLEMS
The function does not account for:
- neighbouring GoIs
- GoIs in the end and the beggining of gbk file

In these cases it might produce incomplete results or raise errors.

**Example of use**
```
select_genes_from_gbk_to_fasta(input_gbk="example_gbk.gbk", genes="iucA",
								n_before=3, n_after=2,
								output_fasta="result_fasta.fasta")
								
select_genes_from_gbk_to_fasta(input_gbk="example_gbk.gbk", 
								genes=["iucA", "kdpD"],
								n_before=3, n_after=2,
								output_fasta="result_fasta.fasta")
```
### parse_blast_output()
Parses through a BLAST results .txt file and extracts best matches.

Arguments:
- input_file: path to the .txt file with BLAST results
- output file: a string containing the name of output .txt file

Returns a list of descriptions with best matches for each query. List is sorted alphabettically.
Writes the resulting list to an output .txt file in the "/processor_output" directory.

**Example of use**
```
parse_blast_output(input_file="example_blast_result.txt",
					output_file="desc_list.txt")
```
## Contacts
Would be glad to hear any suggestions! Especially how to deal with 7 intendation levels...
TG: @small_party


![](Python/Homeworks/HW4/bimboinformatics_tools/pics/image.webp)
![]()
