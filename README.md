# BiMbOiNfOrMaTiCs tools 👁️ 👄 👁️ 
This is my study project where I am supposed to create a python script that performs various manipulations with nucleic acid sequences and filter fastq sequences. Unfortunately for my educators the task didn't say anything specific about the naming of the repository 💀 

The tools here don't require installation, you are free to use them from your IDE. Consider the fact that all modules and scripts here were written with Python 3.14.3.
## bimbo_tools.py

This module contains classes for handling biological sequences and filtering `.fastq` files.
### BiologicalSequences handling
The module allows the user to work with Nucleic Acid Sequences and Amino Acid Sequences. 
Among nucleic acids there are 2 main classes:
- `DNASequence()`
- `RNASequence()`

Both of these classes have these methods:
- `.complement()`: returns complement sequence of original class
- `.reverse()`: returns a reversed sequence of original class
- `.reverse_complement()`: a combination of 2 methods above
- `.check_alphabet()`: checks is the sequence is truly DNA or RNA

Also `DNASequence` Class has `.transcribe()` method: returns `RNASequence` object.

AminoAcidSequence class has these methods:
- `.check_alphabet()`: confirm that given sequence is an amino acid sequence
- `.molecular_weight()`: calculate molecular weight of a given peptide/protein in Da
- `.hydrophobicity_score()`:  calculate the proportion of hydrophobic AAs in a peptide/protein

**Examples of use**
```python
dna = DNASequence('ATGTAGTATCTCATCGT')
rna = RNASequence('CUAGCUAGCUA')
prot = AminoAcidSequence('RNDLIHKMFPS')

print(dna) # ATGTAGTATCTCATCGT
print(rna[1:5]) # UAGC
print(dna.transcrive()) # AUGUAGUAUCUCAUCGU
print(rna.reverse_complement()) # UAGCUAGCUAG
print(prot.molecular_weight()) # 1357.73
print(prot.hydrophobicity_score()) # 0.45454545454545453
```
### filter_fastq()
This function takes in a fastq file. For now the function is capable of taking these parameters for filtering as its arguments:
- gc_bounds: a tuple with GC percentage boundaries (integer or float). Default is (0, 100)
- len_bounds: a tuple with length boundaries (only integer) Default is (0, $2^{32}$)
- quality_threshold: an integer or float number, lower boundary for mean quality. Default is 0.
- save_result: a bool value. If True the results of filtering are saved to a `.fastq` file. Default is True.
- output_file: a string with output file name. Default is "output_fastq.fastq"

The function returns a new, filtered list of SeqRecord objects, prints numbers of taken and filtered sequences, saves the filtered result to a file in a `./filtered` directory. 

**Example of use**
```python
filter_fastq(
	'./test_data/SRR1705851.fastq',
	gc_bounds=(20, 80.5),
	len_bounds=(20, 120),
	quality_threshold=10)
```
Output:
```python
Received 358265 sequences
Saved 0 sequences
Deleted 358265 sequences
Saving result to output_fastq.fastq
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
```python
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

А для проверяющих пасхалка. Если вы знакомы с магичкой - пишите любимого мужчину/женщину оттуда (я в восторге от Хигурумы что после манги, что после последних серий 3 сезона)

![](pics/Gojo.jpg)
