# Bimbo tools 👁️ 👄 👁️ 
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
### fastq_filtrator.py
It is a CLI-style tool for filtering fastq files, capable of parsing through these arguments:
- `-i`, `--input-file` : a path to the input fastq file
- `--gc-lower`, `--gc-upper` : 2 arguments with GC perventage boundaries for filtering, can receive only one of them (the second boundary will be set to default, 0 for lowet and 100 for upper)
- `--len-lower`, `--len-upper` : 2 arguments for length boundaries, only one can be passed and the second will be set to default. Default range is (0, $2^{21}$)
- `--qual` : phred score threshold for mean quality filtering. Default is 0.
- `-o`, `--output-file` : a string with output file path. Default is None, result is printed to stdout in `fastq` format

The tool can filter fastq files based on GC contents percentage, sequence lengths and mean phred33 quality. Filtering statistics are written to `.log` file located in `logs/` directory that is created after the use of the tool. 

**Examples of use**
```bash
python \
	./fastq_filtrator.py \
	--input-file ./test_data/SRR1705851.fastq \
	--gc-lower 25 \
	--gc-upper 75 \
	--len-lower 50 \
	-q 30 \
	-o ./test_data/SRR1705851_filtered.fastq
```
Output (in a `.log` file):
```bash
INFO | 2026-04-22 00:01:15 --> Received 358265 sequences
INFO | 2026-04-22 00:01:15 --> Saved 342092 sequences
INFO | 2026-04-22 00:01:15 --> Deleted 16173 sequences
```

```bash
head ./test_data/SRR1705851_filtered.fastq

@SRR1705851.1 1/1
TTCGTGATTGTTTTCACTATCGTTCCGTTTGGCACTGCATGGTGCCCAAGGCACAGCGTTGCCGTGCTGTTGTCATTTCCAGGAAGTTTTTGAGCGAAAACCAGACATAGAATGTAGCTCAAAGCAATGATAGTCTTCATGGTTAATAG
+
,<==<<<<A@@@@@@@EEE;CEE+AC>EC;>EFFDC@=A@AE999DDD>>@E777EE75C>EF>EDEEFFFF--AE>EDEEEED=C-58AE=<D=<<DD=D9CDD@EEDED@DEDDE*9;@DDED@@@7@E*;*888@*8;@8@;;@@E
@SRR1705851.2 2/1
NATTAACCATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAACTTCCTGGAAATGACAACAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTGCCAAACGGAACGATAGTGAAAACAATCACGAATGA
+
#5<???BBEEEDEDDDGGGGGGIIIIIIIIIIIIIIIIIIIIIHIIIIFHHIIHHHHHIIIIHIIIIIIIHIIIIIIIIIIIIIIHHHHHHHHHHEHHHHHFFHHHHHHFFHHGFGGGGGGGGGGGGGEEEGCEEGGGGGEEGGGGCGEGG
@SRR1705851.4 4/1
GTGCCCAAGGCACAGCGTTGCCGTGCTGTTGTCATTTCCAGGAAGTTTTTGAGCGAAAACCAGACATAGAATGTAGCTCAAAGCAATGATAGTCTTCATGGTTAATAG
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

![](pics/Alice.jpg)
