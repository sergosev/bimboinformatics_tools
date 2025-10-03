# BiMbOiNfOrMaTiCs tools 👁️ 👄 👁️ 
This is my study project where I am supposed to create a python script that performs various manipulations with nucleic acid sequences and filter fastq sequences. Unfortunately for my educators the task didn't say anything specific about the naming of the repository 💀 

The tools here don't require installation, you are free to use them from your IDE. Consider the fact that all modules and scripts here were written with Python 3.12.3.
## main.py
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
This function takes on a dictionary of fastq sequences with theis IDs as keys and items as tuples of sequence itself and phed33 quality scores per each nucleotide. For now the function is capable of taking these parameters for filtering as its arguments:
- gc_bounds: a tuple with GC percentage boundaries (integer or float). Default is (0, 100)
- length_bounds: a tuple with length boundaries (only integer) Default is (0, $2^{32}$)
- quality_threshold: an integer or float number, lower boundary for mean quality. Default is 0.

The function returns a new, filtered dictionary.

**Example of use**
```{python3}
filter_fastq(example_dict,
			gc_bounds=(20, 80.5),
			length_bounds=(50, 5000),
			quality_threshold=99.9)
```

## Contacts
Would be glad to hear any suggestions! 
TG: @small_party

To BI teachers and curators: Please don't expel me for this meme 👺 ![](https://sun9-47.userapi.com/s/v1/if2/v2M8zw2_kcArG0kwZ286hCzTl4CuS5vc8VacQy2gx6SMyxUjEydbNNetLyQaCzGnnIGHUkOPegfpXqdRY-pLzTvH.jpg?quality=95&as=32x29,48x44,72x66,108x99,160x147,240x221,360x331,480x442,540x497,640x589,720x663,828x762&from=bu&u=CS0R8gurS81rusC4tCPCHgewUBC5Sf3pU3mstqEbKeI&cs=828x0)
