import sys
import os

# Setting the path to modules folder
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
modules_path = os.path.join(script_dir, "modules")
sys.path.append(modules_path)

from typing import Union
import processor_tools as pt

def select_genes_from_gbk_to_fasta(input_gbk: str, 
                                   genes: Union[str, list[str]], 
                                   n_before: int = 1, 
                                   n_after: int = 1,
                                   output_fasta: str = "gbk_to_fasta_result.fasta"):
    
    import os
    import sys

    #setting up directories
    work_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    if not os.path.exists(os.path.join(work_dir, 'processor_output')):
        os.mkdir(os.path.join(work_dir, "processor_output"))
    output_path = os.path.join(work_dir, "processor_output", output_fasta)

    buffer = {}
    results = {}

    with open(input_gbk, mode="r") as input_gbk:
        for line in input_gbk:
            line = line.strip()
            if line.startswith('/gene'):
                current_gene = pt.extract_gene_name(line)
                while not line.startswith('/translation'):
                    line = next(input_gbk).strip()
                print(line)
                curr_translation = pt.extract_translation(line) # gotta fix
                buffer[current_gene] = curr_translation
            
                if len(buffer) > n_before + 1:
                    buffer.pop(next(iter(buffer)))
                
                if current_gene in genes:
                    for key in list(buffer.keys()):
                        if key != list(buffer.keys())[-1]:
                            results[key] = buffer.pop(key)
                
                    line = next(input_gbk)
                    i =  0
                    while i < n_after:
                        while not line.startswith('/gene'):
                            line = next(input_gbk).strip()
                        gene = pt.extract_gene_name(line)
                        while not line.startswith('/translation'):
                            line = next(input_gbk).strip()
                        translation = pt.extract_translation(line) # gotta fix
                        results[gene] = translation
                        i+=1

    with open(output_path, mode="w") as output_fasta:
        for gene, translation in results.items():
            output_fasta.write(">" + gene + "\n")
            output_fasta.write(translation + "\n")
            output_fasta.write("\n")

select_genes_from_gbk_to_fasta(input_gbk='/mnt/c/Users/serol/Downloads/example_gbk.gbk',
                               genes = ["kdpB", "kefC_2", "vgrG1_3"],
                               n_after=2, n_before=2)

def parse_blast_output(input_file: str, output_file: str = "parse_output.txt"):
    """
    Parses through a BLAST results .txt file.

    Arguments:
    - input_file: path to the .txt file with BLAST results
    - output file: a string containing the name of output .txt file

    Returns a list of descriptions with best matches for each query. List is sorted alphabettically.
    Writes the resulting list to an output .txt file in the "/processor_output" directory.
    """

    import os
    import sys

    #setting up directories
    work_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    if not os.path.exists(os.path.join(work_dir, 'processor_output')):
        os.mkdir(os.path.join(work_dir, "processor_output"))
    output_path = os.path.join(work_dir, "processor_output", output_file)

    with open(input_file, mode="r") as input_blast:
        proteins = []
        for line in input_blast:
            if line.startswith("Query #"): # found Query
                while not line.startswith("Description"): # Getting to the table
                    line = next(input_blast)
                line = next(input_blast) # now in the first row

                # select exactly Description column value
                symbol_num = 0
                in_desc_col = True
                while in_desc_col:
                    if (line[symbol_num: symbol_num+2] == "  " or
                        line[symbol_num: symbol_num+2] == ". "):
                        in_desc_col = False
                    symbol_num += 1

                # save the name to the resulting list
                proteins.append(line[0:symbol_num])

    proteins.sort(key=str.lower)
    with open(output_path, mode="w") as output_txt:
        # write protein names to an output file
        for protein in proteins:
            output_txt.write(protein + '\n')
        
    return proteins
