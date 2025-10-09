def select_genes_from_gbk_to_fasta():
    pass
    

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

    with (open(input_file, mode="r") as input_blast,
          open(output_path, mode="w") as output_txt):
        
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

        # write names to an output file
        for protein in proteins:
            output_txt.write(protein + '\n')
        
    return proteins
