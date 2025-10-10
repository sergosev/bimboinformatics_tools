def extract_gene_name(line: str):
    """
    Goes though a given line of a gbk file and extrascts gene names.

    Return the gene name from '/gene' annotation.
    """
    if '/gene="' in line:
        start = line.find('/gene="') + 7  # Length of '/gene="'
        end = line.find('"', start)
        if end != -1:
            gene_name = line[start:end]
    return gene_name

def extract_translation(line: str):
    """
    Goes though a given line of a gbk file and extrascts AA sequences.

    Return the aminoacids from '/translation' annotation.
    """
    if '/translation="' in line:
        start = line.find('/translation="') + 13  # Length of '/gene="'
        end = line.find('"', start)
        if end != -1:
            translation = line[start:end]
    return translation