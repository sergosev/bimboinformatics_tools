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

def extract_translation(line: str, file):
    """
    Goes though a given line of a gbk file and extrascts AA sequences.

    Return the aminoacids from '/translation' annotation.
    """
    if '/translation="' in line:
        start = line.find('"') + 1
        if '"' in line[start:-1]:
            end = line.find('"', start)
            translation = line[start:end]
        else:
            translation = line.strip()[start:]
            line = next(file).strip()
            while not line.endswith('"'):
                translation += line.strip()
                line = next(file).strip()
            translation += line[:-1]
            
    return translation