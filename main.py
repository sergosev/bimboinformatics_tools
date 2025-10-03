import sys
import os


# Setting the paths to modules folder
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
modules_path = os.path.join(script_dir, "modules")
sys.path.append(modules_path)


from is_nucleic_acid import is_nucleic_acid
from transcribe import transcribe
from reverse_transcribe import reverse_transcribe
from reverse import reverse
from complement import complement
from reverse_complement import reverse_complement

def run_dna_rna_tools(*seqs: str):
    command = seqs[-1]  #  сохраняем команду, которую нужно выполнить
    sequences = seqs[:-1]  #  сохраняем список строк для работы

    #  создаем список с ранее объявленными функциями
    procedures = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse_transcribe": reverse_transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    #  если была дана одна последовательность - возвращаем ответ для нее,
    #  более 1 - список ответов
    if len(sequences) == 1:
        seq = sequences[0]
        return procedures[command](seq) if is_nucleic_acid(seq) else f'{is_nucleic_acid(seq)}: not a nucleic acid'
    else:
        result = []
        for seq in sequences:
            if is_nucleic_acid(seq):
                result.append(procedures[command](seq))
            else:
                result.append(f'{is_nucleic_acid(seq)}: not a nucleic acid')

        return result