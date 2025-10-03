import sys
import os


# Setting the paths to modules folder
script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
modules_path = os.path.join(script_dir, "modules")
sys.path.append(modules_path)


from is_nucleic_acid import is_nucleic_acid
