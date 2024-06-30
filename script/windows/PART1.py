import os
import sys
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
import ALL_FUNCTIONS as sv

parser = MMCIFParser(QUIET=True)
parserpdb = PDBParser(QUIET=True)

#path is where you want to make the working folders.
#put \\ in windows or / in linux at the end of path to make main_dir
path = r"E:\ese_tesi" # main directory where you are working/saving the resulting structures and matrixes
main_dir = path+"\\" 

data_path = f"{main_dir}PART1\\S2648_renum_pdbnum.mut"
file_type = "pdb"
aasix3 = f"{main_dir}PART1\\aaindex3"
filename = "S2648_renum_pdbnum_5A"

origin_dir, pdb_directory, log_directory, result_dir = sv.make_dir_path(path)

original_stdout = sys.stdout
with open(f"{main_dir}PART1\log\\log {filename}.txt", 'w') as f:
    sys.stdout = f
    print(f"mutation file is in {data_path}")
    pdb_list, file_chain, original_aa, mutation_point, icode, mutated_aa, exp_value = sv.prof_mut_file(data_path)
    #if you need to download pdb files run the next line, otherwise comment it:
    file_ex = os.listdir(pdb_directory)
    if len(file_ex) == 0:
        sv.save_pdbmmcif_files(pdb_directory, pdb_list, file_type)
    near_residues = sv.get_residue_near_aa(pdb_directory, file_type, file_chain, mutation_point, icode)
    gigadict = sv.aasix_dict(aasix3)
    sv.contact_pot(gigadict, file_chain, original_aa, mutation_point, mutated_aa, exp_value, near_residues, result_dir)
    #if you already have run the contact_pot() function, you just have to run the next line using the path to the result dir
    sv.calculate_pearson(result_dir)

sys.stdout = original_stdout