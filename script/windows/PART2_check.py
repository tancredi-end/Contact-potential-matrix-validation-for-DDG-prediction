import os
from Bio.SeqUtils import seq1, seq3
import ALL_FUNCTIONS as sv

main_dir = "E:\\ese_tesi\\"

mut_file_path = f"{main_dir}PART1\\S2648_renum_pdbnum.mut"
modeller_dir = f"{main_dir}PART2\\pymodeller_pdb"
rosetta = f"{main_dir}PART2\\pyrosetta_pdb"
try_dir = f"{main_dir}PART2\\PACK\\try"
original_pdb = f"{main_dir}PART1\\pdblist"


### 1 CHECK IF MUTATIONS IN FILES ARE CORRECT:
#control = sv.get_mutation_info(mut_file_path)
all_data = sv.get_all_info(mut_file_path)
for j in all_data:
    mut_form = ""
    ros_form = ""
    for amin in j[1]:
        mut_form += f"{amin[1]}{amin[2]}"

        #reform = seq3(amin[2]).upper()
        #ros_form += f"{amin[1]}{reform}" #for filenames with 3letter aminoacids

    ori_check = f"{original_pdb}\\{j[0]}.pdb"
    mod_check = f"{modeller_dir}\\{j[0]}_{j[2]}_{mut_form}.pdb"
    ros_check = f"{rosetta}\\{j[0]}_{j[2]}_{mut_form}.pdb"
    if not os.path.exists(ros_check):
        ros_check = f"{rosetta}\\{j[0]}_{j[2]}_{mut_form}_model_1.pdb"


    fact = sv.check_mutation(j, ori_check, mod_check)
    if fact:
        continue
        print(*j, "> Mutation applied correctly with Modeller")
    else:
        print(*j, "> MODELLER WRONG MUTATION")


    #truth = sv.check_mutation(all_data, ori_check, try_dir) # try_dir_check is like ros_check but with PACK absolute path instead of rosetta
    truth = sv.check_mutation(j, ori_check, ros_check)
    if truth:
        continue
        #print(*j, "> Mutation applied correctly with Rosetta")
    else:
        print(*j, "> ROSETTA WRONG MUTATION")