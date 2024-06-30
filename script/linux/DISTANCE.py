from ALL_FUNCTIONS import save_dist_array
import os

##############################################################################################
ESSENTIAL

main_dir = "/wsl.localhost/Ubuntu/home/" # main directory where you are working/saving the resulting structures and matrixes ||| you can delete this line and paste directly the specific folders you want to use. They do not need to be in the same directory

### EXAMPLE: modeller = f"{main_dir}folder1/folder2/here_are_my_things_folder"

pdb_no_h = f"{main_dir}" # abs path folder of original structures
modeller = f"{main_dir}" # abs path folder of mutated structures for modeller

mat_mut = f"{main_dir}" #abs path of distance array to save for modeller structures
mat_ori = f"{main_dir}" #abs path of distance array to save for original structures

##############################################################################################
##############################################################################################
##############################################################################################

list1 = os.listdir(pdb_no_h)
for i in list1:
    filapath1 = os.path.join(pdb_no_h,i)
    save_dist_array(filapath1, mat_ori)

list2 = os.listdir(modeller)
for i in list2:
    filapath2 = os.path.join(modeller,i)
    save_dist_array(filapath2, mat_mut)
    
list3 = os.listdir(modeller)
for i in list3:
    filapath3 = os.path.join(rosetta,i)
    save_dist_array(filapath3, mat_ros)