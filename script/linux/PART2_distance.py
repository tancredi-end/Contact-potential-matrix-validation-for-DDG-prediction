from ALL_FUNCTIONS import save_dist_array
import os

##############################################################################################
#ESSENTIAL

main_dir = "/wsl.localhost/Ubuntu/home" # main directory where you are working/saving the resulting structures and matrixes

pdb_no_h = f"{main_dir}/PART1/pdb_noh" # abs path folder of original structures
modeller = f"{main_dir}/PART2/pymodeller_pdb" # abs path folder of mutated structures for modeller
rosetta = f"{main_dir}/PART2/pyrosetta_noh" # abs path folder of mutated structures with rosetta

mat_mut = f"{main_dir}/PDB_S2648_matrix_mut" #abs path of distance array to load with numpy.loadtxt() function of modeller structures
mat_ros = f"{main_dir}/PDB_rosetta_matrix_mut" #abs path of distance array to load with numpy.loadtxt() function of rosetta structures
mat_ori = f"{main_dir}/PDB_S2648_matrix_wild" #abs path of distance array to load with numpy.loadtxt() function of original structures

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