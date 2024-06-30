from ALL_FUNCTIONS import save_dist_array
import os

##############################################################################################
#ESSENTIAL

main_dir = "/wsl.localhost/Ubuntu/home"

pdb_no_h = f"{main_dir}/PART3/pdb_noh"
modeller = f"{main_dir}/PART3/pymodeller_pdb"

mat_ori = f"{main_dir}/PDB_multi_matrix_wild"
mat_mut = f"{main_dir}/PDB_multi_matrix_mut"

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
    