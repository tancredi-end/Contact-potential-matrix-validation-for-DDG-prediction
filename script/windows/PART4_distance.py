from ALL_FUNCTIONS import save_dist_array
import os

##############################################################################################
#ESSENTIAL

main_dir = "E:\\ese_tesi"

pdb_no_h = f"{main_dir}\\PART3\\pymodeller_pdb" #here the "original" structure is the predicted structure of PART3
modeller = f"{main_dir}\\PART3\\pymodeller_reverse"

mat_ori = f"{main_dir}\\PDB_multi_matrix_mut" #here the "original" structure is the predicted structure of PART3
mat_mut = f"{main_dir}\\PDB_multi_matrix_rev" 


##############################################################################################
##############################################################################################
##############################################################################################

### list1 should not be needed since you already have generated the matrixes in PART3
#list1 = os.listdir(pdb_no_h)
#for i in list1:
    #filapath1 = os.path.join(pdb_no_h,i)
    #save_dist_array(filapath1, mat_ori)

list2 = os.listdir(modeller)
for i in list2:
    filapath2 = os.path.join(modeller,i)
    save_dist_array(filapath2, mat_mut)
    