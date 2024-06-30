import ALL_FUNCTIONS as af
import os
import numpy as np

##############################################################################################
#ESSENTIAL

main_dir = "E:\\ese_tesi"

pdb_no_h = f"{main_dir}\\PART6\\pymodeller_pdb_first" # folder where first applied mutated structures are
modeller = f"{main_dir}\\PART6\\pymodeller_pdb_second" # folder where second consecutive mutated structures are

mat_ori = f"{main_dir}\\PDB_transitivity_matrix_0" #folder with distance matrix of original structure in .txt file that can be loaded with the numpy function loadtxt as an array
mat_mut = f"{main_dir}\\PDB_transitivity_matrix_1" #same but for mutated structure


#where to save
result = f"{main_dir}\\PART5\\end_potential_s669" #where to save ddg results for each mutant

aasix3 = f"{main_dir}\\PART1\\aaindex3"
key = ["THOP960101","BETM990101", "MIYS990106", "BASU010101"]

# cutoff
cutoff = 5 # put here the desired cut-off for ddg calculation. the default value is 5

##############################################################################################
##############################################################################################
##############################################################################################
#CALULATE DDG

af.calculate_all_ddg(modeller,pdb_no_h,mat_mut,mat_ori,aasix3,key, result, choice=cutoff)