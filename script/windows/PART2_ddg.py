import ALL_FUNCTIONS as af
import os
import numpy as np

##############################################################################################
#ESSENTIAL

main_dir = "E:\\ese_tesi" # main directory where you are working/saving the resulting structures and matrixes

pdb_no_h = f"{main_dir}\\PART1\\pdb_noh" # folder where original structure are
modeller = f"{main_dir}\\PART2\\pymodeller_pdb" # folder where mutated structure are generated with modeller
rosetta = f"{main_dir}\\PART2\\pyrosetta_noh" # folder where mutated structure are generated with rosetta

mat_ori = f"{main_dir}\\PDB_S2648_matrix_wild" # folder with distance matrix of original structure in .txt file that can be loaded with the numpy function loadtxt as an array
mat_mut = f"{main_dir}\\PDB_S2648_matrix_mut" # same but for distance matrixes of modeller structure
mat_ros = f"{main_dir}\\PDB_rosetta_matrix_mut" # same but for rosetta structure

#where to save
result = f"{main_dir}\\PART2\\end_potential_s2648" #where to save ddg results for each mutant

aasix3 = f"{main_dir}\\PART1\\aaindex3" # absolute path for the file where are the matrixes of contact potentials of aaindex3
key = ["THOP960101","BETM990101", "MIYS990106", "BASU010101"] # the name of matrixes you want to use from dictionary of dictionaries for contact potential matrixes retrieved from aaindex3. You can use 1 key but it still needs to be a list as a type

# cutoff
cutoff = 5 # put here the desired cut-off for ddg calculation. # arbitrary value. Maximum distance considered to consider two aa in contact and use the contact potential in the calculation. The default value is 5


##############################################################################################
##############################################################################################
##############################################################################################
#CALCULATE DDG


af.calculate_all_ddg(modeller,pdb_no_h,mat_mut,mat_ori,aasix3,key, result, choice=cutoff)

af.calculate_all_ddg(rosetta,pdb_no_h,mat_ros,mat_ori,aasix3,key, result, choice=cutoff)