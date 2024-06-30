import os



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
### ESSENTIAL - PUT HERE WORKING FOLDER PATH TO GENERATE


main_dir = ""
if main_dir == "":
    main_dir = input("insert the absolute path where you want to generate the working folder")


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
if not os.path.exists(main_dir):
    os.makedirs(main_dir)

os.makedirs(f"{main_dir}/PDB_S2648_matrix_mut", exist_ok=True)
os.makedirs(f"{main_dir}/PDB_S2648_matrix_wild", exist_ok=True)
os.makedirs(f"{main_dir}/PDB_rosetta_matrix_mut", exist_ok=True)
os.makedirs(f"{main_dir}/PDB_multi_matrix_mut", exist_ok=True)
os.makedirs(f"{main_dir}/PDB_multi_matrix_wild", exist_ok=True)
os.makedirs(f"{main_dir}/PDB_multi_matrix_rev", exist_ok=True)
os.makedirs(f"{main_dir}/PDB_S669_matrix_mut", exist_ok=True)
os.makedirs(f"{main_dir}/PDB_S669_matrix_wild", exist_ok=True)
os.makedirs(f"{main_dir}/PDB_transitivity_matrix_0", exist_ok=True)
os.makedirs(f"{main_dir}/PDB_transitivity_matrix_1", exist_ok=True)
os.makedirs(f"{main_dir}/PART1", exist_ok=True)
os.makedirs(f"{main_dir}/PART1/log", exist_ok=True)
os.makedirs(f"{main_dir}/PART1/pdblist", exist_ok=True)
os.makedirs(f"{main_dir}/PART1/result", exist_ok=True)
os.makedirs(f"{main_dir}/PART2", exist_ok=True)
os.makedirs(f"{main_dir}/PART2/pymodeller_pdb", exist_ok=True)
os.makedirs(f"{main_dir}/PART2/pyrosetta_pdb", exist_ok=True)
os.makedirs(f"{main_dir}/PART2/pyrosetta_pdb/PACK", exist_ok=True)
os.makedirs(f"{main_dir}/PART2/end_potential_s2648", exist_ok=True)
os.makedirs(f"{main_dir}/PART3", exist_ok=True)
os.makedirs(f"{main_dir}/PART3/pdblist", exist_ok=True)
os.makedirs(f"{main_dir}/PART3/pymodeller_pdb", exist_ok=True)
os.makedirs(f"{main_dir}/PART3/end_potential_multi", exist_ok=True)
os.makedirs(f"{main_dir}/PART4", exist_ok=True)
os.makedirs(f"{main_dir}/PART4/pymodeller_reverse", exist_ok=True)
os.makedirs(f"{main_dir}/PART4/end_potential_reverse", exist_ok=True)
os.makedirs(f"{main_dir}/PART5", exist_ok=True)
os.makedirs(f"{main_dir}/PART5/pdblist", exist_ok=True)
os.makedirs(f"{main_dir}/PART5/pymodeller_pdb", exist_ok=True)
os.makedirs(f"{main_dir}/PART5/end_potential_s669", exist_ok=True)
os.makedirs(f"{main_dir}/PART6", exist_ok=True)
os.makedirs(f"{main_dir}/PART6/pdblist", exist_ok=True)
os.makedirs(f"{main_dir}/PART6/pymodeller_pdb_first", exist_ok=True)
os.makedirs(f"{main_dir}/PART6/pymodeller_pdb_second", exist_ok=True)