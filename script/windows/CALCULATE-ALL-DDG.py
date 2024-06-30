import os
import numpy as np
import ALL_FUNCTIONS as af

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# ESSENTIAL

main_dir = "E:\\ese_tesi" # main directory where you are working/saving the resulting structures and matrixes

pdb_no_h = f"{main_dir}\\PART1\\pdblist_noh" # folder with pdb wild-type with hydrogens removed
modeller_folder = f"{main_dir}\\PART2\\pymodeller_pdb" # folder with pdb mutated

mat_ori = f"{main_dir}\\PDB_S2648_matrix_wild" # folder with matrix of wild-type structure
mat_mut = f"{main_dir}\\PDB_S2648_matrix_mut" # folder with matrix of mutated structures


#where to save
result = f"{main_dir}\\PART2\\end_potential_s2648" #where to save contact potential


### calculate ddg
aasix3 = f"{main_dir}\\PART1\\aaindex3" # absolute path of file where data from contact potential matrixes from aaindex are
key = ["THOP960101","BETM990101", "MIYS990106", "BASU010101"] # CHOOSE THE KEYS YOU WANT TO USE. Only one key is fine but must be as list

choice = 5 # put here the desired cut-off for ddg calculation. the default value is 5


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

for k in key:

    #print(k)
    dictionary = af.aasix_dict(aasix3)

    true_name = os.path.basename(modeller_folder)    

    saving = f"{result}\\{true_name}_{choice}A_{k}.txt"
    #######################mettere script che genera file txt qui se non esiste nella cartella anche in PART2 e PART4
    final_result = []
    to_list1 = os.listdir(modeller_folder)
    for keke in to_list1:
        file_mut = os.path.join(modeller_folder, keke)
        name = keke.split("_")[0]


        ########## >>> RETRIEVE ORIGINAL FILE NAMES FROM MUTATED ONE <<<
        ######## v change wildname depending on situation
        wildname = f"{name}.pdb" # usually use this for names of original PDB files from RCSB. See below for other uses
        ######## ^
        #wildname = keke # this is for reverse mutation that have files with the same names in this project

        ###
        #prepwork = keke.split("_")[:3] # this is to use PDB files with different file names during the
        #joining1 = "_".join(prepwork[:3])
        #wildname = f"{joining1}.pdb"

        file_wild = os.path.join(pdb_no_h,wildname)


        mut_path = f"{os.path.splitext(keke)[0]}.txt"
        matrix_mut = os.path.join(mat_mut, mut_path)

        ######## v change wild_path depending on situation
        wild_path = f"{name}.txt" # usually use this for names for matrix of original PDB files from RCSB. See below for other uses
        ######## ^
        #wild_path = f"{os.path.splitext(keke)[0]}.txt"# the same apply for the matrix file name

        ###
        #wild_path = f"{joining1}.txt" # this is to use matrixes with different file names during the

        matrix_wild = os.path.join(mat_ori, wild_path)

        array_wild = np.loadtxt(matrix_wild)
        array_mut = np.loadtxt(matrix_mut)

        seq1,cc1,seq2,cc2 = af.chain_sequence(file_wild, file_mut)
        print(seq1, len(seq1), cc1, wildname)
        print(seq2, len(seq2), cc2, file_mut)

        checkpoint = None
        black_white = None #to change to cc1 if the mutated sequence is shorter in chained_compare_pdb_arrays
        if len(seq1) != len(seq2):
            print(">>> Difference length in sequences. Trying to rectify... ")
            checkpoint = af.final_ignore_incomplete_residues(file_wild,file_mut)
        if checkpoint:
            if len(seq1) > len(seq2):
                black_white = True
                to_use_seq = list(seq1)
                checkpoint.sort(reverse=True)
                for gh in checkpoint:
                    # print(to_use_seq)
                    # print(gh)
                    # print(checkpoint)
                    del to_use_seq[gh]
                    # print(f"> to delete position {h}")
                    array_wild = np.delete(array_wild, gh, axis=0)
                    array_wild = np.delete(array_wild, gh, axis=1)
            elif len(seq1) < len(seq2):
                to_use_seq = list(seq2)
                checkpoint.sort(reverse=True)
                for gh in checkpoint:
                    # print(to_use_seq)
                    # print(gh)
                    # print(checkpoint)
                    del to_use_seq[gh]
                    # print(f"> to delete position {h}")
                    array_mut = np.delete(array_mut, gh, axis=0)
                    array_mut = np.delete(array_mut, gh, axis=1)

            to_use_seq = ''.join(to_use_seq)
        else:
            to_use_seq = seq1

        if not black_white:
            ddg = af.chained_compare_pdb_arrays(cc1, to_use_seq, seq2, array_wild, array_mut, dictionary, k, cut_off=choice)
        else:
            ddg = af.chained_compare_pdb_arrays(cc2, to_use_seq, seq2, array_wild, array_mut, dictionary, k, cut_off=choice)
        new_line= f'{os.path.splitext(keke)[0]} {ddg}'
        print(new_line)
        final_result.append(new_line)
        #if new_line in final_result:
            #print(f"--- {keke} appended correctly to the list")
    with open(saving, "w") as res:
        print("--- --- ---", len(final_result))
        for koko in final_result:
            res.write(koko+ "\n")