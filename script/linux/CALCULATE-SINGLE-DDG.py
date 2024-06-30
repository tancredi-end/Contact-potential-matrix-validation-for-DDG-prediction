import ALL_FUNCTIONS as af
import math
import numpy as np


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# ESSENTIAL

index_key = "THOP960101" #index key you want to use from the dictionary of dictionaries from aaindex3

mmm = "E:/ese_tesi/PART2/pyrosetta_pdb" ### folder where mutated structure is
ooo = "E:/ese_tesi/PART1/pdblist_noH" ### folder where wild-type structure is

matmmm = "E:/ese_tesi/PDB_rosetta_matrix_mut" ### folder where distance matrix to load with numpy loadtxt() function for mutated structure is
matooo = "E:/ese_tesi/PDB_S2648_matrix_wild" ### folder where distance matrix to load with numpy loadtxt() function for wild-type structure is

###

filename = "1h7m" #id name of file pdb. If there are other part inside the filename, specify them
format = ".pdb" # file extension (please use file pdb to avoid errors)

pdbo = f"{ooo}/{filename}{format}" # absolute path of wild-type structure
pdbm = f"{mmm}/1h7m_A_12ALA.pdb" # absolute path of selected mutated structure FILENAME TO PUT MANUALLY !!!
mao = np.loadtxt(f"{matooo}/{filename}.txt") # absolute path of distance matrix of wild-type structure
mam = np.loadtxt(f"{matmmm}/1h7m_A_12ALA.txt") # absolute path of distance matrix of selected mutated structure FILENAME TO PUT MANUALLY !!!

aafile = "E:/ese_tesi/PART1/aaindex3" # absolute path of file where data from contact potential matrixes from aaindex are


cutoff = 5 #maximum distance to consider to sum contact potential between aminoacids. It is 5 by default

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


def chained_compare_pdb_arrays(counter, seq1, seq2, array1, array2, index_dictionary, index_key, cut_off=cutoff):
    '''
    compare the DDG of an original pdb file with a mutated pdb file.
    counter are the start and end point to retrieve the array position for the sequence in a target chain (retrieved with chain_sequence function)
    seq1 and seq2 = the aminoacid sequences of the 2 structures you want to compare. Single aa nomenclature (type: str)
    array1 and array2 = are arrays generated with _array_pdb() function of the files you want to compare (type: array)
    index_dictionary = the dictionary generated with aasix_dict() function (type: dict of dictionaries)
    index_key = the name(aka: key) of the matrix you want to use to calculate the DDG (type: str)
    checkpoint = (_ignore_incomplete_residues()) list of int of positions where two sequences are different because
                 of a non-equal sequence length
    cut_off = the distance cut_off between a residue and the others. Set normally to 5 (type: int)
    '''
    ddg = []

    start, end = counter[0], counter[1]

    for i in range(start, end+1):
        #the two neighbour residues lists
        totoro = []
        mononoke = []

        for j in range(start, end+1):
            residue_ori = (array1[i][j])
            residue_mut = (array2[i][j])

            #print(i, j)
            # by this equation, around1 and around2 will be around 0.5 if they are in the range of cut-off. If they are out, the value drops near 0
            around1 = 1/(1+math.exp(residue_ori - cut_off))
            around2 = 1/(1+math.exp(residue_mut - cut_off))
            if around1 >= 0.5 and residue_ori != 0: # if residue_ori value is 0 it means it is comparing with itself. So skip it
                aa_key1 = seq1[i-start]+seq1[j-start]
                contact1 = float(index_dictionary[index_key][aa_key1])
                totoro.append(contact1)
            if around2 >= 0.5 and residue_mut != 0:
                aa_key2 = seq2[i-start]+seq2[j-start]
                contact2 = float(index_dictionary[index_key][aa_key2])
                mononoke.append(contact2)
        # totoro are the dg of the original pdb, mononoke are the dg of the mutated pdb
        energy_diff = sum(totoro)-sum(mononoke)
        ddg.append(energy_diff)

    result = sum(ddg)

    return result

#######################################################################################################################

seqori, cc1, seqmut, cc2 = af.chain_sequence(pdbo, pdbm)
print(seqori,cc1,seqmut,cc2)
#print(len(seqori), len(seqmut))
#print(len(mao))
#print(len(mam))


dictionary = af.aasix_dict(aafile)

checkpoint = None
if len(seqori) != len(seqmut):
    print(">>> Difference length in sequences. Trying to rectify... ")
    checkpoint = af.final_ignore_incomplete_residues(pdbo, pdbm)
if checkpoint:
    #print(checkpoint)
    #print(len(mao))
    if len(seqori) > len(seqmut):
        to_use_seq = list(seqori)
        checkpoint.sort(reverse=True)
        for gh in checkpoint:
            # print(to_use_seq)
            # print(gh)
            # print(checkpoint)
            del to_use_seq[gh]
            # print(f"> to delete position {h}")
            mao = np.delete(mao, gh, axis=0)
            mao = np.delete(mao, gh, axis=1)

    to_use_seq = ''.join(to_use_seq)
else:
    to_use_seq = seqori

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

if len(seqori) > len(seqmut): ## cc1 and cc2 can be different when comparing sequences of different lenghts so they must be exchanged to avoid index out of range
    eskere = chained_compare_pdb_arrays(cc2,to_use_seq, seqmut, mao, mam, dictionary, index_key) # calculate ddg
    print("predicted DDG =", eskere)
else:
    askere = chained_compare_pdb_arrays(cc1,to_use_seq, seqmut, mao, mam, dictionary, index_key) #calculate ddg
    print("predicted DDG =",askere)


#vvvv = af.chained_compare_pdb_arrays(cc1,to_use_seq, seqmut, mao, mam, dictionary, index_key)
#print(vvvv)