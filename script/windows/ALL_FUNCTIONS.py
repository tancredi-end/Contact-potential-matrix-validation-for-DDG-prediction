import os
import requests
import Bio
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import NeighborSearch
from Bio.SeqUtils import seq1
from scipy.stats import pearsonr
# |||                          |||
# VVV NEEDED FOR SOME FUNCTION VVV
parsercif = MMCIFParser(QUIET=True)
parserpdb = PDBParser(QUIET=True)
parser = PDBParser(QUIET=True)

###

def make_skeleton(ab_path):
    '''generate the working directories needed in the desired path'''
    print("work in progress")

def make_dir_path(path):
    '''download aaindex3 and returns 3 directories: mother directory, pdb files directory and a log directory'''
    directory_path = os.path.normpath(path)
    pdb_directory = os.path.join(directory_path, "PART1", "pdblist")
    log_directory = os.path.join(directory_path, "PART1","log")
    result_dir = os.path.join(directory_path, "PART1","result")
    origin_dir = directory_path+"\\PART1"
    os.makedirs(pdb_directory, exist_ok=True), os.makedirs(log_directory, exist_ok=True), os.makedirs(result_dir, exist_ok=True)

    ### >>> aaindex3 url management section. Check here if there is any change related to download aaindex3 file or anything related
    index_response = requests.get("https://www.genome.jp/ftp/db/community/aaindex/aaindex3")
    if index_response.status_code == 200:
        if not os.path.exists(origin_dir + "\\aaindex3"):
            with open(origin_dir + "\\aaindex3", "wb") as index_file:
                index_file.write(index_response.content)
                print(f">>> aaindex3 succesfully saved in {origin_dir}+\\aaindex3")
        elif os.path.exists(origin_dir + "\\aaindex3"):
            print(f">>> aaindex3 file is already in the directory")
    elif index_response == 404:
        print(f">>> ERROR: Status code: {index_response.status_code}")
    return origin_dir, pdb_directory, log_directory, result_dir

### 2 >>> File mutation retrieve infos section
def prof_mut_file(data_path):
    '''Retrieve for each line of data_path file and put in different lists: filename, filename+chain, original aa, position of mutation,
     the new aa after mutation and contact potential experimental values'''
    file_chain = []
    pdb_list = []
    original_aa = []
    mutation_point = []
    icode = []
    mutated_aa = []
    exp_value = []
    with open(data_path, "r") as mut_file:
        #retrieve pdb id, chain, original aminoacid, numeric position of the mutation, mutated aminoacid, experimental ddg value
        for r in mut_file:
            get = r.split()
            pdb_id = "".join(lower for lower in get[0] if lower.islower() or lower.isdigit())
            point = ""
            letters = ""
            for n in get[1]:
                if n.isdigit():
                    point += n
                else:
                    letters += n
            # if there is a icode to use to indicate a specific pdb residue, the lenght will be 3 letters. Otherwise only 2
            if len(letters) == 3:
                original_aa.append(letters[0])
                icode.append(letters[1])
                mutated_aa.append(letters[2])
            elif len(letters) == 2:
                original_aa.append(letters[0])
                empty = " "
                icode.append(empty)
                mutated_aa.append(letters[1])
            else:
                print("something is wrong:", r)

            file_chain.append(get[0])
            mutation_point.append(point)
            exp_value.append(get[2])
            #pdb ids are composed of 4 letters right now. If there are more, check if there's a problem
            if len(pdb_id) > 4:
                print(f">>> {pdb_id}, ATTENTION. Odd PDB file nomenclature. Check")
                pdb_id = pdb_id[0:4]
                if pdb_id not in pdb_list:
                    pdb_list.append(pdb_id.strip())
            else:
                pdb_list.append(pdb_id.strip())
    #check if you retrieved correctly each component for each line of the file
    if len(file_chain) == len(pdb_list) == len(exp_value) == len(original_aa) == len(mutation_point) == len(mutated_aa) == len(exp_value):
        print(f">>> Mutations file; all items retrieved have the same LENGTH: {len(mutation_point)}")
    else:
        print(f"{len(file_chain)}, {len(pdb_list)},{len(original_aa)},{len(mutation_point)},{len(mutated_aa)},{len(exp_value)},>>> ERROR: Mutation_point and exp_value DO NOT HAVE THE SAME LENGTH")
    return pdb_list, file_chain, original_aa, mutation_point, icode, mutated_aa, exp_value

### 3 >>> PDB download file section. Here follows the part where you download mmcif IDs you retrieved from the selected file
###   >>> REMEMBER: it uses the items names of the precedent sections. Pay attention if you need to modify them = directory where you save files and the pdb filename.
def save_pdbmmcif_files(pdb_directory, pdb_list, file_type):
    '''To do only once. Retrieve from RCSB.org the desired files from the pdb ids list'''
    for id in pdb_list:
        base_url = "https://files.rcsb.org/download/"
        file_path_pdb = os.path.join(pdb_directory, f"{id}.{file_type}")
        try:
            # Download the PDB file directly using the files.rcsb.org/download/ endpoint
            url_pdb = f"{base_url}{id}.{file_type}"
            response = requests.get(url_pdb)
            if os.path.exists(file_path_pdb):
                continue
            elif response.status_code == 200:
                with open(file_path_pdb, "wb") as file:
                    file.write(response.content)
                print(f"Data saved to: {file_path_pdb}")
            else:
                print(f"Error:{id}.{file_type} Unable to retrieve data. Status code: {response.status_code}")
        except Exception as e:
            print(f"Error: {e}")

### 4 >>> MANAGE PDB FILES: extract infos of residue and sequence
def get_residue_near_aa(pdb_directory, file_type, file_chain, mutation_point, icode):
    '''Returns a list with a sequence of aa near the selected residue for each mutation and with the desired cutoff'''
    count = 0
    angstrom_cutoff = 5.0 # change it to have a different cutoff area
    #the final list will have a sequence of aa near the selected residue with the desired cutoff
    near_residues = []
    for t in file_chain:
        file = t[0:4]
        res_list = []
        #parser iterating through each file
        try:
            the_way = pdb_directory+f"\\{file}.{file_type}"
            if "cif" in file_type:
                structure = parsercif.get_structure(f"{file}.{file_type}", the_way)
            elif "pdb" in file_type:
                structure = parserpdb.get_structure(f"{file}.{file_type}", the_way)
            #get the data you need
            chain = t[4]
            get_atoms = structure.get_atoms()
            target_residue = structure[0][chain][" ",int(mutation_point[count]), icode[count]] # " " because the target residue should always be an ATOM residue not HETATM
            around_atom = NeighborSearch(list(get_atoms))
            # search each residue that has at least 1 atom in the cutoff range | search around each atom of the target residue
            for atom in target_residue:
                # the search command is fast
                neighbors = around_atom.search(atom.coord, angstrom_cutoff)
                for ids in neighbors:
                    residue = ids.get_parent()
                    if residue not in res_list: # the object id of the residue is the same for each atom of its residue
                        res_list.append(residue)
            #Make the residues neighbor sequence
            seqq = ""
            for it in range(len(res_list)):
                res_list[it] = res_list[it].get_resname()
                seqq += (seq1(res_list[it]))
            near_residues.append(seqq)
                #res_list = []
        except Exception as problem:
            print(f">>> ATTENTION: There was a problem with the file {file}: {problem}")
        count += 1
    return near_residues

### 5 >>> Section to make a dictionary of aaindex3. The key is the H of each matrix and it has aa:contact potential as its values
def aasix_dict(aafile):
    combo = [] # aa1+aa1, aa1+aa2 etc...
    triangle = [] # contact potential values (triangle or box shape) of matrixes...
    hotkey = "" # will be the keys of gigadictionary
    dictionary = {} # combo : triangle | it is the specific matrix aa+aa : contact potential value
    gigadictionary = {} # hotkey : dictionary | dictionary of multiple matrixes | hotkey is the name of the matrix

    '''
    from https://www.genome.jp/ftp/db/community/aaindex/aaindex3
    aafile is the aaindex file
    make a dictionary of the dictionaries (the matrixes on the website page)
    '''
    with open(aafile, "r") as aasix:
        file = aasix.readlines()
        for i in file:
            if i[0] == "H":
                # H is where a new matrix start. Here you find the matrix name
                hotkey = i.strip().replace("H ", "")
                check = False #a new matrix start, check must become false to avoid iterating thw wrong lines
            elif i[0] == "M":
                # M is where you find the aminoacids letters used in the matrix. Retrieve them to make the combo aa+aa
                check = True # after the M line, there are all the contact potentiial values. We want to iterate only them
                bepo = i.split()
                prot1 = bepo[3].replace(",", "")
                prot2 = bepo[6]
                for col in prot2:
                    for rows in prot1:
                        add = rows+col
                        if add not in combo:
                            combo.append(add)
            elif check and i[0] != "/" and i[0] != "M":
                #we are at M and check is true. Each line between M and / is iterated and contact potentials are retrieved
                lista = i.strip()
                num = lista.split()
                for k in num:
                    triangle.append(k)
            elif i[0] == "/":
                # / Is the end of the matrix
                if len(triangle) == 400:
                    for t in range(len(triangle)):
                        dictionary[combo[t]] = triangle[t]
                elif len(triangle) == 210:
                    # specular matrixes have only 210 values for the 20x20 aminoacids combinations since the reverse has the same numeric value
                    # we have to iterate well to give the reverse aa+aa key the same value
                    count=0
                    skip = 20 # the aminoacids are 20
                    for j in range(210):
                        test = combo[count]
                        rev_aa = test[::-1]
                        if test not in dictionary or rev_aa not in dictionary:
                            if test[0] == test[1]: # keys such as AA have no reverse and the dictionary would just overwrite them
                                dictionary[test] = triangle[j]
                                count = count+skip
                                skip = skip-1 # each time we skip, we are also moving through the iteration changing the count value so we have to change the skip value
                            else:
                                dictionary[test] = triangle[j]
                                dictionary[rev_aa] = triangle[j]
                                count +=1
                #make the dictionary of the iterated matrix, then reset everything and go to the next matrix
                gigadictionary[hotkey] = dictionary.copy()
                hotkey = ""
                combo = []
                dictionary = {}
                triangle = []
    return gigadictionary

### 6 >>> Section where it calculates the contact potential of original and mutated residues with near residues and saves them in files
def contact_pot(aaindex, file_chain, original_aa, mutation_point, mutated_aa, exp_value, near_residues, result_dir):
    for key,value in aaindex.items():
        results = []
        for j in range(len(file_chain)):
            to_sum = []
            try:
                for jojo in near_residues[j]:
                    if jojo == "X": # ignore hetatm residues
                        continue
                    else:
                        make_ko = original_aa[j]+jojo
                        make_km = mutated_aa[j]+jojo
                        app = float(value[make_ko])-float(value[make_km])
                        to_sum.append(app)
                dg = (sum(to_sum)) #the sum of all the contact potential is the dg of the area
                results.append(f"{file_chain[j]} {original_aa[j]+mutation_point[j]+mutated_aa[j]} {exp_value[j]} {dg}")
            except Exception as key_problem:
                print(f">>> Problem: {key} a key could not exist: {key_problem} in {file_chain[j]} file")
                if key_problem in aaindex[key]:
                    print(f">>> But the key {key_problem} is in the dictionary")
        with open(f"{result_dir}\\results {key}.txt", "w") as save:
            for line in results:
                line_to_save = "".join(map(str, line))
                save.write(str(line_to_save)+ "\n")
        print(f"{result_dir}\\{file_chain[j]} {key}.txt {file_chain[j]} saved succesfully.")

### 7 >>> Pearson calculation section
def calculate_pearson(directory_path):
    file_list = os.listdir(directory_path)
    for file_name in file_list:
        file_path = os.path.join(directory_path, file_name)
        if os.path.isfile(file_path):
            with open(file_path, "r") as file:
                x_values = []
                y_values = []
                try:
                    for e in file:
                        temporary = e.split()
                        x_values.append(float(temporary[2].strip()))
                        y_values.append(float(temporary[3].strip()))
                    correlation_coefficient, p_value = pearsonr(x_values, y_values)
                    print(f"{file_name} ||| Pearson:{correlation_coefficient}, p value:{p_value}")
                except Exception as wrong:
                    print(f"ATTENTION: {file_name} {wrong}")

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

import os
import glob
import shutil
from Bio.PDB import PDBParser

def count_shared_names(directory):
    '''
    make a dictionary where key is each single mutation and values is the number of models for that specific protein mutation.
    mutations generated with the pyrosetta script.
    directory = directory path where the pdb files of mutations with multiple models are
    '''
    # Get all files in the directory
    files = glob.glob(os.path.join(directory, "*"))

    # Create a dictionary to store the counts of shared names
    counts = {}

    # Iterate through each file
    for file in files:
        # Extract the base name of the file (without extension)
        base_name = os.path.splitext(os.path.basename(file))[0]

        # Extract the shared name (name1, name2, etc.) based on the last underscore
        shared_name = base_name.rsplit("_", 1)[0]

        # Increment the count for the shared name in the dictionary
        counts[shared_name] = counts.get(shared_name, 0) + 1

    return counts


def copy_multiple_models(pdb_file,output_directory):
    model_n = _get_number_of_models(pdb_file)
    if model_n >1:
        source_file = pdb_file

        destination_directory = output_directory

        # Create the destination directory if it doesn't exist
        os.makedirs(destination_directory, exist_ok=True)

        shutil.copy(source_file, destination_directory)

        print("File copied successfully.")

def generate_output_files_prefix(file_name):
    '''
    from multiple_model_PDB_file_splitter.py by Wayne Decatur

    Takes a file name as an argument and returns string for the prefix of the
    individual output files. The generated name is based on the original file
    name.

    Specific example
    ================
    Calling function with
        ("reductase_ensemble.pdb")
    returns
        "reductase_ensemble"
    '''
    main_part_of_name, file_extension = os.path.splitext(
        file_name) #from http://stackoverflow.com/questions/541390/extracting-extension-from-filename-in-python
    if '.' in file_name:  #I don't know if this is needed with the os.path.splitext method but I had it before so left it
        return main_part_of_name
    else:
        return file_name

def extract_models(multi_model_PDB_file):
    '''
    I used multiple_model_PDB_file_splitter.py by Wayne Decatur as base

    This function takes a file containing several PDB-formatted structure models
    and extracts each individual model. Saving each individual model to a
    new file based on the name of the original file.

    Arguments for the function are as follows:
        * the file with PDB-formatted models. Requires the PDB file include
        both MODEL and ENDMDL

    The function returns the following:
        * number of models_extracted
        * the first part of the name of the created output_files
    '''
    # in preparation use the file name to generate a prefix for the
    # name of the output files.
    generated_output_files_prefix = generate_output_files_prefix(
        multi_model_PDB_file)

    #initialize values
    with open(multi_model_PDB_file, "r") as the_multi_file_stream:
        model_n = _get_number_of_models(the_multi_file_stream)



    with open(multi_model_PDB_file, "r") as reading:
        header = ""
        header_check = True
        model = ""
        model_check = False
        tail = ""
        only_model = False
        all_lines = reading.readlines()
        for i in range(model_n):
            print(f"Max number models:{model_n}, currently iterating model: {i+1}")
            for line in all_lines:
                line = line.strip()  # for better control of ends of lines
                if not line.startswith("MODEL") and header_check == True:
                    header += line + '\n'
                elif f"MODEL        {i+1}" in line and only_model == False:
                    only_model = True
                    header_check = False
                    model_check = True
                    model += line + '\n'
                elif f"MODEL       {i+1}" in line and only_model == False:
                    only_model = True
                    header_check = False
                    model_check = True
                    model += line + '\n'
                elif model_check == True and not line.startswith("ENDMDL"):
                    model += line + '\n'
                elif line.startswith("ENDMDL") and model_check == True:
                    model += line + '\n'
                    model_check = False
                elif model_check == False and header_check == False:
                    if line.startswith("MODEL") or line.startswith("ENDMDL") or line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                        continue
                    elif line.startswith("END"):
                        tail += line + '\n'
                        if i == 0:
                            new_file_text = header + model + tail
                            with open(generated_output_files_prefix + "_model_" + str(i+1) + ".pdb", "w") as output_file:
                                output_file.write(new_file_text.rstrip('\r\n'))

                            header = ""
                            header_check = True
                            model = ""
                            model_check = False
                            tail = ""
                            only_model = False
                        else:
                            new_file_text = header + model + tail
                            with open(generated_output_files_prefix + "_model_" + str(i+1) + ".pdb", "w") as output_file:
                                output_file.write(new_file_text.rstrip('\r\n'))

                            header = ""
                            header_check = True
                            model = ""
                            model_check = False
                            tail = ""
                            only_model = False
                    else:
                        tail += line + '\n'

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

import os
from Bio.SeqUtils import seq1, seq3
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser

def _get_number_of_models(pdb_file):
    '''get number of models in a pdb_file.pdb'''
    structure = parserpdb.get_structure("pdb_structure", pdb_file)

    model_count = len(structure)

    return model_count


def get_mutation_info(mut_path):
    '''
    mut_file format is: Pdbname-Chain aa_original-Residue-aa_mutated
    without the "-"
    '''
    with open(mut_path, "r") as file:
        modeller_format = []
        for line in file:
            listing = line.strip().split()
            pdb = listing[0][0:4]
            chained = listing[0][4]
            mut = []

            manyres = listing[1].strip().split(",")
            for d in manyres:
                aminoacids = ""
                residue = ""

                for i in d:
                    if i.isdigit():
                        residue = residue + i
                    else:
                        aminoacids = aminoacids + i
                # check for icode
                if len(aminoacids) == 3:
                    residue = residue + aminoacids[1]
                    aa = seq3(aminoacids[2]).upper()
                else:
                    aa = seq3(aminoacids[1]).upper()
                multi = [residue, aa]
                mut.append(multi)
            dataset = [pdb, mut, chained]
            modeller_format.append(dataset)
    return modeller_format

def get_all_info(mut_path):
    '''
    mut_file format is: Pdbname-Chain aa_original-Residue-aa_mutated
    without the "-"

    get name chain original_aa residue_position mutated_aa
    '''
    with open(mut_path, "r") as file:
        info_format = []
        for line in file:
            listing = line.strip().split()
            pdb = listing[0][0:4]
            chained = listing[0][4]
            mut = []

            manyres = listing[1].strip().split(",")
            for d in manyres:
                aminoacids = ""
                residue = ""

                for i in d:
                    if i.isdigit():
                        residue = residue + i
                    else:
                        aminoacids = aminoacids + i
                # check for icode
                if len(aminoacids) == 3:
                    residue = residue + aminoacids[1]
                    aa = aminoacids[2].upper()
                    ori = aminoacids[0].upper()
                else:
                    aa = aminoacids[1].upper()
                    ori = aminoacids[0].upper()
                multi = [ori, residue, aa]
                mut.append(multi)
            dataset = [pdb, mut, chained]
            info_format.append(dataset)
    return info_format

def check_mutation(mut_data, original_pdb, mutated_pdb):
    '''
    check a pdb file with another given a set aminoacid and mutation to see if the mutation is in the same residue position.
    If it is, the mutation was applied correctly by a previous mutational tool.

    PAY ATTENTION TO HOW THE FILES ARE NAMED OR THE SCRIPT WON'T WORK!

    mut_data = list containing pdb_name, residue position, original_aa, mutated_aa, chain (they are all str - the aa must be uppercase in 3-letters format)
    original_pdb = non mutated pdb file (ex file name: 1a43.pdb)
    mutated_pdb = mutated pdb file (ex file name: 1a43_A_167ALA)
    model_number = if it is a multi-models pdb file put the number of models as int (ex multimodel file name: 1a43_A_167ALA_model_1.pdb)
    '''
    file_typeA = os.path.splitext(original_pdb)[1]
    file_typeB = os.path.splitext(mutated_pdb)[1]
    pdb, mut, chain = mut_data
    #print(mut)
    aa_o, residue, aa_m = mut[0][0], mut[0][1], mut[0][2] #original aa, position, mutated aa

    mut_name = ""
    for amin in mut:
        mut_name += f"{amin[1]}{amin[2]}"
    #print(aa_o, residue, aa_m)
    structureA = None
    structureB = None

    check = False
    final = False

    variant_name = f"{pdb}_{chain}_{mut_name}{file_typeB}"
    #print(mutated_pdb)
    #print(variant_name)
    try:
        if ".cif" in file_typeA:
            structureA = parsercif.get_structure(f"{pdb}{file_typeA}", original_pdb)
        elif ".pdb" in file_typeA:
            structureA = parserpdb.get_structure(f"{pdb}{file_typeA}", original_pdb)
        if ".cif" in file_typeB:
            structureB = parsercif.get_structure(variant_name, mutated_pdb)
        elif ".pdb" in file_typeB:
            structureB = parserpdb.get_structure(variant_name, mutated_pdb)
    except FileNotFoundError as file_not:
        print(file_not, "not found")

    position = ""
    hetero = " "
    for e in residue:
        if e.isdigit():
            position = position+e
        else:
            hetero = e

    #just to be sure that the position is the correct one
    for model in structureA:
        aminoacid = model[chain][" ",int(position), hetero]
        #print(aminoacid.get_resname(), seq3(aa_o).upper())
        if aminoacid.get_resname() == seq3(aa_o).upper():
            trial = True
        else:
            print(f">>> Error: original residue does not match. The trial check is not true. {pdb} {chain} {aa_o}{residue}{aa_m}")

    try:
        for model2 in structureB:
            mutation = model2[chain][" ",int(position), hetero]
            #print(mutation.get_resname(), seq3(aa_m).upper())
            if mutation.get_resname() == seq3(aa_m).upper():
                check = True
                #print(f"Mutation applied correctly to the mutated file")
                #print("CORRECTOMUNDO",pdb, residue, aa_m)
            else:
                check = False
                print(f">>> Error: {mutated_pdb} mutated residue does not match. {pdb} {chain} {aa_m} {position}")
                print(f">>> {mutation.get_resname()}")
                #os.remove(mutated_file)
    #except KeyError as key:
        #print(f"||| Residue not found {key} | {pdb}, {chain}, {aa_m}, {residue}, Model: {model_number}")
        #os.remove(mutated_file)
    except FileNotFoundError as notfound:
        print(notfound, mutated_pdb)
    return check


def check_all_pdbs(data, original_dir, mutated_dir, multi_mutated_dir = None):
    '''
    check if a mutation is applied correctly to a new pdb file from an original pdb file for each file in a directory, given a file with all the mutation informations necessary

    data = file path to the text file where the mutation info are (pdb filename, original aa, residue position, mutated aa, chain)
            note: the format used as template to retrieve info is the following: 1a43A C218S
            the get_mutation_info function won't work if you have a different formatting
    original_dir = directory path of the pdb files
    mutated_path = directory path of the SINGLE-model pdb files
    multi_mutated_dir = directory path of MULTI-model pdb files

    it uses _get_number_of_models(), get_mutation_info() and check_mutation() functions to work. Both available in this python script.
    '''
    retrieved = get_mutation_info(data)
    count = 0

    modelname = ""
    for pdb, residue, ori_aa, mut_aa, chain in retrieved:
        looping_data = pdb, residue, ori_aa, mut_aa, chain

        pdb_path = f"{original_dir}\\{pdb}.pdb"
        single_mut_pdb = f"{mutated_dir}\\{pdb}_{chain}_{residue}{mut_aa}.pdb"

        if modelname != pdb:
            numbers = _get_number_of_models(pdb_path)
            modelname = pdb

        if numbers == 1:
            a = check_mutation(looping_data, pdb_path, single_mut_pdb)
            if a == True:
                count+=1

        elif numbers > 1:
            if os.path.exists(single_mut_pdb):
                b = check_mutation(looping_data, pdb_path, single_mut_pdb)
                if b == True:
                    count += 1
            if not os.path.exists(single_mut_pdb) and multi_mutated_dir != None:
                for i in range(numbers):
                    mod = i+1
                    multi_mut_pdb = f"{multi_mutated_dir}\\{pdb}_{chain}_{residue}{mut_aa}_model_{mod}.pdb"
                    c = check_mutation(looping_data, pdb_path, multi_mut_pdb)
                    if c == True:
                        count += 1

    print(count)

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

import numpy as np
import Bio.PDB
from Bio.SeqUtils import seq1
import math


def sequence_pdb(file, model=0):
    structure = parser.get_structure("structure", file)

    residues = []
    sequence = ""

    modelA = structure[model]
    for chain1 in modelA:
        for residue1 in chain1:
            if residue1.get_id()[0] == " ":
                residues.append(residue1)
                res = residue1.get_resname()
                sequence += seq1(res)

    return sequence

def del_hydro(file1):
    '''delete hydrogen atoms from a pdb file. It overwrites the selected file'''
    parser = Bio.PDB.PDBParser(QUIET=True)

    structure = parser.get_structure("structure", file1)

    for model in structure:
        for chain in model:
            for residue in chain:
                atoms_to_remove = [atom for atom in residue if atom.element == "H"]
                for atom in atoms_to_remove:
                    residue.detach_child(atom.id)

    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(file1)

def check_for_hydro(file1):
    '''check if there are still hydrogens in the pdb file structure'''
    parser = Bio.PDB.PDBParser(QUIET=True)

    structure = parser.get_structure("structure", file1)

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == "H":
                        print(f"{file1} has still hydrogens in its structure")
                        break


    ########## residue_length and del_hetatm are used to check if mutated pdb files have more residues than original pdb
def residue_length(dir1, dir2):
    '''see if the length of the original pdb and mutated pdb is the same.
    If mutated is > than original, prints out the file. See what is wrong with it'''
    ori = os.listdir(dir1)
    mut = os.listdir(dir2)
    for i in range(len(mut)):
        name = f"{dir2}\\{mut[i]}"
        sp = mut[i].split("_")

        name1 = f"{dir1}\\{sp[0]}.pdb"
        #print(name, name1)
        seqA = sequence_pdb(name1)
        seqB = sequence_pdb(name)
        seqA = seqA.replace("X", "")
        seqB = seqB.replace("X", "")
        if len(seqA) != len(seqB):
            print(f">>> {name1}, {mut[i]} |||  {len(seqA)} - {len(seqB)}")

def del_hetatm(file1, file2):
    '''delete hetatm aminoacids turned into ATOM residues in rosetta mutations'''
    structure = parserpdb.get_structure("file", file1)
    residues = []

    for model in structure:
        for chain in model:
            for residue1 in chain:
                if residue1.get_id()[0] != " ":
                    res = residue1.get_resname()
                    aa = seq1(res)
                    if aa != "X":
                        position = [chain.id, ( " ", int(residue1.get_id()[1]), residue1.get_id()[2])]
                        residues.append(position)

    structure2 = parserpdb.get_structure("hetatm", file2)
    for r, i in residues:
        delete = structure2[0][r][i]
        for model in structure2:
            for chain in model:
                if r == chain.id:
                    for residue in chain:
                        if residue.id == i:
                            chain.detach_child(delete.id)
                            print("HETATM aminoacid DELETED")
    io = Bio.PDB.PDBIO()
    io.set_structure(structure2)
    io.save(file2)

def _ignore_incomplete_residues(file1, file2):
    '''if the mutated file lacks some residue because it was incomplete the original pdb file, check the position
    of said residue and return each position where the mutated pdb has not said residue.

    the mutated residue must have a nomenclature of name_chain_residue_everythingelse.pdb
    where residue is the numbered position where the mutation is applied

    file1 and file2 are the paths to the 2 pdb files you want to compare'''

    seqA = sequence_pdb(file1)
    seqB = sequence_pdb(file2)

    name = os.path.basename(file2)
    name2 = os.path.basename(file1)
    try:
        split = name.split("_")
        type_chain = split[1]
    except IndexError as index:
        split = name2.split("_")
        type_chain = split[1]
    residue = ""
    icode = ""
    counting = 0

    numerator = []

    for i in split[2]:
        if i.isdigit():
            residue +=i
        else:
            icode += i
    if len(icode) == 4:
        icode = icode[0]
    else:
        icode = " "

    residue = int(residue)

    seqA = seqA.replace("X", "")
    seqB= seqB.replace("X", "")
    structure_ex = parserpdb.get_structure("ex", file1)

    for module in structure_ex:
        if module.id > 0:
            break
        aa = structure_ex[module.id][type_chain][" ", residue, icode]
        for chain in module:
            for res in chain:
                if res == aa:
                    residue = counting
                    break
                if res.get_id()[0] == " ":
                    counting += 1

    count = 0

    true_pos = residue
    for d in range(len(seqA)):
        try:
            if d == true_pos:
                continue
            elif seqA[d] != seqB[d-count] and d != true_pos:
                numerator.append(d)
                count +=1
        except IndexError:
            numerator.append(d)


    return numerator

def _aarray_pdb(file1, model=0):
    '''From PDB file, return array of minimal distances from each residue with each others (no HETATM)'''
    parser = Bio.PDB.PDBParser(QUIET=True)

    structureA = parser.get_structure("structure", file1)

    residuesA = []
    sequence = ""
    #combinations = []

    modelA = structureA[model] 
    for chain1 in modelA:
        for residue1 in chain1:
            #if residue1.id[0] == 'W' or residue1.id[0] == 'H':
            if residue1.id[0] != ' ':
                continue
            else:
                residuesA.append(residue1)
                res = residue1.get_resname()
                sequence += seq1(res)

    array = np.zeros((len(sequence), len(sequence)))

    for i in range(len(array)):
        for j in range(len(array)):
            controller = None


            for t in residuesA[i].get_atoms():
                xyz_zeus = t.get_coord()
                for u in residuesA[j].get_atoms():
                    xyz_mech = u.get_coord()

                    atom_dist = np.linalg.norm(xyz_zeus - xyz_mech)  # Euclidean distance
                    if controller == None or atom_dist < controller:
                        controller = atom_dist
            array[i][j] = controller

            #combinations.append((residuesA[i], residuesA[j]))
    print(f"array generated {file1}")
    return array

def compare_pdb_arrays(seq1, seq2, array1, array2, index_dictionary, index_key, cut_off=5):
    '''
    compare the DDG of an original pdb file with a mutated pdb file.
    seq1 and seq2 = the aminoacid sequences of the 2 structures you want to compare. Single aa nomenclature (type: str)
    array1 and array2 = are arrays generated with _array_pdb() function of the files you want to compare (type: array)
    index_dictionary = the dictionary generated with aasix_dict() function (type: dict of dictionaries)
    index_key = the name(aka: key) of the matrix you want to use to calculate the DDG (type: str)
    checkpoint = (_ignore_incomplete_residues()) list of int of positions where two sequences are different because
                 of a non-equal sequence length
    cut_off = the distance cut_off between a residue and the others. Set normally to 5 (type: int)
    '''
    ddg = []

    for i in range(len(array1)):
        #the two neighbour residues lists
        totoro = []
        mononoke = []

        for j in range(len(array2)):
            residue_ori = array1[i][j]
            residue_mut = array2[i][j]

            # by this equation, around1 and around2 will be around 0.5 if they are in the range of cut-off. If they are out, the value drops near 0
            around1 = 1/(1+math.exp(residue_ori - cut_off))
            around2 = 1/(1+math.exp(residue_mut - cut_off))
            if around1 >= 0.5 and residue_ori != 0: # if residue_ori value is 0 it means it is comparing with itself. So skip it
                aa_key1 = seq1[i]+seq1[j]
                contact1 = float(index_dictionary[index_key][aa_key1])
                totoro.append(contact1)
            if around2 >= 0.5 and residue_mut != 0:
                aa_key2 = seq2[i]+seq2[j]
                contact2 = float(index_dictionary[index_key][aa_key2])
                mononoke.append(contact2)
        # totoro are the dg of the original pdb, mononoke are the dg of the mutated pdb
        energy_diff = sum(totoro)-sum(mononoke)
        ddg.append(energy_diff)

    result = sum(ddg)

    return result


def calculate_directory(calculate_dir, original_pdb, aasix3, key_index, output_file):
    '''
    compare_pdb_arrays in two different directories
    calculate_dir = is the directory of predicted mutation
    original_pdb = is the directory of the original pdb proteins before mutation
    aasaix3 = the path to the aasix 3 matrix file
    key_index = the name of the matrix you want to call from the dictionary generated from aasix3
    output_file = file where you save the results
    '''
    file_list = os.listdir(calculate_dir)
    matrix = aasix_dict(aasix3)
    #create output file if it does not exist
    if not os.path.exists(output_file):
        with open(output_file, "w"):
            pass
    # check if the output file already has calculated some file (ex: if the program was interrupted and you want to restart it)
    # starts from the last non-calculated file
    with open(output_file, "r") as checking:
        reading = checking.readlines()
        for l in reading:
            for rem in file_list:
                if l.startswith(rem):
                    file_list.remove(rem)

    with open(output_file, "a") as file:
        for t in range(len(file_list)):
            pdb_mut_path = f"{calculate_dir}\\{file_list[t]}"

            name = file_list[t].split("_")
            pdb_path = f"{original_pdb}\\{name[0]}.pdb"

            if t == 0:
                seq1 = sequence_pdb(pdb_path)
                seq2 = sequence_pdb(pdb_mut_path)
                array1 = _aarray_pdb(pdb_path)
                array2 = _aarray_pdb(pdb_mut_path)

                arrow = True

            else:
                old_name = file_list[t - 1].split("_")
                if old_name[0] != name[0]:
                    seq1 = sequence_pdb(pdb_path)
                    seq2 = sequence_pdb(pdb_mut_path)
                    array1 = _aarray_pdb(pdb_path)
                    array2 = _aarray_pdb(pdb_mut_path)

                    arrow = True
                elif old_name[0] == name[0]:
                    seq2 = sequence_pdb(pdb_mut_path)
                    array2 = _aarray_pdb(pdb_mut_path)

            seq1 = seq1.replace("X","")
            seq2 = seq2.replace("X","")

            if len(seq1) > len(seq2):
                print("> Different lengths. Modifying array...")
                checkpoint = _ignore_incomplete_residues(pdb_path, pdb_mut_path)
            else:
                checkpoint = None

            if checkpoint != None:
                to_use_seq = list(seq1)
                checkpoint.sort(reverse=True)
                for h in checkpoint:
                    del to_use_seq[h]
                    if arrow == True:
                        #print(f"> to delete position {h}")
                        array1 = np.delete(array1, h, axis=0)
                        array1 = np.delete(array1, h, axis=1)

                to_use_seq = ''.join(to_use_seq)
                arrow = False

            elif checkpoint == None:
                to_use_seq = seq1

            #print(seq1)
            #print(to_use_seq)
            #print(seq2)


            energy = compare_pdb_arrays(to_use_seq, seq2, array1, array2, matrix, key_index)
            line = f"{file_list[t]}, {energy}\n"
            file.write(line)
            print(f"||| Contact energy calculated for {name[0]} /// {file_list[t]} /// {energy}")


### to calculate ddg directly from minimal distance matrix in txt files
def calculate_from_matrix(calculate_dir, original_pdb, matrix_mut, matrix_ori, aasix3, key_index, output_file):
    '''
    compare_pdb_arrays in two different directories using the minimal distance matrix file
    calculate_dir = is the directory of predicted mutation pdbs
    original_pdb = is the directory of the original pdb proteins before mutation pdbs
    matrix_ori = directory absolute path where the minimal distance matrixes are for the original pdbs
    matrix_mut = directory absolute path where minimal distance matrixes are for mutated variants pdbs
    aasaix3 = the path to the aasix3 matrix file
    key_index = the name of the matrix you want to call from the dictionary generated from aasix3
    output_file = file where you save the results
    '''
    file_list = os.listdir(calculate_dir)
    matrix = aasix_dict(aasix3)
    #create output file if it does not exist
    if not os.path.exists(output_file):
        with open(output_file, "w"):
            pass
    # check if the output file already has calculated some file (ex: if the program was interrupted and you want to restart it)
    # starts from the last non-calculated file
    with open(output_file, "r") as checking:
        reading = checking.readlines()
        for l in reading:
            for rem in file_list:
                if l.startswith(rem):
                    file_list.remove(rem)

    with open(output_file, "a") as file:
        for t in range(len(file_list)):
            pdb_mut_path = f"{calculate_dir}\\{file_list[t]}"

            name = file_list[t].split("_")
            pdb_path = f"{original_pdb}\\{name[0]}.pdb"
            #print(name[0], file_list[t])

            matrix1 = f"{matrix_ori}\\{name[0]}.txt" #retrieve array from txt file
            matrix2 = f"{matrix_mut}\\{os.path.splitext(file_list[t])[0]}.txt"

            if t == 0:
                seq1 = sequence_pdb(pdb_path)
                seq2 = sequence_pdb(pdb_mut_path)
                array1 = np.loadtxt(matrix1)
                array2 = np.loadtxt(matrix2)

                arrow = True

            else:
                old_name = file_list[t - 1].split("_")
                if old_name[0] != name[0]:
                    seq1 = sequence_pdb(pdb_path)
                    seq2 = sequence_pdb(pdb_mut_path)
                    array1 = np.loadtxt(matrix1)
                    array2 = np.loadtxt(matrix2)

                    arrow = True
                elif old_name[0] == name[0]:
                    seq2 = sequence_pdb(pdb_mut_path)
                    array2 = np.loadtxt(matrix2)

                    arrow = True

            seq1 = seq1.replace("X","")
            seq2 = seq2.replace("X","")

            if len(seq1) > len(seq2):
                print("> Different lengths. Modifying array...")
                checkpoint = _ignore_incomplete_residues(pdb_path, pdb_mut_path)
            elif len(seq2) > len(seq1):
                checkpoint = _ignore_incomplete_residues(pdb_mut_path,pdb_path)
            else:
                checkpoint = None

            if checkpoint != None:
                #print("ECCOCI")
                if len(seq1) > len(seq2):
                    to_use_seq = list(seq1)
                    checkpoint.sort(reverse=True)
                    for h in checkpoint:
                        del to_use_seq[h]
                        if arrow == True:
                            #print(f"> to delete position {h}")
                            array1 = np.delete(array1, h, axis=0)
                            array1 = np.delete(array1, h, axis=1)
                    to_use_seq = ''.join(to_use_seq)
                    arrow = False
                elif len(seq1) < len(seq2):
                    to_use_seq = list(seq2)
                    checkpoint.sort(reverse=True)
                    for h in checkpoint:
                        del to_use_seq[h]
                        if arrow == True:
                            #print(f"> to delete position {h}")
                            array2 = np.delete(array2, h, axis=0)
                            array2 = np.delete(array2, h, axis=1)
                    to_use_seq = ''.join(to_use_seq)
                    arrow = False

            elif checkpoint == None:
                to_use_seq = seq1

            #print(seq1)
            #print(to_use_seq)
            #print(seq2)


            energy = compare_pdb_arrays(to_use_seq, seq2, array1, array2, matrix, key_index)

            nameid = os.path.splitext(file_list[t])[0]
            line = f"{nameid} {energy}\n"
            file.write(line)
            print(f"||| Contact energy calculated for {name[0]} /// {nameid} /// {energy}")

def save_dist_array(pdb_file, output_directory):
    '''
    save the array of the minimal distance between residues
    pdb_file = the pdb to generate the array
    output_directory =  directory where you save the results
    '''

    filename = os.path.basename(pdb_file)
    name = filename.replace(".pdb","")

    save_output = f"{output_directory}\\{name}.txt"
    if not os.path.exists(save_output):
        array = _aarray_pdb(pdb_file)
        np.savetxt(save_output, array)
        print(f"--- New file generated at {save_output} ---")

    else:
        print(f">>> Minimal distance matrix for {name} already exists")




############### CALCULATE NEW PEARSON
def exp_predict(file_exp, file_predict):
    '''make a file with both predicted and experimental data'''
    file = []
    with open(file_exp, "r") as exp:
        #retrieve from the file name the mutation information
        for line in exp:
            trueline = line.strip()
            elements = trueline.split()
            part1 = elements[0]
            part1 = part1[:4] + "_" + part1[4:]

            l_num = []
            l_icode = []
            l_aa = []
            part2 = elements[1]
            multi_mut = part2.split(",")
            for s in multi_mut:
                number = ""
                icode = ""
                aa = ""
                for o in s:
                    if o.isdigit():
                        number += o
                    else:
                        aa += o
                if len(aa) == 3:
                    icode = aa[1]
                    aa = aa[2].upper()
                else:
                    aa = aa[1].upper()
                l_num.append(number)
                l_icode.append(icode)
                l_aa.append(aa.upper())

            #generate the file name to check if it is the correct file line
            check = f"{part1}"
            for d in range(len(l_num)):
                check += f"_{l_num[d]}{l_icode[d]}{l_aa[d]}"
            print(check)

            #append the new line the experimental ddg value if check is true
            with open(file_predict, "a+") as compare:
                compare.seek(0)
                reading = compare.readlines()

                for li in reading:
                    if check in li:
                        namae = li.replace(",", "")
                        namae = namae.replace(".pdb", "")
                        new_line = namae.strip() + " " + elements[2] + "\n"
                        file.append(new_line)
                        break
    #rewrite the file with all the new lines appende in the list
    with open(file_predict, "w") as rewrite:
        full_data = "".join(file)
        rewrite.write(full_data)

def pearson(file_path):
    #calculate pearcson correlation | experimental ddg (y) | predicted ddg (x)
    with open(file_path, "r") as file:
        x_values = []
        y_values = []
        for e in file:
            temporary = e.split()
            x_values.append(float(temporary[2].strip()))
            y_values.append(float(temporary[1].strip()))
        correlation_coefficient, p_value = pearsonr(x_values, y_values)
        print(f"{os.path.basename(file_path)} ||| Pearson:{correlation_coefficient}, p value:{p_value}")


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

def mut_data(mut_path):
    '''
    retrieve mutation data from file with following format:
    1a43A R167A -4.55
    or for multiple mutations
    1amqA C82A,C191A,C192A,C270A,C401A -1.0
    '''
    with open(mut_path, "r") as file:
        modeller_format = []
        for line in file:
            listing = line.strip().split()
            pdb = listing[0][0:4]
            chained = listing[0][4]
            mut = []

            manyres = listing[1].strip().split(",")
            for d in manyres:
                aminoacids = ""
                residue = ""

                for i in d:
                    if i.isdigit():
                        residue = residue + i
                    else:
                        aminoacids = aminoacids + i
                # check for icode
                if len(aminoacids) == 3:
                    residue = residue + aminoacids[1]
                    aa = seq3(aminoacids[2]).upper()
                else:
                    aa = seq3(aminoacids[1]).upper()
                multi = [residue, aa]
                mut.append(multi)
            dataset = [pdb, mut, chained]
            modeller_format.append(dataset)
        return modeller_format


def s669_data(mut_path):
    '''
    retrieve mutation data from file with following format:
    1A7VA;D3H;D3H;-1.36
    '''
    with open(mut_path, "r") as file:
        s669_format = []
        for line in file:
            if line == "Protein;PDB_Mut;Mut_seq;DDG_checked_dir\n":
                pass
            else:
                listing = line.strip().split(";")
                pdb = listing[0][0:4]
                chained = listing[0][4]
                mut = []

                manyres = listing[1].strip().split(",")
                for d in manyres:
                    aminoacids = ""
                    residue = ""

                    for i in d:
                        if i.isdigit():
                            residue = residue + i
                        else:
                            aminoacids = aminoacids + i
                    # check for icode
                    if len(aminoacids) == 3:
                        residue = residue + aminoacids[1]
                        aa = seq3(aminoacids[2]).upper()
                    else:
                        aa = seq3(aminoacids[1]).upper()
                    multi = [residue, aa]
                    mut.append(multi)
                dataset = [pdb, mut, chained]
                s669_format.append(dataset)
        return s669_format


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


# to calculate reverse ddg

def rev_mut_data(mut_path):
    '''retrieve data for reverse mutation'''
    with open(mut_path, "r") as file:
        modeller_format = []
        for line in file:
            listing = line.strip().split()
            pdb = listing[0][0:4]
            chained = listing[0][4]
            mut = []

            manyres = listing[1].strip().split(",")
            for d in manyres:
                aminoacids = ""
                residue = ""

                for i in d:
                    if i.isdigit():
                        residue = residue + i
                    else:
                        aminoacids = aminoacids + i
                # check for icode
                if len(aminoacids) == 3:
                    residue = residue + aminoacids[1]
                    aa = seq3(aminoacids[0]).upper()
                    ori = seq3(aminoacids[2]).upper()
                else:
                    aa = seq3(aminoacids[0]).upper()
                    ori = seq3(aminoacids[1]).upper()
                multi = [ori, residue, aa]
                mut.append(multi)
            dataset = [pdb, mut, chained]
            modeller_format.append(dataset)
        return modeller_format

def calculate_reverse_directory(calculate_dir, original_pdb, aasix3, key_index, output_file):
    '''
    compare_pdb_arrays in two different directories
    calculate_dir = is the directory of predicted mutation
    original_pdb = is the directory of the original pdb proteins before mutation
    aasaix3 = the path to the aasix 3 matrix file
    key_index = the name of the matrix you want to call from the dictionary generated from aasix3
    output_file = file where you save the results
    '''
    file_list = os.listdir(calculate_dir)
    new_list = os.listdir(original_pdb)
    matrix = aasix_dict(aasix3)
    # check if the output file already has calculated some file (ex: if the program was interrupted and you want to restart it)
    # starts from the last non-calculated file
    with open(output_file, "r") as checking:
        reading = checking.readlines()
        for l in reading:
            for rem in file_list:
                if l.startswith(rem):
                    file_list.remove(rem)

    with open(output_file, "a") as file:
        for t in range(len(file_list)):
            pdb_mut_path = f"{calculate_dir}\\{file_list[t]}"

            name = file_list[t].split("_")

            pdb_path = f"{original_pdb}\\{file_list[t]}"


            seq1 = sequence_pdb(pdb_path)
            seq2 = sequence_pdb(pdb_mut_path)
            array1 = _aarray_pdb(pdb_path)
            array2 = _aarray_pdb(pdb_mut_path)

            arrow = True

            seq1 = seq1.replace("X","")
            seq2 = seq2.replace("X","")

            if len(seq1) > len(seq2):
                print("> Different lengths. Modifying array...")
                checkpoint = _ignore_incomplete_residues(pdb_path, pdb_mut_path)
            else:
                checkpoint = None

            if checkpoint != None:
                to_use_seq = list(seq1)
                checkpoint.sort(reverse=True)
                for h in checkpoint:
                    del to_use_seq[h]
                    if arrow == True:
                        print(f"> position to delete {h}")
                        array1 = np.delete(array1, h, axis=0)
                        array1 = np.delete(array1, h, axis=1)

                to_use_seq = ''.join(to_use_seq)
                arrow = False

            elif checkpoint == None:
                to_use_seq = seq1

            #print(seq1)
            #print(to_use_seq)
            #print(seq2)

            energy = compare_pdb_arrays(seq2, to_use_seq, array2, array1, matrix, key_index)
            line = f"{file_list[t]}, {energy}\n"
            file.write(line)
            print(f"||| Contact energy calculated for {os.path.basename(pdb_mut_path)} /// {os.path.basename(pdb_path)} /// {energy}")

def reverse_ddg(file1, file2):
    '''make a single file with the ddg from file 1 and file 2. The ddg must be in the 2nd column'''

    new_file = []

    with open(file1, "r") as original:
        ori = original.readlines()
        with open(file2, "r") as mutation:
            mut = mutation.readlines()
            for line in ori:
                parts = line.strip().split()
                for name in mut:
                    retrieve = name.strip().split()
                    check = retrieve[0].replace(".pdb,", "")
                    if check == parts[0]:
                        new_line = f"{check} {retrieve[1]} {parts[1]} \n"
                        new_file.append(new_line)

####

def retrieve_transitivity_data(file1, file2):
    '''make a single file with the ddg from file 1 and file 2. The ddg must be in the 2nd column in both files'''
    new_file = []

    with open(file1, "r") as original:
        ori = original.readlines()
        with open(file2, "r") as mutation:
            mut = mutation.readlines()
            for line in ori:
                parts = line.strip().split()
                for name in mut:
                    retrieve = name.strip().split()
                    check = retrieve[0].replace(".pdb,", "")
                    if check == parts[0]:
                        new_line = f"{check} {retrieve[1]} {parts[1]} \n"
                        new_file.append(new_line)

    #print(new_file)
    #rewrite the file with all the new lines appended in the list
    with open(file2, "w") as rewrite:
        full_data = "".join(new_file)
        rewrite.write(full_data)



#######################################################################################################################

# function used in part5 to delete hetatm residue
def remove_hetatm(input_pdb, output_pdb):
    # Create a PDB parser
    parser = PDBParser(QUIET=True)

    # Parse the structure
    structure = parser.get_structure('structure', input_pdb)

    # Create a Select class to filter out HETATM residues
    class NonHetSelect(Bio.PDB.Select):
        def accept_residue(self, residue):
            if residue.id[0] == ' ':
                return True
            else:
                return False

    # Create a PDBIO object
    io = Bio.PDB.PDBIO()

    # Set the structure to save
    io.set_structure(structure)

    # Save the structure while excluding HETATM residues
    io.save(output_pdb, NonHetSelect())


#######################################################################################################################

def final_ignore_incomplete_residues(file1, file2):
    '''if the mutated file lacks some residue because it was incomplete the original pdb file, check the position
    of said residue and return each position where the mutated pdb has not said residue.

    the mutated residue must have a nomenclature of name_chain_residue_everythingelse.pdb
    where residue is the numbered position where the mutation is applied

    file1 and file2 are the paths to the 2 pdb files you want to compare'''

    seqA, cc1, seqB, cc2 = chain_sequence(file1, file2)

    name = os.path.basename(file2)
    name2 = os.path.basename(file1)
    try:
        split = name.split("_")
        type_chain = split[1]
    except IndexError as index:
        split = name2.split("_")
        type_chain = split[1]
    residue = ""
    icode = ""
    counting = 0

    numerator = []

    for i in split[2]:
        if i.isdigit():
            residue +=i
        else:
            icode += i
    if len(icode) == 4:
        icode = icode[0]
    else:
        icode = " "

    residue = int(residue)

    seqA = seqA.replace("X", "")
    seqB= seqB.replace("X", "")

    structure_ex = parserpdb.get_structure("ex", file1)


    for module in structure_ex:
        if module.id > 0:
            break
        aa = structure_ex[module.id][type_chain][" ", residue, icode]
        for chain in module:
            if chain.id == type_chain:
                for res in chain:
                    if res == aa:
                        residue = counting
                    if res.get_id()[0] == " ":
                        counting += 1

    count = 0

    true_pos = residue
    for d in range(len(seqA)):
        try:
            if d == true_pos:
                continue
            elif seqA[d] != seqB[d-count] and d != true_pos:
                numerator.append(d)
                count +=1
        except IndexError:
            numerator.append(d)


    return numerator


def chain_sequence(file1, file2, model=0):
    '''
    file1 is the wild type pdb file
    file2 is the mutated pdb file
    retrieve chain sequence  (using mutated file name to get the chain)
    return seq of original seq chain, mutated seq chain and their starting and ending positions in the overall sequence(as 2 lists)
    '''
    structure1 = parserpdb.get_structure("original", file1)
    structure2 = parserpdb.get_structure("mutated", file2)

    model1 =structure1[model]
    model2 = structure2[model]

    filename2 = os.path.basename(file2)
    try:
        chain_mut = filename2.split("_")[1]
    except IndexError:
        filename1 = os.path.basename(file1)
        chain_mut = filename1.split("_")[1]

    convention = 1
    damn_chain = "A"

    counting = 0
    new_counting = 0
    new_counting1=0
    start=None
    end=None
    lance = True #check if we enter the target chain
    shield = True # check when we exit the target chain
    arrow_ori = False
    counter1 = []
    sequence1 = ""

    initial_d = "A" #check if a pdb file starts with a chain different from A

    for chain1 in model1:
        if counting == 0 and chain1.id != damn_chain:
            sel = chain1.id
            if sel.isdigit():
                damn_chain = convention
            else:
                initial_d = chain1.id
        for residue1 in chain1:
            hellcome = residue1.get_resname()
            if residue1.get_id()[0] == " " and seq1(hellcome) != "":
                #print(residue1, counting)
                counting +=1
                if chain1.id == chain_mut:
                    #print(residue1, counting)
                    if lance:
                        start = counting
                        #print(residue1)
                        if chain_mut == damn_chain:
                            start = 0
                        lance = False
                    if chain_mut == chain1.id and initial_d != damn_chain:
                        arrow_ori = True
                        new_counting1 +=1
                    if residue1.get_id()[0] == " ":
                        #print(residue1.get_id)
                        res = residue1.get_resname()
                        sequence1 += seq1(res)
                elif chain1.id != chain_mut:
                    if shield and not lance:
                        end = counting-1
                        shield = False
                        break
            #end_mut = counting_mut
        if chain_mut == damn_chain and initial_d == damn_chain:
            end = counting - 1
            break
        if not shield and chain_mut != damn_chain:
            end = counting - 1
            break
        elif shield and chain_mut != damn_chain and not lance:
            end = counting - 1
            break

        if arrow_ori == True:
            end = new_counting1 - 1


    counter1.append(start)
    counter1.append(end)

    counting_mut = 0
    start_mut=None
    end_mut=None
    lance_mut = True #check if we enter the target chain
    shield_mut = True # check when we exit the target chain
    arrow = False # to odd chain order
    counter2 = []
    sequence2 = ""


    for chain2 in model2:
        for residue2 in chain2:
            welcome = residue2.get_resname()
            #print(welcome, seq1(welcome), "XXX")
            if residue2.get_id()[0] == " " and seq1(welcome) != "":
                counting_mut +=1
                if chain2.id == chain_mut:
                    #print(residue1, counting)
                    if lance_mut:
                        start_mut = counting_mut
                        if chain_mut == damn_chain:
                            start_mut = 0
                        lance_mut = False
                    if chain_mut == chain2.id and initial_d != damn_chain:
                        arrow = True
                        new_counting +=1
                    if residue2.get_id()[0] == " ":
                        #print(residue1.get_id)
                        res = residue2.get_resname()
                        sequence2 += seq1(res)
                elif chain2.id != chain_mut:
                    if shield_mut and not lance_mut:
                        end_mut = counting_mut-1
                        shield_mut = False
                        break
            #end_mut = counting_mut
        if chain_mut == damn_chain and initial_d == damn_chain:
            end_mut = counting_mut - 1
            break
        if not shield_mut and chain_mut != damn_chain:
            end_mut = counting_mut-1
            break
        elif shield_mut and chain_mut != damn_chain and not lance_mut:
            end_mut = counting_mut-1
            break

        if arrow == True:
            end_mut = new_counting-1




    counter2.append(start_mut)
    counter2.append(end_mut)

    return sequence1, counter1, sequence2, counter2

def chained_compare_pdb_arrays(counter, seq1, seq2, array1, array2, index_dictionary, index_key, cut_off=5):
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
                #print(seq1[i-start])
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


def calculate_all_ddg(mutant_folder, pdb_no_h, mat_mut, mat_ori, aasix3, key, result, choice=5):
    '''
    calculate the ddg of each pdb file in a folder

    mutant_folder = absolute path of folder where the mutant are (.pdb)
    pdb_no_h = absolute path of folder where original structures are (.pdb)
    mat_mut = array of minimal distances of mutant files retrieved from a .txt file using numpy loadtxt()
    mat_mut = array of minimal distances of original files retrieved from a .txt file using numpy loadtxt()
    aasix3 = absolute path of the file with the contact potentials matrixes from aaindex3
    key = list with the name of the matrix from the dictionary of dictionaries with the matrixes of contact potentials retrieved from aaindex3 (ex: BASU010101)
            key must be a list with the name of the matrixes as items (strings type), even if you use a single matrix
    result = absolute path of the folder where you want to save the results

    choice = cutoff to consider the max distance from a residue and its neighbours. It is 5 by default

    '''
    for k in key:

        #print(k)
        dictionary = aasix_dict(aasix3)
        foldname = os.path.basename(mutant_folder)
        saving = f"{result}\\{foldname}_{choice}A_{k}.txt"

        final_result = []
        to_list1 = os.listdir(mutant_folder)
        for keke in to_list1:
            file_mut = os.path.join(mutant_folder, keke)
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

            seq1,cc1,seq2,cc2 = chain_sequence(file_wild, file_mut)
            print(seq1, len(seq1), cc1, wildname)
            print(seq2, len(seq2), cc2, file_mut)

            checkpoint = None
            black_white = None #to change to cc1 if the mutated sequence is shorter in chained_compare_pdb_arrays
            if len(seq1) != len(seq2):
                print(">>> Difference length in sequences. Trying to rectify... ")
                checkpoint = final_ignore_incomplete_residues(file_wild,file_mut)
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
                ddg = chained_compare_pdb_arrays(cc1, to_use_seq, seq2, array_wild, array_mut, dictionary, k, cut_off=choice)
            else:
                ddg = chained_compare_pdb_arrays(cc2, to_use_seq, seq2, array_wild, array_mut, dictionary, k, cut_off=choice)
            new_line= f'{os.path.splitext(keke)[0]} {ddg}'
            print(new_line)
            final_result.append(new_line)
            #if new_line in final_result:
                #print(f"--- {keke} appended correctly to the list")
        with open(saving, "w") as res:
            print("--- --- ---", len(final_result))
            for koko in final_result:
                res.write(koko+ "\n")