import os
import glob
import shutil
from Bio.PDB import PDBParser

main_dir = "home\\dede\\linux\\test"

multi_save = f"{main_dir}\\PART2\\pyrosetta_pdb\\PACK"
pdb_list = f"{main_dir}\\PART1\\pdblist"

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

def _get_number_of_models(pdb_file):
    '''get number of models in a pdb file'''
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb_structure", pdb_file)

    model_count = len(structure)

    return model_count

def copy_multiple_models(pdb_file, output_directory):
    model_n = _get_number_of_models(pdb_file)
    if model_n > 1:
        source_file = pdb_file

        destination_directory = output_directory

        # Create the destination directory if it doesn't exist
        os.makedirs(destination_directory, exist_ok=True)

        shutil.copy(source_file, destination_directory)

        print(f"File copied successfully: {pdb_file}")

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
    main_part_of_name, file_extension = os.path.splitext(file_name) 
    return main_part_of_name

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
    generated_output_files_prefix = generate_output_files_prefix(multi_model_PDB_file)
    
    # Initialize values
    with open(multi_model_PDB_file, "r") as file:
        all_lines = file.readlines()
    
    model_index = 0
    model_content = []
    writing_model = False

    for line in all_lines:
        if line.startswith("MODEL"):
            writing_model = True
            model_index += 1
            model_content = [line]
        elif line.startswith("ENDMDL") and writing_model:
            model_content.append(line)
            writing_model = False
            
            # Write the model to a new file
            output_file_name = f"{generated_output_files_prefix}_model_{model_index}.pdb"
            with open(output_file_name, "w") as output_file:
                output_file.writelines(model_content)
            model_content = []
        elif writing_model:
            model_content.append(line)

##### --------------- script examples: --------------- #####
for file in os.listdir(pdb_list):
    file_path = os.path.join(pdb_list, file)
    if os.path.isfile(file_path):
        # Process the file here
        copy_multiple_models(file_path, multi_save)

# Run extraction on the copied files
for file in os.listdir(multi_save):
    file_path = os.path.join(multi_save, file)
    extract_models(file_path)





##### multiple_model_PDB_file_splitter.py by Wayne Decatur #####

### Then apply mutation using the rosetta script to each model file

### then unite the models in a single pdb file if you want