import os
import ALL_FUNCTIONS as sv
import shutil

##############################################################################################
#ESSENTIAL

main_dir = "/wsl.localhost/Ubuntu/home" # main directory where you are working/saving the resulting structures and matrixes

mut_file_path = f"{main_dir}/PART1/S2648_renum_pdbnum.mut"
modeller_dir = f"{main_dir}/PART2/pymodeller_pdb"
rosetta = f"{main_dir}/PART2/pyrosetta_pdb"
try_dir = f"{main_dir}/PART2/pyrosetta_pdb/PACK/try"
pack_dir = f"{main_dir}/PART2/pyrosetta_pdb/PACK"
original_pdb = f"{main_dir}/PART1/pdblist"

part_dir = f"{main_dir}/PART2"


##############################################################################################
##############################################################################################
##############################################################################################

### 1 FOR ROSETTA: MANUALLY COPY AND PASTE THE MODEL_1 OF PDB FILES FROM THE MULTI-MODELS MUTATED FILES
###              INTO THE FOLDER WITH THE SINGLE-MODEL PDB MUTATED FILES
get_model1 = os.listdir(try_dir)
for mod in get_model1:
    if "_model_1.pdb" in mod:
        source_file = os.path.join(try_dir,mod)
        final_file = os.path.join(rosetta,mod)

        shutil.copy(source_file,final_file)

shutil.move(pack_dir, part_dir)


### 2 GENERATE A COPY OF STRUCTURE FOLDERS TO MAKE PREPARATIONS (DELETE H AND HETATM) FOR NEXT STEP
part1 = f"{main_dir}/PART1/pdb_noh"
part2 = f"{main_dir}/PART2/pyrosetta_noh"
if not os.path.exists(part1):
    shutil.copytree(original_pdb, part1)
if not os.path.exists(part2):
    shutil.copytree(rosetta, part2)


### 3 DELETE HYDROGEN ATOMS FROM ALL FILES TO CALCULATE DDG
li = os.listdir(part1)
for i in li:
    filename = f"{part1}/{i}"
    sv.del_hydro(filename)

li = os.listdir(part2)
for i in li:
    filename = f"{part2}/{i}"
    sv.del_hydro(filename)



### 4 SEEK ERROR RESIDUES (delete hetatm aminoacidic residues in the rosetta file that became ATOM residues during the elaboration of structure)
for gg in range(len(li)):
    wild_name = li[gg].split("_")[0]
    hetatm = f"{original_pdb}/{wild_name}.pdb"
    sv.del_hetatm(hetatm, li[gg]) ### example: 3eca pdb file and its mutations