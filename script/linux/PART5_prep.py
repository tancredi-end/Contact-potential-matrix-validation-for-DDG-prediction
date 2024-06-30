import os
import shutil
import ALL_FUNCTIONS as sv
from Bio.SeqUtils import seq3


main_dir = "/wsl.localhost/Ubuntu/home"

pdbexp_path = f"{main_dir}/PART5/pdblist"
data_path = f"{main_dir}/PART1/Data_s669_with_predictions_sign_corrected_csv.csv"

pdb_no_h = f"{main_dir}/PART5/pdb_noh"
modeller_folder = f"{main_dir}/PART5/pymodeller_pdb"


### Save the new pdb files using the pdb file name from the new mutation file.
os.makedirs(pdbexp_path, exist_ok=True)
rei = sv.s669_data(data_path)
pdbnames = []
for i in range(len(rei)):
    pdbfile = rei[i][0]
    pdbnames.append(pdbfile)

fill = os.listdir(pdbexp_path)
if len(fill) == 0:
    sv.save_pdbmmcif_files(pdbexp_path, pdbnames,"pdb")


### make mutation files

#make a copy of original_files to modify them and keep a clean copy
if not os.path.exists(pdb_no_h):
    os.makedirs(pdb_no_h)
    shutil.copytree(pdbexp_path, pdb_no_h, dirs_exist_ok=True)

### delete H atoms from original pdb files and hetatm residues from original structures
li = os.listdir(pdb_no_h)
for i in li:
    filename = f"{pdb_no_h}/{i}"
    sv.del_hydro(filename)
    sv.del_hetatm(filename, filename)

#bi = os.listdir(modeller_folder) # it should not be needed
#for ff in bi:
    #filename = f"{modeller_folder}/{ff}"
    #sv.del_hydro(filename)

