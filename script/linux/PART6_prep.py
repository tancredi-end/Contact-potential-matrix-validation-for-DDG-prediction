import os
import ALL_FUNCTIONS as sv
from Bio.SeqUtils import seq3
import shutil

main_dir = "/wsl.localhost/Ubuntu/home"

pdbexp_path = f"{main_dir}/PART6/pdblist/"
multi_path = f"{main_dir}/PART1/DA_th1.0_s669.mut"

pdb_no_h = f"{main_dir}/PART6/pdb_noh/"

### Save the new pdb files using the pdb file name from the new mutation file.
rei = sv.get_mutation_info(multi_path)
pdbnames = []
for i in range(len(rei)):
    pdbfile = rei[i][0]
    pdbnames.append(pdbfile)
sv.save_pdbmmcif_files(pdbexp_path, pdbnames,"pdb")

if not os.path.exists(pdb_no_h):
    shutil.copytree(pdbexp_path,pdb_no_h)

### delete H atoms from original pdb files and hetatm residues
li = os.listdir(pdb_no_h)
for i in li:
    filename = f"{pdb_no_h}{i}"
    sv.del_hydro(filename)
    sv.remove_hetatm(filename, filename)