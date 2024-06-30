
import ALL_FUNCTIONS as af

main_dir = "E:\\ese_tesi"

pdb_no_h = f"{main_dir}\\PART1\\pdblist_noH"
modeller_folder = f"{main_dir}\\PART2\\pymodeller_pdb"

mat_mut = f"{main_dir}\\PDB_S2648_matrix_mut"
mat_ori = f"{main_dir}\\PDB_S2648_matrix_wild"


### calculate ddg
aasix3 = f"{main_dir}\\PART1\\aaindex3"
key = ["THOP960101","BETM990101", "MIYS990106", "BASU010101"]

dictionary = af.aasix_dict(aasix3)

s2648_path = f"{main_dir}\\PART1\\S2648_renum_pdbnum.mut"

result = f"{main_dir}\\PART2\\end_potential_cutoff"
tried = [4,4.5,6,8,10,12]

### make a file with predicted and exp ddg and calculate pearson coefficient

for ang in tried:
    for n in key:
        modified = f"{result}\\S2648_MODELLER_{ang}A_{n}.txt"
        #af.exp_predict(s2648_path, modified)


for s in key:
    print(f"> Pearson different cutoff {s}:")
    for ang in tried:
        modified = f"{result}\\S2648_MODELLER_{ang}A_{s}.txt"
        af.pearson(modified)


