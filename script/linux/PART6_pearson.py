import os
import ALL_FUNCTIONS as sv

##############################################################################################
#ESSENTIAL

main_dir = "/wsl.localhost/Ubuntu/home"

pdb_no_h = f"{main_dir}/PART6/pymodeller_pdb_first"
modeller = f"{main_dir}/PART6/pymodeller_pdb_second"

#where to save
result = f"{main_dir}/PART5/end_potential_s669" #where to save ddg results for each mutant

aasix3 = f"{main_dir}/PART1/aaindex3"
key = ["THOP960101","BETM990101", "MIYS990106", "BASU010101"]

true_name = os.path.basename(modeller)

cutoff = 5

##############################################################################################
### calculate ddg
##############################################################################################


### make files with predicted ddg and experimental ddg
experimental = f"{main_dir}/PART1/DA_th1.0_s669_compare.mut"
for n in key:
    modi = f"{result}/{true_name}_{cutoff}A_{n}.txt"

    checker = False
    if os.path.exists(modi):
        with open(modi, "r") as checking:
            for line in checking.readlines():
                clean = line.strip()
                items = line.split()
                try:
                    control = items[2]
                    break
                except IndexError:
                    print("> Generating file with predicted adn experimental ddg...")
                    checker = True
                    break

    if checker:
        sv.retrieve_transitivity_data(experimental, modi)

### calculate pearson coefficient
iterate = os.listdir(result)
for f in iterate:
    path = result + "/" + f
    sv.pearson(path)