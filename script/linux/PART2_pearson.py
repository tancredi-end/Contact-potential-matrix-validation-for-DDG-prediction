import os
import ALL_FUNCTIONS as sv

##############################################################################################
#ESSENTIAL

main_dir = "/wsl.localhost/Ubuntu/home" # main directory where you are working/saving the resulting structures and matrixes

pdb_no_h = f"{main_dir}/PART1/pdb_noh" # folder where original structure are
modeller = f"{main_dir}/PART2/pymodeller_pdb" # folder where mutated structure are generated with modeller
rosetta = f"{main_dir}/PART2/pyrosetta_noh" # folder where mutated structure are generated with rosetta

#where to save
result = f"{main_dir}/PART2/end_potential_s2648" # where to save ddg results for each mutant

aasix3 = f"{main_dir}/PART1/aaindex3" # absolute path for the file where are the matrixes of contact potentials of aaindex3
key = ["THOP960101","BETM990101", "MIYS990106", "BASU010101"] # the name of matrixes you want to use from dictionary of dictionaries for contact potential matrixes retrieved from aaindex3. You can use 1 key but it still needs to be a list as a type

cutoff = 5 # arbitrary value. Maximum distance considered to consider two aa in contact and use the contact potential in the calculation. It is 5 by default

experimental = f"{main_dir}/PART1/S2648_renum_pdbnum.mut" # FILE PATH OF EXPERIMENTAL MUTATIONS

##############################################################################################
### make files with predicted ddg and experimental ddg

true_name1 = os.path.basename(modeller)
true_name2 = os.path.basename(rosetta)

for n in key:
    modi = f"{result}/{true_name1}_{cutoff}A_{n}.txt"
    rosi = f"{result}/{true_name2}_{cutoff}A_{n}.txt"

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
        sv.exp_predict(experimental, modi)
        sv.exp_predict(experimental, rosi)

### calculate pearson coefficient
iterate = os.listdir(result)
for f in iterate:
    path = result + "/" + f
    sv.pearson(path)