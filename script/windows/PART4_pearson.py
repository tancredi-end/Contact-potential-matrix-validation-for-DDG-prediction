import os
import ALL_FUNCTIONS as sv

##############################################################################################
#ESSENTIAL

main_dir = "E:\\ese_tesi"

pdb_no_h = f"{main_dir}\\PART3\\pymodeller_pdb"
modeller = f"{main_dir}\\PART4\\pymodeller_reverse"

#where to save
result = f"{main_dir}\\PART3\\end_potential_reverse" #where to save ddg results for each mutant

aasix3 = f"{main_dir}\\PART1\\aaindex3"
key = ["THOP960101","BETM990101", "MIYS990106", "BASU010101"]

true_name = os.path.basename(modeller)
primal_name = os.path.basename(pdb_no_h)

cutoff = 5

##############################################################################################
### calculate ddg
##############################################################################################


### make files with predicted ddg and experimental ddg

for n in key:
    experimental = f"{main_dir}PART3\\{primal_name}_{cutoff}A_{n}.mut"
    modi = f"{result}\\{true_name}_{cutoff}A_{n}.txt"

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
        sv.reverse_ddg(experimental, modi)

### calculate pearson coefficient
iterate = os.listdir(result)
for f in iterate:
    path = result + "\\" + f
    sv.pearson(path)