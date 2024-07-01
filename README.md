Here is the code for my thesis.
The objective of the work was a validation to check the ability of contact potential matrixes from AAindx3 to predict the DDG.

The Pearson linear correlation coefficient was used to compare the experimental values and the corrispective predicted one for mutations of 3 different datasets.

Both MODELLER and ROSETTA were used to generate structures of the mutated protein to the sum the contact potential of each residue with other aminoacids considering a selected cut-off (5 Angstrom is the default).

The "datasets" folder contains the id_name of PDB files, the mutation/s (considering the numbering convention of PDB positions) and the experimental DDG value of the specific mutation/s.

The scripts are divided in the "scripts" folder for windows and linux, mainly for absolute path formatting difference.

The scripts named "PART" should automatically reproduce the passages to generate each result. 
More information on the correct order and on what each script does into the GUIDE.txt file.

The main scripts for generic use (outside the 3 used dataset if needed) are:

Calculate-distance.py

calculate-single-ddg.py

calculate-all-ddg.py

PART3_modeller.py

PART2_pyrosetta.py

ALL_FUNCTIONS.py contains all function used in the different scripts with the exception of "calculate-single-ddg.py" which is not a function and it's used to calculate the ddg of a single mutation and not the entire folder.

REMEMBER: by generating new structures, it is expected that results will differ from the ones in the "data" folder.

NOTE: due to the size of the files with arrays of minimal distances between each residue, it was not possible to load them into this repository. You'll need to remake them using the calculate-distance.py script
  Beware that it is a time-consuming activity.
