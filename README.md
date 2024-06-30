Here is the code for my thesis.
The objective of the work was a validation to check the ability of contact potential matrixes from AAindx3 to predict the DDG.

The Pearson linear correlation coefficient was used to compare the experimental values and the corrispective predicted one for mutations of 3 different datasets.

Both MODELLER and ROSETTA were used to generate structures of the mutated protein to the sum the contact potential of each residue with other aminoacids considering a selected cut-off (5 Angstrom is the default).

The "datasets" folder contains the id_name of PDB files, the mutation (considering the numbering convention of PDB possitions) and the experimental DDG value of the specific mutation.

The scripts are divided in the "scripts" folder for windows and linux, mainly for absolute path formatting difference.

The scripts named "PART" should automatically reproduce the passages to generate each result. 
REMEMBER: by generating new structures, it is expected that results will differ from the ones in the data folder.
