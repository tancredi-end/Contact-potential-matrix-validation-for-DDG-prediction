#import pyrosettacolabsetup; pyrosettacolabsetup.install_pyrosetta()
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq3
import Bio.SeqUtils
import Bio.Data.IUPACData as IUPACData
import pyrosetta; pyrosetta.init()
import pyrosetta
from pyrosetta import *
from pyrosetta import rosetta
from pyrosetta.rosetta.core.pose import Pose
from pyrosetta.rosetta.core.pose import PDBInfo
import logging
logging.basicConfig(level=logging.INFO)
import pandas
import seaborn
import matplotlib
import pyrosetta
import pyrosetta.distributed.io as io
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.io.pdb import dump_multimodel_pdb, add_to_multimodel_pdb
import pyrosetta.distributed.packed_pose as packed_pose
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
import pyrosetta.distributed.tasks.score as score
import pyrosetta.rosetta.core.pose as pose
import os,sys,platform
from pyrosetta.distributed.packed_pose.core import to_pose

# ESSENTIAL

main_dir = "/home/dede/linux/test"

pdbexp_path = f"{main_dir}/PART1/pdblist"

mut_path = f"{main_dir}/PART1/S2648_renum_pdbnum.mut"
packing_path = f"{main_dir}/PART2/pyrosetta_pdb/"
multi_save = f"{main_dir}/PART2/pyrosetta_pdb/PACK"
pack_mutated = f"{main_dir}/PART2/pyrosetta_pdb/PACK/try"


###############################################



if not os.path.exists(pack_mutated):
    os.makedirs(pack_mutated, exist_ok=True)
    print("FOLDER GENERATED ------------------------------------------------------------------------------------------")

file_list = os.listdir(multi_save)
file_list.remove("try")
for i in range(len(file_list)):
    filename = os.path.basename(file_list[i])
    name_parts = filename.split("_")
    file_list[i] = name_parts[0]


with open(mut_path, "r") as file:
    rosetta_format = []
    for line in file:
        listing = line.strip().split()
        pdb = listing[0][0:4]
        chained = listing[0][4]
        mut = listing[1]
        residue = ""
        aminoacids = ""
        icode = " "
        for i in mut:
            if i.isdigit():
                residue = residue + i
            else:
                aminoacids = aminoacids + i
        #if len(aminoacids) vedere se lungo 2 o 3
        if len(aminoacids) == 3:
            residue = int(residue)
            aa = aminoacids[2].upper()
            icode = aminoacids[1].upper()
        else:
            aa = aminoacids[1].upper()

        to_append = [pdb, chained, residue, aa, icode]
        rosetta_format.append(to_append)



input_protocol = """
<ROSETTASCRIPTS>
  <TASKOPERATIONS>
    <RestrictToRepacking name="only_pack"/>
  </TASKOPERATIONS>

  <MOVERS>
    <PackRotamersMover name="pack" task_operations="only_pack" />
  </MOVERS>

  <PROTOCOLS>
    <Add mover="pack"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>
"""

refine = """
<ROSETTASCRIPTS>

  <RESIDUE_SELECTORS>
    <ResiduePDBInfoHasLabel name="mutation" property="mutation" />
    <Not name="not_neighbor">
      <Neighborhood selector="mutation" distance="12.0" />
    </Not>
  </RESIDUE_SELECTORS>

  <TASKOPERATIONS>
    <RestrictToRepacking name="only_pack"/>
    <OperateOnResidueSubset name="only_repack_neighbors" selector="not_neighbor" >
      <PreventRepackingRLT/>
    </OperateOnResidueSubset>
  </TASKOPERATIONS>

  <MOVERS>
    <PackRotamersMover name="pack_area" task_operations="only_pack,only_repack_neighbors" />
  </MOVERS>

  <PROTOCOLS>
    <Add mover="pack_area"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>
    """

def get_number_of_models(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("pdb_structure", pdb_file)

    model_count = len(structure)

    return model_count

def retrieve_number_res(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb_structure", pdb_file)

    results = []

    for model in structure:
        residue_count = 0
        hetatm_count = 0
        for chain in model:
            for residue in chain:
                residue_count += 1
            for hetatm in chain.get_list():
                if hetatm.id[0] == "H":
                    hetatm_count += 1
        added = residue_count + hetatm_count
        results.append(added)

    return results

def mutate_residue(input_pose, res_index, new_aa, res_label=None):

    work_pose = packed_pose.to_pose(input_pose)

    # Annotate strucure with reslabel, for use in downstream protocol
    # Add parameters as score, for use in downstream analysis
    if res_label:
        work_pose.pdb_info().add_reslabel(res_index, res_label)
        pose.setPoseExtraScore(work_pose, "mutation_index", res_index)
        pose.setPoseExtraScore(work_pose, "mutation_aa", new_aa)

    if len(new_aa) == 1:
        new_aa = str.upper(Bio.SeqUtils.seq3(new_aa))
    assert new_aa in map(str.upper, IUPACData.protein_letters_3to1)

    protocol = """
<ROSETTASCRIPTS>
    <MOVERS>
        <MutateResidue name="mutate" new_res="%(new_aa)s" target="%(res_index)i" />
    </MOVERS>
    <PROTOCOLS>
        <Add mover_name="mutate"/>
    </PROTOCOLS>
</ROSETTASCRIPTS>
    """ % locals()

    return rosetta_scripts.SingleoutputRosettaScriptsTask(protocol)(work_pose)

# Initialize PyRosetta
pyrosetta.distributed.init

with open("rosetta_error.txt", "a") as file:
    for pdbfile, chain, residue, aminoacids, icode in rosetta_format:
        filepath = f"{packing_path}{pdbfile}_{chain}_{residue}{aminoacids.upper()}.pdb"
        check = get_number_of_models(f"{pdbexp_path}/{pdbfile}.pdb")

        if os.path.exists(filepath):
            continue
        elif not os.path.exists(filepath):
            try:
                print(">>> >>> >>> >>> >>> GENERATING NEW MUTATION:",pdbfile, residue, aminoacids)
                platform.python_version()
                input_relax = rosetta_scripts.SingleoutputRosettaScriptsTask(input_protocol)
                # Syntax check via setup
                input_relax.setup()

                if check > 1:
                    modelname = os.listdir(multi_save)
                    modelname.remove("try")
                    for t in range(check):
                        filename = f"{multi_save}/{pdbfile}_model_{t+1}.pdb"
                        if os.path.exists(filename):
                            basic = os.path.basename(filename)
                            get_parts = basic.split("_")
                            print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", get_parts)
                            filepath = f"{pack_mutated}/{pdbfile}_{chain}_{residue}{aminoacids.upper()}_model_{get_parts[2]}"
                        else: print("||||| ||||| ||||| ||||| ||||| ||||| SOMETHING WRONG ||||| ||||| ||||| ||||| ||||| |||||")

                        if os.path.exists(filepath):
                            continue
                        elif not os.path.exists(filepath):
                            try:
                                print(">>> >>> >>> >>> >>> GENERATING NEW MUTATION:", pdbfile, residue,
                                        aminoacids)
                                platform.python_version()
                                input_relax = rosetta_scripts.SingleoutputRosettaScriptsTask(input_protocol)
                                # Syntax check via setup
                                input_relax.setup()

                                print(
                                    ">>> MULTIMODEL PDB FILE - SINGLE MODEL RETRIEVED LOADING... ... ... ... ...")

                                print(
                                    f">>> >>> >>> >>> >>> PRINTING {filename} PRINTING <<< <<< <<< <<< <<<")
                                raw_input_pose = pose_from_pdb(filename)
                                score_task = score.ScorePoseTask()
                                scored_pose = score_task(raw_input_pose)
                                input_pose = input_relax(raw_input_pose)

                                # does not work with the notebook env python 3.7
                                pose_posi = raw_input_pose.pdb_info().pdb2pose(chain, int(residue), icode)
                                print(pose_posi)

                                mutation = mutate_residue(input_pose, pose_posi, aminoacids)

                                refine_mutation = rosetta_scripts.SingleoutputRosettaScriptsTask(refine)
                                output_packed_pose = refine_mutation(mutation)

                                # Convert packed pose back to regular pose
                                output_pose = to_pose(output_packed_pose)
                                if icode == " ":
                                    output_pose.dump_pdb(filepath)
                                else:
                                    output_pose.dump_pdb(f"{pack_mutated}/{pdbfile}_{chain}_{residue}{icode}{aminoacids.upper()}_model_{get_parts[2]}")
                            except RuntimeError as key:
                                error = f"> ERROR: {pdbfile} {get_parts[1]}_{get_parts[2]} {chain} {residue} {aminoacids} ||| {key} ||| {RuntimeError}\n"
                                file.write(error)
                        else:
                            print(f"{pdbfile} ||| {get_parts[0]} are not the same. Iterating to the next item.")
                            continue

                elif check == 1:
                    print(">>> SINGLE MODEL PDB FILE LOADING... ... ... ... ...")
                    raw_input_pose = pose_from_pdb(f"{pdbexp_path}/{pdbfile}.pdb")
                    score_task = score.ScorePoseTask()
                    scored_pose = score_task(raw_input_pose)
                    input_pose = input_relax(raw_input_pose)

                    #does not work with the notebook env python 3.7
                    pose_posi = raw_input_pose.pdb_info().pdb2pose(chain, int(residue), icode)
                    #print(pose_posi, icode)

                    mutation = mutate_residue(input_pose, pose_posi, aminoacids)

                    refine_mutation = rosetta_scripts.SingleoutputRosettaScriptsTask(refine)
                    output_packed_pose = refine_mutation(mutation)

                    # Convert packed pose back to regular pose
                    output_pose = to_pose(output_packed_pose)
                    if icode == " ":
                        output_pose.dump_pdb(f"{packing_path}{pdbfile}_{chain}_{residue}{aminoacids.upper()}.pdb")
                    else:
                        output_pose.dump_pdb(f"{packing_path}{pdbfile}_{chain}_{residue}{icode}{aminoacids.upper()}.pdb")
            except RuntimeError as key:
                error = f"> ERROR: {pdbfile} {chain} {residue} {aminoacids} ||| {key} ||| {RuntimeError}\n"
                file.write(error)