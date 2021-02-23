import os
import argparse
import shutil


def cleanPrevRunData(families_dir, free_energy, numPerInternalNode, maxGroups, delete):
    families = os.listdir(families_dir)
    for family in families:
        full_family_path = os.path.join(families_dir, family)
        results_dir = os.path.join(full_family_path, "Node_"+numPerInternalNode+"_MaxSubgroups_"+ maxGroups+
                                   "_energy_"+str(free_energy))
        if not delete:
            if not os.path.exists(results_dir):
                os.makedirs(results_dir)
        files = os.listdir(full_family_path)
        for file in files:
            if os.path.isdir(os.path.join(full_family_path, file)):
                continue
            if file.endswith(".pkl"):
                if file == "geneIsoform.pkl":
                    continue
                elif file == "subgroups.pkl":
                    continue
                else:
                    if delete:
                        os.remove(os.path.join(full_family_path, file))
                    else:
                        shutil.move(os.path.join(full_family_path, file), os.path.join(results_dir, file))
            elif file == "AmiRNAs"+ free_energy+".csv":
                if delete:
                    os.remove(os.path.join(full_family_path, file))
                else:
                    shutil.move(os.path.join(full_family_path, file), os.path.join(results_dir, file))
            elif file == "AmiRNAs_final_"+ free_energy+".csv":
                if delete:
                    os.remove(os.path.join(full_family_path, file))
                else:
                    shutil.move(os.path.join(full_family_path, file), os.path.join(results_dir, file))

    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input_dir", "-i", help = "families directory")
    parser.add_argument("--numInternalNodes", "-n", help = "number of AmiRNAs per internal node")
    parser.add_argument("--maxNumSubgroups", "-g",  help= "max number of subgroups")
    parser.add_argument("--free_energy", "-f", help = "free energy threshold")
    parser.add_argument("--delete", "-d", type= int, help = "delete or not")
    args = parser.parse_args()
    input_file_path = args.input_dir
    free_energy = args.free_energy
    numPerInternalNode = args.numInternalNodes
    maxGroups = args.maxNumSubgroups
    delete = args.delete

    cleanPrevRunData(input_file_path, free_energy, numPerInternalNode, maxGroups, delete)
