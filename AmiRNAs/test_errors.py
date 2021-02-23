import os
import pickle
import argparse
def summarizeResults(dir_of_families, outputPath):
    dic_final = {}
    families = os.listdir(dir_of_families)
    for family in families:
        full_family_dir = os.path.join(dir_of_families, family)
        path_for_Initial_AmiRNAs = os.path.join(full_family_dir, "initial_AmiRNAs.pkl")
        if not os.path.exists(path_for_Initial_AmiRNAs):
            continue
        path_res_per_subgroups = os.path.join(full_family_dir, "subgroups_res.pkl")
        res_subgroups_path_final = os.path.join(full_family_dir, "subgroups_res_final.pkl")
        family_initial = getAmiRNAs(path_for_Initial_AmiRNAs)
        if not family_initial:
            continue

        if os.path.exists(res_subgroups_path_final):
            path_final = res_subgroups_path_final
        else:
            path_final = path_res_per_subgroups
        getFromSubgroupsOutput(path_final, dic_final, family)
    with open(outputPath, 'wb') as handle:
        pickle.dump(dic_final, handle)




def getAmiRNAs(path):
    setAmi = set()
    with open(path, 'rb') as handle:
        amiRNAs = pickle.load(handle)
    for amiRNA in amiRNAs:
        setAmi.add(amiRNA)
    return len(setAmi)


def getFromSubgroupsOutput(path, dic, family):
    with open(path, 'rb') as handle:
        subgroups_dic = pickle.load(handle)
    for subgroup in subgroups_dic:
        amiRNAs = subgroups_dic[subgroup]
        for amiRNA in amiRNAs:
            dic[amiRNA] = family
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="creates sequence files")
    parser.add_argument("--families_dir", "-i", help = "output directory")
    args = parser.parse_args()
    output = args.families_dir
    res_path = "/groups/itay_mayrose/anatshafir1/AmiRNA_lib/data/RESULTS_PLAZA_PARTITION/internal_nodes_6_subgroups_10/dict_final.pkl"
    summarizeResults(output, res_path)








