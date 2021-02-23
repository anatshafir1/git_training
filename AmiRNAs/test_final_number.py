import os
import pickle
import argparse


def summarizeResults(dir_of_families, second_families_dir):
    families_first = os.listdir(dir_of_families)
    families_second = os.listdir(second_families_dir)

    amiRNAs_f = getAmiRNAs(families_first, dir_of_families)
    amiRNAs_s = getAmiRNAs(families_second, second_families_dir)
    union = amiRNAs_f.union(amiRNAs_s)
    print("The final number: ", str(len(union)))

def getAmiRNAs(families, dir_of_families):
    setOfAmiRNAs = set()

    for family in families:
        full_family_dir = os.path.join(dir_of_families, family)
        path_res_per_subgroups = os.path.join(full_family_dir, "subgroups_res.pkl")
        res_subgroups_path_final = os.path.join(full_family_dir, "subgroups_res_final.pkl")
        path_for_Initial_AmiRNAs = os.path.join(full_family_dir, "initial_AmiRNAs.pkl")
        if not os.path.exists(path_for_Initial_AmiRNAs):
            continue

        family_initial = getInitial(path_for_Initial_AmiRNAs)
        if not family_initial:
            print(family)
            continue
        if not os.path.exists(res_subgroups_path_final):
            res_path = path_res_per_subgroups
        else:
            res_path = res_subgroups_path_final
        family_set_of_Amis = getFromSubgroupsOutput(res_path)
        setOfAmiRNAs = setOfAmiRNAs.union(family_set_of_Amis)
    return setOfAmiRNAs


def getInitial(path):
    setAmi = set()
    with open(path, 'rb') as handle:
        amiRNAs = pickle.load(handle)
    for amiRNA in amiRNAs:
        setAmi.add(amiRNA)
    return len(setAmi)


def getFromSubgroupsOutput(path):
    setAmi = set()
    with open(path, 'rb') as handle:
        subgroups_dic = pickle.load(handle)
    for subgroup in subgroups_dic:
        amiRNAs = subgroups_dic[subgroup]
        for amiRNA in amiRNAs:
            setAmi.add(amiRNA)
    return setAmi

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="creates sequence files")
    parser.add_argument("--families_dir", "-i", help = "output directory")
    parser.add_argument("--families_second_dir", "-o", help = "second output directory")
    args = parser.parse_args()
    families_dir = args.families_dir
    second_families_dir = args.families_second_dir
    summarizeResults(families_dir, second_families_dir)
