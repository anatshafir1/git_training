import os
import pickle
import argparse
def summarizeResults(dir_of_families):
    families = os.listdir(dir_of_families)
    numberOfInitialAmiRNAs = 0
    numberOfEnergyFilteredAmiRNAs = 0
    numberOfSimilarityInternalNodeFiltering = 0
    numberOfMaxGroupFiltering = 0
    numberOfFinal = 0
    numberOfFamilies = 0
    numberOfFamiliesForAdded = 0
    for family in families:

        full_family_dir = os.path.join(dir_of_families, family)
        path_for_Initial_AmiRNAs = os.path.join(full_family_dir, "initial_AmiRNAs.pkl")
        if not os.path.exists(path_for_Initial_AmiRNAs):
            continue
        numberOfFamilies += 1
        pathForAmiRNAsAfterEnergyFiltering = os.path.join(full_family_dir, "subgroups_energy_filtered.pkl")
        path_res_per_subgroups = os.path.join(full_family_dir, "subgroups_res.pkl")
        path_filtered_maxGroups = os.path.join(full_family_dir, "filteredMaxNumGroups.pkl")
        res_subgroups_path_final = os.path.join(full_family_dir, "subgroups_res_final.pkl")


        family_initial = getAmiRNAs(path_for_Initial_AmiRNAs)
        family_energyFilter = getFromSubgroupsOutput(pathForAmiRNAsAfterEnergyFiltering)
        family_maxGroupsFilter = getFromSubgroupsOutput(path_filtered_maxGroups)
        family_internalNodeFilter = getFromSubgroupsOutput(path_res_per_subgroups)
        if os.path.exists(res_subgroups_path_final):
            numberOfFamiliesForAdded += 1
            family_farSubgroups = getFromSubgroupsOutput(res_subgroups_path_final)

        numberOfInitialAmiRNAs += family_initial
        numberOfEnergyFilteredAmiRNAs += family_energyFilter
        numberOfMaxGroupFiltering += family_maxGroupsFilter
        numberOfSimilarityInternalNodeFiltering += family_internalNodeFilter
        if os.path.exists(res_subgroups_path_final):
            numberOfFinal += family_farSubgroups
        else:
            numberOfFinal += family_internalNodeFilter


    print("number of Initial: ", str(numberOfInitialAmiRNAs))
    print("AmiRNAs after energy filtering: ", str(numberOfEnergyFilteredAmiRNAs))
    print("AmiRNAs after maxGroup: ", str(numberOfMaxGroupFiltering))
    print("AmiRNAs after internal node filtering: ", str(numberOfSimilarityInternalNodeFiltering))
    print("AmiRNAs after adding far subgroups: ", str(numberOfFinal))

    print("*************************************\n")


def getAmiRNAs(path):
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
    return len(setAmi)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="creates sequence files")
    parser.add_argument("--families_dir", "-i", help = "output directory")
    args = parser.parse_args()
    output = args.families_dir
    summarizeResults(output)








