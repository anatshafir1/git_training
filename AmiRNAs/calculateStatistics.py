import os
import pickle
import argparse
import pandas as pd

def showStatistics(dir_of_families, numPerInterNode, results_dir):
    csv_file = os.path.join(results_dir, "statistics.csv")
    families = os.listdir(dir_of_families)
    dic = {"family":[], "numOfFoundAmiRNAs":[], "numOfGenes":[], "notFoundGenesInitialSearch":[],
           "NumOfNotFoundGenesAfterFiltering0.8":[], "NumOfNotFoundGenesAfterFiltering0.7":[],
           "NumOfSubgroupsWithLessThanThresholdAmiRNAs":[], "NumberOfEmptySubgroups":[],
           "InitialNumberOfSubgroups": [], "InitialNumberOfAmiRNAs": [], "NumberOfFilteredInternalNodeCriteria":[], "NumberOfFilteredMaxNumGroups":[],
           "NumberOfFilteredSimilarityCriterion":[], "NumberOfFilteredAmiRNAFreeEnergyCriteion":[]}
    for family in families:
        full_family_dir = os.path.join(dir_of_families, family)
        if not os.path.isdir(full_family_dir):
            continue
        dicGeneIsoform_path = os.path.join(full_family_dir, "geneIsoform.pkl")
        with open (dicGeneIsoform_path, 'rb') as handle_gene:
            dicGeneIsoform = pickle.load(handle_gene)
        if len(dicGeneIsoform) == 1:
            continue
        calculateStatistics(full_family_dir, family, numPerInterNode, dic)
    df = pd.DataFrame(dic)
    df.to_csv(csv_file)


def calculateStatistics(full_family_dir, family, numPerInterNode, dataFrameDic):
    print("****************************************")
    print(family)
    path_found_genes_path = os.path.join(full_family_dir, "found_genes_AmiRNAs.pkl")
    path_res_amiRNAs_path = os.path.join(full_family_dir, "subgroups_res.pkl")
    path_subgroups = os.path.join(full_family_dir, "subgroups.pkl")
    path_for_Initial_AmiRNAs = os.path.join(full_family_dir, "initial_AmiRNAs.pkl")
    path_for_found_genes = os.path.join(full_family_dir, "found_genes_AmiRNAs.pkl")
    pathForAmiRNAsAfterEnergyFiltering = os.path.join(full_family_dir, "subgroups_energy_filtered.pkl")
    pathForAmiRNAsAfterSimilarityFiltering = os.path.join(full_family_dir, "subgroups_similarity_removed.pkl")
    pathForAmiRNAsAdded = os.path.join(full_family_dir, "subgroups_added_from_filtered.pkl")
    path_filtered_maxGroups = os.path.join(full_family_dir, "filteredMaxNumGroups.pkl")
    path_filtered_maxPerInternal_node = os.path.join(full_family_dir, "filteredMaxPerInternalNode.pkl")

    with open(path_found_genes_path, 'rb') as handle:
        found_genes_initial = pickle.load(handle)
    with open(path_res_amiRNAs_path, 'rb') as handle:
        res_amiRNAs = pickle.load(handle)
    with open(path_subgroups, 'rb') as handle:
        subgroups_intial = pickle.load(handle)
    #print("Number of found AmiRNAs:", getNumberOfAmiRNAs(res_amiRNAs))
    dataFrameDic["family"].append(family)
    dataFrameDic["numOfFoundAmiRNAs"].append(getNumberOfAmiRNAs(res_amiRNAs))
    #print("Number of genes:", getNumberOfGenes(subgroups_intial))
    dataFrameDic["numOfGenes"].append(getNumberOfGenes(subgroups_intial))
    #print("Number of not found genes from initial search:", getNumberOfNotFoundGenes(found_genes_initial, subgroups_intial))
    dataFrameDic["notFoundGenesInitialSearch"].append(getNumberOfNotFoundGenes(found_genes_initial, subgroups_intial))
    aboveThr, aboveAndBelow = getNumberOfNotFoundGenesAfterFiltering(res_amiRNAs, subgroups_intial)
    #print("Number of not found genes after filtering > 0.8:", aboveThr)
    dataFrameDic["NumOfNotFoundGenesAfterFiltering0.8"].append(aboveThr)
    #print("Number of not found genes after filtering > 0.7:",aboveAndBelow)
    dataFrameDic["NumOfNotFoundGenesAfterFiltering0.7"].append(aboveAndBelow)
    #print("Number of subgroups with less than 4 amiRNAs:", getNumberOfBelowThrNodes(res_amiRNAs, numPerInterNode))
    dataFrameDic["NumOfSubgroupsWithLessThanThresholdAmiRNAs"].append(getNumberOfBelowThrNodes(res_amiRNAs, numPerInterNode))
    #print("Number of EMPTY subgroups:", getNumberOfEmptyGroups(res_amiRNAs))
    dataFrameDic["NumberOfEmptySubgroups"].append(getNumberOfEmptyGroups(res_amiRNAs))
    #print("Initial number of subgroups:", len(subgroups_intial))
    dataFrameDic["InitialNumberOfSubgroups"].append(len(subgroups_intial))
    numberOfInitialAmiRNAs = getInitialNumOfAmiRNAs(path_for_Initial_AmiRNAs)
    #print("Initial number of AmiRNAs is:", str(numberOfInitialAmiRNAs))
    dataFrameDic["InitialNumberOfAmiRNAs"].append(numberOfInitialAmiRNAs)
    filteredFreeEnergy = getNumberOfFiltered(pathForAmiRNAsAfterEnergyFiltering, pathForAmiRNAsAdded)
    filteredInternalNodes = getNumberOfFiltered(path_filtered_maxPerInternal_node, pathForAmiRNAsAdded, 1)
    filteredMaxGroups = getNumberOfFiltered(path_filtered_maxGroups, pathForAmiRNAsAdded, 1)
    filteredSimilarity = getNumberOfFiltered(pathForAmiRNAsAfterSimilarityFiltering, pathForAmiRNAsAdded)
    #print("Number of filtered internal node criteria:", str(filteredInternalNodes))
    dataFrameDic["NumberOfFilteredInternalNodeCriteria"].append(filteredInternalNodes)
    #print("Number of filtered max number of groups:", str(filteredMaxGroups))
    dataFrameDic["NumberOfFilteredMaxNumGroups"].append(filteredMaxGroups)
    #print("Number of filtered similarity criterion:", str(filteredSimilarity))
    dataFrameDic["NumberOfFilteredSimilarityCriterion"].append(filteredSimilarity)
    #print("Number of filtered AmiRNA free energy criteion:", str(filteredFreeEnergy))
    dataFrameDic["NumberOfFilteredAmiRNAFreeEnergyCriteion"].append(filteredFreeEnergy)
    


#############################################################################################
def getNumberOfFiltered(removedAmiRNAsFiltered_path, pathForAmiRNAsAdded, setOfMirs = 0):
    with open(removedAmiRNAsFiltered_path, 'rb') as handle:
        removedAmiRNAsFiltered = pickle.load(handle)
    with open(pathForAmiRNAsAdded, 'rb') as handle:
        setAdded = pickle.load(handle)
    removed = set()
    if setOfMirs:
        for amirna in removedAmiRNAsFiltered:
            if amirna in setAdded:
                continue
            else:
                removed.add(amirna)
    else:

        for subgroup in removedAmiRNAsFiltered:
            for amirna in removedAmiRNAsFiltered[subgroup]:
                if amirna in setAdded:
                    continue
                else:
                    removed.add(amirna)
    return len(removed)






############################################################################################
def getInitialNumOfAmiRNAs(path_for_Initial_AmiRNAs):
    with open(path_for_Initial_AmiRNAs, 'rb') as handle:
        amiRNAs = pickle.load(handle)
    return len(amiRNAs)

#############################################################################################
def getNumberOfAmiRNAs(res_amiRNAs):
    setAmiRNAs = set()
    for subgroup in res_amiRNAs:
        for amiRNA in res_amiRNAs[subgroup]:
            setAmiRNAs.add(amiRNA)
    return len(setAmiRNAs)
################################################################################################
def getNumberOfGenes(subgroups_intial):
    return len (subgroups_intial[-1])
################################################################################################
def getNumberOfNotFoundGenes(found_genes_initial, subgroups_intial):
    counter = 0
    allGenes = set(list(subgroups_intial[-1]))
    for gene in allGenes:
        if not gene in found_genes_initial:
            counter += 1
    return counter
###############################################################################################
def getNumberOfNotFoundGenesAfterFiltering(res_amiRNAs, subgroups_intial):
    allGenes = set(list(subgroups_intial[-1]))
    genesAboveThreshold = set()
    genesAboveAndBelow = set()
    for subgroup in res_amiRNAs:
        for amiRNA in res_amiRNAs[subgroup]:
            targetsAboveThreshold = set([target[0] for target in res_amiRNAs[subgroup][amiRNA]["targets"]])
            genesAboveThreshold = genesAboveThreshold.union(targetsAboveThreshold)
            genesAboveAndBelow = genesAboveAndBelow.union(targetsAboveThreshold)
            belowThr = set([target[0] for target in res_amiRNAs[subgroup][amiRNA]["belowThresholdTargets"]])
            genesAboveAndBelow = genesAboveAndBelow.union(belowThr)
    countAbove = 0
    countBelowAndAbove = 0
    for gene in allGenes:
        if not gene in genesAboveThreshold:
            countAbove += 1
        if not gene in genesAboveAndBelow:
            countBelowAndAbove += 1
    return countAbove, countBelowAndAbove
###############################################################################################
def getNumberOfBelowThrNodes(res_amiRNAs, numPerInterNode):
    count = 0
    for subgroup in res_amiRNAs:
        if len(res_amiRNAs[subgroup]) < numPerInterNode:
            count += 1
    return count
##############################################################################################
def getNumberOfEmptyGroups(res_amiRNAs):
    count = 0
    for subgroup in res_amiRNAs:
        if len(res_amiRNAs[subgroup]) == 0:
            count += 1
    return count
##############################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input_file_path", "-i", help = "path for the WMD3 output file")
    parser.add_argument("--numPerInterNode", "-n", type = int, help = "how many AmiRNAs per internal node?")
    parser.add_argument("--results_dir", "-r", help = "results directory")
    args = parser.parse_args()
    input_file_path = args.input_file_path
    numPerInterNode = args.numPerInterNode
    results_dir = args.results_dir
    showStatistics(input_file_path, numPerInterNode, results_dir)

