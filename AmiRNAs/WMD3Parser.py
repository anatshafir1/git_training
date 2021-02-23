import os
import argparse
import pickle
import pandas as pd
########################################################################################################################
def parseWMD3Results(dir_of_families, free_energy_thr, numPerInterNode, maxNumOfGeneGroups, numOfAllowedMismatches):
    families = os.listdir(dir_of_families)
    for family in families:
        full_family_dir = os.path.join(dir_of_families, family)
        if not os.path.isdir(full_family_dir):
            continue
        uniteAmiRNAsOfSubgroups(full_family_dir, free_energy_thr, family, numPerInterNode, maxNumOfGeneGroups,
                                numOfAllowedMismatches)
########################################################################################################################
def uniteAmiRNAsOfSubgroups(familyPath, free_energy_thr, family, numPerInterNode, maxNumOfGeneGroups,
                            numOfAllowedMismatches):
    subgroup_pickle_path = os.path.join(familyPath, "subgroups.pkl")
    path_for_Initial_AmiRNAs = os.path.join(familyPath, "initial_AmiRNAs.pkl")
    path_for_found_genes = os.path.join(familyPath, "found_genes_AmiRNAs.pkl")
    dicGeneIsoform_path = os.path.join(familyPath, "geneIsoform.pkl")
    pathForAmiRNAsAfterEnergyFiltering = os.path.join(familyPath, "subgroups_energy_filtered.pkl")
    path_res_per_subgroups = os.path.join(familyPath, "subgroups_res.pkl")
    path_filtered_maxGroups = os.path.join(familyPath, "filteredMaxNumGroups.pkl")
    path_filtered_maxPerInternal_node = os.path.join(familyPath, "filteredMaxPerInternalNode.pkl")
    Initial_AmiRNAs = set()
    found_genes = set()

    with open(dicGeneIsoform_path, 'rb') as handleIsoform:
        geneIsoDic = pickle.load(handleIsoform)
    if len(geneIsoDic) < 2: # no need for families with a single gene
        print(family)
        return
    elif len(geneIsoDic) == 2:  # create dummy subgroups of genes for families with two genes that do not have a tree
        if not os.path.exists(subgroup_pickle_path):
            createDummySubgroups(geneIsoDic, subgroup_pickle_path)


    with open(subgroup_pickle_path, 'rb') as handle:
        lst_of_subgroups = pickle.load(handle)

    allSubgroupsAmiRNAs = {}
    numOfSubgroups = len(lst_of_subgroups)
    removeLowEnergy(allSubgroupsAmiRNAs, numOfSubgroups, familyPath, free_energy_thr,
                    family, lst_of_subgroups, Initial_AmiRNAs, found_genes)


    root = tuple(sorted(list(lst_of_subgroups[len(lst_of_subgroups)-1])))   # it's a postorder -> the root is the last
    sonsSubgroups, fatherSubgroups = createSubgroupsTree(lst_of_subgroups)
    subtrees = splitToSubTrees(sonsSubgroups, root, maxNumOfGeneGroups)

    # remove duplicates
    removeDuplicates(allSubgroupsAmiRNAs, subtrees, sonsSubgroups, fatherSubgroups)

    with open(pathForAmiRNAsAfterEnergyFiltering, 'wb') as handle:
        pickle.dump(allSubgroupsAmiRNAs, handle)
    # delete unneeded subgroups
    remained_subgroups = getListOfRemainedInternalNodes(subtrees, sonsSubgroups)
    internal_nodes_to_del = set(allSubgroupsAmiRNAs.keys()).difference(set(remained_subgroups))
    print("#############################################################################\n")
    print(family)
    for internal_node in internal_nodes_to_del:
        print("\tInternal node of length:", str(len(internal_node)))
        print("\tremoved AmiRNAs:")
        amiRNAs_to_del = allSubgroupsAmiRNAs[internal_node]
        for mi in amiRNAs_to_del:
            print("\t"+mi)
        print("\tThe number of deleted AmiRNAs:", str(len(amiRNAs_to_del)))
        del allSubgroupsAmiRNAs[internal_node]

    # Store the results after using max number of internal nodes factor
    with open(path_filtered_maxGroups, 'wb') as handle:
        pickle.dump(allSubgroupsAmiRNAs, handle)
    updatedAmiRNAsPerInternalNode = {}
    for subtree in subtrees:
        newAmiRNAs = filterAccordingToInternalNodeAndSimilarity(sonsSubgroups, numPerInterNode, allSubgroupsAmiRNAs,
                                                   numOfAllowedMismatches,
                                                   subtree)
        updatedAmiRNAsPerInternalNode.update(newAmiRNAs)

    #print("initial number:", str(len(Initial_AmiRNAs)))
    with open(path_for_Initial_AmiRNAs, 'wb') as handle:
        pickle.dump(Initial_AmiRNAs, handle)
    with open(path_for_found_genes, 'wb') as handle:
        pickle.dump(found_genes, handle)
    with open(path_filtered_maxPerInternal_node, 'wb') as handle:
        pickle.dump(updatedAmiRNAsPerInternalNode, handle)
    with open(path_res_per_subgroups, 'wb') as handle:
        pickle.dump(updatedAmiRNAsPerInternalNode, handle)
    amiRNAs_dict = {}
    for interNode in updatedAmiRNAsPerInternalNode:
        amiRNAs_per_internalNode = updatedAmiRNAsPerInternalNode[interNode]
        for amiRNA in amiRNAs_per_internalNode:
            amiRNAs_dict[amiRNA] = amiRNAs_per_internalNode[amiRNA]
            amiRNAs_dict[amiRNA]['internal_node'] = interNode

    listOfAmiRNAs = sorted(list(amiRNAs_dict), key=lambda x: (-len(amiRNAs_dict[x]['internal_node']),
                                                              amiRNAs_dict[x]['internal_node'],
                                                              -len(amiRNAs_dict[x]["targets"]),
                                                              amiRNAs_dict[x]["score"]))
    if listOfAmiRNAs == []:
        return
    else:
        df = pd.DataFrame(columns= list(amiRNAs_dict[listOfAmiRNAs[0]].keys()))
        for i in listOfAmiRNAs:
            df.loc[i] = [amiRNAs_dict[i][df.columns[k]] for k in range(len(df.columns))]
        df.to_csv(os.path.join(familyPath, "AmiRNAs" +str(free_energy_thr)+ ".csv"))
    return
########################################################################################################################
def createDummySubgroups(geneIsoDic, subgroup_pickle_path):
    subgroups_to_create = []
    setOfGenes = set()
    for geneIso in geneIsoDic:
        setOfGenes.add(geneIso)
    subgroups_to_create.append(setOfGenes)
    with open(subgroup_pickle_path, 'wb') as handle_subgr:
        pickle.dump(subgroups_to_create, handle_subgr)
    return
########################################################################################################################
def removeLowEnergy(allSubgroupsAmiRNAs, numOfSubgroups, familyPath, free_energy_thr,
                    family, family_dic_of_subgroups, Initial_AmiRNAs, found_genes,
                    dicAmiRNAWithCorrespondingSetOfGenes = None):
    amiRNAs = {}
    removedAmiRNAs = {}
    for i in range(numOfSubgroups):
        #print(i)
        subgroup_path = os.path.join(familyPath, str(i))
        dictAmiRNASubgroup = parseWMD3File(os.path.join(subgroup_path, "WMD3.OU"),
                                           free_energy_thr,
                                           family,
                                           removedAmiRNAs,
                                           Initial_AmiRNAs,
                                           found_genes,
                                           dicAmiRNAWithCorrespondingSetOfGenes)
        allSubgroupsAmiRNAs[tuple(sorted(list(family_dic_of_subgroups[i])))] = dictAmiRNASubgroup

        amiRNAs.update(allSubgroupsAmiRNAs[tuple(sorted(list(family_dic_of_subgroups[i])))])
    return removedAmiRNAs
########################################################################################################################
def parseWMD3File(path, free_energy_thr, family, removedAmiRNAs, Initial_AmiRNAs, found_genes,
                  dicAmiRNAWithCorrespondingSetOfGenes = None):
    dictOfAmiRNAs = {}
    file = open (path, 'r')
    content = file.readlines()
    file.close()
    #print(path)
    count = 0
    for line in content:
        count += 1
        listOfInfo = line.split('\t')
        AmiRNA = listOfInfo[0]
        Initial_AmiRNAs.add(AmiRNA)
        count += 1
        dictOfAmiRNAs[AmiRNA] = getAmiRNAAttributes(listOfInfo, free_energy_thr)
        dictOfAmiRNAs[AmiRNA]["family"] = family
        subgroup = [target[0] for target in dictOfAmiRNAs[AmiRNA]["targets"]] + [target[0] for target in
                                                                                 dictOfAmiRNAs[AmiRNA][
                                                                                     "belowThresholdTargets"]]
        subgroup = tuple(sorted(subgroup))

        for gene in subgroup:
            found_genes.add(gene)
        if len(dictOfAmiRNAs[AmiRNA]["targets"]) <= 1:
            if not subgroup in removedAmiRNAs:
                removedAmiRNAs[subgroup] = {}
            removedAmiRNAs[subgroup][AmiRNA] = dictOfAmiRNAs[AmiRNA]
            del dictOfAmiRNAs[AmiRNA]
        else:
            if not dicAmiRNAWithCorrespondingSetOfGenes is None:
                dicAmiRNAWithCorrespondingSetOfGenes[AmiRNA] = subgroup
        #addAmiRNA(line, dictOfAmiRNAs, free_energy_thr, family, removedAmiRNAs, Initial_AmiRNAs, found_genes,
        # dicAmiRNAWithCorrespondingSetOfGenes)
    #print(count)
    return dictOfAmiRNAs
########################################################################################################################
def getAmiRNAAttributes(listOfInfo, free_energy_thr):
    #check it: is 'AT is at the beginning of all the genes??
    dictOfAttributes = {}
    dictOfAttributes["score"] = float(listOfInfo[-1])
    dictOfAttributes["perfectMatch"] = float(listOfInfo[6])
    dictOfAttributes["targets"] = []
    dictOfAttributes["belowThresholdTargets"]= []
    dictOfAttributes["color"] = getColorOfScore(dictOfAttributes["score"])
    for i in range(7, len(listOfInfo)-2, 2):
        free_energy = float(listOfInfo[i+1])
        free_energy_frac = free_energy/dictOfAttributes["perfectMatch"]
        if free_energy_frac < free_energy_thr:
            dictOfAttributes["belowThresholdTargets"].append(tuple((listOfInfo[i], float(listOfInfo[i+1]))))
        else:
            dictOfAttributes["targets"].append(tuple((listOfInfo[i], float(listOfInfo[i+1]))))
    return dictOfAttributes
########################################################################################################################
def getColorOfScore(score):
    if score <= 15:
        return "green"
    elif score <= 20:
        return "yellow"
    elif score <= 25:
        return "orange"
    return "red"
########################################################################################################################
def createSubgroupsTree(lstOfSubgroups):
    fatherSubgroups = {}
    subgroupsTree = {}
    setOfPossibleSons = []
    for i in range(len(lstOfSubgroups)):
        currNode = lstOfSubgroups[i]
        subgroupsTree[tuple(sorted(list(currNode)))] = []
        sons = []
        for sonCandidate in setOfPossibleSons:
            if sonCandidate.issubset(currNode):
                subgroupsTree[tuple(sorted(list(currNode)))].append(tuple(sorted(list(sonCandidate))))
                fatherSubgroups[tuple(sorted(list(sonCandidate)))] = tuple(sorted(list(currNode)))
                sons.append(sonCandidate)
        for son in sons:
            setOfPossibleSons.remove(son)
        setOfPossibleSons.append(currNode)
        if i == len(lstOfSubgroups)-1:
            fatherSubgroups[tuple(sorted(list(currNode)))] = None
    return subgroupsTree, fatherSubgroups
########################################################################################################################
def splitToSubTrees(sonsSubgroups, root, maxNumInternalNodes):
    subtrees = set()
    #root = getRoot(fatherSubgroups)
    dict_num_internal_nodes = {}
    createDictNumInternalNodes(root, sonsSubgroups, dict_num_internal_nodes, subtrees, maxNumInternalNodes)
    #print(dict_num_internal_nodes)
    #for subtree in subtrees:
        #fatherSubgroups[subtree] = None
    return subtrees
########################################################################################################################
def createDictNumInternalNodes(root, sonsSubgroups, dict_num_internal_nodes, subtrees, maxNumInternalNodes):
    sons = sonsSubgroups[root]
    if sons == []:
        dict_num_internal_nodes[root] = 1
        subtrees.add(root)
        return
    for son in sons:
        createDictNumInternalNodes(son, sonsSubgroups, dict_num_internal_nodes, subtrees, maxNumInternalNodes)
        if not root in dict_num_internal_nodes:
            dict_num_internal_nodes[root] = dict_num_internal_nodes[son] + 1
        else:
            dict_num_internal_nodes[root] += dict_num_internal_nodes[son]
    if dict_num_internal_nodes[root] <= maxNumInternalNodes:
        for son in sons:
            subtrees.remove(son)
        subtrees.add(root)
########################################################################################################################
def removeDuplicates(allSubgroupsAmiRNAs, subtrees, sonsSubgroups, fatherSubgroups):
    setOfUsedAmiRNAsInDifferentSubtrees = {}
    setOfAmiRNAs = set()
    for subtree in subtrees:
        dict_amiRNA_subgroup = {}
        setOfUsedAmiRNAs = set(allSubgroupsAmiRNAs[subtree].keys())
        removeDuplicatesRec(allSubgroupsAmiRNAs, subtree, sonsSubgroups, setOfUsedAmiRNAs, dict_amiRNA_subgroup)
        setOfUsedAmiRNAsInDifferentSubtrees[subtree] = dict_amiRNA_subgroup
        setOfAmiRNAs.update(set(dict_amiRNA_subgroup))
    for amiRNA in setOfAmiRNAs:
        lstScoreSubtree = []
        for subtree in subtrees:
            if amiRNA in setOfUsedAmiRNAsInDifferentSubtrees[subtree]:
                subgroupOfMiR = setOfUsedAmiRNAsInDifferentSubtrees[subtree][amiRNA]
                score = allSubgroupsAmiRNAs[subgroupOfMiR][amiRNA]['score']
                lstScoreSubtree.append((score, subtree, subgroupOfMiR))
        lstScoreSubtree.sort(key = lambda x: x[0])
        if len(lstScoreSubtree) > 1:
            for i in range(1, len(lstScoreSubtree)):
                del allSubgroupsAmiRNAs[lstScoreSubtree[i][2]][amiRNA]
    deleteBackwards(subtrees, allSubgroupsAmiRNAs, fatherSubgroups, setOfAmiRNAs)
    for subtree in subtrees:
        fatherSubgroups[subtree] = None
########################################################################################################################
def deleteBackwards(subtrees, allSubgroupsAmiRNAs, fatherSubgroups, setOfAmiRNAs):
    for subtree in subtrees:
        currNode = fatherSubgroups[subtree]
        while not currNode is None:
            for amiRNA in setOfAmiRNAs:
                if amiRNA in allSubgroupsAmiRNAs[currNode]:
                    del allSubgroupsAmiRNAs[currNode][amiRNA]
            currNode = fatherSubgroups[currNode]

########################################################################################################################
def removeDuplicatesRec(allSubgroupsAmiRNAs, subtree, sonsSubgroups, setOfUsedAmiRNAs, dict_amiRNA_subgroup):
    sons = sonsSubgroups[subtree]
    amiRNAs = allSubgroupsAmiRNAs[subtree]
    for miR in amiRNAs:
        dict_amiRNA_subgroup[miR] = subtree
    if sons == []:
        return

    setOfUsedAmiRNAs = setOfUsedAmiRNAs.union(set(amiRNAs.keys()))
    for son in sons:
        amiRNAs_sons = allSubgroupsAmiRNAs[son]
        for amiRNA in setOfUsedAmiRNAs:
            if amiRNA in amiRNAs_sons:
                del allSubgroupsAmiRNAs[son][amiRNA]
    for son in sons:
        removeDuplicatesRec(allSubgroupsAmiRNAs, son, sonsSubgroups, setOfUsedAmiRNAs, dict_amiRNA_subgroup)
########################################################################################################################
def getListOfRemainedInternalNodes(subtrees, sonsSubgroups):
    internalNodes = []
    for subtree in subtrees:
        addInternalNodesRec(subtree, sonsSubgroups, internalNodes)
    return internalNodes
########################################################################################################################
def addInternalNodesRec(subtree, sonsSubgroups, internalNodes):
    internalNodes.append(subtree)
    sons = sonsSubgroups[subtree]
    if sons == []:
        return
    for son in sons:
        addInternalNodesRec(son, sonsSubgroups, internalNodes)
########################################################################################################################
def filterAccordingToInternalNodeAndSimilarity(sonsSubgroups, numPerInterNode, allSubgroupsAmiRNAs,
                                               numOfAllowedMismatches,
                                               rootSubgroup, existedResultsSubgroups = None):
    if existedResultsSubgroups is None:
        updatedAmiRNAs = {}
    else:
        updatedAmiRNAs = existedResultsSubgroups
    addedAmiRNAsForAncestors = {rootSubgroup:{}}
    fillTreeWithAmiRNAs(rootSubgroup, sonsSubgroups, allSubgroupsAmiRNAs, updatedAmiRNAs, addedAmiRNAsForAncestors,
                        numOfAllowedMismatches, numPerInterNode, existedResultsSubgroups)
    return updatedAmiRNAs
########################################################################################################################
def fillTreeWithAmiRNAs(internalNode, sonsSubgroups, allSubgroupsAmiRNAs, updatedAmiRNAs, addedAmiRNAsForAncestors,
                        numOfAllowedMismatches, numPerInterNode, existedResultsSubgroups):
    #print("internal node:", internalNode)
    #print("subgroups:", addedAmiRNAsForAncestors.keys())
    #for key in addedAmiRNAsForAncestors:
        #print("\t", addedAmiRNAsForAncestors[key].keys())
    #print("###################################################")
    if internalNode in allSubgroupsAmiRNAs:
        amiRNAs = allSubgroupsAmiRNAs[internalNode]  # to be tested
        amiRNAs_candidates = sorted(list(amiRNAs), key=lambda x: (-len(amiRNAs[x]["targets"]), amiRNAs[x]["score"]))
        alreadyAdded = addedAmiRNAsForAncestors[internalNode]
        if existedResultsSubgroups:
            alreadyAdded.update(existedResultsSubgroups[internalNode])
        if not internalNode in updatedAmiRNAs:
            updatedAmiRNAs[internalNode] = {}  # AmiRNAs added for current node
        for candidate_amiRNA in amiRNAs_candidates:
            if len(updatedAmiRNAs[internalNode]) >= numPerInterNode:
                break
            filtered = False
            for existedAmiRNA in alreadyAdded:
                if isSimilar(candidate_amiRNA, existedAmiRNA, numOfAllowedMismatches):
                    if amiRNAs[candidate_amiRNA]['score'] >= alreadyAdded[existedAmiRNA]['score']:
                        filtered = True
                        break
            if (not filtered) and (len(updatedAmiRNAs[internalNode]) < numPerInterNode):
                # add AmiRNA
                updatedAmiRNAs[internalNode][candidate_amiRNA] = allSubgroupsAmiRNAs[internalNode][candidate_amiRNA]
                alreadyAdded[candidate_amiRNA] = allSubgroupsAmiRNAs[internalNode][candidate_amiRNA]
    else:
        #if existedResultsSubgroups:
        alreadyAdded = addedAmiRNAsForAncestors[internalNode]
        alreadyAdded.update(existedResultsSubgroups[internalNode])
        if existedResultsSubgroups is None:
            raise Exception("fillTreeWithAmiRNAs(): Key error!")

    sons = sonsSubgroups[internalNode]
    for son in sons:
        addedAmiRNAsForAncestors[son] = alreadyAdded

    if sonsSubgroups[internalNode] == []:
        return
    else:
        for i in range(len(sonsSubgroups[internalNode])):
            fillTreeWithAmiRNAs(sonsSubgroups[internalNode][i], sonsSubgroups, allSubgroupsAmiRNAs, updatedAmiRNAs,
                                addedAmiRNAsForAncestors, numOfAllowedMismatches, numPerInterNode, existedResultsSubgroups)

        del addedAmiRNAsForAncestors[internalNode]
########################################################################################################################
def isSimilar(amiRNA1, amiRNA2, MAX_DIVERGENCE):
    divergence = 0
    for i in range(len(amiRNA1)):
        if amiRNA1[i] != amiRNA2[i]:
            divergence += 1
    if divergence <= MAX_DIVERGENCE:
        return True
    return False
########################################################################################################################
if __name__ == "__main__":
    # full_family_dir = "D:\\ItayMNB9\\Documents\\amiRNA_lib\\families\\ABC_A"
    # free_energy_thr = 0.8
    # numPerInterNode = 3
    # maxNumOfGeneGroups = 4
    # numOfAllowedMismatches = 1

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input_file_path", "-i", help = "path for the WMD3 output file")
    parser.add_argument("--free_energy_thr", "-f", type = float, help = "the threshold for energey")
    parser.add_argument("--numPerInterNode", "-n", type = int, help = "how many AmiRNAs per internal node?")
    parser.add_argument("--maxNumOfGeneGroups", "-g", type = int, help = "The max number of internal nodes to consider")
    parser.add_argument("--numOfMismatches", "-m", type= int, default= 1, help = "number of allowed mismatches")
    args = parser.parse_args()
    input_file_path = args.input_file_path
    free_energy = args.free_energy_thr
    numPerInterNode = args.numPerInterNode
    maxNumOfGeneGroups = args.maxNumOfGeneGroups
    numOfAllowedMismatches = args.numOfMismatches
    parseWMD3Results(input_file_path, free_energy, numPerInterNode, maxNumOfGeneGroups, numOfAllowedMismatches)
