from WMD3Parser import *


#####################################################################
def create_tree(subgroup_lst):
    sons_dic, fathers_dic = createSubgroupsTree(subgroup_lst)
    for node in sons_dic:
        if fathers_dic[node] is None:
            del fathers_dic[node]
        if len(sons_dic[node]) == 0:
            singletons = [{sub} for sub in node]
            sons_dic[node] = [tuple(singleton) for singleton in singletons]
            for singleton in singletons:
                fathers_dic[tuple(sorted(list(singleton)))] = node
        elif len(sons_dic[node]) == 1:
            son_as_set = set(list(sons_dic[node][0]))
            node_as_set = set(list(node))
            son_to_add = node_as_set.difference(son_as_set)
            sons_dic[node].append(tuple(list(son_to_add)))
            fathers_dic[tuple(list(son_to_add))] = node
    return sons_dic, fathers_dic

##################################################################################################
def calculateDistBetweenGenes(sons_dic, father_dic, subgroups):
    dic_dist = {}
    genes = list(subgroups[-1])
    for i in range(len(genes)-1):
        gene1 = genes[i]
        for j in range(i+1, len(genes)):
            subset = {gene1}
            subset_tuple = tuple([gene1])
            gene2 = genes[j]
            distance = 0
            while not (gene2 in subset):
                subset_tuple = father_dic[subset_tuple]
                subset = set(list(subset_tuple))
                distance += 1
            while subset != {gene2}:
                son1 = sons_dic[subset_tuple][0]
                setSon1 = set(list(son1))
                son2 = sons_dic[subset_tuple][1]
                setSon2 = set(list(son2))
                if gene2 in setSon1:
                    subset = setSon1
                    subset_tuple = son1
                elif gene2 in setSon2:
                    subset = setSon2
                    subset_tuple = son2
                else:
                    raise Exception("Error in calculate distance!\n")
                distance += 1

            dic_dist[tuple(sorted([gene1, gene2]))] = distance-1
    return dic_dist
##############################################################################################################
def addFarSubgroupsOfGenes(dir_of_families, free_energy_thr, numPerInterNode, numOfMismatches, maxNumOfGeneGroups):
    families = os.listdir(dir_of_families)
    for family in families:
        full_family_dir = os.path.join(dir_of_families, family)
        dicGeneIsoform_path = os.path.join(full_family_dir, "geneIsoform.pkl")
        with open (dicGeneIsoform_path, 'rb') as handle_gene:
            dicGeneIsoform = pickle.load(handle_gene)
        if len(dicGeneIsoform) <= 2:
            continue
        if not os.path.isdir(full_family_dir):
            continue
        addNonAdjacentGenesAmiRNAs(full_family_dir, free_energy_thr, family, numPerInterNode,
                                   numOfMismatches, maxNumOfGeneGroups)
#############################################################################################################
def addNonAdjacentGenesAmiRNAs(full_family_dir, free_energy_thr, family, numPerInterNode,
                               numOfMismatches, maxNumOfGeneGroups):
    path_initial_miRs = os.path.join(full_family_dir, "nonAdjacentGeneSubgroups.pkl")
    path_MiRs_after_primary_remove = os.path.join(full_family_dir, "subgroups_gene_distance_remained.pkl")
    path_MiRs_after_energy_remove = os.path.join(full_family_dir, "second_round_energy_filtered.pkl")
    path_filtered_maxGroups = os.path.join(full_family_dir, "filteredMaxNumGroupsSecondRound.pkl")
    subgroups_path = os.path.join(full_family_dir, "subgroups.pkl")
    res_subgroups_path = os.path.join(full_family_dir, "subgroups_res.pkl")
    res_subgroups_path_final = os.path.join(full_family_dir, "subgroups_res_final.pkl")
    with open(subgroups_path, 'rb') as handle:
        subgroups = pickle.load(handle)
    largest_subgroup = str(len(subgroups)-1)

    with open(res_subgroups_path, 'rb') as handle:
        subgroups_res_prev_round =  pickle.load(handle)

    dirNewAmiRNAs = os.path.join(full_family_dir, largest_subgroup, "WMD3_1.OU")
    dic_amiRNAs = getNewAmiRNAs(dirNewAmiRNAs, free_energy_thr, family, subgroups)
    with open(path_initial_miRs, 'wb') as handle:
        pickle.dump(dic_amiRNAs, handle)
    res_amiRNAs = removeFarGenesGroups(dic_amiRNAs, subgroups, res_subgroups_path, numPerInterNode)
    with open(path_MiRs_after_primary_remove, 'wb') as handle:
        pickle.dump(res_amiRNAs, handle)
    removeAccordingEnergyThreshold(res_amiRNAs, dic_amiRNAs)
    with open (path_MiRs_after_energy_remove, 'wb') as handle:
        pickle.dump(res_amiRNAs, handle)


    root = tuple(sorted(list(subgroups[len(subgroups)-1])))   # it's a postorder -> the root is the last
    sonsSubgroups, fatherSubgroups = createSubgroupsTree(subgroups)
    subtrees = splitToSubTrees(sonsSubgroups, root, maxNumOfGeneGroups)
    remained_subgroups = getListOfRemainedInternalNodes(subtrees, sonsSubgroups)
    internal_nodes_to_del = set(res_amiRNAs.keys()).difference(set(remained_subgroups))
    for internal_node in internal_nodes_to_del:
        del res_amiRNAs[internal_node]

    # Store the results after using max number of internal nodes factor
    with open(path_filtered_maxGroups, 'wb') as handle:
        pickle.dump(res_amiRNAs, handle)
    updatedAmiRNAsPerInternalNode = {}
    for subtree in subtrees:
        newAmiRNAs = filterAccordingToInternalNodeAndSimilarity(sonsSubgroups, numPerInterNode, res_amiRNAs,
                                                   numOfMismatches,
                                                   subtree, subgroups_res_prev_round)
        updatedAmiRNAsPerInternalNode.update(newAmiRNAs)
    with open(res_subgroups_path_final, 'wb') as handle:
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
        df.to_csv(os.path.join(full_family_dir, "AmiRNAs_final_" +str(free_energy_thr)+ ".csv"))
    return


###################################################################################################################
def fillFinalDic(res_amiRNAs, res_subgroups_path):
    with open(res_subgroups_path, 'rb') as handle:
        subgroups_res = pickle.load(handle)
    for subgroup in subgroups_res:
        if len(subgroups_res[subgroup]) == 0:
            subgroups_res[subgroup] = res_amiRNAs[subgroup]
    return subgroups_res


################################################################################################################
def removeAccordingEnergyThreshold(res_amiRNAs, dic_amiRNAs):
    removed = set()
    for subgroup in res_amiRNAs:
        for amiRNA in res_amiRNAs[subgroup]:
            if len(res_amiRNAs[subgroup][amiRNA]["targets"]) <= 1:
                removed.add(amiRNA)
                del dic_amiRNAs[amiRNA]
    for amiRNA in removed:
        deleted = False
        for subgroup in res_amiRNAs:
            if amiRNA in res_amiRNAs[subgroup]:
                del res_amiRNAs[subgroup][amiRNA]
                deleted = True
                break
        if deleted:
            continue
    return removed
###################################################################################################################
def getGeneClusters(tree_subgroup_tuple, subgroup_amiRNA, sons_dic):
    lst_of_clusters = []
    initial_subgroup = subgroup_amiRNA
    son1 = sons_dic[tree_subgroup_tuple][0]
    son2 = sons_dic[tree_subgroup_tuple][1]

    getGeneClustersRec(lst_of_clusters, initial_subgroup, son1, son2, sons_dic)
    return lst_of_clusters
###################################################################################################################
def getGeneClustersRec(lst_of_clusters, initial_subgroup, son1, son2, sons_dic):
    notMonoSubgrpups = [set(), set()]
    for gene in initial_subgroup:
        if gene in son1:
            notMonoSubgrpups[0].add(gene)
        elif gene in son2:
            notMonoSubgrpups[1].add(gene)
        else:
            #print("gene", gene)
            #print("son1:", son1)
            #print("son2:", son2)
            raise Exception("ERROR!!! getGeneClustersRec(): gene is not found in any subgroup!")

    if len(initial_subgroup) == 0:
        return
    if notMonoSubgrpups[0] == set(son1):
        #print(notMonoSubgrpups[0], "son1", son1)
        lst_of_clusters.append(notMonoSubgrpups[0])

    else:
        if len(notMonoSubgrpups[0]) != 0:
            son1_for_first = sons_dic[son1][0]
            son2_for_first = sons_dic[son1][1]
            getGeneClustersRec(lst_of_clusters, notMonoSubgrpups[0], son1_for_first, son2_for_first, sons_dic)
    if notMonoSubgrpups[1] == set(son2):
        #print(notMonoSubgrpups[1], "son2", son2)
        lst_of_clusters.append(notMonoSubgrpups[1])
    else:
        if len(notMonoSubgrpups[1]) != 0:
            son1_for_second = sons_dic[son2][0]
            son2_for_second = sons_dic[son2][1]
            getGeneClustersRec(lst_of_clusters, notMonoSubgrpups[1], son1_for_second, son2_for_second, sons_dic)

###################################################################################################################
def areAllClustersNear(genes_clusters, dic_genes_distance):
    for i in range(len(genes_clusters)-1):
        for j in range(i+1, len(genes_clusters)):
            near = False
            cluster1 = genes_clusters[i]
            cluster2 = genes_clusters[j]
            #print(cluster1, cluster2)
            for gene_cluster1 in cluster1:
                for gene_cluster2 in cluster2:
                    #print("genes:", gene_cluster1, gene_cluster2, "distance", dic_genes_distance [tuple(sorted([gene_cluster1, gene_cluster2]))])
                    if dic_genes_distance [tuple(sorted([gene_cluster1, gene_cluster2]))] <= 3:
                        near = True
                        break
                if near:
                    break
            if not near:
                return False
    return True

###################################################################################################################
def removeFarGenesGroups(dic_amiRNAs, subgroups, res_subgroups_path, numPerInterNode):
    dic_subgroup_AmiRNA = {}
    removed = set()
    with open(res_subgroups_path, 'rb') as handle:
        res_subgroups = pickle.load(handle)
    tree_subgroups = [subgroup for subgroup in res_subgroups if len(res_subgroups[subgroup]) < numPerInterNode]
    for subgroup_in_tree in tree_subgroups:
        dic_subgroup_AmiRNA[subgroup_in_tree] = {}
    sons_dic, fathers_dic = create_tree(subgroups)
    dic_genes_distance = calculateDistBetweenGenes(sons_dic, fathers_dic, subgroups)

    subgroups_small_to_large = sorted (tree_subgroups, key= lambda x: len(x))
    subgroups_small_to_large = [set(tree_subgroup) for tree_subgroup in subgroups_small_to_large]
    for amiRNA in dic_amiRNAs:
        subgroup_amiRNA = getSubgroupOfAmiRNA(dic_amiRNAs, amiRNA)
        for tree_subgroup in subgroups_small_to_large:
            if subgroup_amiRNA.issubset(tree_subgroup):
                tree_subgroup_tuple = tuple(sorted(list(tree_subgroup)))
                genes_clusters = getGeneClusters(tree_subgroup_tuple, subgroup_amiRNA, sons_dic)
                near = areAllClustersNear(genes_clusters, dic_genes_distance)
                if not near:
                    removed.add(amiRNA)
                    break
                dic_subgroup_AmiRNA[tree_subgroup_tuple][amiRNA] = dic_amiRNAs[amiRNA]
                break
    remained_amiRNAs = set()
    for subgroup in dic_subgroup_AmiRNA:
        for amiRNA in dic_subgroup_AmiRNA[subgroup]:
            remained_amiRNAs.add(amiRNA)

    ## remove
    for amiRNA in dic_amiRNAs:
        if not (amiRNA in remained_amiRNAs):
            removed.add((amiRNA))
    for miR in removed:
        del dic_amiRNAs[miR]
    return dic_subgroup_AmiRNA
#####################################################################################################################
def getSubgroupOfAmiRNA(dic_amiRNAs, amiRNA):
    subgroup = [target[0] for target in dic_amiRNAs[amiRNA]["targets"]] + [target[0] for target in
                                                                                 dic_amiRNAs[amiRNA][
                                                                                     "belowThresholdTargets"]]
    subgroup = set(subgroup)
    return subgroup

#####################################################################################################################################
def getNewAmiRNAs(path, free_energy_thr, family, subgroups):
    set_of_tree_subgroups = set([tuple(sorted(list(subgroup))) for subgroup in subgroups])
    dictOfAmiRNAs = {}
    file = open (path, 'r')
    content = file.readlines()
    file.close()
    #print(path)
    for line in content:
        listOfInfo = line.split('\t')
        AmiRNA = listOfInfo[0]
        dictOfAmiRNAs[AmiRNA] = getAmiRNAAttributes(listOfInfo, free_energy_thr)
        dictOfAmiRNAs[AmiRNA]["family"] = family
        subgroup = [target[0] for target in dictOfAmiRNAs[AmiRNA]["targets"]] + [target[0] for target in
                                                                                 dictOfAmiRNAs[AmiRNA][
                                                                                     "belowThresholdTargets"]]
        subgroup = tuple(sorted(subgroup))
        if (subgroup in set_of_tree_subgroups) or (len(subgroup) == 1):
            del dictOfAmiRNAs[AmiRNA]

        #addAmiRNA(line, dictOfAmiRNAs, free_energy_thr, family, removedAmiRNAs, Initial_AmiRNAs, found_genes, dicAmiRNAWithCorrespondingSetOfGenes)
    return dictOfAmiRNAs

#######################################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input_file_path", "-i", help = "path for the WMD3 output file")
    parser.add_argument("--free_energy_thr", "-f", type = float, help = "the threshold for energey")
    parser.add_argument("--numPerInterNode", "-n", type = int, help = "how many AmiRNAs per internal node?")
    parser.add_argument("--maxNumOfGeneGroups", "-g", type=int, help="The max number of internal nodes to consider")
    parser.add_argument("--numOfMismatches", "-m", type= int, default= 1, help = "number of allowed mismatches")


    args = parser.parse_args()
    input_file_path = args.input_file_path
    free_energy = args.free_energy_thr
    numPerInterNode = args.numPerInterNode
    maxNumOfGeneGroups = args.maxNumOfGeneGroups
    numOfMismatches = args.numOfMismatches
    addFarSubgroupsOfGenes(input_file_path, free_energy, numPerInterNode, numOfMismatches, maxNumOfGeneGroups)











