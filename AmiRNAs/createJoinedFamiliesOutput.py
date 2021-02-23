import pandas as pd
import pickle
import re
from WMD3Parser import isSimilar
import os
import argparse

def parseResultsFilesOfDiffPartitions(path_first_partition, path_plaza_partition,
                                      numberPerInternalNode, numberOfMismatches, final_output_path):
    new_dic = {}
    df_first = pd.read_csv(path_first_partition, index_col=0)
    df_plaza = pd.read_csv(path_plaza_partition, index_col=0)
    dict_results_first, dic_family_node_first = getDictOfResults(df_first)
    dict_results_plaza, dic_family_node_plaza = getDictOfResults(df_plaza)
    constructNewAmiRNAsDict(dict_results_first, dict_results_plaza, new_dic)
    finalResults = {}
    filterFinalOutputAccordingToSimilarity(new_dic, numberPerInternalNode, numberOfMismatches, finalResults)

    final_results_for_output = {}
    for internalNode in finalResults:
        amiRNAs = finalResults[internalNode]
        for amiRNA in amiRNAs:
            final_results_for_output[amiRNA] = {}
            final_results_for_output[amiRNA]["family"] = finalResults[internalNode][amiRNA][0]
            final_results_for_output[amiRNA]["internal_node"] = internalNode
            final_results_for_output[amiRNA]["score"] = finalResults[internalNode][amiRNA][1]
            final_results_for_output[amiRNA]['perfectMatch'] = finalResults[internalNode][amiRNA][2]
            final_results_for_output[amiRNA]["targets"] = finalResults[internalNode][amiRNA][3]
            final_results_for_output[amiRNA]['belowThresholdTargets'] = finalResults[internalNode][amiRNA][4]
            final_results_for_output[amiRNA]['color'] = finalResults[internalNode][amiRNA][5]

    listOfAmiRNAs = sorted(list(final_results_for_output), key=lambda x: (final_results_for_output[x]["family"],
                                                              -len(final_results_for_output[x]['internal_node']),
                                                               final_results_for_output[x]['internal_node'],
                                                              -len(final_results_for_output[x]['targets']),
                                                              final_results_for_output[x]['score']))
    df = pd.DataFrame(columns=list(final_results_for_output[listOfAmiRNAs[0]].keys()))
    for i in listOfAmiRNAs:
        df.loc[i] = [final_results_for_output[i][df.columns[k]] for k in range(len(df.columns))]
    df.to_csv(final_output_path)



######################################################################################
def constructNewAmiRNAsDict(dict_results_first, dict_results_plaza, new_dic):
    all_internal_nodes = set(list(dict_results_first.keys()) + list(dict_results_plaza.keys()))
    for internalNode in all_internal_nodes:
        if (internalNode in dict_results_plaza) and (internalNode in dict_results_first):
            amiRNAs_first = list(dict_results_first[internalNode].keys())
            amiRNAs_plaza = list(dict_results_plaza[internalNode].keys())
            family = dict_results_first[internalNode][amiRNAs_first[0]][0] +\
                     "+"+ dict_results_plaza[internalNode][amiRNAs_plaza[0]][0]
            for amiRNA in amiRNAs_first:
                dict_results_first[internalNode][amiRNA][0] = family
            for amiRNA in amiRNAs_plaza:
                dict_results_plaza[internalNode][amiRNA][0] = family
            amiRNAs = set(amiRNAs_first + amiRNAs_plaza)
            new_dic[internalNode] = {}
            for amiRNA in amiRNAs:
                if (amiRNA in dict_results_first[internalNode]) and (not amiRNA in dict_results_plaza[internalNode]):
                    new_dic[internalNode][amiRNA] = dict_results_first[internalNode][amiRNA]
                elif (amiRNA in dict_results_plaza[internalNode]) and (not amiRNA in dict_results_first[internalNode]):
                    new_dic[internalNode][amiRNA] = dict_results_plaza[internalNode][amiRNA]
                else:
                    # intersection
                    if len(dict_results_first[internalNode][amiRNA][3]) > len(dict_results_plaza[internalNode][amiRNA][3]):
                        new_dic[internalNode][amiRNA] = dict_results_first[internalNode][amiRNA]

                    else:
                        new_dic[internalNode][amiRNA] = dict_results_plaza[internalNode][amiRNA]
        else:
            if internalNode in dict_results_plaza:
                new_dic[internalNode] = dict_results_plaza[internalNode]
            else:
                new_dic[internalNode] = dict_results_first[internalNode]

#####################################################################################################################
def filterFinalOutputAccordingToSimilarity(new_dic, numberPerInternalNode, numberOfMismatches, finalResults):
    allPossibleInternalNodes = list(new_dic.keys())
    allPossibleInternalNodes.sort(key = lambda x: len(x), reverse=True)
    ancestral_results = {}
    numOfInternalNodes = len(allPossibleInternalNodes)
    for i in range(numOfInternalNodes):
        internalNode = allPossibleInternalNodes[i]
        finalResults[internalNode] = {}
        if internalNode in ancestral_results:
            ancetors_MiRs = ancestral_results[internalNode]
        else:
            ancetors_MiRs = {}
        miRs_to_add = new_dic[internalNode]
        amiRNAs_candidates = sorted(list(miRs_to_add), key=lambda x: (-len(miRs_to_add[x][3]), miRs_to_add[x][1]))
        counterOfAdded = 0
        for new_miR in amiRNAs_candidates:
            if counterOfAdded >= numberPerInternalNode:
                break
            filtered = False
            for added_miR in ancetors_MiRs:
                if (isSimilar(new_miR, added_miR, numberOfMismatches)):
                    if new_dic[internalNode][new_miR][1] >= ancetors_MiRs[added_miR][1]:
                        filtered = True
                        break
            if not filtered:
                finalResults[internalNode][new_miR] = new_dic[internalNode][new_miR]
                if not internalNode in ancestral_results:
                    ancestral_results[internalNode] = {}
                ancestral_results[internalNode][new_miR] = new_dic[internalNode][new_miR]
                counterOfAdded += 1

        for j in range(i, numOfInternalNodes):
            if set(allPossibleInternalNodes[j]).issubset(set(internalNode)):
                if not allPossibleInternalNodes[j] in ancestral_results:
                    ancestral_results[allPossibleInternalNodes[j]] = ancestral_results[internalNode]
                else:
                    ancestral_results[allPossibleInternalNodes[j]].update(ancestral_results[internalNode])
        del ancestral_results[internalNode]


##################################################################################
def getDictOfResults(data_frame):
    dic = {}
    dic_family = {}
    amiRNAs = list(data_frame.index)
    for amiRNA in amiRNAs:
        row = data_frame.loc[amiRNA]
        score = row['score']
        perfectMatch = row['perfectMatch']
        targets = parseTargets(row['targets'])
        belowThresholdTargets = parseTargets(row['belowThresholdTargets'])
        color = row['color']
        family = row['family']
        internal_node = parseInternalNode(row['internal_node'])
        amiRNA_attr = [family, score, perfectMatch, targets, belowThresholdTargets, color]
        if internal_node in dic:
            dic[internal_node][amiRNA] = amiRNA_attr
        else:
            dic[internal_node] = {amiRNA: amiRNA_attr}
            dic_family[family] = internal_node
    return dic, dic_family
#####################################################################################################################
def parseTargets(targets_str):
    targets = []
    targets_str = targets_str[1:-1]
    if targets_str == '':
        return []
    p = re.compile("\(.*?\)")
    pattern_target = re.compile('([\w]+).*?,[\s]*([\S]+)?\)')
    lstOfTargets = p.findall(targets_str)
    for target_str in lstOfTargets:
        target_score = pattern_target.findall(target_str)[0]
        target = target_score[0]
        score = float(target_score[1])
        targets.append((target, score))
    return targets
####################################################################################################################
def parseInternalNode(internal_node_str):
    p = re.compile('([\w]+)')
    internal_node = tuple(p.findall(internal_node_str))
    return internal_node
##################################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="creates sequence files")
    parser.add_argument("--first_partition", "-f", help = "output of first partition")
    parser.add_argument("--second_partition", "-s", help = "output of second partition")
    parser.add_argument("--numPerInternalNode", "-n", type = int, help = "number of AmiRNAs per internal node")
    parser.add_argument("--numOfMismatches", "-m", type=int, help="number of allowed mismatches")
    parser.add_argument("--output_file_path", "-o", help="output file path")
    args = parser.parse_args()
    first_partition_path = args.first_partition
    second_partition_path = args.second_partition
    numPerInternalNode = args.numPerInternalNode
    numOfMismatches = args.numOfMismatches
    output_file_path = args.output_file_path
    parseResultsFilesOfDiffPartitions(first_partition_path, second_partition_path,
                                      numPerInternalNode, numOfMismatches, output_file_path)
