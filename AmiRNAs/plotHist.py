import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
import argparse



def createDataForHist(dir_of_families, pathOfFigure):
    families = os.listdir(dir_of_families)
    familiesNumPerInternalNode = []
    numOfFamilies = 0

    for family in families:
        full_family_dir = os.path.join(dir_of_families, family)
        path_res_per_subgroups = os.path.join(full_family_dir, "subgroups_res.pkl")
        res_subgroups_path_final = os.path.join(full_family_dir, "subgroups_res_final.pkl")
        dicGeneIsoform_path = os.path.join(full_family_dir, "geneIsoform.pkl")
        if not os.path.isdir(full_family_dir):
            continue
        with open(dicGeneIsoform_path, 'rb') as handleIsoform:
            geneIsoDic = pickle.load(handleIsoform)
        if len(geneIsoDic) < 2:  # no need for families with a single gene
            print(family)
            continue
        res_path = None
        numOfFamilies += 1
        if os.path.exists(res_subgroups_path_final):
            res_path = res_subgroups_path_final
        else:
            res_path = path_res_per_subgroups
        perFamilyStat = getCounts(res_path)
        familiesNumPerInternalNode.append(perFamilyStat)
    finalStats = np.array(familiesNumPerInternalNode)
    n, bins, patches = plt.hist(x=finalStats, bins='auto', color='#0504aa',
                                alpha=0.7, rwidth=0.85)
    plt.title('Number of AmiRNAs per internal node for each family')
    plt.xlabel('Number of AmiRNAs')
    plt.ylabel('Number of families')
    plt.savefig(pathOfFigure)
    print("Number of families is:", str(numOfFamilies))

def getCounts(res_path):
    with open(res_path, 'rb') as handle:
        subgroupRes = pickle.load(handle)
        lstCounts = []
    for subgroup in subgroupRes:
        AmiRNAs = subgroupRes[subgroup]
        count = 0
        for miR in AmiRNAs:
            count += 1
        lstCounts.append(count)
    if lstCounts == []:
        return 0
    return int(np.array(lstCounts).mean())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--input_file_path", "-i", help = "input directory")
    parser.add_argument("--output_file_path", "-o", help = "path for the figure")


    args = parser.parse_args()
    input_file_path = args.input_file_path
    output_file_path = args.output_file_path

    createDataForHist(input_file_path, output_file_path)

