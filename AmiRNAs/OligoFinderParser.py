import pandas as pd
import os
import re
import pickle
import argparse

def parseOligoResults(dirOfResults, dictAmiRNAOligosPath):
    files = os.listdir(dirOfResults)
    seqs = ["miR-s", "miR-a", "miR*s", "miR*a"]
    dictAmiOligo = {}
    for file in files:
        if not file.endswith(".OU"):
            continue
        AmiRNA_pattern = re.compile("[A|T|C|G]+")
        AmiRNA_lst = AmiRNA_pattern.findall(file)
        if len(AmiRNA_lst) == 0:
            raise Exception("ERROR!!! parseOligoResults: Undefined output job file!\n")
        AmiRNA = AmiRNA_lst[0]
        if len(AmiRNA) != 21:
            raise Exception("ERROR!!! parseOligoResults: Not AmiRNA!\n")
        fullPathOligoOutput = os.path.join(dirOfResults, file)
        oligo_output_file = open(fullPathOligoOutput, 'r')
        content = oligo_output_file.read()
        oligo_output_file.close()
        #delete files#
        os.remove(fullPathOligoOutput)
        os.remove(os.path.join(dirOfResults, AmiRNA+".ER"))
        os.remove(os.path.join(dirOfResults, "Oligo_"+ AmiRNA+".sh"))
        # extract info to dictionary
        lines = re.split("[\s]+", content)
        lines = list(filter(lambda x: x != '', lines))
        lines = lines[1:]
        dictOligo = {}
        for i in range(len(lines)):
            dictOligo[seqs[i]] = lines[i]
        dictAmiOligo[AmiRNA] = dictOligo
    if not os.path.exists(dictAmiRNAOligosPath):
        with open (dictAmiRNAOligosPath, 'wb') as handle:
            pickle.dump(dictAmiOligo, handle)
    else:
        with open(dictAmiRNAOligosPath, 'rb') as handle:
            soFarOligos = pickle.load(handle)
        soFarOligos.update(dictAmiOligo)
        with open(dictAmiRNAOligosPath, 'wb') as handle:
            pickle.dump(soFarOligos, handle)


###########################################################################3
def getOligos(pathOfMiRs, pathOfDictOligos):
    with open(pathOfDictOligos, 'rb') as handle:
        dictOligos = pickle.load(handle)
    df = pd.read_csv(pathOfMiRs, index_col=0)
    AmiRNAs = list(df['AmiRNAs'])
    oligos = []
    pattern_antiSense = re.compile("[A|T|C|G]+")
    for AmiRNA in AmiRNAs:
        antiSenseWithBackbone = dictOligos[AmiRNA]["miR*s"]
        antiSense = pattern_antiSense.findall(antiSenseWithBackbone)[0]
        oligos.append(antiSense)
    df['miR*s'] = oligos
    df.to_csv(pathOfMiRs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run analyses")
    parser.add_argument("--directory_of_results", "-d",  help = "directory where to store the results")
    parser.add_argument("--inputMiRPath", "-i", help="MiR file")
    parser.add_argument("--oligo_dict", '-o', help = "oligo dictionary file")
    parser.add_argument("--function", '-f', type=int, help = "function")
    args = parser.parse_args()
    input = args.inputMiRPath
    res_dir = args.directory_of_results
    oligo_dict = args.oligo_dict
    func = args.function
    if func:
        parseOligoResults(res_dir, oligo_dict)
    else:
        getOligos(input, oligo_dict)






