import pickle
import os
import argparse
import re

def checkProteinSeq(dataSetPath, tester_path):
    listOfDirs = os.listdir(dataSetPath)
    tester_file = open(tester_path, 'w')
    text_Log = ""
    for family in listOfDirs:
        pathOfFamily = os.path.join(dataSetPath, family)
        if not os.path.isdir(pathOfFamily):
            continue
        path_prot = os.path.join(pathOfFamily, "aaSequences.fasta")
        path_nuc = os.path.join(pathOfFamily, "nucSequences.fasta")
        toWrite = compareNucToProt(path_prot, path_nuc, family)
        text_Log += toWrite
    tester_file.write(text_Log)

#####################################################################3
def compareNucToProt(path_prot, path_nuc, family):
    log = ""
    dic_name_nuc = getDictNameSeq(path_nuc)
    dic_name_prot = getDictNameSeq(path_prot)
    STOP_CODONS = ["TAG", "UAG", "TAA", "UAA", "TGA", "UGA"]
    START_CODONS = ["ATG", "AUG"]
    print("*****\n")
    log += "*****\n"
    for gene in dic_name_nuc:
        nuc_seq = dic_name_nuc[gene]
        prot_seq = dic_name_prot[gene]
        if len(nuc_seq)%3 != 0:
            print(family, gene, "ERROR: Not divisable by zero!")
            log += (family+" "+ gene + " "+ "ERROR: Not divisable by zero!")
        elif not (nuc_seq[len(nuc_seq)-3:] in STOP_CODONS):
            print(family, gene, "ERROR: the last codon is not a stop codon!")
            log += (family +" "+ gene+" "+ "ERROR: the last codon is not a stop codon!")
        elif not (nuc_seq[:3] in START_CODONS):
            print(family, gene, "ERROR: the start codon is not M!")
            log += (family + " "+ gene+ " "+ "ERROR: the start codon is not M!")
        elif int((len(nuc_seq)/3)-1) != len(prot_seq):
            print(family, gene, "ERROR: length of the protein doesn't match the length of the nuc seq!")
            log += (family+ " "+ gene+" "+"ERROR: length of the protein doesn't match the length of the nuc seq!")
    return log




#######################################################################################################
def getDictNameSeq(path):
    file = open(path, 'r')
    content = file.read()
    file.close()
    dic_name_seq = {}
    reg = re.compile(">([\S]+)[\s]+([\S]+)", re.MULTILINE)
    geneNames_seq = reg.findall(content)
    for gene, seq in geneNames_seq:
        dic_name_seq[gene] = seq
    return dic_name_seq
#######################################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run analyses")
    parser.add_argument("--inputFolder", "-i", help = "families folder")
    parser.add_argument("--tester_file", "-t", help = "output of tester file")
    args = parser.parse_args()
    families_dir = args.inputFolder
    tester_file = args.tester_file
    checkProteinSeq(families_dir, tester_file)




