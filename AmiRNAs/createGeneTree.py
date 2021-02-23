from ete3 import Tree
import re
import pickle
from Bio import AlignIO
import argparse


def createSubgroupsOfGenes(tree_path, pathForSubgroups):
    #setOfGenes = set()
    setOfGenes = []
    tree = Tree(tree_path)
    dicNodesLeaves = tree.get_cached_content()
    for node in tree.traverse("postorder"):
        #setOfGenes.add(dicNodesLeaves[node])
        currSet = dicNodesLeaves[node]
        setOfNames = set()
        for n in currSet:
            setOfNames.add(n.name)
        if (len(setOfNames) == 1):
            continue
        setOfGenes.append(setOfNames)
    with open(pathForSubgroups, 'wb') as handle:
        pickle.dump(setOfGenes, handle)
    return

#####################################################################################
def convert_alignment_format(input_msa_clustal, output_msa_phylip):
    AlignIO.convert(input_msa_clustal, "clustal", output_msa_phylip, "phylip-sequential")
    return
######################################################################################

def runCommandMAFFT(fasta_input, msa_output):
    cmd = "module load mafft/mafft7149\n"
    #cmd = "mafft --nuc --auto --maxiterate 1000 --clustalout --preservecase "+ fasta_input+" > "+ msa_output+"\n"
    cmd += "mafft --amino --auto --maxiterate 1000 --clustalout --preservecase "+ fasta_input+" > "+ msa_output+"\n"
    return cmd

########################################################################################
def runCommandTreeCreation(msa_file):
    # Sequential of interleaved? Need to check- I don't remember
    #cmd = "PhyML -i "+msa_file +" --sequential -d nt -b 0 -m 012345 --run_id 1 -f m -c 4 -a e -v e -t e -o tlr --no_memory_check\n"
    cmd = "PhyML -i " + msa_file + " --sequential -d aa -b 0 --run_id 1 -f m -c 4 -a e -v e -o tlr --no_memory_check\n"
    return cmd

##########################################################################################3
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run analyses")
    parser.add_argument("--input", "-i", help="input file")
    parser.add_argument("--output", "-o", help="output_file")
    parser.add_argument("--function", "-f", type = int, help="0- mafft, 1- msa conversion, 2- phyML, 3- gene subgroups")
    args = parser.parse_args()
    input = args.input
    output = args.output
    function = args.function
    if function == 0:
        runCommandMAFFT(input, output)
    elif function == 1:
        convert_alignment_format(input, output)
    elif function == 2:
        runCommandTreeCreation(input)
    elif function == 3:
        createSubgroupsOfGenes(input, output)
    else:
        print("No chosen function!")

