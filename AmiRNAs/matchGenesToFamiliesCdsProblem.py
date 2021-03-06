import re
import os
import pickle
import argparse
from Bio import Seq

def getFamilyGeneSequences(outputDir, all_families_short_dir, first_partition):
    homeDir = "/groups/itay_mayrose/anatshafir1/AmiRNA_lib/data"
    all_families_dir = os.path.join(homeDir, all_families_short_dir)
    if first_partition:
        allSeqPath = os.path.join(homeDir, "allSequences.fasta")
        family_dic_path = os.path.join(homeDir, "familyAndGenes.pkl")
    else:
        allSeqPath = os.path.join(homeDir, "all_nuc_sequences_completion.fasta")
        family_dic_path = os.path.join(homeDir, "families_genes_for_completion.pkl")
    output_dir = os.path.join(homeDir, outputDir)
    dicOfGenesSeq = parseAllSeqFile(allSeqPath)
    counter = 0
    with open(family_dic_path, 'rb') as handle:
        family_dic = pickle.load(handle)
    families = list(family_dic.keys())
    #families = ['ABC-A group', 'ABC-B group', 'ABC-C group', 'FAX-type transporter', 'ABC-I group']
    p = re.compile("[\s]+|-|\/")
    for i in range(len(families)):
        counter += 1
        family = families[i]
        genes = family_dic[family]
        family = re.sub(p, "_", family)
        if not os.path.exists(os.path.join(all_families_dir, family)):
            os.makedirs(os.path.join(output_dir, family))
        else:
            continue
        print(family, "Number:", str(i))
        seqPerFamilyPath = os.path.join(output_dir, family,  "nucSequences.fasta")
        seqPerFamilyPathAA = os.path.join(output_dir, family,  "aaSequences.fasta")
        dictOfGeneAndIsoformPath = os.path.join(output_dir, family, "geneIsoform.pkl")
        dicOfGeneIsoforms = {}
        seqPerFamilyFile = open(seqPerFamilyPath, 'w')
        #seqPerFamilyPathAAFile = open(seqPerFamilyPathAA, 'w')
        for gene in genes:
            for geneVariant in dicOfGenesSeq:
                if geneVariant.startswith(gene):
                    dicOfGeneIsoforms[gene] = geneVariant
                    seqPerFamilyFile.write(">"+ gene + "\n")
                    #seqPerFamilyPathAAFile.write(">"+ gene + "\n")
                    seqPerFamilyFile.write(dicOfGenesSeq[geneVariant]+"\n")
                    if len(dicOfGenesSeq[geneVariant]) % 3 != 0:
                        print(family, geneVariant)
                    #seqPerFamilyPathAAFile.write(Seq.translate(dicOfGenesSeq[geneVariant],  cds=True) + "\n")
        seqPerFamilyFile.close()
        #seqPerFamilyPathAAFile.close()
        with open(dictOfGeneAndIsoformPath, 'wb') as handle:
            pickle.dump(dicOfGeneIsoforms, handle)
    print("Number of families is ", str(counter))


def parseAllSeqFile(allSeqPath):
    dic = {}
    file = open(allSeqPath, 'r')
    lines = file.readlines()
    file.close()
    geneNamePattern = re.compile(">[\S]+")
    seqPattern = re.compile("[\S]+")
    sequence = ""
    gene = ""
    for line in lines:
        if ">" in line:
            if gene != "":
                dic[gene] = sequence
            gene = geneNamePattern.findall(line)[0][1:]
            sequence = ""
        else:
            seqExpression = seqPattern.findall(line)
            if seqExpression != []:
                seqPart = seqExpression[0]
                sequence += seqPart
    # should add the last gene
    dic[gene] = sequence
    return dic



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="creates sequence files")
    parser.add_argument("--output", "-o", help = "output directory")
    parser.add_argument("--source", "-s", help =  "source for correct cds")
    parser.add_argument("--first_partition", "-f", type = int, help = "is it the first partition?")

    args = parser.parse_args()
    output = args.output
    source_cds = args.source
    first_partition = args.first_partition
    getFamilyGeneSequences(output, source_cds, first_partition)



