import os
from subprocess import Popen
import argparse
import createGeneTree

def prepareForWMD3Run(familyFullDir):
    listOfFamilies = os.listdir(familyFullDir)
    for family in listOfFamilies:
        fullFamilyDir = os.path.join(familyFullDir, family)
        createJobGeneTree(fullFamilyDir)

################################################################3
def createJobGeneTree(fullFamilyDir):
    scripts_dir = "/groups/itay_mayrose/anatshafir1/AmiRNA_lib/scripts"
    createGeneTreeScript = os.path.join(scripts_dir, "createGeneTree.py")
    cmd = createHeaderForJob(fullFamilyDir, "createGeneTree")
    cmd += "module load R/3.5.1\n"
    cmd += "cd "+ fullFamilyDir +"\n"
    cmd += createGeneTree.runCommandMAFFT(os.path.join(fullFamilyDir, 'aaSequences.fasta'), os.path.join(fullFamilyDir, "msa.clustal"))
    cmd += "python "+ createGeneTreeScript + " -i "+os.path.join(fullFamilyDir,  "msa.clustal") + " -o "+ os.path.join(fullFamilyDir, "msa.phy") + " -f 1\n"
    cmd += createGeneTree.runCommandTreeCreation(os.path.join(fullFamilyDir,  "msa.phy"))
    cmd += "R CMD BATCH '--args input="+ "\""+ os.path.join(fullFamilyDir, "msa.phy_phyml_tree_1.txt") + "\" output=\""+ os.path.join(fullFamilyDir, "tree.newick") + "\"' "+ os.path.join(scripts_dir, "rootPhyMLTree.R") + " " + os.path.join(fullFamilyDir, "R_output.out") +"\n"
    cmd += "python "+ createGeneTreeScript + " -i "+os.path.join(fullFamilyDir,  "tree.newick") + " -o "+ os.path.join(fullFamilyDir, "subgroups.pkl")+" -f 3\n"
    jobPath = os.path.join(fullFamilyDir, "GeneTree.sh")
    new_job_file = open(jobPath, 'w')
    new_job_file.write(cmd)
    new_job_file.close()
    Popen(["qsub", jobPath])


###############################################################
def createHeaderForJob(path, job_name, ncpu=1):
    text = ""
    text += "#!/bin/bash\n\n"
    text += "#PBS -S /bin/bash\n"
    text += "#PBS -r y\n"
    text += "#PBS -q itay_25_2\n"
    text += "#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n"
    text += "#PBS -N " + job_name + "\n"
    text += "#PBS -e " + path + "/" + job_name + ".ER" + "\n"
    text += "#PBS -o " + path + "/" + job_name + ".OU" + "\n"
    text += "#PBS -l select=ncpus=" + str(ncpu) + ":mem=2gb\n"
    return text






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="creates sequence files")
    parser.add_argument("--output", "-o", help = "output directory")
    args = parser.parse_args()
    output = args.output
    prepareForWMD3Run(output)