import pickle
import os
import argparse
from subprocess import Popen


def runWMD3AmiRNA_all(dataSetPath, restart):
    counter = 0
    listOfDirs = os.listdir(dataSetPath)
    for family in listOfDirs:
        pathOfFamily = os.path.join(dataSetPath, family)
        if not os.path.isdir(pathOfFamily):
            continue
        if not os.path.exists(os.path.join(pathOfFamily, "tree.newick")):
            print(pathOfFamily)
            continue
        runWMD3(pathOfFamily, restart)
        counter += 1
        print("number of families is:", str(counter))
#####################################################################3
def runWMD3OneTwoGenesAmiRNAAll(dataSetPath):
    counter = 0
    listOfDirs = os.listdir(dataSetPath)
    for family in listOfDirs:
        pathOfFamily = os.path.join(dataSetPath, family)
        if not os.path.isdir(pathOfFamily):
            continue
        if os.path.exists(os.path.join(pathOfFamily, "tree.newick")):
            continue
        runWMD3OneTwoGenes(pathOfFamily)
        counter += 1
    print("number of families is:", str(counter))
#######################################################################
def runWMD3OneTwoGenes(familyPath):
    dicGeneIsoform_path = os.path.join(familyPath, "geneIsoform.pkl")
    with open(dicGeneIsoform_path, 'rb') as handleIsoform:
        geneIsoDic = pickle.load(handleIsoform)
    if len(geneIsoDic) > 2:
        print(familyPath)
        file_err = open(os.path.join(familyPath, "ERROR.txt"), 'w')
        file_err.write("No tree but more than two genes!\n")
        file_err.close()
    subgroup_folder = os.path.join(familyPath, str(0))
    if not os.path.exists(subgroup_folder):
        os.makedirs(subgroup_folder)
    if not os.path.exists(os.path.join(subgroup_folder, "WMD3.sh")):
        runWMD3PerSubgroup(subgroup_folder, list(geneIsoDic.keys()), geneIsoDic)

#######################################################################
def runWMD3(familyPath, restart):

    subgroups_path = os.path.join(familyPath, "subgroups.pkl")
    dicGeneIsoform_path = os.path.join(familyPath, "geneIsoform.pkl")
    if not os.path.exists(subgroups_path):
        print("file does not exist", subgroups_path)
        return
    with open(subgroups_path, 'rb') as handle:
        subgroups = pickle.load(handle)
    with open(dicGeneIsoform_path, 'rb') as handleIsoform:
        geneIsoDic = pickle.load(handleIsoform)
    for i in range(len(subgroups)):
        subgroup_folder = os.path.join(familyPath, str(i))
        if not restart:
            if os.path.exists(os.path.join(subgroup_folder, "WMD3.sh")):
                continue
        if not os.path.exists(subgroup_folder):
            os.makedirs(subgroup_folder)
        if not restart:
            runWMD3PerSubgroup(subgroup_folder, subgroups[i], geneIsoDic)
        else:
            if i == len(subgroups)-1:
                runWMD3PerSubgroup(subgroup_folder, subgroups[i], geneIsoDic, restart)

    return

#####################################################################################

def runWMD3PerSubgroup(subgroup_folder, subgroup, dicGeneIsoform, restart = 0):
    input_file_path = os.path.join(subgroup_folder, "genes.txt")
    if not restart:
        file = open(input_file_path, 'w')
        for gene in subgroup:
            file.write(dicGeneIsoform[gene]+'\n')
        file.close()
    jobName = "WMD3"
    if restart:
        jobName += "_1"
    cmd = createHeaderForJob(subgroup_folder, jobName)
    cmd += "cd "+ subgroup_folder +"\n"
    if not restart:
        cmd += "singularity exec --bind /groups/itay_mayrose/anatshafir1/AmiRNA_lib:/mnt --bind /tmp:/var/opt/amirna/tmp /powerapps/singularity/sing_files/AmiRNA_WMD3/AmiRNA.simg perl /usr/local/sbin/batch_Designer.pl -c /etc/AmiRNA/AmiRNA.xml -g Araport11_genes.201606.cdna.fasta -f "+input_file_path.replace('/groups/itay_mayrose/anatshafir1/AmiRNA_lib', '/mnt')+"\n"
    else:
        cmd+= "singularity exec --bind /groups/itay_mayrose/anatshafir1/AmiRNA_lib:/mnt --bind /tmp:/var/opt/amirna/tmp /powerapps/singularity/sing_files/AmiRNA_WMD3/AmiRNA.simg perl /groups/itay_mayrose/anatshafir1/AmiRNA_lib/scripts/batch_Designer_all_subgroups.pl -c /etc/AmiRNA/AmiRNA.xml -g Araport11_genes.201606.cdna.fasta -f "+input_file_path.replace('/groups/itay_mayrose/anatshafir1/AmiRNA_lib', '/mnt')+"\n"
    if restart:
        job_file = os.path.join(subgroup_folder, "WMD3_1.sh")
    else:
        job_file = os.path.join(subgroup_folder, "WMD3.sh")
    file = open(job_file, 'w')
    file.write(cmd)
    file.close()
    Popen(["qsub", job_file])

def createHeaderForJob(path, job_name, ncpu = 1):
    text = ""
    text += "#!/bin/bash\n\n"
    text += "#PBS -S /bin/bash\n"
    text += "#PBS -r y\n"
    text += "#PBS -q itay_25_2\n"
    text += "#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n"
    text += "#PBS -N "+ job_name+"\n"
    text += "#PBS -e " + path + "/"+job_name+".ER" + "\n"
    text += "#PBS -o " + path + "/" + job_name +".OU"+ "\n"
    text += "#PBS -l select=ncpus="+ str(ncpu)+ ":mem=2gb\n"
    return text

#######################################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run analyses")
    parser.add_argument("--inputFolder", "-i", help = "families folder")
    parser.add_argument("--restart", "-r", type = int, default= 0, help = "restart run")
    parser.add_argument("--oneTwoGenes", '-g', type = int, default=0, help = "1 if one or two genes")
    args = parser.parse_args()
    families_dir = args.inputFolder
    restart = args.restart
    oneTwoGenes = args.oneTwoGenes
    if not oneTwoGenes:
        runWMD3AmiRNA_all(families_dir, restart)
    else:
        runWMD3OneTwoGenesAmiRNAAll(families_dir)



