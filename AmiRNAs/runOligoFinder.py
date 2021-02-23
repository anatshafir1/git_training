import argparse
import os
from subprocess import Popen
import pandas as pd
from createJob import createHeaderForJob



def runOligo(pathAmiRNAs, pathResDir, queue, start, stop):
    if not os.path.exists(pathResDir):
        os.makedirs(pathResDir)
    df = pd.read_csv(pathAmiRNAs, index_col=0)
    AmiRNAs = list(df["AmiRNAs"])
    if not stop:
        stop = len(AmiRNAs)
    for i in range(start, stop):
        AmiRNA = AmiRNAs[i]
        job_path = createOligoJob(AmiRNA, pathResDir, queue)
        Popen(['qsub', job_path])


def createOligoJob(AmiRNA, pathResDir, queue):
    COMMAND = "singularity exec --bind /groups/itay_mayrose/anatshafir1/AmiRNA_lib:/mnt --bind /tmp:/var/opt/amirna/tmp /powerapps/singularity/sing_files/AmiRNA_WMD3/AmiRNA.simg perl /groups/itay_mayrose/anatshafir1/AmiRNA_lib/scripts/testOligoDesign.pl -s "
    header = createHeaderForJob(pathResDir, "Oligo", AmiRNA, queue, 1, ncpu=1)
    command = "cd "+ pathResDir +"\n"
    COMMAND += AmiRNA +"\n"
    path =  os.path.join(pathResDir, "Oligo_"+ AmiRNA+".sh")
    file = open(path, 'w')
    job_content = header + command + COMMAND
    file.write(job_content)
    file.close()
    return path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="run analyses")
    parser.add_argument("--inputFile", "-i", help = "input file with AmiRNAs")
    parser.add_argument("--directory_of_results", "-d",  help = "directory where to store the results")
    parser.add_argument("--queue", '-q', default="itaym", help = "queue to use")
    parser.add_argument("--start", '-s', type = int, help= "start from this AmiRNA")
    parser.add_argument("--end", '-e', type = int, default = 0, help = "stop at this AmiRNA")
    args = parser.parse_args()
    input = args.inputFile
    res_dir = args.directory_of_results
    queue = args.queue
    start = args.start
    stop = args.end
    runOligo(input, res_dir, queue, start, stop)
