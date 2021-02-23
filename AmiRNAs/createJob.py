def createHeaderForJob(path, job_name, identifier, queue, memory, ncpu = 1):
    text = ""
    text += "#!/bin/bash\n\n"
    text += "#PBS -S /bin/bash\n"
    text += "#PBS -r y\n"
    text += "#PBS -q "+ queue+ "\n"
    text += "#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH\n"
    text += "#PBS -N "+ job_name+"\n"
    text += "#PBS -e " + path + "/"+identifier+".ER" + "\n"
    text += "#PBS -o " + path + "/" + identifier +".OU"+ "\n"
    text += "#PBS -l select=ncpus="+ str(ncpu)+ ":mem="+ str(memory) +"gb\n"
    return text