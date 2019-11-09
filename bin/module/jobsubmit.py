#!/usr/bin/env python
docstring='''module to submit and query jobs using PBS job scheduling system
'''
import sys,os
from configure import sbatch,qstat,jobID_prefix
import subprocess
from string import Template
import time

#PBS_template=Template("""#!/bin/bash
##PBS -q urgent
##PBS -l nodes=1:ppn=1
##PBS -l walltime=48:00:00,mem=4000mb
##PBS -N $TAG
##PBS -o $JOBNAME.out
##PBS -e $JOBNAME.err
#$CMD
#""") # $TAG - PBS job name; $JOBNAME - PBS script full path; 
#     # $CMD - full command to execute

PBS_template=Template("""#!/bin/bash
#SBATCH --partition=batch
#SBATCH --account=zhanglab
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4000mb
#SBATCH --time=48:00:00
#SBATCH --job-name=$TAG
#SBATCH --output=$JOBNAME.out
#SBATCH --error=$JOBNAME.err
$CMD
""") # $TAG - PBS job name; $JOBNAME - PBS script full path; 
     # $CMD - full command to execute



def submit_job(jobname,cmd):
    '''write command "cmd" to PBS script "jobname", submit the script'''
    fp=open(jobname,'w')
    fp.write(PBS_template.substitute(dict(
        TAG=os.path.basename(jobname),
        JOBNAME=os.path.abspath(jobname),
        CMD=cmd,
    )))
    fp.close()
    os.chmod(jobname, os.stat(jobname).st_mode|0111)
    p=subprocess.Popen(sbatch+' '+jobname,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr=p.communicate()
    if stdout.strip():
        #sys.stdout.write(jobname+" submitted\n")
        print(jobname+" submitted")
    else:
        print("something wrong with sbatch, sleep 10s")
        time.sleep(10) # something wrong with sbatch
    return stdout.strip() # return sbatch job ID


def showq():
    '''return job queue status'''
    #cmd="qstat -r|grep zhanglabs|grep SSIPE"
    p=subprocess.Popen(qstat,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr=p.communicate()
    return stdout

if __name__=="__main__":
    sys.stderr.write(docstring)
