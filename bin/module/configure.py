#!/usr/bin/env python
'''configuration file for webserver'''
import os

# root path for SSIPe backend server: /nfs/amino-home/zhanglabs/SSIPe
root_dir=os.path.dirname(os.path.dirname(
    os.path.abspath(os.path.dirname(__file__))))
bindir=os.path.dirname(os.path.abspath(os.path.dirname(__file__)))

#### webpage display ####
# HTTP URL to SSIPe web portal http://zhanglab.ccmb.med.umich.edu/SSIPe
http_link="http://zhanglab.ccmb.med.umich.edu/"+os.path.basename(root_dir)
# citation to paper
citation="Xiaoqiang Huang, Wei Zheng, Robin Pearce, Yang Zhang. SSIPe: Accurately estimating protein-protein binding affinity change upon SNP mutations using evolutionary profiles in combination with an optimized physical energy function. (In preparation)"
# javascripts for JSmol
javascript_list=["/jmol/JSmol.min.js", "/jmol/Jmol2.js",
    "/3Dmol/3Dmol-min.js"]

#### job scheduling system ####
# excutable to query job statistics
#qstat="qstat -f|grep Job_Name"
qstat="/sw/dcmb/zhng/bin/qstat -r|grep SSIPE"
# excutable to query job statistics
#qsub="qsub "
sbatch="/sw/dcmb/zhng/bin/sbatch "
# maximum number of jobs to run
max_job_num=300

#### webportal frontend server ####
frontend_server="zhanglab.ccmb.med.umich.edu"
# user name at frontend server
frontend_username="zhanglabs"
# directory of data folder for individual job at frontend server
frontend_output_dir="/www/html/SSIPe/output"

#### directory of executables ####
#bindir=os.path.join(root_dir,"bin")
# where to find executables for calculating profile conservation score
run_align=os.path.join(bindir,"ssip","run_ialign.pl")
run_seqalign=os.path.join(bindir,"ssip","run_seqalign.pl")
get_ssip_score=os.path.join(bindir,"ssip","get_ssip_score.pl")
# where to find executables for calculating physics potential
run_evoef=os.path.join(bindir,"evoef","run_evoef.pl")
# where to find executable to get final overall score
get_final_score=os.path.join(bindir,"get_final_score.py")
get_final_score2=os.path.join(bindir,"get_final_score2.py")
get_final_score3=os.path.join(bindir,"get_final_score3.py")
# where to find executable to extract interface
extint=os.path.join(bindir,"extint","extint")

#### user submission log ####
log_dir=os.path.join(root_dir,"log")
# list of unfinished jobs
# jobID, email, IP, date, chainA length, chainB length, force_field
pending_txt=os.path.join(log_dir,"list_pending.txt")
# log file only present when there are too many pending jobs
log_txt=os.path.join(log_dir,"log.txt")

#### data output ####
output_dir=os.path.join(root_dir,"output") 
# prefix of job ID
jobID_prefix="SSIPE"
# list of output files
output_list=[
    "structure_align.out",  # iAlign alignment output
    "sequence_align.out",   # PSI-BLAST alignment output
    "result.txt", # final score
    "index.html", # result webpage
]
output_list3=[
    "result.txt", # final score
    "index.html", # result webpage
]

#### flock ####
# lock target file whom crontab scripts try to put flock to 
# avoid multiple instance of crontab scripts
lock_target=os.path.join(log_dir,"locker")
