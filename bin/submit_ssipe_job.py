#!/usr/bin/env python
'''crontab script to submit SSIPe jobs by qsub'''
import os,sys
import re
import datetime
import fcntl

from module.configure import max_job_num
from module.configure import frontend_server,frontend_username,frontend_output_dir
from module.configure import pending_txt,log_txt
from module.configure import output_dir,jobID_prefix,output_list,output_list3
from module.configure import get_final_score,get_final_score2,get_final_score3
from module.jobsubmit import showq,submit_job
from module.configure import lock_target

def check_job_overload(showq_txt):
  '''check if there are too many jobs in the queue'''
  #job_pattern=re.compile("\\b"+jobID_prefix+"\d+")
  #current_job_num=len([line for line in showq_txt.splitlines() \
  #  if job_pattern.findall(line)])
  current_job_num=len([line for line in showq_txt.splitlines()])
  if current_job_num>max_job_num:
    fp=open(log_txt,'w')
    fp.write("Too many jobs: %d>%d\n"%(current_job_num,max_job_num))
    fp.close()
    exit()
  else:
    if os.path.isfile(log_txt):
      os.unlink(log_txt)
  return

def get_pending_list():
  '''get a list of pending/running jobs'''
  fp=open(pending_txt)
  pending_list=[line.split()[0] for line in fp.read().splitlines() if line.strip()]
  fp.close()
  #for debug
  print "pending jobs:"
  print '\n'.join(pending_list)
  return pending_list

def get_unsubmitted_list(pending_list,showq_txt):
  '''get a list of pending jobs that are neither submitted nor finished'''
  unsubmitted_list=[]
  for jobID in pending_list:
    if re.findall(r"\s+"+jobID+"\s+",showq_txt):
    #if re.search(jobID,showq_txt):
      #print("job %s found in running list, continue"%jobID)
      continue # submitted
    data_dir=os.path.join(output_dir,jobID)
    force_field_file=os.path.join(data_dir,"force_field")
    fp=open(force_field_file,'rU')
    force_field=fp.read().strip()
    fp.close()
    if force_field=="SSIPE" or force_field=="SSIP":
      for item in output_list:
        if not os.path.isfile(os.path.join(output_dir,jobID,item)):
          unsubmitted_list.append(jobID)
          break # some output is not generated
    else:
      for item in output_list3:
        if not os.path.isfile(os.path.join(output_dir,jobID,item)):
          unsubmitted_list.append(jobID)
          break # some output is not generated
  print "unsubmitted jobs:"
  print '\n'.join(unsubmitted_list)
  return unsubmitted_list

def scp_data_folder_from_frontend_server(pending_list):
  '''scp pending jobs data folder to output_dir if the data folder is absent'''
  tranfered_list=os.listdir(output_dir)
  untransfered_list=list(set(pending_list)-set(tranfered_list))
  #for debug
  print "untransfered list:"
  print '\n'.join(untransfered_list)
  if not untransfered_list:
    return []

  frontend_dir="%s@%s:%s"%(frontend_username,frontend_server,frontend_output_dir)
  for jobID in untransfered_list:
    sys.stderr.write("scp %s from %s\n"%(jobID,frontend_server))
    os.system("scp -r "+os.path.join(frontend_dir,jobID)+" "+output_dir)
  return untransfered_list

def submit_pending_list(pending_list):
  '''submit prediction jobs specified by jobID list "pending_list"'''
  for jobID in pending_list:
    data_dir=os.path.join(output_dir,jobID)
    force_field_file=os.path.join(data_dir,"force_field")
    fp=open(force_field_file,'rU')
    force_field=fp.read().strip()
    fp.close()
    if force_field=="SSIPE":
      submit_job(
        jobname=os.path.join(data_dir,jobID), # PBS script name
        cmd=get_final_score+' '+data_dir,    # full command
      )
    else:
      if force_field=="SSIP":
        submit_job(
          jobname=os.path.join(data_dir,jobID), # PBS script name
          cmd=get_final_score2+' '+data_dir,    # full command
        )
      else: # EVOEF
        submit_job(
          jobname=os.path.join(data_dir,jobID), # PBS script name
          cmd=get_final_score3+' '+data_dir,    # full command
        )
    date=str(datetime.datetime.now())
    fp=open(os.path.join(data_dir,"status"),'a')
    #fp.write("qsub\t"+date+'\n')
    fp.write("sbatch\t"+date+'\n')
    fp.close()
  return

if __name__=="__main__":
  #use flock to avoid multiple instance of crontab
  #fp=open(lock_target)
  #try:
  #  fcntl.flock(fp,fcntl.LOCK_EX|fcntl.LOCK_NB)
  #except IOError:
  #  sys.stdout.write("Unable to lock %s\n"%lock_target)
  #  exit()

  #query job status
  showq_txt=showq() # list running jobs
  #print("running jobs:")
  #print(showq_txt)
  # die if there are too many jobs
  check_job_overload(showq_txt)

  # get a list of pending/running jobs
  pending_list=get_pending_list()
  
  # if absent, scp data folder from web portal
  scp_data_folder_from_frontend_server(pending_list)
  #scp_data_folder_from_frontend_server(unsubmitted_list)
  
  # get a list of pending jobs that are neither submitted nor finished
  #note previously there's a bug
  unsubmitted_list=get_unsubmitted_list(pending_list,showq_txt)

  # submit jobs
  submit_pending_list(unsubmitted_list)
