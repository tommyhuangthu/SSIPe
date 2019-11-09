#!/usr/bin/env python
'''crontab script to transfer finished BindProfX jobs to web portal'''
import os,sys
import datetime
import fcntl

from module.configure import pending_txt
from module.configure import output_dir,output_list,output_list3
from module.configure import frontend_server,frontend_username,frontend_output_dir
from module.configure import lock_target

def get_pending_list():
  '''get a list of pending/running jobs'''
  fp=open(pending_txt)
  pending_list=[line.split()[0] for line in fp.read().splitlines() if line.strip()]
  fp.close()
  print "pending jobs:"
  print '\n'.join(pending_list)
  return pending_list

def get_finish_list(pending_list):
  '''get a list of pending jobs that are finished'''
  finish_list=[]
  for jobID in pending_list:
    data_dir=os.path.join(output_dir,jobID)
    force_field_file=os.path.join(data_dir,"force_field")
    fp=open(force_field_file,'rU')
    force_field=fp.read().strip()
    fp.close()
    if force_field=="SSIPE" or force_field=="SSIP":
      if len(output_list)==len([f for f in output_list if \
      os.path.isfile(os.path.join(output_dir,jobID,f))]):
        finish_list.append(jobID) # some output is not generated
    else:
      if len(output_list3)==len([f for f in output_list3 if \
      os.path.isfile(os.path.join(output_dir,jobID,f))]):
        finish_list.append(jobID) # some output is not generated
  print "finished jobs:"
  print '\n'.join(finish_list)
  return finish_list

def scp_data_folder_to_frontend_server(finish_list):
  '''scp job data folder to web portal'''
  frontend_dir="%s@%s:%s"%(frontend_username,frontend_server,frontend_output_dir)
  for jobID in finish_list:
    data_dir=os.path.join(output_dir,jobID)
    date=str(datetime.datetime.now())
    fp=open(os.path.join(data_dir,"status"),'a')
    fp.write("send_to_web_portal\t"+date+'\n')
    fp.close()
    os.system("scp -r "+data_dir+' '+frontend_dir)
  return

if __name__=="__main__":
  # use flock to avoid multiple instance of crontab
  #fp=open(lock_target)
  #try:
  #  fcntl.flock(fp,fcntl.LOCK_EX|fcntl.LOCK_NB)
  #except IOError:
  #  sys.stdout.write("Unable to lock %s\n"%lock_target)
  #  exit()
    
  # get a list of pending jobs
  pending_list=get_pending_list()
  # get a list of pending jobs that are finished 
  # TODO send warning to webmaster if too many failed qsub
  finish_list=get_finish_list(pending_list)
  # transfer finished jobs to web portal
  scp_data_folder_to_frontend_server(finish_list)
