#!/usr/bin/env python
docstring='''
extint.py complex.pdb
    list all interface residues within 4.5 A between chain A and chain B
'''
import sys,os
import re
import subprocess

extint_exe=os.path.join(os.path.abspath(os.path.dirname(__file__)),"extint")
dist_cutoff=4.5
interface_residue_pattern=re.compile("\n\s+(\w\s+\d+)\s+\w{3}\s+\w+\s+(\w\s+\d+)\s+\w{3}\s+\w+\s+[.\d]+")

def get_interface_residue_list(PDBtxt,chain_list=["A","B"]):
    '''get a list of interface residues'''
    cmd=' '.join([extint_exe,
        "-r",chain_list[0],"-l",chain_list[1],"-d",str(dist_cutoff)])
    p=subprocess.Popen(cmd,shell=True,stdin=subprocess.PIPE,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout,stderr=p.communicate(input=PDBtxt)
    interface_list=[]
    for residue_A,residue_B in set(interface_residue_pattern.findall(stdout)):
        interface_list.append(residue_A.replace(' ',''))
        interface_list.append(residue_B.replace(' ',''))
    interface_list=sorted(set(interface_list))
    return interface_list

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    fp=open(sys.argv[1],'rU')
    PDBtxt=fp.read()
    fp.close()
    interface_list=get_interface_residue_list(PDBtxt)
    txt='\n'.join(interface_list)+'\n'
    if len(sys.argv)<3:
        sys.stdout.write(txt)
    else:
        fp=open(sys.argv[2],'w')
        fp.write(txt)
        fp.close()
