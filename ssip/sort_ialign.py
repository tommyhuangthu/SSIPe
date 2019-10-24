#!/usr/bin/env python
docstring='''sort_ialign_profile.py align.out align.out.sort
    sort ialign profile so that ialign match is sorted by IS-score in
    descending order
'''
import sys

def sort_ialign_profile(align_txt=''):
    '''sort ialign output profile'''
    align_lines=align_txt.splitlines()
    isscore,sort_align_lines=zip(*sorted([(float(line.split()[5]),line
        ) for line in align_lines[7:]],reverse=True))
    txt=''.join([line+'\n' for line in align_lines[:7]+list(sort_align_lines)])
    return txt

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    fp=open(sys.argv[1],'rU')
    align_txt=fp.read()
    fp.close()

    txt=sort_ialign_profile(align_txt)
    if len(sys.argv)<3:
        sys.stdout.write(txt)
    else:
        fp=open(sys.argv[2],'w')
        fp.write(txt)
        fp.close()
