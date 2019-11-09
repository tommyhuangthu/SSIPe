#!/usr/bin/env python
docstring='''
get_final_score.py data_dir
    main script to perform SSIPe prediction

Input files:
    complex.pdb (complex structure)
    mutList.txt (list of mutations to make)
    force_field (optional: EVOEF - evoef score only)
    isscore     (optional: IS-score cutoff for interface similarity.
                 default is 0.5)

Output files:
    structure_align.out (output iAlign alignment from NIL)
    sequence_align.out  (output PSI-BLAST alignment from STRING)   
    ssip_score.txt      (profile score)
    evoef_score.txt     (evoef score)
    result.txt          (final prediction)

    seq.fasta       (fasta format sequence)
    interface.pdb   (interface structure)
    contact.list    (list of inter-chain contact residues)
    index.html      (webpage output)
'''
import os,sys
import shutil
from string import Template
import textwrap
import re

from module.configure import run_align      # iAlign interface alignment
from module.configure import run_seqalign   # psiblast to build sequence interface alignment
from module.configure import get_ssip_score # get ssip profile score
from module.configure import run_evoef      # get evoef physcial score
from module.configure import extint         # extract interface
from module.configure import http_link,javascript_list,citation # make html web page

#from module.reindex_complex import reindex_complex # renumber residue index
from module.pdb2fasta import pdb2seq # convert pdb to fasta

#def make_webpage3(inv_resi_dict,data_dir,isscore,force_field,result_list):
def make_webpage3(data_dir,isscore,force_field,result_list):
    '''make index.html for data folder "data_dir".'''
    html_template=Template("""<html>
<head><meta http-equiv="content-type" content="text/html; charset=UTF-8">
<title>SSIPe job finished</title></head>
<body>
[<a href="$HTTP_LINK">back to SSIPe server</a>]<br/>
<h1 align="center">SSIPe job id $JOBID</h1>

<div style="background:#6599FF;width:350px;">
<b>User input parameters and structure</b></div>
<ul>

<li>Interface similarity cutoff<br>
<font face="Courier">IS-score = $ISSCORE</font>
<br>

<li>Force field used for &Delta;&Delta;G prediction<br>
<font face="Courier">$FF</font>

<li>Sequence [<a href=seq.fasta>download</a>]<br>
<font face="Courier">
&gt;$CHAIN0<br>$SEQUENCE0<br>&gt;$CHAIN1<br>$SEQUENCE1<br></font>
<br>

<li>View of 3D structure structure
<script type="text/javascript" src="$JSMOL0"></script>
<script type="text/javascript" src="$JSMOL1"></script>
<table>
    <td align ="center">
    complex structure [<a href=complex.pdb>download</a>]
    <script type="text/javascript">
    jmolInitialize("$JSMOL_HOME/");
    jmolSetAppletColor("#000000");
    jmolApplet(300,"load complex.pdb; set frank off; calculate structure rama; select all; color chain; spacefill off; wireframe off; cartoons; set ambient 40; spin off; $JSMOL_CMD");
    </script>
    </td>

    <td align ="center">
    interface structure [<a href=interface.pdb>download</a>]
    <script type="text/javascript">
    jmolInitialize("$JSMOL_HOME/");
    jmolSetAppletColor("#000000");
    jmolApplet(300,"load interface.pdb; select all; color chain; set ambient 40; spin off; $JSMOL_CMD");
    </script>
    </td>
</table>
</ul>

<div style="background:#6599FF;width:150px;">
<b>Results</b></div>
<ul>
<li>Binding affinity change upon mutations [<a href=result.txt>download</a>] [<a href=$HTTP_LINK/help.html#ddG>Explanation</a>]<br>
<font face="Courier">
DDG&nbsp;&nbsp;&nbsp;&nbsp; EvoEFscore mutations<br>
$RESULTS<br>
</font>
<br>

<li>Download structure models of all mutants [<a href=mutation_structures.tar.bz2>download</a>] [<a href=$HTTP_LINK/help.html#forcefield>Explanation</a>]
<br><br>

</ul>
<hr>
<b>Reference:</b>
<ul>$CITATION</ul>
[<a href="$HTTP_LINK">back to SSIPe server</a>]<br/>
</body></html>""")
    # $HTTP_LINK  - HTTP url to SSIPe web portal
    # $JSMOL0     - first javascript for JSmol
    # $JSMOL1     - second javascript for JSmol
    # $JSMOL_HOME - JSmol javascript home
    # $JSMOL_CMD  - JSmol commands for highlighting interface residues
    # $CITATION   - citation for SSIPe paper
    # $JOBID      - webserver jobID
    # $ALIGN      - interface alignment
    # $RESULTS    - result list
    # $ISSCORE    - IS-score cutoff for interface simarity
    # $FF         - interface profile (and physics potential)
    # $CHAIN0     - chain 0
    # $SEQUENCE0  - sequence 0
    # $CHAIN1     - chain 1
    # $SEQUENCE1  - sequence 1

    # template for JSMOL_CMD
    jsmol_cmd_template=Template("select $RESI:$CHAIN; wireframe 75; spacefill 100; color $COLOR; select $RESI:$CHAIN and *.ca; label $RESN$CHAIN$RESI; color label magenta;")
    # $RESI       - residue index
    # $CHAIN      - chain ID
    # $COLOR      - residue color
    # $RESN       - residue amino acid type

    #### read seq.fasta ####
    chain_list=[]
    sequence_list=[]
    fasta_file=os.path.join(data_dir,"seq.fasta")
    if os.path.isfile(fasta_file):
        fp=open(fasta_file,'rU')
        txt=fp.read()
        fp.close()
        for block in txt.split('>'):
            if block.strip():
                block=block.splitlines()
                chain_list.append(block[0])
                sequence_list.append(''.join(block[1:]))
    else:
        header_list,sequence_list=pdb2seq(
            os.path.join(data_dir,"complex.pdb"))
        chain_list=[h.split(':')[-1] for h in header_list]
        fp=open(fasta_file,'w')
        fp.write( ''.join(['>'+h+'\n'+s+'\n' \
            for h,s in zip(header_list,sequence_list)]))
        fp.close()

    #### get all mutated residues ####
    mutation_residue_list=[]
    for mutation_set in result_list:
        for mutation in mutation_set.split()[1].rstrip(';').split(','):
            resi=mutation[2:-1]
            chain=mutation[1]
            mutation_residue_list.append(chain+resi)

    jsmol_cmd=''
    #### write index.html ####
    txt=html_template.substitute(dict(
        HTTP_LINK=http_link,
        JSMOL0=javascript_list[0],
        JSMOL1=javascript_list[1],
        JSMOL_HOME=os.path.dirname(javascript_list[1]),
        JSMOL_CMD=jsmol_cmd,
        CITATION=citation,
        
        JOBID=os.path.basename(data_dir),
        RESULTS="<br>".join(result_list),
        ISSCORE=isscore,
        FF=force_field,
        CHAIN0=chain_list[0],
        SEQUENCE0="<br>".join(textwrap.wrap(sequence_list[0],60)),
        CHAIN1=chain_list[1],
        SEQUENCE1="<br>".join(textwrap.wrap(sequence_list[1],60)),
    ))
    fp=open(os.path.join(data_dir,"index.html"),'w')
    fp.write(txt)
    fp.close()


def get_final_score3(data_dir):
    '''
    4. run EvoEF,      -> evoef_score.txt
    5. get final score -> result.txt
    '''
    #### make temporary directory ####
    tmp_dir="/tmp/"+os.getenv("USER")+"/"+os.path.basename(data_dir)
    if os.path.isdir(tmp_dir):
        shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir)

    #### parse input files ####
    isscore="0.5"
    force_field="EVOEF"

    isscore_file=os.path.join(data_dir,"isscore")
    force_field_file=os.path.join(data_dir,"force_field")
    pdb_file=os.path.join(data_dir,"complex.pdb")
    seq_file=os.path.join(data_dir,"seq.fasta")
    mutlist_file=os.path.join(data_dir,"mutList.txt")
    align_file1=os.path.join(data_dir,"structure_align.out")
    align_file2=os.path.join(data_dir,"sequence_align.out")
    ssipscore_file=os.path.join(data_dir,"ssip_score.txt")
    evoefscore_file=os.path.join(data_dir,"evoef_score.txt")
    result_file=os.path.join(data_dir,"result.txt")

    if os.path.isfile(isscore_file):
        fp=open(isscore_file,'rU')
        isscore=fp.read().strip()
        fp.close()
    if os.path.isfile(force_field_file):
        fp=open(force_field_file,'rU')
        force_field=fp.read().strip()
        fp.close()

    #### define temporary files ####
    pdb_tmp=os.path.join(tmp_dir,"complex.pdb")
    seq_tmp=os.path.join(tmp_dir,"seq.fasta")
    align_tmp1=os.path.join(tmp_dir,"structure_align.out")
    align_tmp2=os.path.join(tmp_dir,"sequence_align.out")
    mutlist_tmp=os.path.join(tmp_dir,"mutList.txt")
    ssipscore_tmp=os.path.join(tmp_dir,"ssip_score.txt")
    evoefscore_tmp=os.path.join(tmp_dir,"evoef_score.txt")

    #### renumber residue index and output the files to temporary folder ####
    #resi_dict,inv_resi_dict=reindex_complex(mutlist_file,pdb_file,mutlist_tmp,pdb_tmp)
    shutil.copy(mutlist_file,mutlist_tmp)
    shutil.copy(pdb_file,pdb_tmp)
    
    #### get evoef score ####
    sys.stdout.write("step 4: get_evoef_score\n")
    if force_field=="EVOEF":
        if os.path.isfile(evoefscore_file):
            shutil.copy(evoefscore_file,evoefscore_tmp)
        else:
            #shutil.copy(pdb_file,pdb_tmp)
            cmd=' '.join(["cd",tmp_dir,';',
                run_evoef,pdb_tmp,mutlist_tmp,evoefscore_tmp])
            os.system(cmd)
            shutil.copy(evoefscore_tmp,evoefscore_file)

    #### get final score ####
    sys.stdout.write("step 5: get_final_score\n")
    fp=open(mutlist_file,'rU')
    mutation_list=fp.read().strip().splitlines()
    fp.close()
    fp=open(evoefscore_tmp,'rU')
    evoefscore_list=[float(f) for f in fp.read().splitlines()]
    fp.close()
    finalscore_list=evoefscore_list
    result_list=["%-7.3f %-10.3f %s"%(f,e,m) for f,e,m in zip(finalscore_list,evoefscore_list,mutation_list)]
    fp=open(result_file,'w')
    fp.write('DDG     EvoEFscore mutations\n')
    fp.write('\n'.join(result_list)+'\n')
    fp.close()

    #### copy mutant PDB structure if available ####
    for mutation in mutation_list:
        mutant_file=os.path.join(data_dir,"complex_"+mutation.rstrip(';').replace(',','_')+".pdb")
        mutant_tmp =os.path.join(tmp_dir,"complex_"+mutation.rstrip(';').replace(',','_')+".pdb")
        if os.path.isfile(mutant_tmp):
            shutil.copy(mutant_tmp,mutant_file)
    cmd=' '.join(["cd", tmp_dir, ";", "tar -jcvf mutation_structures.tar.bz2 complex_*.pdb"])
    os.system(cmd)
    shutil.copy(os.path.join(tmp_dir,"mutation_structures.tar.bz2"),os.path.join(data_dir,"mutation_structures.tar.bz2"))

    #get chain list
    header_list,sequence_list=pdb2seq(pdb_file)
    chain_list=[h.split(':')[-1] for h in header_list]
    #### extract interface PDB file ####
    cmd=' '.join([extint, 
        "-s",pdb_file,
        "-r",chain_list[0],"-l",chain_list[1],
        "-c",os.path.join(data_dir,"contact.list"),
        "-i",os.path.join(data_dir,"interface.pdb")])
    os.system(cmd)
    
    #### make html output ####
    sys.stdout.write("step 6: index.html\n")
    #make_webpage3(inv_resi_dict,data_dir,isscore,force_field,result_list)
    make_webpage3(data_dir,isscore,force_field,result_list)

    #### clean up temporary directory ####
    shutil.rmtree(tmp_dir)
    return

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    for arg in sys.argv[1:]:
        get_final_score3(os.path.abspath(arg))
