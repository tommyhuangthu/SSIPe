#!/usr/bin/env python
docstring='''
get_final_score.py data_dir
    main script to perform SSIPe prediction

Input files:
    complex.pdb (complex structure)
    mutList.txt (list of mutations to make)
    force_field (optional: SSIPE - (default) profile score + evoef
                           SSIP  - profile score only
                           EVOEF - evoef score only)
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

#def make_webpage(inv_resi_dict,data_dir,isscore,force_field,result_list):
def make_webpage(data_dir,isscore,force_field,result_list):
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
    jmolApplet(300,"load complex.pdb; set frank off; calculate structure rama; select all; color chain; spacefill off; wireframe off; cartoons; set ambient 40; set spin x 10; set spin y 10; spin off; $JSMOL_CMD");
    </script>
    </td>

    <td align ="center">
    interface structure [<a href=interface.pdb>download</a>]
    <script type="text/javascript">
    jmolInitialize("$JSMOL_HOME/");
    jmolSetAppletColor("#000000");
    jmolApplet(300,"load interface.pdb; select all; color chain; set ambient 40; set spin x 10; set spin y 10; spin off; $JSMOL_CMD");
    </script>
    </td>
</table>
</ul>

<div style="background:#6599FF;width:150px;">
<b>Results</b></div>
<ul>
<li>Binding affinity change upon mutations [<a href=result.txt>download</a>] [<a href=$HTTP_LINK/help.html#ddG>Explanation</a>]<br>
<font face="Courier">
DDG&nbsp;&nbsp;&nbsp;&nbsp; SSIPscore EvoEFscore mutations<br>
$RESULTS<br>
</font>
<br>

<li>Download structure models of all mutants [<a href=mutation_structures.tar.bz2>download</a>]
[<a href=$HTTP_LINK/help.html#forcefield>Explanation</a>]<br>
<br><br>

<li>Structure-based interface MSA by iAlign from NIL [<a href=structure_align.out>download</a>]
[<a href=$HTTP_LINK/help.html#structure_alignment>Explanation</a>]<br>
<font face="Courier">$ALIGN</font><br><br>

<li>Sequence-based interface MSA by PSI-BLAST from STRING [<a href=sequence_align.out>download</a>]
[<a href=$HTTP_LINK/help.html#sequence_alignment>Explanation</a>]<br>
<font face="Courier">$ALIGN2</font><br><br>

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

    #### read structure_align.out ####
    fp=open(os.path.join(data_dir,"structure_align.out"),'rU')
    txt=fp.read().replace('query','target')
    fp.close()
    align_txt=''
    jsmol_cmd=''
    for line in txt.splitlines():
        align_txt+=line+'<br>'
    #don't label the mutated residues to make webpage faster
    '''for line in txt.splitlines():
        if line.startswith("chain_ID_A:"):
            chain0=line[len("chain_ID_A:"):].strip()
        elif line.startswith("chain_ID_B:"):
            chain1=line[len("chain_ID_B:"):].strip()
        elif line.startswith("residue_A:"):
            txt="residue_A: "
            for i in line[len("residue_A:"):].strip().split():
                resi=inv_resi_dict[chain0][i] # residue index
                txt+=' '+resi
                if chain0+resi in mutation_residue_list:
                    resn=sequence_list[0][int(i)-1]
                    jsmol_cmd+=jsmol_cmd_template.substitute(dict(
                        RESI=resi,CHAIN=chain0,COLOR="blue",RESN=resn))
            line=txt
        elif line.startswith("residue_B:"):
            txt="residue_B: "
            for i in line[len("residue_B:"):].strip().split():
                resi=inv_resi_dict[chain1][i] # residue index
                txt+=' '+resi
                if chain1+resi in mutation_residue_list:
                    resn=sequence_list[1][int(i)-1]
                    jsmol_cmd+=jsmol_cmd_template.substitute(dict(
                        RESI=resi,CHAIN=chain1,COLOR="yellow",RESN=resn))
            line=txt
        align_txt+=line+'<br>'
    '''
    #### read sequence_align.out ####
    fp=open(os.path.join(data_dir,"sequence_align.out"),'rU')
    txt=fp.read().replace('query','target')
    fp.close()
    align_txt2=''
    jsmol_cmd=''
    for line in txt.splitlines():
        align_txt2+=line+'<br>'

    #### write index.html ####
    txt=html_template.substitute(dict(
        HTTP_LINK=http_link,
        JSMOL0=javascript_list[0],
        JSMOL1=javascript_list[1],
        JSMOL_HOME=os.path.dirname(javascript_list[1]),
        JSMOL_CMD=jsmol_cmd,
        CITATION=citation,
        
        JOBID=os.path.basename(data_dir),
        ALIGN=align_txt,
        ALIGN2=align_txt2,
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


def get_final_score(data_dir):
    '''
    1. run iAlign,  -> structure_align.out 
    2. run PSI-BLAST,  -> sequence_align.out
    3. get SSIP score, -> ssip_score.txt
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
    force_field="SSIPE" # SSIPe force field by default

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

    #### interface alignment ####
    sys.stdout.write("step 1: run iAlign\n")
    if os.path.isfile(align_file1):
        shutil.copy(align_file1,align_tmp1)
    else:
        newscore = max(float(isscore),0.45)
        newscore = min(newscore, 0.55)
        isscore = str(newscore)
        cmd=' '.join(["cd",tmp_dir,';',
            run_align,pdb_tmp,str(newscore),align_tmp1])
        os.system(cmd)
        align_file1=os.path.join(data_dir,"structure_align.out")
        shutil.copy(align_tmp1,align_file1)

    #### fix insertion code in align.out ####
    fp=open(align_tmp1,'rU')
    lines=fp.read().splitlines()
    fp.close()
    txt=''
    insert_code_pattern=re.compile("\d+[A-Za-z]")
    for line in lines:
        if not line.startswith("residue_") or not insert_code_pattern.search(line):
            txt+=line+'\n'
        else:
            for resi in line.split(' '):
                if resi=='':
                    txt+=' '
                elif not insert_code_pattern.match(resi):
                    txt+=resi
                else:
                    txt+=' '+resi[:-1]
            txt+='\n'
    fp=open(align_tmp1,'w')
    fp.write(txt)
    fp.close()

    #### run psiblast to get sequence interface alignment ####
    sys.stdout.write("step 2: run psiblast\n")
    if(os.path.isfile(align_tmp1)):
      if(os.path.isfile(align_file2)):
        shutil.copy(align_file2,align_tmp2)
      else:
        cmd=' '.join(["cd", tmp_dir, ";", run_seqalign, pdb_tmp,align_tmp1, align_tmp2])
        os.system(cmd)
        align_file2=os.path.join(data_dir,"sequence_align.out")
        shutil.copy(align_tmp2,align_file2)
    else:
      sys.stdout.write("step 2: skip psiblast because iAlign alignment file doesn't exist")

    #### get SSIP score ####
    sys.stdout.write("step 3: get_ssip_score\n")
    if(os.path.isfile(align_tmp1) or os.path.isfile(align_tmp2)):
      cmd=' '.join(["cd",tmp_dir,';',get_ssip_score,mutlist_tmp,align_tmp1,align_tmp2,isscore,ssipscore_tmp])
      os.system(cmd)
      shutil.copy(ssipscore_tmp,ssipscore_file)
      #copy the sorted file back to the data_dir
      shutil.copy(align_tmp1,align_file1)
      shutil.copy(align_tmp2,align_file2)
    else:
      sys.stdout.write("step 3: skip get_ssip_score because structure_align.out and sequence_align.out don't exist")

    #### get evoef score ####
    sys.stdout.write("step 4: get_evoef_score\n")
    if force_field=="SSIPE":
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
    fp=open(ssipscore_tmp,'rU')
    ssipscore_list=[float(l) for l in fp.read().splitlines()]
    fp.close()
    finalscore_list = []
    if os.path.isfile(evoefscore_tmp):
        fp=open(evoefscore_tmp,'rU')
        evoefscore_list=[float(f) for f in fp.read().splitlines()]
        fp.close()
        for x,f,m in zip(ssipscore_list,evoefscore_list,mutation_list):
            if x==0. and f==0.:
                finalscore_list.append(0.)
            else:
                finalscore_list.append(0.734*x+0.341*f+0.205)
        #finalscore_list=[0.74*x+0.33*f+0.26 for x,f,m in \
        #    zip(ssipscore_list,evoefscore_list,mutation_list)]
    else:
        finalscore_list=ssipscore_list
    result_list=["%-7.3f %-9.3f %-10.3f %s"%(f,s,e,m) for f,s,e,m in zip(finalscore_list,ssipscore_list,evoefscore_list,mutation_list)]
    fp=open(result_file,'w')
    fp.write('DDG     SSIPscore EvoEFscore mutations\n')
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
    #make_webpage(inv_resi_dict,data_dir,isscore,force_field,result_list)
    make_webpage(data_dir,isscore,force_field,result_list)

    #### clean up temporary directory ####
    shutil.rmtree(tmp_dir)
    return

if __name__=="__main__":
    if len(sys.argv)<2:
        sys.stderr.write(docstring)
        exit()

    for arg in sys.argv[1:]:
        get_final_score(os.path.abspath(arg))
