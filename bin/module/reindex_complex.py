#!/usr/bin/env python
docstring='''
reindex_complex.py mutList.txt complex.pdb reindex.txt reindex.pdb
    Reindex residue number in BindProf format mutation file "mutList.txt"
    and dimer PDB file "complex.pdb" so that the first residue start from 
    one. if the structure contains any missing residue, residue number gaps
    will be removed
'''
import Bio.PDB
import sys,os
import re

class NonHetSelect(Bio.PDB.Select): # class to select ATOM entries
    def accept_residue(self,residue): # remove heteroatoms
        return 1 if residue.id[0]==' ' else 0
    def accept_atom(self, atom): # remove alternative location
        return not atom.is_disordered() or atom.get_altloc() == 'A'

def reindex_complex(mutList_in,pdb_in,mutList_out,pdb_out):
    '''
    Read PDB file "pdb_in", reindex residue number of both chains so that
    residue number starts from 1. Make corresponding change to BindProf format
    mutation file "mutList_in" and make corresponding change.

    Output reindex mutation file and PDB file to "mutList_out" and "pdb_out"

    return two dict, first has old residue index as key & new residue index as
    value. second has new residue index as key & old residue index as value
    '''
    #### reindex complex pdb ####
    resi_dict=dict()
    inv_resi_dict=dict()
    struct = Bio.PDB.PDBParser(PERMISSIVE=1).get_structure("complex",pdb_in)
    model=struct[0]
    for chain in model:
        resi_dict[chain.id]=dict()
        inv_resi_dict[chain.id]=dict()
        index=1
        for residue in chain:
            if residue.id[0]==' ':
                resi_dict[chain.id][str(residue.id[1])]=str(index)
                inv_resi_dict[chain.id][str(index)]=str(residue.id[1])
                residue.id=(residue.id[0],index,residue.id[-1])
                index+=1

    io=Bio.PDB.PDBIO()
    io.set_structure(model)
    io.save(pdb_out,NonHetSelect())

    #### reindex mutation file ####
    fp=open(mutList_in,'rU')
    mutation_list=[l.strip() for l in \
        fp.read().replace(';','').splitlines() if l.strip()]
    fp.close()

    mutation_pattern=re.compile("^(\w)(\w)(\d+)(\w)")
    mutation_txt=''
    for mutation_set in mutation_list:
        for mutation in mutation_set.split(','):
            wt_res,chain,resi,mut_res=mutation_pattern.findall(mutation)[0]
            mutation_txt+=wt_res+chain+resi_dict[chain][resi]+mut_res+','
        mutation_txt=mutation_txt[:-1]+';\n' # remove trailing comma
    fp=open(mutList_out,'w')
    fp.write(mutation_txt)
    fp.close()
    return resi_dict,inv_resi_dict
    
if __name__=="__main__":
    if len(sys.argv)<5:
        print >>sys.stderr,docstring
        exit()

    reindex_complex(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
