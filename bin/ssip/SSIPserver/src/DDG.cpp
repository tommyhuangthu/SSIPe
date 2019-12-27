/*******************************************************************************************************************************
This file is a part of project SSIPe

Copyright (c) 2019 Xiaoqiang Huang (tommyhuangthu@foxmail.com, xiaoqiah@umich.edu)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE 
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
********************************************************************************************************************************/

#pragma warning(disable:4996)
#include "DDG.h"
#include "ErrorHandling.h"
#include "AminoacidOrder.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

SIP::DDG::DDG(){
  ;
}

SIP::DDG::~DDG(){
  ;
}

int SIP::DDG::ddg_initialize(){
  ms.mutant_set_initialize();
  ia_structure.interface_alignment_initialize();
  ia_sequence.interface_alignment_initialize();
  return Success;
}

int SIP::DDG::find_index_for_single_mutation(SIP::SingleMutant *sm){
  int i = -1;
  if(ia_structure.lig_chn_name[0] == sm->chain_name){
    for(int i = 0; i < ia_structure.lig_res_count; i++){
      if(sm->pos == ia_structure.lig_res_index[i]){
        return i;
      }
    }
  }
  else if(ia_structure.rec_chn_name[0] == sm->chain_name){
    for(int i = 0; i < ia_structure.rec_res_count; i++){
      if(sm->pos == ia_structure.rec_res_index[i]){
        return i;
      }
    }
  }

  else if(ia_sequence.lig_chn_name[0] == sm->chain_name){
    for(int i = 0; i < ia_sequence.lig_res_count; i++){
      if(sm->pos == ia_sequence.lig_res_index[i]){
        return i;
      }
    }
  }
  else if(ia_sequence.rec_chn_name[0] == sm->chain_name){
    for(int i = 0; i < ia_sequence.rec_res_count; i++){
      if(sm->pos == ia_sequence.rec_res_index[i]){
        return i;
      }
    }
  }

  return i;
}

///////////////////////////////////////////////////////////////////////////////////////////
// calculate three kinds of pseudo-counts
//////////////////////////////////////////////////////////////////////////////////////////

double SIP::DDG::calc_aat_pseudo_count_for_single_mutation(char aatype, char chain_name, int alnindex, SIP::Parameter *par){
  SIP::AminoacidOrder aao;
  aao.get_aminoacid_order(aatype);
  return par->pc.aat_pseudo_count[aao.aaOrder];
}

double SIP::DDG::calc_gap_pseudo_count_for_single_mutation(char aatype, char chain_name, int alnindex, SIP::Parameter *par){
  double count = 0.0;
  SIP::AminoacidOrder aao;
  aao.get_aminoacid_order(aatype);

  // cannot find the mutation, return 0
  if(alnindex < 0){
    return 0.0;
  }

  if(chain_name == ia_structure.lig_chn_name[0]){
    count = par->pc.gap_pseudo_count[aao.aaOrder] * ia_structure.lig_aligntable[alnindex][20];
  }
  else if(chain_name == ia_structure.rec_chn_name[0]){
    count = par->pc.gap_pseudo_count[aao.aaOrder] * ia_structure.rec_aligntable[alnindex][20];
  }

  return count;
}

double SIP::DDG::calc_evo_pseudo_count_for_single_mutation(char aatype, char chainname, int alnindex, SIP::Parameter *par){
  double count = 0.0;
  SIP::AminoacidOrder aao;
  aao.get_aminoacid_order(aatype);
  
  if(alnindex < 0){
    return 0.0;
  }

  double ntotal = 0.0;
  if(chainname == ia_structure.lig_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_structure.lig_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      return 0.0;
    }
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_structure.lig_aligntable[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
    }
  }
  else if(chainname == ia_structure.rec_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_structure.rec_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      return 0.0;
    }
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_structure.rec_aligntable[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
    }
  }

  return count;
}

///////////////////////////////////////////////////////////////////////////////////////////
// calculate DDG with observed and pseudo counts
//////////////////////////////////////////////////////////////////////////////////////////
int SIP::DDG::calc_ddg_for_single_mutation(SIP::SingleMutant *sm, SIP::Parameter *par){
  SIP::AminoacidOrder nat, mut;
  nat.get_aminoacid_order(sm->native_aa);
  mut.get_aminoacid_order(sm->mutant_aa);
  int aln_index = SIP::DDG::find_index_for_single_mutation(sm);
  if(aln_index < 0){
    sm->smddg = 0.0;
    return Success;
  }
 
  double nat_count = 0.0, mut_count = 0.0;
  nat_count += par->fixed_count;
  nat_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->native_aa, sm->chain_name, aln_index, par);
  nat_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation(sm->native_aa, sm->chain_name, aln_index, par);
  nat_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation(sm->native_aa, sm->chain_name, aln_index, par);
  mut_count += par->fixed_count;
  mut_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->mutant_aa, sm->chain_name, aln_index, par);
  mut_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation(sm->mutant_aa, sm->chain_name, aln_index, par);
  mut_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation(sm->mutant_aa, sm->chain_name, aln_index, par);
  if(nat_count < 0.0 || mut_count < 0.0){
    printf("pseudo counts for mutation ");
    sm->single_mutant_print();
    printf(" is less than zero, that's impossible, please check\n");
    exit(ValueError);
  }
  // add observed counts
  if(ia_structure.lig_chn_name[0] == sm->chain_name){
    double obs_nat=ia_structure.lig_aligntable[aln_index][nat.aaOrder];
    double obs_mut=ia_structure.lig_aligntable[aln_index][mut.aaOrder];
    if(obs_nat-obs_mut>40.0){
      obs_nat=obs_mut+5.0;
    }
    else if(obs_mut-obs_nat>40.0){
      obs_mut=obs_nat+5.0;
    }
    nat_count += obs_nat;
    mut_count += obs_mut;
    //nat_count += ia_from_bindprofx.lig_aligntable[aln_index][nat.aaOrder];
    //mut_count += ia_from_bindprofx.lig_aligntable[aln_index][mut.aaOrder];
  }
  else if(ia_structure.rec_chn_name[0] == sm->chain_name){
    double obs_nat=ia_structure.rec_aligntable[aln_index][nat.aaOrder];
    double obs_mut=ia_structure.rec_aligntable[aln_index][mut.aaOrder];
    if(obs_nat-obs_mut>40.0){
      obs_nat=obs_mut+5.0;
    }
    else if(obs_mut-obs_nat>40.0){
      obs_mut=obs_nat+5.0;
    }
    nat_count += obs_nat;
    mut_count += obs_mut;
    //nat_count += ia_from_bindprofx.rec_aligntable[aln_index][nat.aaOrder];
    //mut_count += ia_from_bindprofx.rec_aligntable[aln_index][mut.aaOrder];
  }
  else{
    sm->smddg = 0.0;
    return IndexError;
  }

  if(nat_count > 1e-3){
    sm->smddg = par->lamda0*(-1.0)*log(mut_count/nat_count);
    //sm->smddg = par->lamda01[nat.aaOrder][mut.aaOrder]*(-1.0)*log(mut_count/nat_count) + par->lamda23[nat.aaOrder][mut.aaOrder];
  }
  return Success;
}

int SIP::DDG::calc_ddg_for_multiple_mutation(SIP::MultipleMutant *mm, SIP::Parameter *par){
  mm->mmddg = 0.0;
  for(int i = 0; i < mm->nmut; i++){
    SIP::SingleMutant *sm = &mm->mutants[i];
    SIP::DDG::calc_ddg_for_single_mutation(sm, par);
    mm->mmddg += sm->smddg;
  }

  return Success;
}


int SIP::DDG::calc_ddg_for_mutation_set(SIP::Parameter *par){
  for(int i = 0; i < ms.nset; i++){
    SIP::MultipleMutant *mm = &ms.mutset[i];
    SIP::DDG::calc_ddg_for_multiple_mutation(mm, par);
  }
  
  return Success;
}

////////////////////////////////////////////////////////////////////////////////////////////
// calculte DDG of multiple-point mutation in a different method
// base on my test, this method is not as good as the decomposition method
////////////////////////////////////////////////////////////////////////////////////////////

int SIP::DDG::calc_ddg_for_multiple_mutation_new(SIP::MultipleMutant *mm, SIP::Parameter *par){
  mm->mmddg = 0.0;
  double nat_count = 0.0, mut_count = 0.0;
  for(int i = 0; i < mm->nmut; i++){
    SIP::SingleMutant *sm = &mm->mutants[i];
    SIP::AminoacidOrder nat, mut;
    nat.get_aminoacid_order(sm->native_aa);
    mut.get_aminoacid_order(sm->mutant_aa);
    int aln_index = SIP::DDG::find_index_for_single_mutation(sm);
    if(aln_index < 0){
      mm->mmddg += 0.0;
      continue;
    }
    
    // calculate pseudo-count for native amino acid
    nat_count += par->fixed_count;
    nat_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->native_aa, sm->chain_name, aln_index, par);
    nat_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation(sm->native_aa, sm->chain_name, aln_index, par);
    nat_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation(sm->native_aa, sm->chain_name, aln_index, par);
    // calculate pseudo-count for mutant amino acid
    mut_count += par->fixed_count;
    mut_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->mutant_aa, sm->chain_name, aln_index, par);
    mut_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation(sm->mutant_aa, sm->chain_name, aln_index, par);
    mut_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation(sm->mutant_aa, sm->chain_name, aln_index, par);
    if(nat_count < 0 || mut_count < 0){
      printf("pseudo counts for mutation ");
      sm->single_mutant_print();
      printf(" is less than zero, that's impossible, please check\n");
      exit(ValueError);
    }
    // add observed counts
    // single mutation is from ligand protein
    if(ia_structure.lig_chn_name[0] == (sm->chain_name)){
      nat_count += ia_structure.lig_aligntable[aln_index][nat.aaOrder];
      mut_count += ia_structure.lig_aligntable[aln_index][mut.aaOrder];
    }
    // single mutation is from receptor protein
    else if(ia_structure.rec_chn_name[0] == sm->chain_name){
      nat_count += ia_structure.rec_aligntable[aln_index][nat.aaOrder];
      mut_count += ia_structure.rec_aligntable[aln_index][mut.aaOrder];
    }
    // chain name for the mutation is wrong
    else{
      printf("Chain name %c is not identified for mutation ", sm->chain_name);
      sm->single_mutant_print();
      printf("please check!\n");
      exit(FormatError);
    }
  }
  if(nat_count > 1e-2){
    mm->mmddg += par->lamda0*(-1.0)*log((double)mut_count/nat_count);
  }

  return Success;
}

int SIP::DDG::calc_ddg_for_mutation_set_new(SIP::Parameter *par){
  for(int i = 0; i < ms.nset; i++){
    SIP::MultipleMutant *mm = &ms.mutset[i];
    SIP::DDG::calc_ddg_for_multiple_mutation_new(mm, par);
  }

  return Success;
}

///////////////////////////////////////////////////////////////////////////////////////////
// write the DDG values to file
///////////////////////////////////////////////////////////////////////////////////////////

int SIP::DDG::write_ddg_for_mutation_set(char *filepath){
  FILE* file = fopen(filepath, "w");
  if(file == NULL){
    printf("**** cannot open file %s for writing, output to screen instead ****\n", filepath);
    file = stdout;
  }

  for(int i = 0; i < ms.nset; i++){
    SIP::MultipleMutant *mm = &ms.mutset[i];
    fprintf(file, "%f\n", mm->mmddg);
  }

  if(file != NULL){
    fclose(file);
  }

  return Success;
}


///////////////////////////////////////////////////////////////////////////////////////////
// the following are new methods to calculate DDG with henikoff weighted AA frequencies
//////////////////////////////////////////////////////////////////////////////////////////
double SIP::DDG::calc_gap_pseudo_count_for_single_mutation_henikoff(char aatype, char chain_name, int alnindex, SIP::Parameter *par){
  double count = 0.0;
  SIP::AminoacidOrder aao;
  aao.get_aminoacid_order(aatype);
  double ntotal = 0.0;

  // cannot find the mutation, return 0
  if(alnindex < 0){
    return count;
  }

  if(chain_name == ia_structure.lig_chn_name[0]){
    count = par->pc.gap_pseudo_count[aao.aaOrder] * ia_structure.eff_ligalign[alnindex][20];
    return count;
  }
  else if(chain_name == ia_structure.rec_chn_name[0]){
    count = par->pc.gap_pseudo_count[aao.aaOrder] * ia_structure.eff_recalign[alnindex][20];
    return count;
  }
  else{
    printf("Chain name %s is not identified for mutation, ", chain_name);
    printf("please check!\n");
    exit(FormatError);
    return count;
  }

  return count;
}

double SIP::DDG::calc_evo_pseudo_count_for_single_mutation_henikoff(char aatype, char chainname, int alnindex, SIP::Parameter *par){
  double count = 0.0;
  SIP::AminoacidOrder aao;
  aao.get_aminoacid_order(aatype);
  double ntotal = 0.0;

  // cannot find the mutation, return 0
  if(alnindex < 0){
    return count;
  }

  if(chainname == ia_structure.lig_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_structure.eff_ligalign[alnindex][i];
    }
    // if ntotal == 0, return 0;
    if(ntotal < 1e-3){
      return 0.0;
    }
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      count += ia_structure.eff_ligalign[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
    }
    count *= par->pc.evo_pseudo_count[aao.aaOrder];
    return count;
  }
  else if(chainname == ia_structure.rec_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_structure.eff_recalign[alnindex][i];
    }
    // if ntotal == 0, return 0;
    if(ntotal < 1e-3){
      return 0.0;
    }
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      count += ia_structure.eff_recalign[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
    }
    count *= par->pc.evo_pseudo_count[aao.aaOrder];
    return count;
  }
  else{
    printf("Chain name %s is not identified for mutation, ", chainname);
    printf("please check!\n");
    exit(FormatError);
    return count;
  }

  return count;
}

int SIP::DDG::calc_ddg_for_single_mutation_henikoff(SIP::SingleMutant *sm, SIP::Parameter *par){
  SIP::AminoacidOrder nat, mut;
  nat.get_aminoacid_order(sm->native_aa);
  mut.get_aminoacid_order(sm->mutant_aa);
  int aln_index = SIP::DDG::find_index_for_single_mutation(sm);
  if(aln_index < 0){
    sm->smddg = 0.0;
    return Success;
  }

  double nat_count = 0.0, mut_count = 0.0;
  // calculate pseudo-count for native amino acid
  nat_count += par->fixed_count;
  nat_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->native_aa, sm->chain_name, aln_index, par);
  nat_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation_henikoff(sm->native_aa, sm->chain_name, aln_index, par);
  nat_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation_henikoff(sm->native_aa, sm->chain_name, aln_index, par);
  // calculate pseudo-count for mutant amino acid
  mut_count += par->fixed_count;
  mut_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->mutant_aa, sm->chain_name, aln_index, par);
  mut_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation_henikoff(sm->mutant_aa, sm->chain_name, aln_index, par);
  mut_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation_henikoff(sm->mutant_aa, sm->chain_name, aln_index, par);
  if(nat_count < 0 || mut_count < 0){
    printf("pseudo counts for mutation ");
    sm->single_mutant_print();
    printf(" is less than zero, that's impossible, please check\n");
    exit(ValueError);
  }
  // add observed counts
  // single mutation is from ligand protein
  if(ia_structure.lig_chn_name[0] == sm->chain_name){
    nat_count += ia_structure.eff_ligalign[aln_index][nat.aaOrder];
    mut_count += ia_structure.eff_ligalign[aln_index][mut.aaOrder];
    if(nat_count > 1e-3){
      sm->smddg = par->lamda0*(-1.0)*log(mut_count/nat_count);
    }
  }
  // single mutation is from receptor protein
  else if(ia_structure.rec_chn_name[0] == sm->chain_name){
    nat_count += ia_structure.eff_recalign[aln_index][nat.aaOrder];
    mut_count += ia_structure.eff_recalign[aln_index][mut.aaOrder];
    if(nat_count > 1e-3){
      sm->smddg = par->lamda0*(-1.0)*log(mut_count/nat_count);
    }
  }
  // chain name for the mutation is wrong
  else{
    printf("Chain name %c is not identified for mutation ", sm->chain_name);
    sm->single_mutant_print();
    printf("please check!\n");
    exit(FormatError);
  }
  return Success;
}

int SIP::DDG::calc_ddg_for_multiple_mutation_henikoff(SIP::MultipleMutant *mm, SIP::Parameter *par){
  mm->mmddg = 0.0;
  for(int i = 0; i < mm->nmut; i++){
    SIP::SingleMutant *sm = &mm->mutants[i];
    SIP::DDG::calc_ddg_for_single_mutation_henikoff(sm, par);
    mm->mmddg += sm->smddg;
  }

  return Success;
}


int SIP::DDG::calc_ddg_for_mutation_set_henikoff(SIP::Parameter *par){
  for(int i = 0; i < ms.nset; i++){
    SIP::MultipleMutant *mm = &ms.mutset[i];
    SIP::DDG::calc_ddg_for_multiple_mutation_henikoff(mm, par);
  }

  return Success;
}

/////////////////////////////////////////////////////////////////////////////////////////////
// calculate DDG using bindprofx and sip combined interface alignment
/////////////////////////////////////////////////////////////////////////////////////////////

double SIP::DDG::calc_gap_pseudo_count_for_single_mutation12(char aatype, char chain_name, int alnindex, SIP::Parameter *par){
  double count = 0.0;
  SIP::AminoacidOrder aao;
  aao.get_aminoacid_order(aatype);

  // cannot find the mutation, return 0
  if(alnindex < 0){
    return 0.0;
  }

  if(chain_name == ia_structure.lig_chn_name[0]){
    count = par->pc.gap_pseudo_count[aao.aaOrder] * ia_structure.lig_aligntable[alnindex][20];
  }
  else if(chain_name == ia_structure.rec_chn_name[0]){
    count = par->pc.gap_pseudo_count[aao.aaOrder] * ia_structure.rec_aligntable[alnindex][20];
  }
  else{
    count = 0.0;
  }

  if(chain_name == ia_sequence.lig_chn_name[0]){
    count += par->pc.gap_pseudo_count[aao.aaOrder] * ia_sequence.lig_aligntable[alnindex][20];
  }
  else if(chain_name == ia_sequence.rec_chn_name[0]){
    count += par->pc.gap_pseudo_count[aao.aaOrder] * ia_sequence.rec_aligntable[alnindex][20];
  }
  else{
    count += 0.0;
  }

  return count;
}

double SIP::DDG::calc_evo_pseudo_count_for_single_mutation12(char aatype, char chainname, int alnindex, SIP::Parameter *par){
  double count = 0.0;
  SIP::AminoacidOrder aao;
  aao.get_aminoacid_order(aatype);
  if(alnindex < 0){
    return 0.0;
  }

  double ntotal = 0.0;
  if(chainname == ia_structure.lig_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_structure.lig_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      count = 0.0;
    }
    else{
      for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
        count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_structure.lig_aligntable[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
      }
    }
  }
  else if(chainname == ia_structure.rec_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_structure.rec_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      count = 0.0;
    }
    else{
      for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
        count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_structure.rec_aligntable[alnindex][i] * par->pc.iPTM[aao.aaOrder][i]/ntotal;
      }
    }
  }

  ntotal = 0.0;
  if(chainname == ia_sequence.lig_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_sequence.lig_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      count += 0.0;
    }
    else{
      for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
        count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_sequence.lig_aligntable[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
      }
    }
  }
  else if(chainname == ia_sequence.rec_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_sequence.rec_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      count += 0.0;
    }
    else{
      for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
        count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_sequence.rec_aligntable[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
      }
    }
  }

  return count;
}

///////////////////////////////////////////////////////////////////////////////////////////
// calculate DDG with observed and pseudo counts
//////////////////////////////////////////////////////////////////////////////////////////

int SIP::DDG::calc_ddg_for_single_mutation12(SIP::SingleMutant *sm, SIP::Parameter *par){
  sm->smddg = 0.0;

  SIP::AminoacidOrder nat, mut;
  nat.get_aminoacid_order(sm->native_aa);
  mut.get_aminoacid_order(sm->mutant_aa);
  int aln_index = SIP::DDG::find_index_for_single_mutation(sm);
  if(aln_index < 0){
    return IndexError;
  }

  double nat_count = 0.0, mut_count = 0.0;
  // calculate pseudo-count for native amino acid
  nat_count += par->fixed_count;
  nat_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->native_aa, sm->chain_name, aln_index, par);
  nat_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation12(sm->native_aa, sm->chain_name, aln_index, par);
  nat_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation12(sm->native_aa, sm->chain_name, aln_index, par);
  mut_count += par->fixed_count;
  mut_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->mutant_aa, sm->chain_name, aln_index, par);
  mut_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation12(sm->mutant_aa, sm->chain_name, aln_index, par);
  mut_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation12(sm->mutant_aa, sm->chain_name, aln_index, par);
  if(nat_count < 0 || mut_count < 0){
    printf("pseudo counts for mutation ");
    sm->single_mutant_print();
    printf(" is less than zero, that's impossible, please check\n");
    exit(ValueError);
  }
  // add observed counts
  if(ia_structure.lig_chn_name[0] == sm->chain_name){
    double obs_nat=ia_structure.lig_aligntable[aln_index][nat.aaOrder];
    double obs_mut=ia_structure.lig_aligntable[aln_index][mut.aaOrder];
    if(obs_nat-obs_mut>40.0){
      obs_nat=obs_mut+5.0;
    }
    else if(obs_mut-obs_nat>40.0){
      obs_mut=obs_nat+5.0;
    }
    nat_count += obs_nat;
    mut_count += obs_mut;
    //nat_count += ia_from_bindprofx.lig_aligntable[aln_index][nat.aaOrder];
    //mut_count += ia_from_bindprofx.lig_aligntable[aln_index][mut.aaOrder];
  }
  else if(ia_structure.rec_chn_name[0] == sm->chain_name){
    double obs_nat=ia_structure.rec_aligntable[aln_index][nat.aaOrder];
    double obs_mut=ia_structure.rec_aligntable[aln_index][mut.aaOrder];
    if(obs_nat-obs_mut>40.0){
      obs_nat=obs_mut+5.0;
    }
    else if(obs_mut-obs_nat>40.0){
      obs_mut=obs_nat+5.0;
    }
    nat_count += obs_nat;
    mut_count += obs_mut;
    //nat_count += ia_from_bindprofx.rec_aligntable[aln_index][nat.aaOrder];
    //mut_count += ia_from_bindprofx.rec_aligntable[aln_index][mut.aaOrder];
  }
  else{
    sm->smddg = 0.0;
    return IndexError;
  }
  //printf("before sip, nat_count & mut_count: %.12f, %.12f\n", nat_count, mut_count);

  if(ia_sequence.lig_chn_name[0] == sm->chain_name){
    nat_count += ia_sequence.lig_aligntable[aln_index][nat.aaOrder];
    mut_count += ia_sequence.lig_aligntable[aln_index][mut.aaOrder];
  }
  else if(ia_sequence.rec_chn_name[0] == sm->chain_name){
    nat_count += ia_sequence.rec_aligntable[aln_index][nat.aaOrder];
    mut_count += ia_sequence.rec_aligntable[aln_index][mut.aaOrder];
  }
  else{
    sm->smddg = 0.0;
    return IndexError;
  }
  //printf("after sip, nat_count & mut_count: %.12f, %.12f\n\n", nat_count, mut_count);

  if(nat_count > 1e-3){
    sm->smddg = par->lamda0*(-1.0)*log(mut_count/nat_count);
  }

  return Success;
}

int SIP::DDG::calc_ddg_for_multiple_mutation12(SIP::MultipleMutant *mm, SIP::Parameter *par){
  mm->mmddg = 0.0;
  for(int i = 0; i < mm->nmut; i++){
    SIP::SingleMutant *sm = &mm->mutants[i];
    SIP::DDG::calc_ddg_for_single_mutation12(sm, par);
    mm->mmddg += sm->smddg;
  }

  return Success;
}


int SIP::DDG::calc_ddg_for_mutation_set12(SIP::Parameter *par){
  for(int i = 0; i < ms.nset; i++){
    SIP::MultipleMutant *mm = &ms.mutset[i];
    SIP::DDG::calc_ddg_for_multiple_mutation12(mm, par);
  }

  return Success;
}



double SIP::DDG::calc_gap_pseudo_count_for_single_mutation1234(char aatype, char chain_name, int alnindex, SIP::Parameter *par){
  double count = 0.0;
  SIP::AminoacidOrder aao;
  aao.get_aminoacid_order(aatype);

  // cannot find the mutation, return 0
  if(alnindex < 0){
    return 0.0;
  }

  if(chain_name == ia_structure.lig_chn_name[0]){
    count = par->pc.gap_pseudo_count[aao.aaOrder] * ia_structure.lig_aligntable[alnindex][20];
  }
  else if(chain_name == ia_structure.rec_chn_name[0]){
    count = par->pc.gap_pseudo_count[aao.aaOrder] * ia_structure.rec_aligntable[alnindex][20];
  }
  else{
    count = 0.0;
  }

  if(chain_name == ia_sequence.lig_chn_name[0]){
    count += par->pc.gap_pseudo_count[aao.aaOrder] * ia_sequence.lig_aligntable[alnindex][20];
  }
  else if(chain_name == ia_sequence.rec_chn_name[0]){
    count += par->pc.gap_pseudo_count[aao.aaOrder] * ia_sequence.rec_aligntable[alnindex][20];
  }
  else{
    count += 0.0;
  }

  return count;
}

double SIP::DDG::calc_evo_pseudo_count_for_single_mutation1234(char aatype, char chainname, int alnindex, SIP::Parameter *par){
  double count = 0.0;
  SIP::AminoacidOrder aao;
  aao.get_aminoacid_order(aatype);
  if(alnindex < 0){
    return 0.0;
  }

  double ntotal = 0.0;
  if(chainname == ia_structure.lig_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_structure.lig_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      count = 0.0;
    }
    else{
      for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
        count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_structure.lig_aligntable[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
      }
    }
  }
  else if(chainname == ia_structure.rec_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_structure.rec_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      count = 0.0;
    }
    else{
      for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
        count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_structure.rec_aligntable[alnindex][i] * par->pc.iPTM[aao.aaOrder][i]/ntotal;
      }
    }
  }

  ntotal = 0.0;
  if(chainname == ia_sequence.lig_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_sequence.lig_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      count += 0.0;
    }
    else{
      for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
        count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_sequence.lig_aligntable[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
      }
    }
  }
  else if(chainname == ia_sequence.rec_chn_name[0]){
    for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
      ntotal += ia_sequence.rec_aligntable[alnindex][i];
    }
    if(ntotal < 1e-3){
      count += 0.0;
    }
    else{
      for(int i = 0; i < MAX_NUM_AA_TYPE-1; i++){
        count += par->pc.evo_pseudo_count[aao.aaOrder] * ia_sequence.rec_aligntable[alnindex][i]*par->pc.iPTM[aao.aaOrder][i]/ntotal;
      }
    }
  }

  return count;
}






int SIP::DDG::calc_ddg_for_single_mutation1234(SIP::SingleMutant *sm, SIP::Parameter *par){
  sm->smddg = 0.0;

  SIP::AminoacidOrder nat, mut;
  nat.get_aminoacid_order(sm->native_aa);
  mut.get_aminoacid_order(sm->mutant_aa);
  int aln_index = SIP::DDG::find_index_for_single_mutation(sm);
  if(aln_index < 0){
    return IndexError;
  }

  double nat_count = 0.0, mut_count = 0.0;
  // calculate pseudo-count for native amino acid
  nat_count += par->fixed_count;
  nat_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->native_aa, sm->chain_name, aln_index, par);
  nat_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation1234(sm->native_aa, sm->chain_name, aln_index, par);
  nat_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation1234(sm->native_aa, sm->chain_name, aln_index, par);
  mut_count += par->fixed_count;
  mut_count += SIP::DDG::calc_aat_pseudo_count_for_single_mutation(sm->mutant_aa, sm->chain_name, aln_index, par);
  mut_count += SIP::DDG::calc_gap_pseudo_count_for_single_mutation1234(sm->mutant_aa, sm->chain_name, aln_index, par);
  mut_count += SIP::DDG::calc_evo_pseudo_count_for_single_mutation1234(sm->mutant_aa, sm->chain_name, aln_index, par);
  if(nat_count < 0 || mut_count < 0){
    printf("pseudo counts for mutation ");
    sm->single_mutant_print();
    printf(" is less than zero, that's impossible, please check\n");
    exit(ValueError);
  }
  // add observed counts
  double MAX_COUNT_DIFF = 40;
  if(ia_structure.lig_chn_name[0] == sm->chain_name){
    double obs_nat=ia_structure.lig_aligntable[aln_index][nat.aaOrder]*par->ial_weight;
    double obs_mut=ia_structure.lig_aligntable[aln_index][mut.aaOrder]*par->ial_weight;
    if(obs_nat-obs_mut>MAX_COUNT_DIFF){
      //printf("In single mutant ");
      //sm->single_mutant_print();
      //printf(", observed count for native is: %.0f, for mutant is: %.0f\n", obs_nat, obs_mut);
      obs_nat=obs_mut+5.0;
    }
    else if(obs_mut-obs_nat>MAX_COUNT_DIFF){
      //printf("In single mutant ");
      //sm->single_mutant_print();
      //printf(", observed count for native is: %.0f, for mutant is: %.0f\n", obs_nat, obs_mut);
      obs_mut=obs_nat+5.0;
    }
    nat_count += obs_nat;
    mut_count += obs_mut;
    //nat_count += ia_from_bindprofx.lig_aligntable[aln_index][nat.aaOrder];
    //mut_count += ia_from_bindprofx.lig_aligntable[aln_index][mut.aaOrder];
  }
  else if(ia_structure.rec_chn_name[0] == sm->chain_name){
    double obs_nat=ia_structure.rec_aligntable[aln_index][nat.aaOrder]*par->ial_weight;
    double obs_mut=ia_structure.rec_aligntable[aln_index][mut.aaOrder]*par->ial_weight;
    if(obs_nat-obs_mut>MAX_COUNT_DIFF){
      //printf("In single mutant ");
      //sm->single_mutant_print();
      //printf(", observed count for native is: %.0f, for mutant is: %.0f\n", obs_nat, obs_mut);
      obs_nat=obs_mut+5.0;
    }
    else if(obs_mut-obs_nat>MAX_COUNT_DIFF){
      //printf("In single mutant ");
      //sm->single_mutant_print();
      //printf(", observed count for native is: %.0f, for mutant is: %.0f\n", obs_nat, obs_mut);
      obs_mut=obs_nat+5.0;
    }
    nat_count += obs_nat;
    mut_count += obs_mut;
    //nat_count += ia_from_bindprofx.rec_aligntable[aln_index][nat.aaOrder];
    //mut_count += ia_from_bindprofx.rec_aligntable[aln_index][mut.aaOrder];
  }
  else{
    sm->smddg = 0.0;
    return IndexError;
  }

  if(ia_sequence.lig_chn_name[0] == sm->chain_name){
    nat_count += ia_sequence.lig_aligntable[aln_index][nat.aaOrder]*par->psi_weight;
    mut_count += ia_sequence.lig_aligntable[aln_index][mut.aaOrder]*par->psi_weight;
  }
  else if(ia_sequence.rec_chn_name[0] == sm->chain_name){
    nat_count += ia_sequence.rec_aligntable[aln_index][nat.aaOrder]*par->psi_weight;
    mut_count += ia_sequence.rec_aligntable[aln_index][mut.aaOrder]*par->psi_weight;
  }
  else{
    nat_count += 0.0;
    mut_count += 0.0;
  }

  if(nat_count > 1e-3){
    sm->smddg = par->lamda0*(-1.0)*log(mut_count/nat_count);
  }

  return Success;
}


int SIP::DDG::calc_ddg_for_multiple_mutation1234(SIP::MultipleMutant *mm, SIP::Parameter *par){
  mm->mmddg = 0.0;
  for(int i = 0; i < mm->nmut; i++){
    SIP::SingleMutant *sm = &mm->mutants[i];
    SIP::DDG::calc_ddg_for_single_mutation1234(sm, par);
    mm->mmddg += sm->smddg;
  }

  return Success;
}


int SIP::DDG::calc_ddg_for_mutation_set1234(SIP::Parameter *par){
  for(int i = 0; i < ms.nset; i++){
    SIP::MultipleMutant *mm = &ms.mutset[i];
    SIP::DDG::calc_ddg_for_multiple_mutation1234(mm, par);
  }

  return Success;
}