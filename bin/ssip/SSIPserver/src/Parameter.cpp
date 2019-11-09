#pragma warning(disable:4996)
#include "Parameter.h"
#include "ErrorHandling.h"
#include "DDG.h"

#include <stdio.h>
#include <stdlib.h>
//#include <io.h>
#include <string.h>
#include <ctype.h>

SIP::Parameter::Parameter(){
  ;
}

SIP::Parameter::~Parameter(){
  ;
}


int SIP::Parameter::parameter_print(){
  printf("amino-acid_types                 A     C     D     E     F     G     H     I     K     L     M     N     P     Q     R     S     T     V     W     Y\n");
  printf("fix_dependent_pseudo-counts  %5d\n", fixed_count);
  printf("aas_dependent_pseudo-counts ");
  for(int i = 0; i < 20; i++){
    printf(" %5d", pc.aat_pseudo_count[i]);
  }
  printf("\n");
  printf("gap_dependent_pseudo-counts ");
  for(int i = 0; i < 20; i++){
    printf(" %5d", pc.gap_pseudo_count[i]);
  }
  printf("\n");
  printf("evo_dependent_pseudo-counts ");
  for(int i = 0; i < 20; i++){
    printf(" %5d", pc.evo_pseudo_count[i]);
  }
  printf("\n");
  printf("aas_dependent_lamda0        ");
  printf(" %5.2f", lamda0);
  printf("\n");
  printf("ial_weight           %.2f\n", ial_weight);
  printf("ial_isscore_cutoff   %.2f\n", isscore_cutoff);
  printf("ial_count_cutoff     %d\n",ial_count_cutoff);
  printf("psi_weight           %.2f\n",psi_weight);
  printf("psi_linkscore_cutoff %.2f\n",psi_linkscore_cutoff);
  printf("psi_seqidhigh_cutoff %.2f\n",seqid_high_cutoff);
  printf("psi_seqidlow_cutoff  %.2f\n",seqid_low_cutoff);
  printf("psi_count_cutoff     %d\n",psi_count_cutoff);
  
  return Success;
}

int SIP::Parameter::parameter_copy(Parameter *newpar){
  fixed_count = newpar->fixed_count;
  pc.pseudo_counts_copy(&newpar->pc);
  lamda0 = newpar->lamda0;
  
  //ialign
  ial_weight = newpar->ial_weight;
  isscore_cutoff = newpar->isscore_cutoff;
  ial_count_cutoff = newpar->ial_count_cutoff;
  //psiblast
  psi_weight = newpar->psi_weight;
  psi_linkscore_cutoff = newpar->psi_linkscore_cutoff;
  seqid_high_cutoff = newpar->seqid_high_cutoff;
  seqid_low_cutoff = newpar->seqid_low_cutoff;
  psi_count_cutoff = newpar->psi_count_cutoff;

  return Success;
}

int SIP::Parameter::parameter_randomize_pseudo_count(int nmut){
  int MAX_COUNT = 50;
  for(int i = 0; i < nmut; i++){
    int index = rand()%60+1;
    // 20 parameters
    if(index >= 1 && index <= 20){
      pc.set_aat_pseudo_count(index-1, rand()%(MAX_COUNT+1));
    }
    // 20 parameters
    else if(index >= 21 && index <= 40){
      pc.set_gap_pseudo_count(index-21, rand()%(MAX_COUNT+1));
    }
    // 20 parameters
    else if(index >= 41 && index <= 60){
      pc.set_evo_pseudo_count(index-41, rand()%(MAX_COUNT+1));
    }
  }  

  return Success;
}


int SIP::Parameter::parameter_randomize_lamda(int nmut){
  int MAX_COUNT = 50;
  for(int i = 0; i < nmut; i++){
    int index = rand()%2;
    lamda0 = (rand()+1.0)/(RAND_MAX+1.0)*20.0;
  }

  return Success;
}


int SIP::Parameter::parameter_initialize(){
  fixed_count = 1;
  for(int i = 0; i < 20; i++){
    pc.aat_pseudo_count[i] = 24;
    pc.gap_pseudo_count[i] = 15;
    pc.evo_pseudo_count[i] = 5;
  }
  lamda0 = 15.0;
  ial_weight = 1.0;
  psi_weight = 0.0;
  return Success;
}


int SIP::Parameter::parameter_initialize_random(){
  fixed_count = 1;
  pc.set_random_pseudo_count();
  lamda0 = 1.0;
  ial_weight = 1.0;
  psi_weight = 0.0;
  return Success;
}


int SIP::Parameter::parameter_add(Parameter *newpar){
  for(int i = 0; i < 20; i++){
    pc.aat_pseudo_count[i] += newpar->pc.aat_pseudo_count[i];
    pc.evo_pseudo_count[i] += newpar->pc.evo_pseudo_count[i];
    pc.gap_pseudo_count[i] += newpar->pc.gap_pseudo_count[i];
  }
  lamda0 += newpar->lamda0;
  fixed_count += newpar->fixed_count;
  return Success;
}

int SIP::Parameter::parameter_average(int count){
  if(count <= 0){
    printf("the number used to average parameters cannot be <= 0, please check\n");
    exit(ValueError);
  }

  for(int i = 0; i < 20; i++){
    pc.aat_pseudo_count[i] /= count;
    pc.evo_pseudo_count[i] /= count;
    pc.gap_pseudo_count[i] /= count;
  }
  lamda0 /= count;
  fixed_count /= count;
  return Success;
}

int SIP::Parameter::parameter_read(char* parameterfile){
  FileReader fr;
  FileReaderCreate(&fr, parameterfile);
  int counter = 0;
  for(int j = 0; j < StringArrayGetCount(&fr.lines); j++){
    char* line = StringArrayGet(&fr.lines, j);
    StringArray columns;
    StringArrayCreate(&columns);
    StringArraySplitString(&columns, line, ' ');
    if(!strcmp(StringArrayGet(&columns, 0), "amino-acid_types")){
      ;
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "fix_dependent_pseudo-counts")){
      fixed_count = atoi(StringArrayGet(&columns, 1));
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "aas_dependent_pseudo-counts")){
      for(int i = 0; i < 20; i++){
        pc.aat_pseudo_count[i] = atoi(StringArrayGet(&columns, i+1));
      }
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "gap_dependent_pseudo-counts")){
      for(int i = 0; i < 20; i++){
        pc.gap_pseudo_count[i] = atoi(StringArrayGet(&columns, i+1));
      }
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "evo_dependent_pseudo-counts")){
      for(int i = 0; i < 20; i++){
        pc.evo_pseudo_count[i] = atoi(StringArrayGet(&columns, i+1));
      }
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "aas_dependent_lamda0")){
      lamda0 = atof(StringArrayGet(&columns, 1));
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "ial_weight")){
      ial_weight = atof(StringArrayGet(&columns, 1));
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "ial_isscore_cutoff")){
      isscore_cutoff = atof(StringArrayGet(&columns, 1));
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "ial_count_cutoff")){
      ial_count_cutoff = atoi(StringArrayGet(&columns, 1));
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "psi_weight")){
      psi_weight = atof(StringArrayGet(&columns, 1));
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "psi_linkscore_cutoff")){
      psi_linkscore_cutoff = atof(StringArrayGet(&columns, 1));
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "psi_seqidhigh_cutoff")){
      seqid_high_cutoff = atof(StringArrayGet(&columns, 1));
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "psi_seqidlow_cutoff")){
      seqid_low_cutoff = atof(StringArrayGet(&columns, 1));
    }
    else if(!strcmp(StringArrayGet(&columns, 0), "psi_count_cutoff")){
      psi_count_cutoff = atoi(StringArrayGet(&columns, 1));
    }
    else{
      ;
    }
    StringArrayDestroy(&columns);
  }
  return Success;
}


int SIP::Parameter::parameter_randomize_21par_model_aat(int nmut){
  if(nmut==0) return Success;
  int MAX_COUNT = 50;
  for(int i=0; i<20; i++){
    pc.set_gap_pseudo_count(i, 0);
    pc.set_evo_pseudo_count(i, 0);
  }
  for(int i = 0; i < nmut; i++){
    int index = rand()%25;
    // 20 parameters
    if(index >= 0 && index <= 19){
      pc.set_aat_pseudo_count(index, rand()%(MAX_COUNT+1));
    }
    // 5 parameters
    else if(index >= 20 && index <= 24){
      lamda0 = (rand()+1.0)/(RAND_MAX+1.0)*20.0;
    }
  }

  return Success;
}

int SIP::Parameter::parameter_randomize_41par_model_gap(int nmut){
  if(nmut==0) return Success;
  int MAX_COUNT = 50;
  for(int i=0; i<20; i++){
    pc.set_evo_pseudo_count(i, 0);
  }
  for(int i = 0; i < nmut; i++){
    int index = rand()%45;
    // 20 parameters
    if(index >= 0 && index <= 19){
      pc.set_aat_pseudo_count(index, rand()%(MAX_COUNT+1));
    }
    // 20 parameters
    else if(index >= 20 && index <= 39){
      pc.set_gap_pseudo_count(index-20, rand()%(MAX_COUNT+1));
    }
    // 10 parameters
    else if(index >= 40 && index <= 44){
      lamda0 = (rand()+1.0)/(RAND_MAX+1.0)*20.0;
    }
  }

  return Success;
}

int SIP::Parameter::parameter_randomize_41par_model_evo(int nmut){
  if(nmut==0) return Success;
  int MAX_COUNT = 50;
  for(int i=0; i<20; i++){
    pc.set_gap_pseudo_count(i, 0);
  }
  for(int i = 0; i < nmut; i++){
    int index = rand()%45;
    if(index >= 0 && index <= 19){
      pc.set_aat_pseudo_count(index, rand()%(MAX_COUNT+1));
    }
    else if(index >= 20 && index <= 39){
      pc.set_evo_pseudo_count(index-20, rand()%(MAX_COUNT+1));
    }
    else if(index >= 40 && index <= 44){
      lamda0 = (rand()+1.0)/(RAND_MAX+1.0)*20.0;
    }
  }

  return Success;
}

int SIP::Parameter::parameter_randomize_61par_model(int nmut){
  if(nmut==0) return Success;
  int MAX_COUNT = 50;
  for(int i = 0; i < nmut; i++){
    int index = rand()%65;
    // 20 parameters
    if(index >= 0 && index <= 19){
      pc.set_aat_pseudo_count(index, rand()%(MAX_COUNT+1));
    }
    // 20 parameters
    else if(index >= 20 && index <= 39){
      pc.set_gap_pseudo_count(index-20, rand()%(MAX_COUNT+1));
    }
    // 20 parameters
    else if(index >= 40 && index <= 59){
      pc.set_evo_pseudo_count(index-40, rand()%(MAX_COUNT+1));
    }
    // 10 parameters
    else if(index >= 60 && index <= 64){
      lamda0 = (rand()+1.0)/(RAND_MAX+1.0)*20.0;
    }
  }

  return Success;
}


