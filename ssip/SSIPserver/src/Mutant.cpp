#pragma warning(disable:4996)
#include "Mutant.h"
#include "Utility.h"
#include "ErrorHandling.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

SIP::SingleMutant::SingleMutant(){
  ;
}

SIP::SingleMutant::~SingleMutant(){
  ;
}

int SIP::SingleMutant::single_mutant_initialize(){
  native_aa = 'A';
  chain_name = 'A';
  pos = -1;
  mutant_aa = 'A';
  index_in_alignment = -1;
  smddg = 0.0;
  return Success;
}



int SIP::SingleMutant::single_mutant_create(char *mutStr){
  char aa1, chn, aa2;
  int index;
  sscanf(mutStr, "%c%c%d%c", &aa1, &chn, &index, &aa2);
  native_aa = aa1;
  chain_name = chn;
  pos = index;
  mutant_aa = aa2;

  return Success;
}

int SIP::SingleMutant::single_mutant_destroy(){
  return Success;
}

int SIP::SingleMutant::single_mutant_print(){
  printf("%c%c%d%c", native_aa, chain_name, pos, mutant_aa);
  return Success;
}


SIP::MultipleMutant::MultipleMutant(){
  mutants = NULL;
  nmut = Success;
}

SIP::MultipleMutant::~MultipleMutant(){
  if(mutants != NULL){
    free(mutants);
    mutants = NULL;
  }
  
}


int SIP::MultipleMutant::multiple_mutant_initialize(){
  nmut = 0;
  mmddg = 0.0;
  mutants = NULL;
  return Success;
}


int SIP::MultipleMutant::multiple_mutant_create(StringArray *mutStrs){
  nmut = StringArrayGetCount(mutStrs);
  mutants = (SingleMutant*)malloc(sizeof(SingleMutant)*nmut);
  for(int i = 0; i < nmut; i++){
    mutants[i].single_mutant_create(StringArrayGet(mutStrs, i));
  }

  return Success;
}

int SIP::MultipleMutant::multiple_mutant_destroy(){
  if(mutants != NULL){
    free(mutants);
    mutants = NULL;
  }
  return Success;
}

int SIP::MultipleMutant::multiple_mutant_print(){
  for(int i = 0; i < nmut; i++){
    if(i == 0){
      mutants[i].single_mutant_print();
    }
    else{
      printf(",");
      mutants[i].single_mutant_print();
    }
  }
  printf(";\n");
  return Success;
}


SIP::MutantSet::MutantSet(){
  mutset = NULL;
  nset = 0;
}

SIP::MutantSet::~MutantSet(){
  if(mutset != NULL){
    free(mutset);
    mutset = NULL;
  }
}

int SIP::MutantSet::mutant_set_initialize(){
  nset = 0;
  mutset = NULL;
  return Success;
}


int SIP::MutantSet::mutant_set_create(char* mutantfile){
  FileReader fr;
  if(FAILED(FileReaderCreate(&fr, mutantfile))){
    printf("cannot open file %s for reading\n", mutantfile);
    exit(-1);
  }
  nset = StringArrayGetCount(&fr.lines);
  mutset = (MultipleMutant*)malloc(sizeof(MultipleMutant)*nset);
  char line[1000];
  int counter = 0;
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    line[strlen(line)-1] = '\0';
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, line, ',');
    mutset[counter].multiple_mutant_create(&strings);
    counter++;
    StringArrayDestroy(&strings);
  }
  FileReaderDestroy(&fr);
  return Success;
}

int SIP::MutantSet::mutant_set_append(char *multipleMutStr){
  nset++;
  mutset = (MultipleMutant*)realloc(mutset, sizeof(MultipleMutant)*nset);
  StringArray strings;
  StringArrayCreate(&strings);
  StringArraySplitString(&strings, multipleMutStr, ',');
  mutset[nset-1].multiple_mutant_create(&strings);
  StringArrayDestroy(&strings);
  return Success;
}




int SIP::MutantSet::mutant_set_destroy(){
  if(mutset != NULL){
    for(int i = 0; i < nset; i++){
      mutset[i].multiple_mutant_destroy();
    }
    free(mutset);
    mutset = NULL;
  }
  return Success;
}

int SIP::MutantSet::mutant_set_print(){
  printf("The mutants are:\n");
  for(int i = 0; i < nset; i++){
    mutset[i].multiple_mutant_print();
  }

  return Success;
}



int SIP::MutantSet::mutant_set_get_count(){
  return nset;
}

