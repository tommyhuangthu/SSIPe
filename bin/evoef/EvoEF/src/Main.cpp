///////////////////////////////////////////////////////////////////////////////////////
//Copyright (C) 2017-2019 Yang Zhang Lab, University of Michigan - All Rights Reserved
//Unauthorized copying of this file, via any medium is strictly prohibited
//Written by Xiaoqiang Huang <xiaoqiah@umich.edu>,2017 Winter
//////////////////////////////////////////////////////////////////////////////////////

#pragma warning(disable:4996)
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "Getopt.h"
#include "ProgramFunction.h"

clock_t timeStart;
clock_t timeEnd;
clock_t timePassed;

// global variables, file paths
char *atom_param_file = "./data/param_charmm19_lk_ref2015.prm";
char *residue_top_file = "./data/top_polh19_prot.inp";
char *rotamer_lib_file = "./data/rotlib984.txt";
char *pdb_structure_file = "../example/1A22.pdb";

int main(int argc, char* argv[]){
  // show EvoEF interface
  EvoEF_interface();
  timeStart = clock();
  setvbuf(stdout, NULL, _IONBF, 0);

  char usrMsg[MAX_LENGTH_ERR_MSG+1];
  //char *cmdname = "ComputeStability";
  char *cmdname = "ComputeBinding";
  //char *cmdname = "RepairStructure";
  //char *cmdname = "BuildMutant";
  //char *cmdname = "OptimizeHydrogen";
  //char *cmdname = "ComputeResiEnergy";
  //char *cmdname = "ShowResiComposition";

  int opt;
  char* mutant_file = NULL;
  char *output_file = "evoef.log";
  const char *short_opts = "-vhc:i:";
  struct option long_opts[] = {
    {"help",          no_argument,       NULL, 1},
    {"version",       no_argument,       NULL, 2},
    {"command",       required_argument, NULL, 3},
    {"pdb",           required_argument, NULL, 4},
    {"mutant-file",   required_argument, NULL, 6},
    {"output-file",   optional_argument, NULL, 8},
    {"cutoff",        required_argument, NULL, 9},
    {NULL,            no_argument,       NULL, 0}
  };
  
  while(TRUE){
    opt=getopt_long(argc, argv, short_opts, long_opts, NULL);
    if(opt == -1){
      break;
    }
    switch(opt){
      // deal with short options
      case 'h':
        EvoEF_help();
        exit(Success);
      case 'v':
        EvoEF_version();
        exit(Success);
      case 'c':
        cmdname = optarg;
        if(!CheckCommandName(cmdname)){
          printf("Command %s is not supported by EvoEF, program exits.\n", cmdname);
          exit(ValueError);
        }
        else{
          printf("Set computational type to '%s'.\n", optarg);
        }
        break;
      case 'i':
        pdb_structure_file = optarg;
        break;
      
      // deal with long options
      case 1:
        EvoEF_help();
        exit(Success);
      case 2:
        EvoEF_version();
        exit(Success);
      case 3:
        cmdname = optarg;
        if(!CheckCommandName(cmdname)){
          printf("Command %s is not supported by EvoEF, program exits.\n", cmdname);
          exit(ValueError);
        }
        else{
          printf("Command %s works.\n", cmdname);
        }
        break;
      case 4: 
        pdb_structure_file = optarg;
        break;
      case 6:
        mutant_file = optarg;
        break;
      case 8:
        FILE *fout;
        output_file = optarg;
        fout = freopen(output_file, "w", stdout);
        setvbuf(fout, NULL, _IONBF, 0);
        break;
      default:
        sprintf(usrMsg, "in file %s function %s() line %d, unknown option, EvoEF will exit.", __FILE__, __FUNCTION__, __LINE__);
        TraceError(usrMsg, ValueError);
        exit(ValueError);
        break;
    }
  }


  // deal with file name
  char pdbid[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  strcpy(pdbid,pdb_structure_file);
  int i=0;
  for(i=strlen(pdbid); i>=0; i--){
    if(pdbid[i]=='/' || pdbid[i]=='\\'){
      strcpy(pdbid,pdbid+i+1);
      break;
    }
  }
  for(int i=0; i<(int)strlen(pdbid);i++){
    if(pdbid[i]=='.'){
      strncpy(pdbid,pdbid,i);
      pdbid[i]='\0';
    }
  }

  AtomParamsSet atomParam;
  ResiTopoSet resiTopo;
  Structure structure;
  AtomParamsSetCreate(&atomParam);
  ResiTopoSetCreate(&resiTopo);
  AtomParameterRead(&atomParam, atom_param_file);
  AtomparamsSetAssignEEF1Parameters(&atomParam);
  AtomparamsSetAssignFOLDEFParameters(&atomParam);
  ResiTopoSetRead(&resiTopo, residue_top_file);
  StructureCreate(&structure);
  StructureConfig(&structure, pdb_structure_file, &atomParam, &resiTopo);
  printf("pdb file %s.pdb was read by EvoEF.\n", pdbid);

  if(!strcmp(cmdname, "ComputeStability")){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM];
    EvoEF_Stability(&structure,energyTerms);
  }
  else if(!strcmp(cmdname, "ComputeBinding")){
    double energyTerms[MAX_EVOEF_ENERGY_TERM_NUM];
    EvoEF_AnalyseComplex(&structure, energyTerms);
  }
  else if(!strcmp(cmdname, "RepairStructure")){
    RotamerLib rotlib;
    RotamerLibCreate(&rotlib,rotamer_lib_file);
    EvoEF_RepairPDB(&structure, &rotlib, &atomParam, &resiTopo,pdbid);
    RotamerLibDestroy(&rotlib);
  }
  else if(!strcmp(cmdname, "BuildMutant")){
    RotamerLib rotlib;
    RotamerLibCreate(&rotlib,rotamer_lib_file);
    mutant_file = "individual_list.txt";
    EvoEF_BuildModel(&structure, mutant_file, &rotlib, &atomParam, &resiTopo,pdbid);
    RotamerLibDestroy(&rotlib);
  }
  else if(!strcmp(cmdname, "ComputeResiEnergy")){
    for(int i=0; i<StructureGetChainCount(&structure); ++i){
      Chain* pChain=StructureGetChain(&structure,i);
      for(int j=0; j<ChainGetResidueCount(pChain); ++j){
        Residue* pResidue=ChainGetResidue(pChain,j);
        printf("residue %s%d%s energy details:\n",ChainGetName(pChain), ResidueGetPosInChain(pResidue), ResidueGetName(pResidue));
        StructureComputeResidueInteractionWithFixedSurroundingResidues(&structure, i, j);
      }
    }
  }
  else if(!strcmp(cmdname,"ShowResiComposition")){
    int aas[20]={0};
    StructureGetAminoAcidComposition(&structure,aas);
  }
  else if(!strcmp(cmdname,"OptimizeHydrogen")){
    EvoEF_OptimizeHydrogen(&structure,&atomParam, &resiTopo,pdbid);
  }
  else{
    printf("Unknown command name: %s\n, EvoEF will exit.\n", cmdname);
    exit(ValueError);
  }
  StructureDestroy(&structure);
  ResiTopoSetDestroy(&resiTopo);
  AtomParamsSetDestroy(&atomParam);

  timeEnd = clock();
  SpentTimeShow(timeStart, timeEnd);
  
  return Success;
}

