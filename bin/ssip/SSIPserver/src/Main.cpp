#pragma warning(disable:4996)

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <io.h>

#include "Utility.h"
#include "ErrorHandling.h"
#include "Mutant.h"
#include "InterfaceAlignment.h"
#include "DDG.h"
#include "ParameterOptimization.h"
#include "Structure.h"
#include "Getopt.h"
#include "CrossValidation.h"

// global variables
char *file_path_charmm_atom_param   = "./parameter/param_charmm19_lk.prm";
char *file_path_topology_amino_acid = "./parameter/top_polh19_prot.inp";
char *rotamer_lib_file              = "./parameter/rotlib984.txt";
char *file_path_pdb                 = "./example/1A22.pdb";
char *parameterfile                 = "./parameter/parameter_optimal";


int AtomParameterRead(AtomParamsSet* pAtomParam, char* filePath){
  char usrMsg[MAX_LENGTH_ERR_MSG+1];

  if( FAILED(AtomParamsSetAddFromFile(pAtomParam, filePath)) ){
    sprintf(usrMsg, "in file %s function %s() line %d, when opening:\n%s", 
      __FILE__, __FUNCTION__, __LINE__, filePath);
    TraceError(usrMsg, FormatError);
    exit(FormatError);
  }
  return Success;
}

int TopologyRead(ResiTopoSet* pResiTopo, char* filePath){
  if( filePath == NULL  || FAILED(ResiTopoCollectionAddFromFile(pResiTopo, filePath))){
      return FormatError;
  }
  return Success;
}

int StructureInitialize(Structure* pStructure, char* pdbFile, AtomParamsSet* pAtomParams, ResiTopoSet* pTopos, StringArray* pResiduePatches){
  int i, chainCounter = -1;
  BOOL firstResidueInChain = TRUE;
  int status=Success;
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  char initChainID[MAX_LENGTH_CHAIN_NAME+1];
  char initResSeq[MAX_LENGTH_ONE_LINE_IN_FILE+1];

  // read in the protein structure
  FileReader file;
  status = FileReaderCreate(&file, pdbFile);
  strcpy(initChainID, "UNK");
  initChainID[strlen(initChainID)]='\0';
  strcpy(initResSeq, "UNKNOWN");
  initResSeq[strlen(initResSeq)]='\0';
  while(!FAILED(FileReaderGetNextLine(&file, line))){
    char strType[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "";
    char strAtomName[MAX_LENGTH_ATOM_NAME+1] = "";
    char strResName[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "";
    char strChainID[MAX_LENGTH_CHAIN_NAME+1] = "";
    char strResPos[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "";

    ExtractTargetStringFromSourceString(strType, line, 0, 4);
    if(strcmp(strType, "ATOM") != 0/* && strcmp(strType, "HETA") != 0*/){
      continue;
    }
    ExtractTargetStringFromSourceString(strAtomName, line, 12, 4);
    ExtractTargetStringFromSourceString(strResName, line, 17, 4);
    ExtractTargetStringFromSourceString(strChainID, line, 21, 1);
    ExtractTargetStringFromSourceString(strResPos, line, 22, 5);

    /*if(strcmp(strChainID, "") == 0){
      int result = Warning;
      char usrMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(usrMsg, "in file %s function %s() line %d, "
        "no chain ID identified from line %s", 
        __FILE__, __FUNCTION__, __LINE__, line);
      TraceError(usrMsg, Warning);
      strcpy(strChainID, "A");
    }*/
    // if a new chain encountered
    if(strcmp(initChainID, strChainID) != 0){
      Chain newChain;
      Type_Chain chainType;
      ChainCreate(&newChain);

      chainType = ChainTypeIdentifiedFromResidueName(strResName);
      ChainSetType(&newChain, chainType);
      ChainSetName(&newChain, strChainID);
      StructureAddChain(pStructure, &newChain);

      // determine if the chain is protein or small molecule
      ChainDestroy(&newChain);
      strcpy(initChainID, strChainID);
      FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file)-1);
      firstResidueInChain = TRUE;
      chainCounter++;
    }
    // else the chain already existed
    else{
      if(strcmp(initResSeq, strResPos) != 0){
        Residue newResi;

        ResidueCreate(&newResi);
        if(strcmp(strResName, "HIS")==0){
          strcpy(strResName, "HSD");
        }
        ResidueSetName(&newResi, strResName);
        // convert the residue position into integer
        ResidueSetPosInChain(&newResi, atoi(strResPos));

        ResidueAddAtomsFromAtomParams(&newResi, pAtomParams);
        ResidueAddBondsFromResiTopos(&newResi, pTopos);

        if(firstResidueInChain == TRUE){
          ResiduePatchNTERorCTER(&newResi, "NTER", pAtomParams, pTopos);
          newResi.residueTerminalFlag = Type_ResidueIsNter;
          firstResidueInChain = FALSE;
        }
        
        FileReaderSetCurrentPos(&file, FileReaderGetCurrentPos(&file)-1);
        ResidueReadXYZFromPDB(&newResi, &file, pAtomParams, pTopos);
        if(StructureGetChain(pStructure, chainCounter) == NULL){
          int result = ValueError;
          char usrMsg[MAX_LENGTH_ERR_MSG+1];
          sprintf(usrMsg, "in file %s function %s() line %d", 
            __FILE__, __FUNCTION__, __LINE__);
          TraceError(usrMsg, ValueError);
          exit(ValueError);
        }
        else{
          ChainAppendResidue(StructureGetChain(pStructure, chainCounter), &newResi);
        }
        ResidueDestroy(&newResi);
        strcpy(initResSeq, strResPos);
      }
    }
  }

  // make patches to residues if needed and calculate all atom residues
  for(i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChain = pStructure->chains + i;
    Residue *pFirsResidueInChain = pChain->residues + 0;
    Residue *pLastResidueInChain = pChain->residues + pChain->residueNum - 1;
    if(ResidueGetAtomByName(pFirsResidueInChain, "HT1") != NULL || ResidueGetAtomByName(pFirsResidueInChain, "HN1") !=NULL){
      // if the first residue is patched with NTER, 
      // the last residue will be patched CTER to keep charge balance
      if(ResidueGetAtomByName(pLastResidueInChain, "OXT") == NULL){
        ResiduePatchCTER(pLastResidueInChain, "CTER", pAtomParams, pTopos);
        pLastResidueInChain->residueTerminalFlag = Type_ResidueIsCter;
      }
    }
    ChainCalcAllAtomXYZ(pChain, pTopos);
  }

  FileReaderDestroy(&file);
  return Success;
}


int StructureConfig(Structure *pStructure, char* configFile, char* pdbFile){
  AtomParamsSet atomParam;
  ResiTopoSet resiTopo;
  StringArray residuePatches;

  AtomParamsSetCreate(&atomParam);
  ResiTopoCollectionCreate(&resiTopo);
  StringArrayCreate(&residuePatches);

  // read configuration file
  if(configFile != NULL){
    //;
  }

  AtomParameterRead(&atomParam, file_path_charmm_atom_param);
  TopologyRead(&resiTopo, file_path_topology_amino_acid);
  StructureInitialize(pStructure, pdbFile, &atomParam, &resiTopo, &residuePatches);

  //ProteinRotamerGenerate(pStructure, &atomParam, &resiTopo);
  //StructureShowDesignSites(pStructure, stdout);
  //StructureSetEEF1AtomParameter(pStructure);
  //StructureSetHbondDonorAndAcceptor(pStructure);
  //StructureCheckIntraBondType(pStructure);
  //StructureCheckNeighbouringBondType(pStructure);
  //StructureSetFOLDEFAtomParameter(pStructure);
  //StructureCalcAtomOccupancy(pStructure);
  //StructureCalcResiduePosition(pStructure);
  //CalcResidueSelfEnergy(pStructure, "A", 71);

  ResiTopoCollectionDestroy(&resiTopo);
  AtomParamsSetDestroy(&atomParam);
  StringArrayDestroy(&residuePatches);
  return Success;
}


char AminoAcid3To1(const char* three){
  char *list[] = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "HSD", "HSE", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR", "XXX"};
  int index = -1;
  for(int i = 0; i < 23; i++){
    if(!strcmp(list[i], three)){
      index = i;
      break;
    }
  }

  switch(index){
    case 0:
      return 'A';
      break;
    case 1:
      return 'C';
    case 2:
      return 'D';
    case 3:
      return 'E';
    case 4:
      return 'F';
    case 5:
      return 'G';
    case 6:
      return 'H';
    case 7:
      return 'H';
    case 8:
      return 'H';
    case 9:
      return 'I';
    case 10:
      return 'K';
    case 11:
      return 'L';
    case 12:
      return 'M';
    case 13:
      return 'N';
    case 14:
      return 'P';
    case 15:
      return 'Q';
    case 16:
      return 'R';
    case 17:
      return 'S';
    case 18:
      return 'T';
    case 19:
      return 'V';
    case 20:
      return 'W';
    case 21:
      return 'Y';
    case 22:
      return 'X';
    default:
      return 'X';
  }
  return 'X';
}


int StructureFindInterfaceResidues(Structure *pStructure, double cutoff){
  if(pStructure->chainNum < 2){
    printf("There is only one chain in the whole structure, on protein-protein interface found.\n");
    exit(ValueError);
  }

  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChainI = StructureGetChain(pStructure, i);
    for(int k = i+1; k < pStructure->chainNum; k++){
      Chain *pChainK = StructureGetChain(pStructure, k);
      IntArray arrayI, arrayK;
      IntArrayCreate(&arrayI, pChainI->residueNum);
      for(int j = 0; j < pChainI->residueNum; j++){
        IntArraySet(&arrayI, j, -100);
      }
      IntArrayCreate(&arrayK, pChainK->residueNum);
      for(int s = 0; s < pChainK->residueNum; s++){
        IntArraySet(&arrayK, s, -100);
      }
      for(int j = 0; j < pChainI->residueNum; j++){
        Residue* pResiIJ = ChainGetResidue(pChainI, j);
        for(int s = 0; s < pChainK->residueNum; s++){
          Residue* pResiKS = ChainGetResidue(pChainK, s);
          if(AtomArrayCalcMinDistance(&pResiIJ->atoms, &pResiKS->atoms) < cutoff){
            IntArraySet(&arrayI, j, pResiIJ->posInChain);
            IntArraySet(&arrayK, s, pResiKS->posInChain);
          }

        }
      }
      printf("interface residues between Chain %s and %s with cutoff = %.2f by sip\n", pChainI->name, pChainK->name, cutoff);
      for(int j = 0; j < IntArrayGetLength(&arrayI); j++){
        if(IntArrayGet(&arrayI, j) >= 0){
          printf("%c%s%d ", AminoAcid3To1(pChainI->residues[j].name), pChainI->name, IntArrayGet(&arrayI, j));
        }
      }
      for(int s = 0; s < IntArrayGetLength(&arrayK); s++){
        if(IntArrayGet(&arrayK, s) >= 0){
          printf("%c%s%d ", AminoAcid3To1(pChainK->residues[s].name), pChainK->name, IntArrayGet(&arrayK, s));
        }
      }
      printf("\n");
    }
  }

  return Success;
}

Residue* StructureFindMutationResidue(Structure *pStructure, char *mutant){
  char nativeResi[2], mutantResi[2], chainName[2];
  int posInChain;
  sscanf(mutant, "%c%c%d%c", &nativeResi[0], &chainName[0], &posInChain, &mutantResi[0]);
  nativeResi[1] = chainName[1] = mutantResi[1] = '\0';
  for(int i = 0; i < pStructure->chainNum; i++){
    if(pStructure->chains[i].name[0] != chainName[0]) continue;
    for(int j = 0; j < pStructure->chains[i].residueNum; j++){
      if(posInChain == pStructure->chains[i].residues[j].posInChain){
        return &pStructure->chains[i].residues[j];
      }
    }
  }

  return NULL;
}

int ResidueGetNeighbourChainCount(Residue* pThis, Structure *pStructure, double cutoff, StringArray *pChains){
  // the residue must contact with its own chain
  int neighbour_chain_count = 1;
  for(int i = 0; i < pStructure->chainNum; i++){
    Chain *pChain = StructureGetChain(pStructure, i);
    if(pChain->name[0] == pThis->chainName[0]) continue;
    for(int j = 0; j < pChain->residueNum; j++){
      Residue *pResi = ChainGetResidue(pChain, j);
      if(AtomArrayCalcMinDistance(&pThis->atoms, &pResi->atoms) < cutoff){
        // for the first residue in the mutation list, pos is always equal to -1
        // so always add the chain name
        int pos = -1;
        StringArrayFind(pChains, pThis->chainName, &pos);
        if(pos == -1){
          StringArrayAppend(pChains, pThis->chainName);
        }
        pos = -1;
        StringArrayFind(pChains, pResi->chainName, &pos);
        if(pos == -1){
          StringArrayAppend(pChains, pResi->chainName);
        }
        neighbour_chain_count++;
        break;
      }
    }
  }
  return neighbour_chain_count;
}

int FindMutationsInvolvedInOnlyDimerComplexFromSKEMPI()
{
  char *mutants_in_more_than_two_chains = "D:/work/database/SKEMPIv2.0/skempi_v2_part1_mutants_in_more_than_two_chains.txt";
  char *mutants_not_in_interface = "D:/work/database/SKEMPIv2.0/skempi_v2_part2_mutants_not_in_interface.txt";
  char *mutants_in_two_sides = "D:/work/database/SKEMPIv2.0/skempi_v2_part3_mutants_in_two_sides.txt";
  char *mutants_only_in_one_side = "D:/work/database/SKEMPIv2.0/skempi_v2_part4_mutants_only_in_one_side.txt";
  char *complex_has_more_than_two_chains = "D:/work/database/SKEMPIv2.0/skempi_v2_part5_complex_has_more_than_two_chains.txt";
  char *residues_not_point_to_interface="D:/work/database/SKEMPIv2.0/skempi_v2_part6_residue_not_point_to_interface.txt";
  FILE* file1 = fopen(mutants_in_more_than_two_chains, "w");
  FILE* file2 = fopen(mutants_not_in_interface, "w");
  FILE* file3 = fopen(mutants_in_two_sides, "w");
  FILE* file4 = fopen(mutants_only_in_one_side, "w");
  FILE* file5 = fopen(complex_has_more_than_two_chains, "w");
  // process skempiv2 files
  char *skempiv2 = "D:/work/database/SKEMPIv2.0/skempi_v2_corrected.txt";
  FileReader fr;
  FileReaderCreate(&fr, skempiv2);
  char line[1024];
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, line, ';');
    char *first = StringArrayGet(&strings, 0);
    if(first[0] != '#'){
      StringArray pdbs;
      StringArrayCreate(&pdbs);
      StringArraySplitString(&pdbs, first, '_');
      char *pdbid = StringArrayGet(&pdbs, 0);
      //if(strcmp(pdbid, "1FY8")!=0) continue;
      //printf("processing %s ...\n", pdbid);
      char *chain1 = StringArrayGet(&pdbs, 1);
      char *chain2 = StringArrayGet(&pdbs, 2);
      char *mutantStr = StringArrayGet(&strings, 2);
      StringArray mutants;
      StringArrayCreate(&mutants);
      StringArraySplitString(&mutants, mutantStr, ',');
      // discard the complexes that have more than two chains
      if(strlen(chain1) > 1 || strlen(chain2) > 1){
        fprintf(file5, "%s\n", line);
        StringArrayDestroy(&strings);
        StringArrayDestroy(&mutants);
        continue;
      }

      // read PDB into structure;
      Structure structure;
      StructureCreate(&structure);
      char cmd[1024];
      sprintf(cmd, "D:/work/database/SKEMPIv2.0/pdbs_corrected/%s.pdb", pdbid);
      StructureConfig(&structure, "", cmd);
      StringArray chainnames_involved_in_all_mutants;
      StringArrayCreate(&chainnames_involved_in_all_mutants);
      BOOL all_mutation_in_interface = TRUE;
      for(int i = 0; i < StringArrayGetCount(&mutants); i++){
        char *mutant = StringArrayGet(&mutants, i);
        Residue *pCurResi = StructureFindMutationResidue(&structure, mutant);
        
        int residue_contact_chain_count = 1;
        residue_contact_chain_count = ResidueGetNeighbourChainCount(pCurResi, &structure, 5.0, &chainnames_involved_in_all_mutants);
        // the residue is not in the interface
        if(residue_contact_chain_count == 1){
          //printf("%s: Multiple mutant %s has mutation %s not in the interface\n", first, mutantStr, mutant);
          fprintf(file2, "%s\n", line);
          all_mutation_in_interface = FALSE;
          break;
        }
      }
      if(all_mutation_in_interface == TRUE){
        if(StringArrayGetCount(&chainnames_involved_in_all_mutants) >= 3){
          //printf("%s: Multiple mutant %s is involved in 3 or more chains\n", first, mutantStr);
          fprintf(file1, "%s\n", line);
        }
        else if(StringArrayGetCount(&chainnames_involved_in_all_mutants) == 2){
          // check if the mutants involved in only one side or 2+ sides
          char firstChainName;
          int side_num = 0;
          for(int i = 0; i < StringArrayGetCount(&mutants); i++){
            char *mutant = StringArrayGet(&mutants, i);
            if(i == 0){firstChainName = mutant[1]; side_num++;}
            else{
              if(mutant[1] != firstChainName){
                side_num++;
              }
            }
          }
          
          if(side_num >= 2){// the mutants is involved in 2 or more sides
            //printf("%s: Mutant %s is involved in 2 or more sides\n", first, mutantStr);
            if(strlen(chain1) == 1 && strlen(chain2) == 1){
              fprintf(file3, "%s\n", line);
            }
            else{
              char firstStr[1024];
              ExtractFirstStringFromSourceStringNew(firstStr, line, ';');
              char *chn1 = StringArrayGet(&chainnames_involved_in_all_mutants, 0);
              char *chn2 = StringArrayGet(&chainnames_involved_in_all_mutants, 1);
              if(chn1[0] < chn2[0]){
                //printf("%s_%s_%s;", pdbid, chn1, chn2);
                fprintf(file3, "%s_%s_%s;", pdbid, chn1, chn2);
              }
              else{
                //printf("%s_%s_%s;", pdbid, chn2, chn1);
                fprintf(file3, "%s_%s_%s;", pdbid, chn2, chn1);
              }
              //printf("%s\n", line);
              fprintf(file3, "%s\n", line);
            }
          }

          else{ // all_mutation_in_one_side
            if(strlen(chain1) == 1 && strlen(chain2) == 1){
              //printf("%s\n",line);
              fprintf(file4, "%s\n", line);
            }
            else{
              char firstStr[1024];
              ExtractFirstStringFromSourceStringNew(firstStr, line, ';');
              char *chn1 = StringArrayGet(&chainnames_involved_in_all_mutants, 0);
              char *chn2 = StringArrayGet(&chainnames_involved_in_all_mutants, 1);
              if(side_num >= 2){

              }
              else{
                if(chn1[0] < chn2[0]){
                  //printf("%s_%s_%s;", pdbid, chn1, chn2);
                  fprintf(file4, "%s_%s_%s;", pdbid, chn1, chn2);
                }
                else{
                  //printf("%s_%s_%s;", pdbid, chn2, chn1);
                  fprintf(file4, "%s_%s_%s;", pdbid, chn2, chn1);
                }
                //printf("%s\n", line);
                fprintf(file4, "%s\n", line);
              }
            }
          }
        }

      } // if all_mutation_in_interface == TRUE
      
      StringArrayDestroy(&chainnames_involved_in_all_mutants);
      StructureDestroy(&structure);
      StringArrayDestroy(&mutants);
      StringArrayDestroy(&pdbs);
    }
    /*else{
      printf("%s\n", line);
    }*/
    StringArrayDestroy(&strings);
  }
  FileReaderDestroy(&fr);
  fclose(file1);
  fclose(file2);
  fclose(file3);
  fclose(file4);
  fclose(file5);

  exit(Success);
  return Success;
}



void program_help(){
  printf(  
    "SSIP program options:\n\n"
    "Basic OPTIONS:\n"
    "--version\n"
    "--help             produce help message\n"
    "--command=arg      choose command:\n"
    "                   DDG\n"
    "                   FindInterfaceResidue\n" 
    "                   OptimizeParameter\n"
    "                   CrossValidation\n"
    "                   ExtractMutation\n"
    "                   ShowAlignment\n"
    "                   TestDDGPrediction\n"
    );
  return;
}

void program_version(){
  printf("Sequence and Structural Interface Profile potential version 1.0\n");
  printf("used to calculate binding energy change due to mutation at the interface\n");
  return;
}

void program_interface(){
  printf("******************************************************\n");
  printf("*  Sequence-Strcuture based Interface Profile (SSIP) *\n");
  printf("*                                                    *\n");
  printf("*  Copyright (C) 2018-2019 The Yang Zhang Lab        *\n");
  printf("*  Prof. Yang Zhang's Research Group                 *\n");
  printf("*  Dept. of Computational Medicine & Bioinformatics  *\n");
  printf("*  Medical School                                    *\n");
  printf("*  University of Michigan                            *\n");
  printf("*  Ann Arbor, MI 48109-2218, USA                     *\n");
  printf("*  E-mail: zhng@umich.edu                            *\n");
  printf("******************************************************\n");
  return;
}

BOOL CmdNameExisted(char* queryname){
  int MAX_CMD_NUM = 100;
  char *supportedcmd[] = {
    "OptimizeParameter", 
    "DDG", 
    "FindInterfaceResidue", 
    "CrossValidation", 
    "ExtractMutation", 
    "ShowAlignment",
    "TestDDGPrediction",
    NULL};

  BOOL exist = FALSE;
  for(int i = 0; i < MAX_CMD_NUM; i++){
    if(supportedcmd[i] == NULL) break;
    else{
      if(strcmp(queryname, supportedcmd[i]) == 0){
        exist = TRUE;
        break;
      }
    }
  }
  return exist;
}


int main(int argc, char** argv){
  const char *short_opts = "";
  struct option long_opts[] = {
    {"help",                 no_argument,       NULL, 1 },
    {"version",              no_argument,       NULL, 2 },
    {"command",              required_argument, NULL, 3 },
    {"pdb",                  required_argument, NULL, 4 },
    {"working_path",         required_argument, NULL, 5 },
    {"mutantlist_file",      required_argument, NULL, 6 },
    {"isscore",              required_argument, NULL, 7 },
    {"output",               optional_argument, NULL, 8 },
    {"cutoff",               required_argument, NULL, 9 },
    {"mutant_file",          required_argument, NULL, 10},
    {"structure_align",      required_argument, NULL, 11},
    {"sequence_align",       required_argument, NULL, 12},
    {"parameter_file",       required_argument, NULL, 15},
    {"kfold",                required_argument, NULL, 16},
    {"cvcount",              required_argument, NULL, 17},
    {"ssipscore_file",       required_argument, NULL, 18},
    {NULL,                   no_argument,       NULL, 0 }
  };

  //char *cmdname = "DDG";
  //char *cmdname = "OptimizeParameter";
  //char *cmdname = "FindInterfaceResidue";
  char *cmdname = "CrossValidation";
  //char *cmdname = "ExtractMutation";
  //char *cmdname = "ShowAlignment";
  //char *cmdname = "TestDDGPrediction";


  char usrMsg[MAX_LENGTH_ERR_MSG+1];
  double interface_dist_cutoff = 5.0;
  double isscore = 0.55;

  //set parameters for cross validation
  int cvcount = 2;
  int kfold = 5;

  char *output_file      = "sip.out";
  char *ssipscore_file   = "ssip_score.txt";
  char* trainingdatafile = "data/skempiv2/pengx/part34/mutant2204.txt";
  char* workingpath      = "interface_alignment_xq_ev1e-3";
  char* mutantfile       = NULL;

  //for cross validation
  char* structure_align  = "structure_align.out";
  char* sequence_align   = "sequence_align.out";

  setvbuf(stdout, NULL, _IONBF, 0);
  int opt;
  while(TRUE){
    opt=getopt_long(argc, argv, short_opts, long_opts, NULL);
    if(opt == -1){
      break;
    }
    switch(opt){
    case 1:
      program_help();
      exit(Success);
    case 2:
      program_version();
      exit(Success);
    case 3:
      cmdname = optarg;
      if(!CmdNameExisted(cmdname)){
        printf("Command %s is not supported by SIP, program exits.\n", cmdname);
        exit(ValueError);
      }
      else{
        printf("Command %s works.\n", cmdname);
      }
      break;
    case 4: 
      file_path_pdb = optarg;
      break;
    case 5:
      workingpath = optarg;
      break;
    case 6:
      trainingdatafile = optarg;
      break;
    case 7:
      isscore = atof(optarg);
      break;
    case 8:
      FILE *fout;
      output_file = optarg;
      fout = freopen(output_file, "w", stdout);
      setvbuf(fout, NULL, _IONBF, 0);
      break;
    case 9:
      interface_dist_cutoff = atof(optarg);
      break;
    case 10:
      mutantfile = optarg;
      break;
    case 11:
      structure_align = optarg;
      break;
    case 12:
      sequence_align = optarg;
      break;
    case 15:
      parameterfile = optarg;
      break;
    case 16:
      kfold = atoi(optarg);
      break;
    case 17:
      cvcount = atoi(optarg);
      break;
    case 18:
      ssipscore_file = optarg;
      break;
    default:
      sprintf(usrMsg, "in file %s function %s() line %d, unknown option, program will exit.", __FILE__, __FUNCTION__, __LINE__);
      TraceError(usrMsg, ValueError);
      exit(ValueError);
      break;
    }
  }


  if(!strcmp(cmdname, "DDG")){
    SIP::Parameter par;
    par.parameter_read(parameterfile);

    SIP::DDG ddg;
    // read mutant information
    ddg.ms.mutant_set_initialize();
    ddg.ms.mutant_set_create(mutantfile);
    ddg.ms.mutant_set_print();
    // read interface alignment information
    ddg.ia_structure.interface_alignment_read_bindprofx_with_cutoff(structure_align, par.ial_weight, isscore, par.ial_count_cutoff);
    ddg.ia_sequence.interface_alignment_read_sip_with_cutoff(sequence_align, par.psi_weight, 
      par.psi_linkscore_cutoff, par.seqid_high_cutoff, par.seqid_low_cutoff, par.psi_count_cutoff);
    // calculate DDG change
    ddg.calc_ddg_for_mutation_set1234(&par);
    ddg.write_ddg_for_mutation_set(ssipscore_file);
    ddg.ms.mutant_set_destroy();
  }
  else if(!strcmp(cmdname, "TestDDGPrediction")){
    SIP::CrossValidation cv;
    SIP::Parameter par;
    par.parameter_read(parameterfile);
    par.parameter_print();
    double pearson, rmse;
    cv.cross_validation_check_test_performance1234(&par, trainingdatafile, workingpath, structure_align, sequence_align, NULL, NULL, &pearson, &rmse);
    SIP::print_pearson_and_rmse(pearson, rmse);
  }
  else if(!strcmp(cmdname, "OptimizeParameter")){
    for(int i=0; i<cvcount; i++){
      SIP::MCoptimization mcs;
      SIP::Parameter par;
      mcs.optimization_by_sa1234(&par, trainingdatafile, workingpath, structure_align, sequence_align, NULL, NULL, parameterfile);
      mcs.optimization_release_memory();
    }
  }
  else if(!strcmp(cmdname, "CrossValidation")){
    for(int i=0; i<cvcount; i++){
      printf("cross validation cycle %d:\n", i+1);
      SIP::CrossValidation cv;
      //cv.cross_validation_divide_train_and_test_set(trainingdatafile, 0.20);
      cv.cross_validation_divide_train_set_into_kfold(trainingdatafile, kfold);
      cv.cross_validation_run1234(trainingdatafile, kfold, workingpath, structure_align, sequence_align, NULL, NULL, parameterfile);
    }
  }
  else if(!strcmp(cmdname, "FindInterfaceResidue")){
    Structure structure;
    StructureCreate(&structure);
    StructureConfig(&structure, "", file_path_pdb);
    StructureFindInterfaceResidues(&structure, interface_dist_cutoff);
    StructureDestroy(&structure);
  }
  else if(!strcmp(cmdname, "ExtractMutation")){
    FindMutationsInvolvedInOnlyDimerComplexFromSKEMPI();
  }
  else if(!strcmp(cmdname, "ShowAlignment")){
    SIP::MCoptimization mcs;
    SIP::Parameter par;
    par.parameter_read(parameterfile);
    par.parameter_print();
    mcs.optimization_read_data_from_file1234(&par, trainingdatafile, workingpath, structure_align, sequence_align, NULL, NULL);
    int averageNum = 0;
    for(int i = 0; i < mcs.expdata.npdb; i++){
      SIP::InterfaceAlignment *ia = &mcs.ddgs[i].ia_structure;
      printf("structure align file: %s, alignment number: %d\n", ia->alignment_filename, ia->lig_alignment_num);
      averageNum += ia->lig_alignment_num;
    }
    printf("The average number of structure align is: %f\n", (double)averageNum/mcs.expdata.npdb);

    averageNum = 0;
    for(int i = 0; i < mcs.expdata.npdb; i++){
      SIP::InterfaceAlignment *ia = &mcs.ddgs[i].ia_sequence;
      printf("sequence align file: %s, alignment number: %d\n", ia->alignment_filename, ia->lig_alignment_num);
      averageNum += ia->lig_alignment_num;
    }
    printf("The average number of sequence align is: %f\n", (double)averageNum/mcs.expdata.npdb);
    mcs.optimization_release_memory();
  }
  
  return Success;
}
