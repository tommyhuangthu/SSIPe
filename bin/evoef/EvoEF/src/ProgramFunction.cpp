/*******************************************************************************************************************************
This file is a part of the EvoDesign physical Energy Function (EvoEF)

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

#include "ProgramFunction.h"
#include <string.h>
#include <ctype.h>

int EvoEF_help(){
  printf(  
    "EvoEF program options:\n\n"
    "Basic OPTIONS:\n"
    "Usage: EvoEF [OPTIONS] [pdbfile]\n"  
    "   [pdbfile]       input a valid pdb format file\n"
    "   -v              print verion info of EvoEF\n"
    "   -h              print help message of EvoEF\n"  
    "   -c arg          choose your computation type:\n"
    "                   AnalyseComplex\n"
    "                   Stability\n"
    "                   BuildModel\n"
    "   -d arg          debug, produce more output\n"
    "   -f arg          config file location\n"  
    "   -o arg          output file location\n"
    );
  return Success;
}

int EvoEF_version(){
  printf("EvoDesign Force Field version 1.0\n");
  return Success;
}


int EvoEF_interface(){
  printf("******************************************************\n");
  printf("*              EvoDesign Energy Function             *\n");
  printf("*                                                    *\n");
  printf("*  Copyright (c) Xiaoqiang Huang                     *\n");
  printf("*  tommyhuangthu@foxmail.com, xiaoqiah@umich.edu     *\n");
  printf("*  The Yang Zhang Lab                                *\n");
  printf("*  Dept. of Computational Medicine & Bioinformatics  *\n");
  printf("*  Medical School                                    *\n");
  printf("*  University of Michigan                            *\n");
  printf("******************************************************\n");
  return Success;
}


BOOL CheckCommandName(char* queryname){
  int MAX_CMD_NUM = 100;
  char *supportedcmd[] = {
    "RepairStructure", 
    "ComputeStability", 
    "ComputeBinding", 
    "BuildMutant",
    "ComputeResiEnergy",
    "OptimizeHydrogen",
    "ShowResiComposition",
    NULL
  };

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

int EvoEF_Stability(Structure *pStructure, double *energyTerms){
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++) energyTerms[i] = 0.0;
  int aas[20]={0}; //ACDEFGHIKLMNPQRSTVWY, only for regular amino acid
  StructureGetAminoAcidComposition(pStructure, aas);
  // if the structure is composed of several chains, the residue position could be different in the whole structure from that in the separate chain
  StructureComputeResiduePosition(pStructure);
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    Chain *pChainI = StructureGetChain(pStructure,i);
    for(int ir = 0; ir < ChainGetResidueCount(pChainI); ir++){
      Residue *pResIR = ChainGetResidue(pChainI,ir);
      double ratio1 = CalcResidueBuriedRatio(pResIR);
      double refer=0.0;
      ResidueReferenceEnergy(pResIR, energyTerms);
      EVOEF_EnergyResidueSelfEnergy(pResIR,ratio1,energyTerms);
      for(int is = ir+1; is < ChainGetResidueCount(pChainI); is++){
        Residue *pResIS = ChainGetResidue(pChainI,is);
        double ratio2 = CalcResidueBuriedRatio(pResIS);
        double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
        if(is==ir+1) EVOEF_EnergyResidueAndNextResidue(pResIR,pResIS,ratio12,energyTerms);
        else EVOEF_EnergyResidueAndOtherResidueSameChain(pResIR,pResIS,ratio12,energyTerms);
      }
      for(int k = i+1; k < StructureGetChainCount(pStructure); k++){
        Chain *pChainK = StructureGetChain(pStructure,k);
        for(int ks = 0; ks < ChainGetResidueCount(pChainK); ks++){
          Residue *pResKS = ChainGetResidue(pChainK,ks);
          double ratio2 = CalcResidueBuriedRatio(pResKS);
          double ratio12 = CalcAverageBuriedRatio(ratio1, ratio2);
          EVOEF_EnergyResidueAndOtherResidueDifferentChain(pResIR,pResKS,ratio12,energyTerms);
        }
      }
    }
  }

  //total energy: weighted
  EnergyTermWeighting(energyTerms);
  for(int i = 1; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTerms[0] += energyTerms[i];
  }

  //energy term details: not weighted
  printf("\nStructure energy details:\n");
  printf("reference_ALA         =            %8.2f\n", energyTerms[21]);
  printf("reference_CYS         =            %8.2f\n", energyTerms[22]);
  printf("reference_ASP         =            %8.2f\n", energyTerms[23]);
  printf("reference_GLU         =            %8.2f\n", energyTerms[24]);
  printf("reference_PHE         =            %8.2f\n", energyTerms[25]);
  printf("reference_GLY         =            %8.2f\n", energyTerms[26]);
  printf("reference_HIS         =            %8.2f\n", energyTerms[27]);
  printf("reference_ILE         =            %8.2f\n", energyTerms[28]);
  printf("reference_LYS         =            %8.2f\n", energyTerms[29]);
  printf("reference_LEU         =            %8.2f\n", energyTerms[30]);
  printf("reference_MET         =            %8.2f\n", energyTerms[31]);
  printf("reference_ASN         =            %8.2f\n", energyTerms[32]);
  printf("reference_PRO         =            %8.2f\n", energyTerms[33]);
  printf("reference_GLN         =            %8.2f\n", energyTerms[34]);
  printf("reference_ARG         =            %8.2f\n", energyTerms[35]);
  printf("reference_SER         =            %8.2f\n", energyTerms[36]);
  printf("reference_THR         =            %8.2f\n", energyTerms[37]);
  printf("reference_VAL         =            %8.2f\n", energyTerms[38]);
  printf("reference_TRP         =            %8.2f\n", energyTerms[39]);
  printf("reference_TYR         =            %8.2f\n", energyTerms[40]);
  printf("intraR_vdwatt         =            %8.2f\n", energyTerms[6]);
  printf("intraR_vdwrep         =            %8.2f\n", energyTerms[7]);
  printf("intraR_electr         =            %8.2f\n", energyTerms[8]);
  printf("intraR_deslvP         =            %8.2f\n", energyTerms[9]);
  printf("intraR_deslvH         =            %8.2f\n", energyTerms[10]);
  printf("intraR_hbbbbb_dis     =            %8.2f\n", energyTerms[41]);
  printf("intraR_hbbbbb_the     =            %8.2f\n", energyTerms[42]);
  printf("intraR_hbbbbb_phi     =            %8.2f\n", energyTerms[43]);
  printf("intraR_hbscbb_dis     =            %8.2f\n", energyTerms[44]);
  printf("intraR_hbscbb_the     =            %8.2f\n", energyTerms[45]);
  printf("intraR_hbscbb_phi     =            %8.2f\n", energyTerms[46]);
  printf("intraR_hbscsc_dis     =            %8.2f\n", energyTerms[47]);
  printf("intraR_hbscsc_the     =            %8.2f\n", energyTerms[48]);
  printf("intraR_hbscsc_phi     =            %8.2f\n", energyTerms[49]);
  printf("interS_vdwatt         =            %8.2f\n", energyTerms[1]);
  printf("interS_vdwrep         =            %8.2f\n", energyTerms[2]);
  printf("interS_electr         =            %8.2f\n", energyTerms[3]);
  printf("interS_deslvP         =            %8.2f\n", energyTerms[4]);
  printf("interS_deslvH         =            %8.2f\n", energyTerms[5]);
  printf("interS_hbbbbb_dis     =            %8.2f\n", energyTerms[11]);
  printf("interS_hbbbbb_the     =            %8.2f\n", energyTerms[12]);
  printf("interS_hbbbbb_phi     =            %8.2f\n", energyTerms[13]);
  printf("interS_hbscbb_dis     =            %8.2f\n", energyTerms[14]);
  printf("interS_hbscbb_the     =            %8.2f\n", energyTerms[15]);
  printf("interS_hbscbb_phi     =            %8.2f\n", energyTerms[16]);
  printf("interS_hbscsc_dis     =            %8.2f\n", energyTerms[17]);
  printf("interS_hbscsc_the     =            %8.2f\n", energyTerms[18]);
  printf("interS_hbscsc_phi     =            %8.2f\n", energyTerms[19]);
  printf("interD_vdwatt         =            %8.2f\n", energyTerms[51]);
  printf("interD_vdwrep         =            %8.2f\n", energyTerms[52]);
  printf("interD_electr         =            %8.2f\n", energyTerms[53]);
  printf("interD_deslvP         =            %8.2f\n", energyTerms[54]);
  printf("interD_deslvH         =            %8.2f\n", energyTerms[55]);
  printf("interD_hbbbbb_dis     =            %8.2f\n", energyTerms[61]);
  printf("interD_hbbbbb_the     =            %8.2f\n", energyTerms[62]);
  printf("interD_hbbbbb_phi     =            %8.2f\n", energyTerms[63]);
  printf("interD_hbscbb_dis     =            %8.2f\n", energyTerms[64]);
  printf("interD_hbscbb_the     =            %8.2f\n", energyTerms[65]);
  printf("interD_hbscbb_phi     =            %8.2f\n", energyTerms[66]);
  printf("interD_hbscsc_dis     =            %8.2f\n", energyTerms[67]);
  printf("interD_hbscsc_the     =            %8.2f\n", energyTerms[68]);
  printf("interD_hbscsc_phi     =            %8.2f\n", energyTerms[69]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %8.2f\n\n", energyTerms[0]);
  return Success;
}



int EvoEF_AnalyseComplex(Structure *pStructure, double *energyTerms){
  double energyTermsStructure[MAX_EVOEF_ENERGY_TERM_NUM];
  double energyTermsChain[MAX_EVOEF_ENERGY_TERM_NUM];
  double energyTermsChainSum[MAX_EVOEF_ENERGY_TERM_NUM];
  for(int i = 0; i < MAX_EVOEF_ENERGY_TERM_NUM; i++){
    energyTermsChain[i] = 0.0;
    energyTermsStructure[i] = 0.0;
    energyTermsChainSum[i] = 0.0;
  }
  EvoEF_Stability(pStructure, energyTermsStructure);
  for(int i = 0; i < StructureGetChainCount(pStructure); i++){
    EvoEF_ComputeChainStability(pStructure, i, energyTermsChain);
    for(int j = 0; j < MAX_EVOEF_ENERGY_TERM_NUM; j++){
      energyTermsChainSum[j] += energyTermsChain[j];
    }
  }

  // energy terms are weighted during the calculation, don't weight them for the difference
  printf("Binding energy details:\n");
  printf("reference_ALA         =            %8.2f\n", energyTermsStructure[21] - energyTermsChainSum[21]);
  printf("reference_CYS         =            %8.2f\n", energyTermsStructure[22] - energyTermsChainSum[22]);
  printf("reference_ASP         =            %8.2f\n", energyTermsStructure[23] - energyTermsChainSum[23]);
  printf("reference_GLU         =            %8.2f\n", energyTermsStructure[24] - energyTermsChainSum[24]);
  printf("reference_PHE         =            %8.2f\n", energyTermsStructure[25] - energyTermsChainSum[25]);
  printf("reference_GLY         =            %8.2f\n", energyTermsStructure[26] - energyTermsChainSum[26]);
  printf("reference_HIS         =            %8.2f\n", energyTermsStructure[27] - energyTermsChainSum[27]);
  printf("reference_ILE         =            %8.2f\n", energyTermsStructure[28] - energyTermsChainSum[28]);
  printf("reference_LYS         =            %8.2f\n", energyTermsStructure[29] - energyTermsChainSum[29]);
  printf("reference_LEU         =            %8.2f\n", energyTermsStructure[30] - energyTermsChainSum[30]);
  printf("reference_MET         =            %8.2f\n", energyTermsStructure[31] - energyTermsChainSum[31]);
  printf("reference_ASN         =            %8.2f\n", energyTermsStructure[32] - energyTermsChainSum[32]);
  printf("reference_PRO         =            %8.2f\n", energyTermsStructure[33] - energyTermsChainSum[33]);
  printf("reference_GLN         =            %8.2f\n", energyTermsStructure[34] - energyTermsChainSum[34]);
  printf("reference_ARG         =            %8.2f\n", energyTermsStructure[35] - energyTermsChainSum[35]);
  printf("reference_SER         =            %8.2f\n", energyTermsStructure[36] - energyTermsChainSum[36]);
  printf("reference_THR         =            %8.2f\n", energyTermsStructure[37] - energyTermsChainSum[37]);
  printf("reference_VAL         =            %8.2f\n", energyTermsStructure[38] - energyTermsChainSum[38]);
  printf("reference_TRP         =            %8.2f\n", energyTermsStructure[39] - energyTermsChainSum[39]);
  printf("reference_TYR         =            %8.2f\n", energyTermsStructure[40] - energyTermsChainSum[40]);
  printf("intraR_vdwatt         =            %8.2f\n", energyTermsStructure[6] - energyTermsChainSum[6]);
  printf("intraR_vdwrep         =            %8.2f\n", energyTermsStructure[7] - energyTermsChainSum[7]);
  printf("intraR_electr         =            %8.2f\n", energyTermsStructure[8] - energyTermsChainSum[8]);
  printf("intraR_deslvP         =            %8.2f\n", energyTermsStructure[9] - energyTermsChainSum[9]);
  printf("intraR_deslvH         =            %8.2f\n", energyTermsStructure[10] - energyTermsChainSum[10]);
  printf("intraR_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[41] - energyTermsChainSum[41]);
  printf("intraR_hbbbbb_the     =            %8.2f\n", energyTermsStructure[42] - energyTermsChainSum[42]);
  printf("intraR_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[43] - energyTermsChainSum[43]);
  printf("intraR_hbscbb_dis     =            %8.2f\n", energyTermsStructure[44] - energyTermsChainSum[44]);
  printf("intraR_hbscbb_the     =            %8.2f\n", energyTermsStructure[45] - energyTermsChainSum[45]);
  printf("intraR_hbscbb_phi     =            %8.2f\n", energyTermsStructure[46] - energyTermsChainSum[46]);
  printf("intraR_hbscsc_dis     =            %8.2f\n", energyTermsStructure[47] - energyTermsChainSum[47]);
  printf("intraR_hbscsc_the     =            %8.2f\n", energyTermsStructure[48] - energyTermsChainSum[48]);
  printf("intraR_hbscsc_phi     =            %8.2f\n", energyTermsStructure[49] - energyTermsChainSum[49]);
  printf("interS_vdwatt         =            %8.2f\n", energyTermsStructure[1] - energyTermsChainSum[1]);
  printf("interS_vdwrep         =            %8.2f\n", energyTermsStructure[2] - energyTermsChainSum[2]);
  printf("interS_electr         =            %8.2f\n", energyTermsStructure[3] - energyTermsChainSum[3]);
  printf("interS_deslvP         =            %8.2f\n", energyTermsStructure[4] - energyTermsChainSum[4]);
  printf("interS_deslvH         =            %8.2f\n", energyTermsStructure[5] - energyTermsChainSum[5]);
  printf("interS_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[11] - energyTermsChainSum[11]);
  printf("interS_hbbbbb_the     =            %8.2f\n", energyTermsStructure[12] - energyTermsChainSum[12]);
  printf("interS_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[13] - energyTermsChainSum[13]);
  printf("interS_hbscbb_dis     =            %8.2f\n", energyTermsStructure[14] - energyTermsChainSum[14]);
  printf("interS_hbscbb_the     =            %8.2f\n", energyTermsStructure[15] - energyTermsChainSum[15]);
  printf("interS_hbscbb_phi     =            %8.2f\n", energyTermsStructure[16] - energyTermsChainSum[16]);
  printf("interS_hbscsc_dis     =            %8.2f\n", energyTermsStructure[17] - energyTermsChainSum[17]);
  printf("interS_hbscsc_the     =            %8.2f\n", energyTermsStructure[18] - energyTermsChainSum[18]);
  printf("interS_hbscsc_phi     =            %8.2f\n", energyTermsStructure[19] - energyTermsChainSum[19]);
  printf("interD_vdwatt         =            %8.2f\n", energyTermsStructure[51] - energyTermsChainSum[51]);
  printf("interD_vdwrep         =            %8.2f\n", energyTermsStructure[52] - energyTermsChainSum[52]);
  printf("interD_electr         =            %8.2f\n", energyTermsStructure[53] - energyTermsChainSum[53]);
  printf("interD_deslvP         =            %8.2f\n", energyTermsStructure[54] - energyTermsChainSum[54]);
  printf("interD_deslvH         =            %8.2f\n", energyTermsStructure[55] - energyTermsChainSum[55]);
  printf("interD_hbbbbb_dis     =            %8.2f\n", energyTermsStructure[61] - energyTermsChainSum[61]);
  printf("interD_hbbbbb_the     =            %8.2f\n", energyTermsStructure[62] - energyTermsChainSum[62]);
  printf("interD_hbbbbb_phi     =            %8.2f\n", energyTermsStructure[63] - energyTermsChainSum[63]);
  printf("interD_hbscbb_dis     =            %8.2f\n", energyTermsStructure[64] - energyTermsChainSum[64]);
  printf("interD_hbscbb_the     =            %8.2f\n", energyTermsStructure[65] - energyTermsChainSum[65]);
  printf("interD_hbscbb_phi     =            %8.2f\n", energyTermsStructure[66] - energyTermsChainSum[66]);
  printf("interD_hbscsc_dis     =            %8.2f\n", energyTermsStructure[67] - energyTermsChainSum[67]);
  printf("interD_hbscsc_the     =            %8.2f\n", energyTermsStructure[68] - energyTermsChainSum[68]);
  printf("interD_hbscsc_phi     =            %8.2f\n", energyTermsStructure[69] - energyTermsChainSum[69]);
  printf("----------------------------------------------------\n");
  printf("Total                 =            %8.2f\n", energyTermsStructure[0] - energyTermsChainSum[0]);
  return Success;
}


//this function is used to build the structure model of mutations
int EvoEF_BuildModel(Structure* pStructure, char* mutantfile, RotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  FileReader fr;
  FileReaderCreate(&fr, mutantfile);
  int mutantcount = FileReaderGetLineCount(&fr);
  if(mutantcount<=0){
    printf("There is no mutant found in the mutant file\n");
    return DataNotExistError;
  }

  StringArray* mutants = (StringArray*)malloc(sizeof(StringArray)*mutantcount);
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  int mutantIndex=0;
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArrayCreate(&mutants[mutantIndex]);
    StringArraySplitString(&mutants[mutantIndex], line, ',');
    char lastMutant[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    int lastmutindex = StringArrayGetCount(&mutants[mutantIndex])-1;
    strcpy(lastMutant, StringArrayGet(&mutants[mutantIndex], lastmutindex));
    //deal with the last char of the last single mutant
    if((!isdigit(lastMutant[strlen(lastMutant)-1])) && !isalpha(lastMutant[strlen(lastMutant)-1])){
      lastMutant[strlen(lastMutant)-1] = '\0';
    }
    StringArraySet(&mutants[mutantIndex], lastmutindex, lastMutant);
    mutantIndex++;
  }
  FileReaderDestroy(&fr);

  for(int mutantIndex = 0; mutantIndex < mutantcount; mutantIndex++){
    //initialize designsites first
    StructureInitializeDesignSites(pStructure);
    //for each mutant, build the rotamer-tree
    IntArray mutantArray,rotamersArray;
    IntArrayCreate(&mutantArray,0);
    IntArrayCreate(&rotamersArray,0);
    for(int cycle=0; cycle<StringArrayGetCount(&mutants[mutantIndex]); cycle++){
      char mutstr[10];
      char aa1, chn, aa2;
      int posInChain;
      strcpy(mutstr, StringArrayGet(&mutants[mutantIndex], cycle));
      sscanf(mutstr, "%c%c%d%c", &aa1, &chn, &posInChain, &aa2);
      int chainIndex = -1, residueIndex = -1;
      char chainname[MAX_LENGTH_CHAIN_NAME]; chainname[0] = chn; chainname[1] = '\0';
      StructureFindChain(pStructure, chainname, &chainIndex);
      if(chainIndex==-1){
        printf("in file %s function %s() line %d, cannot find mutation %s\n", __FILE__, __FUNCTION__, __LINE__, mutstr);
        exit(ValueError);
      }
      ChainFindResidueByPosInChain(StructureGetChain(pStructure, chainIndex), posInChain, &residueIndex);
      if(residueIndex==-1){
        printf("in file %s function %s() line %d, cannot find mutation %s\n", __FILE__, __FUNCTION__, __LINE__, mutstr);
        exit(ValueError);
      }
      char mutaatype[MAX_LENGTH_RESIDUE_NAME];
      OneLetterAAToThreeLetterAA(aa2, mutaatype);
      StringArray designType, patchType;
      StringArrayCreate(&designType);
      StringArrayCreate(&patchType);
      // for histidine, the default mutaatype is HSD, we need to add HSE
      StringArrayAppend(&designType, mutaatype); StringArrayAppend(&patchType, "");
      if(aa2=='H'){StringArrayAppend(&designType, "HSE"); StringArrayAppend(&patchType, "");}
      ProteinSiteBuildMutatedRotamers(pStructure, chainIndex, residueIndex, rotlib, atomParams, resiTopos, &designType, &patchType);
      IntArrayAppend(&mutantArray, chainIndex);
      IntArrayAppend(&mutantArray, residueIndex);
      IntArrayAppend(&rotamersArray,chainIndex);
      IntArrayAppend(&rotamersArray,residueIndex);
      StringArrayDestroy(&designType);
      StringArrayDestroy(&patchType);
    }

    // for each mutant, find the surrounding residues and build the wild-type rotamer-tree
    for(int ii=0; ii<IntArrayGetLength(&mutantArray); ii+=2){
      int chainIndex = IntArrayGet(&mutantArray,ii);
      int resiIndex = IntArrayGet(&mutantArray,ii+1);
      Residue *pResi1 = ChainGetResidue(StructureGetChain(pStructure, chainIndex), resiIndex);
      for(int j = 0; j < StructureGetChainCount(pStructure); ++j){
        Chain* pChain = StructureGetChain(pStructure,j);
        for(int k=0; k<ChainGetResidueCount(pChain); k++){
          Residue* pResi2 = ChainGetResidue(pChain,k);
          if(AtomArrayCalcMinDistance(&pResi1->atoms,&pResi2->atoms)<VDW_DISTANCE_CUTOFF){
            if(pResi2->designSiteType==Type_ResidueDesignType_Fixed){
              ProteinSiteBuildWildtypeRotamers(pStructure,j,k,rotlib,atomParams,resiTopos);
              ProteinSiteAddCrystalRotamer(pStructure,j,k,resiTopos);
              IntArrayAppend(&rotamersArray,j);
              IntArrayAppend(&rotamersArray,k);
            }
          }
        }
      }
    }

    // optimization rotamers sequentially
    printf("EvoEF Building Mutation Model %d, the following sites will be optimized:\n",mutantIndex+1);
    IntArrayShow(&rotamersArray);
    printf("\n");
    for(int cycle=0; cycle<3; cycle++){
      printf("optimization cycle %d ...\n",cycle+1);
      for(int ii=0; ii<IntArrayGetLength(&rotamersArray); ii+=2){
        int chainIndex = IntArrayGet(&rotamersArray, ii);
        int resiIndex = IntArrayGet(&rotamersArray, ii+1);
        //ProteinSiteOptimizeRotamer(pStructure, chainIndex, resiIndex);
        ProteinSiteOptimizeRotamerLocally(pStructure,chainIndex,resiIndex,1.0);
      }
    }
    IntArrayDestroy(&mutantArray);
    IntArrayDestroy(&rotamersArray);
    //remember to delete rotamers for previous mutant
    StructureDeleteRotamers(pStructure);

    char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    if(pdbid!=NULL)
      sprintf(modelfile,"%s_Model_%d.pdb",pdbid,mutantIndex+1);
    else
      sprintf(modelfile,"EvoEF_Model_%d.pdb",mutantIndex+1);
    FILE* pf=fopen(modelfile,"w");
    fprintf(pf,"REMARK EvoEF generated pdb file\n");
    fprintf(pf,"REMARK Output generated by EvoEF <BuildMutant>\n");
    StructureShowInPDBFormat(pStructure,TRUE,pf);
    fclose(pf);
  }

  return Success;
}


int EvoEF_RepairPDB(Structure* pStructure, RotamerLib* rotlib, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  StructureInitializeDesignSites(pStructure);
  for(int cycle=0; cycle<1; cycle++){
    printf("EvoEF Repairing PDB: optimization cycle %d ...\n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        //skip CYS which may form disulfide bonds
        if(strcmp(ResidueGetName(pResi),"CYS")==0) continue;
        if(strcmp(ResidueGetName(pResi),"ASN")==0||strcmp(ResidueGetName(pResi),"GLN")==0||strcmp(ResidueGetName(pResi),"HSD")==0||strcmp(ResidueGetName(pResi),"HSE")==0){
          printf("We will flip residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteBuildFlippedCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
        }
        else if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
          printf("We will rotate hydroxyl group of residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
        }
        if(TRUE){
          printf("We optimize side chain of residue %s%d%c\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteBuildWildtypeRotamers(pStructure,i,j,rotlib,atomParams,resiTopos);
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          //ProteinSiteOptimizeRotamer(pStructure,i,j);
          ProteinSiteOptimizeRotamerLocally(pStructure,i,j, 1.0);
        }
        ProteinSiteDeleteRotamers(pStructure,i,j);
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_Repair.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_Repair.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <RepairStructure>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);

  return Success;
}


int EvoEF_WriteStructureToFile(Structure* pStructure, char* pdbfile){
  FILE* pf=fopen(pdbfile,"w");
  if(pf!=NULL){
    StructureShowInPDBFormat(pStructure,TRUE,pf);
    fclose(pf);
  }
  else{
    printf("failed to open file for writing structure coordinates\n");
    return IOError;
  }
  return Success;
}

int EvoEF_AddHydrogens(Structure* pStructure, char* pdbid){
  //polar hydrogens are automatically added, so we just output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_PolH.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_PolH.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <AddHydrogens>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);
  return Success;
}


int EvoEF_OptimizeHydrogen(Structure* pStructure, AtomParamsSet* atomParams,ResiTopoSet* resiTopos, char* pdbid){
  StructureInitializeDesignSites(pStructure);
  for(int cycle=0; cycle<1; cycle++){
    printf("EvoEF Repairing PDB: optimization cycle %d ...\n",cycle+1);
    for(int i=0; i<StructureGetChainCount(pStructure); ++i){
      Chain* pChain = StructureGetChain(pStructure, i);
      for(int j=0; j<ChainGetResidueCount(pChain); j++){
        Residue* pResi = ChainGetResidue(pChain, j);
        if(strcmp(ResidueGetName(pResi),"SER")==0 || strcmp(ResidueGetName(pResi),"THR")==0 || strcmp(ResidueGetName(pResi),"TYR")==0){
          printf("We will rotate hydroxyl group of residue %s%d%c to optimize hbond\n", ResidueGetChainName(pResi),ResidueGetPosInChain(pResi),ThreeLetterAAToOneLetterAA(ResidueGetName(pResi)));
          ProteinSiteAddCrystalRotamer(pStructure,i,j,resiTopos);
          ProteinSiteExpandHydroxylRotamers(pStructure,i,j,resiTopos);
          ProteinSiteOptimizeRotamerHBondEnergy(pStructure,i,j);
          ProteinSiteDeleteRotamers(pStructure,i,j);
        }
      }
    }
  }

  //output the repaired structure
  char modelfile[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  if(pdbid!=NULL){sprintf(modelfile,"%s_OptH.pdb",pdbid);}
  else{strcpy(modelfile,"EvoEF_OptH.pdb");}
  FILE* pf=fopen(modelfile,"w");
  fprintf(pf,"REMARK EvoEF generated pdb file\n");
  fprintf(pf,"REMARK Output generated by EvoEF <OptimizeHydrogen>\n");
  StructureShowInPDBFormat(pStructure,TRUE,pf);
  fclose(pf);

  return Success;
}
