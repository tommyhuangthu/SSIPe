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

#include "Structure.h"
//#include "EnergyFunction.h"
#include <string.h>

int StructureCreate(Structure* pThis)
{
  strcpy(pThis->name, "");
  pThis->chainNum = 0;
  pThis->chains = NULL;
  return Success;
}
int StructureDestroy(Structure* pThis)
{
  int i;
  for(i=0;i<pThis->chainNum;i++)
  {
    ChainDestroy(&pThis->chains[i]);
  }
  free(pThis->chains);

  strcpy(pThis->name, "");
  pThis->chainNum = 0;
  pThis->chains = NULL;
  return Success;
}

int StructureCopy(Structure* pThis, Structure* pOther)
{
  return Success;
}

char* StructureGetName(Structure* pThis)
{
  return pThis->name;
}

int StructureSetName(Structure* pThis, char* newName)
{
  if(strlen(newName)>MAX_LENGTH_STRUCTURE_NAME)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->name, newName);
  return Success;
}

int StructureGetChainCount(Structure* pThis)
{
  return pThis->chainNum;
}

Chain* StructureGetChain(Structure* pThis, int index)
{
  if(index<0 || index>=pThis->chainNum)
  {
    return NULL;
  }
  return pThis->chains + index;
}

Chain* StructureGetChainByName(Structure* pThis, char* segName)
{
  int    result;
  int index = -1;
  result = StructureFindChain(pThis, segName, &index);
    if(FAILED(result))
  {
    return NULL;
    }
    else
  {
        return StructureGetChain(pThis, index);
    }
}

int StructureFindChain(Structure* pThis, char* segName, int* index)
{
  int i;
  for(i=0;i<pThis->chainNum;i++)
  {
    if(strcmp(ChainGetName(&pThis->chains[i]), segName)==0)
    {
      *index = i;
      return Success;
    }
  }
  return DataNotExistError;
}
int StructureFindSmallMol(Structure* pThis, Residue** ppSmallMol){
  int i;
  int    result = DataNotExistError;
  for(i=0;i<pThis->chainNum;i++){
    if(pThis->chains[i].type == Type_Chain_SmallMol){
      *ppSmallMol = ChainGetResidue(&pThis->chains[i], 0);
      if(*ppSmallMol != NULL)
      {
         result = Success;
         break;
      } 
    }
  }
  return result;
}

int StructureAddChain(Structure* pThis, Chain* newSeg)
{
  int index = -1;
  int result;
  result = StructureFindChain(pThis, ChainGetName(newSeg), &index);
  if(FAILED(result))
  {
    (pThis->chainNum)++;
    pThis->chains = (Chain*)realloc(pThis->chains, sizeof(Chain)*pThis->chainNum);
    ChainCreate(&pThis->chains[pThis->chainNum-1]);
    return ChainCopy(&pThis->chains[pThis->chainNum-1], newSeg);
  }
  else
  {
    return ChainCopy(&pThis->chains[index], newSeg);
  }
}

int StructureDeleteChain(Structure* pThis, char* segName)
{
  int index, i;
  int    result;
  result = StructureFindChain(pThis, segName, &index);
  if(FAILED(result))
    return result;
  for(i=index;i<pThis->chainNum-1;i++)
  {
    ChainCopy(&pThis->chains[i], &pThis->chains[i+1]);
  }
  ChainDestroy(&pThis->chains[pThis->chainNum-1]);
  (pThis->chainNum)--;
  return Success;
}

int StructureShowInPDBFormat(Structure* pThis, BOOL showHydrogen, FILE* pFile)
{
  Chain* pChain;
  int chainNum, atomNum = 1, i, j;

  chainNum = StructureGetChainCount(pThis);
  for(i=0;i<chainNum;i++)
  {
    int resiNum;
    pChain = StructureGetChain(pThis, i);
    resiNum = pChain->residueNum;
    for(j = 0; j < resiNum; j++)
    {
      Residue *pResiIJ = pChain->residues + j;
      ResidueShowInPDBFormat(pResiIJ, "ATOM", pResiIJ->chainName, atomNum, pResiIJ->posInChain, TRUE, pFile);
      atomNum += pResiIJ->atoms.atomNum;
    }
  }
  return Success;
}



