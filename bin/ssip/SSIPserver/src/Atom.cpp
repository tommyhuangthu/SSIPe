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

#include "Atom.h"
#include <string.h>

int AtomCreate(Atom* pThis)
{
  strcpy(pThis->name, "");
  pThis->xyz.X = pThis->xyz.Y = pThis->xyz.Z = 0.0;
  pThis->isXyzValid = FALSE;
  strcpy(pThis->chainName, "");
  pThis->posInChain = -1;
  pThis->isInHBond = FALSE;
  pThis->isBBAtom = FALSE;

  pThis->FOLDEF_OccCal = 0.0;

  return Success;
}
void AtomDestroy(Atom* pThis)
{
  strcpy(pThis->name, "");
  pThis->xyz.X = pThis->xyz.Y = pThis->xyz.Z = 0.0;
  pThis->isXyzValid = 0;
}

int AtomCopy(Atom* pThis, Atom* pOther)
{
  *pThis = *pOther;
  return Success;
}

char* AtomGetName(Atom* pThis)
{
  return pThis->name;
}

int AtomSetName(Atom* pThis, char* newName)
{
  if(newName == NULL || strlen(newName)>MAX_LENGTH_ATOM_NAME)
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

char* AtomGetType(Atom* pThis)
{
  return pThis->type;
}

Type_AtomHybridType AtomGetHybridType(Atom* pThis)
{
  return pThis->hybridType;
}

char* AtomGetHbHorA(Atom* pThis)
{
  return pThis->hbHorA;
}

char* AtomGetHbDorB(Atom* pThis)
{
  return pThis->hbDorB;
}

char* AtomGetHbB2(Atom* pThis)
{
  return pThis->hbB2;
}

BOOL AtomIsHydrogen(Atom* pThis)
{
  if(pThis->name[0] == 'H')
  {
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

//int AtomSetParamsByStringArray(Atom* pThis, StringArray* pParams)
//{
//  char* atomName;
//  char* atomType;
//  char* isBB;
//  char* polar;
//  double epsilon;
//  double rmin;
//  double charge;
//  char* donor;
//  char* acceptor;
//  char* hybrid;
//
//  if(StringArrayGetCount(pParams)!=11)
//  {
//    return ValueError;
//  }
//
//  atomName = StringArrayGet(pParams, 1);
//  atomType = StringArrayGet(pParams, 2);
//  isBB     = StringArrayGet(pParams, 3);
//  polar    = StringArrayGet(pParams, 4);
//  epsilon  = atof(StringArrayGet(pParams, 5));
//  rmin     = atof(StringArrayGet(pParams, 6));
//  charge   = atof(StringArrayGet(pParams, 7));
//  donor    = StringArrayGet(pParams, 8);
//  acceptor = StringArrayGet(pParams, 9);
//  hybrid   = StringArrayGet(pParams, 10);
//
//  if( strlen(atomName) > MAX_LENGTH_ATOM_NAME ||
//    strlen(atomType) > MAX_LENGTH_ATOM_TYPE ||
//    strlen(donor)    > MAX_LENGTH_ATOM_DONOR ||
//    strlen(acceptor) > MAX_LENGTH_ATOM_ACCEPTOR)
//  {
//    return ValueError;
//  }
//
//  pThis->CHARMM22_epsilon = epsilon;
//  pThis->CHARMM22_radius = rmin;
//  pThis->CHARMM22_charge = charge;
//
//  //AtomSetName(pThis, atomName);
//  strcpy(pThis->name, atomName);
//  strcpy(pThis->type, atomType);
//
//  if(strcmp(donor, "Y") == 0)
//  {
//    pThis->isHBatomH = TRUE;
//  }
//  else
//  {
//    pThis->isHBatomH = FALSE;
//  }
//  
//  if(strcmp(acceptor, "Y") == 0)
//  {
//    pThis->isHBatomA = TRUE;
//  }
//  else
//  {
//    pThis->isHBatomA = FALSE;
//  }
//
//  strcpy(pThis->donor, donor);
//  strcpy(pThis->acceptor, acceptor);
//  //strcpy(pThis->hybridType, hybrid);
//  
//
//  if(strcmp(hybrid,"SP2") == 0)
//  {
//    pThis->hybridType = Type_AtomHybridType_SP2;
//  }
//  else if(strcmp(hybrid,"SP3") == 0)
//  {
//    pThis->hybridType = Type_AtomHybridType_SP2;
//  }
//  else if(strcmp(hybrid,"SP") == 0)
//  {
//    pThis->hybridType = Type_AtomHybridType_SP;
//  }
//
//  switch(isBB[0])
//  {
//      case 'Y':
//        pThis->isBBAtom = TRUE;break;
//      case 'N':
//        pThis->isBBAtom = FALSE;break;
//      default:
//        return ValueError;
//  }
//
//  if(pThis->isBBAtom == TRUE)
//  {
//    if(strcmp(atomName, "HN") == 0
//      || strcmp(atomName, "HT1") == 0
//      || strcmp(atomName, "HT2") == 0
//      || strcmp(atomName, "HT3") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_HN;
//    }
//    else if(strcmp(atomName, "N") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_N;
//    }
//    else if(strcmp(atomName, "HA") == 0
//      || strcmp(atomName, "HA1") == 0
//      || strcmp(atomName, "HA2") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_HA;
//    }
//    else if(strcmp(atomName, "CA") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_CA;
//    }
//    else if(strcmp(atomName, "C") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_C;
//    }
//    else if(strcmp(atomName, "O") == 0
//      || strcmp(atomName, "OXT") == 0)
//    {
//      pThis->bbAtomType = Type_BB_Atom_O;
//    }
//  }
//  else
//  {
//    pThis->bbAtomType = Type_NotBBAtom;
//  }
//
//  if(strcmp(polar, "P")==0)
//  {
//    pThis->polarity = Type_AtomPolarity_P;
//  }
//  else if(strcmp(polar, "C")==0)
//  {
//    pThis->polarity = Type_AtomPolarity_C;
//  }
//  else if(strcmp(polar, "NP1")==0)
//  {
//    pThis->polarity = Type_AtomPolarity_NPAliphatic;
//  }
//  else if(strcmp(polar, "NP2")==0)
//  {
//    pThis->polarity = Type_AtomPolarity_NPAromatic;
//  }
//  else
//  {
//    return ValueError;
//  }
//
//  return Success;
//}

int AtomSetParamsByStringArray(Atom* pThis, StringArray* pParams)
{
  char* atomName;
  char* atomType;
  char* isBB;
  char* polar;
  double epsilon;
  double rmin;
  double charge;
  char* hbHorA;
  char* hbDorB;
  char* hbB2;
  char* hybrid;

  if(StringArrayGetCount(pParams)!=12)
  {
    return ValueError;
  }

  atomName = StringArrayGet(pParams, 1);
  atomType = StringArrayGet(pParams, 2);
  isBB     = StringArrayGet(pParams, 3);
  polar    = StringArrayGet(pParams, 4);
  epsilon  = atof(StringArrayGet(pParams, 5));
  rmin     = atof(StringArrayGet(pParams, 6));
  charge   = atof(StringArrayGet(pParams, 7));
  hbHorA   = StringArrayGet(pParams, 8);
  hbDorB   = StringArrayGet(pParams, 9);
  hbB2     = StringArrayGet(pParams, 10);
  hybrid   = StringArrayGet(pParams, 11);

  if( strlen(atomName) > MAX_LENGTH_ATOM_NAME ||
    strlen(atomType) > MAX_LENGTH_ATOM_TYPE ||
    strlen(hbHorA)    > MAX_LENGTH_ATOM_DONOR ||
    strlen(hbDorB) > MAX_LENGTH_ATOM_ACCEPTOR)
  {
    return ValueError;
  }

  pThis->CHARMM_epsilon = epsilon;
  pThis->CHARMM_radius = rmin;
  pThis->CHARMM_charge = charge;

  //AtomSetName(pThis, atomName);
  strcpy(pThis->name, atomName);
  strcpy(pThis->type, atomType);

  if(strcmp(hbHorA, "H") == 0)
  {
    pThis->isHBatomH = TRUE;
  }
  else
  {
    pThis->isHBatomH = FALSE;
  }

  if(strcmp(hbHorA, "A") == 0)
  {
    pThis->isHBatomA = TRUE;
  }
  else
  {
    pThis->isHBatomA = FALSE;
  }

  strcpy(pThis->hbHorA, hbHorA);
  strcpy(pThis->hbDorB, hbDorB);
  strcpy(pThis->hbB2, hbB2);


  if(strcmp(hybrid,"SP2") == 0)
  {
    pThis->hybridType = Type_AtomHybridType_SP2;
  }
  else if(strcmp(hybrid,"SP3") == 0)
  {
    pThis->hybridType = Type_AtomHybridType_SP3;
  }
  else if(strcmp(hybrid,"SP") == 0)
  {
    pThis->hybridType = Type_AtomHybridType_SP;
  }
  else
  {
    pThis->hybridType = Type_AtomHybridType_None;
  }

  switch(isBB[0])
  {
  case 'Y':
    pThis->isBBAtom = TRUE;break;
  case 'N':
    pThis->isBBAtom = FALSE;break;
  default:
    return ValueError;
  }

  if(pThis->isBBAtom == TRUE)
  {
    if(strcmp(atomName, "HN") == 0
      || strcmp(atomName, "HT1") == 0
      || strcmp(atomName, "HT2") == 0
      || strcmp(atomName, "HT3") == 0)
    {
      pThis->bbAtomType = Type_BB_Atom_HN;
    }
    else if(strcmp(atomName, "N") == 0)
    {
      pThis->bbAtomType = Type_BB_Atom_N;
    }
    else if(strcmp(atomName, "HA") == 0
      || strcmp(atomName, "HA1") == 0
      || strcmp(atomName, "HA2") == 0)
    {
      pThis->bbAtomType = Type_BB_Atom_HA;
    }
    else if(strcmp(atomName, "CA") == 0)
    {
      pThis->bbAtomType = Type_BB_Atom_CA;
    }
    else if(strcmp(atomName, "C") == 0)
    {
      pThis->bbAtomType = Type_BB_Atom_C;
    }
    else if(strcmp(atomName, "O") == 0
      || strcmp(atomName, "OXT") == 0)
    {
      pThis->bbAtomType = Type_BB_Atom_O;
    }
  }
  else
  {
    pThis->bbAtomType = Type_NotBBAtom;
  }

  if(strcmp(polar, "P")==0)
  {
    pThis->polarity = Type_AtomPolarity_P;
  }
  else if(strcmp(polar, "C")==0)
  {
    pThis->polarity = Type_AtomPolarity_C;
  }
  else if(strcmp(polar, "NP1")==0)
  {
    pThis->polarity = Type_AtomPolarity_NPAliphatic;
  }
  else if(strcmp(polar, "NP2")==0)
  {
    pThis->polarity = Type_AtomPolarity_NPAromatic;
  }
  // for debug
  //AtomShowParams(pThis);

  return Success;
}



int AtomShowParams(Atom* pThis)
{
  char polarity[32];
  switch(pThis->polarity)
  {
    case Type_AtomPolarity_P:
      strcpy(polarity, "P  ");break;
    case Type_AtomPolarity_C:
      strcpy(polarity, "C  ");break;
    case Type_AtomPolarity_NPAliphatic:
      strcpy(polarity, "NP1");break;
    case Type_AtomPolarity_NPAromatic:
      strcpy(polarity, "NP2");break;
    default:
      strcpy(polarity, "polarity type Unknown");break;
  }

  printf("%5s %5s %c %3s %8.4f %8.4f %8.4f %5s %5s %5s %d\n", 
    pThis->name, pThis->type, 
    pThis->isBBAtom? 'Y':'N', 
    polarity, pThis->CHARMM_epsilon, pThis->CHARMM_radius, pThis->CHARMM_charge, 
    pThis->hbHorA, pThis->hbDorB, pThis->hbB2, pThis->hybridType);
  return Success;
}

char* AtomGetChainName(Atom* pThis)
{
  return pThis->chainName;
}

int AtomSetChainName(Atom* pThis, char* newSegName)
{
  if(strlen(newSegName)>MAX_LENGTH_CHAIN_NAME)
  {
    return ValueError;
  }
  else
  {
    strcpy(pThis->chainName, newSegName);
    return Success;
  }
}

int AtomGetPosInChain(Atom* pThis)
{
  return pThis->posInChain;
}

int AtomSetPosInChain(Atom* pThis, int newSegPos)
{
  pThis->posInChain = newSegPos;
  return Success;
}

int AtomShowInPDBFormat(Atom* pThis, char* header, 
               char* resiName, char* chainName, 
               int atomIndex, int resiIndex, 
               BOOL showHydrogen, FILE* pFile)
{
  char atomName[MAX_LENGTH_ATOM_NAME+1];
  if(pFile==NULL)
  {
    pFile = stdout;
  }

  if(showHydrogen==FALSE && AtomIsHydrogen(pThis))
  {
    return Success;
  }

  // type, serial, name, altLoc, resName, chainID, resSeq, iCode, X, Y, Z
  // 0,  6,    12, 16,   17,    21,    22,   26,  30, 38, 46
  // 6,  5,    4,  1,    4,     1,     4,    1,   8, 8, 8
  //
  // printf("%-6.6s%5d %-5.5s%-4s%1.1s%4d     %7.3f %7.3f %7.3f  %d                    %c\n";

  strcpy(atomName, pThis->name);
  if(strlen(atomName) >= 4)
  {
    //char newAtomName[MAX_LENGTH_ATOM_NAME+1];
    //newAtomName[0] = atomName[3];
    //newAtomName[1] = atomName[0];
    //newAtomName[2] = atomName[1];
    //newAtomName[3] = atomName[2];
    //newAtomName[4] = '\0';
    //strcpy(atomName, newAtomName);

    fprintf(pFile, "%-6.6s%5d %-4.4s %-4s%1.1s%4d     %7.3f %7.3f %7.3f  %d                    %c\n", 
      header, 
      atomIndex, 
      atomName, 
      resiName, 
      chainName, 
      resiIndex, 
      pThis->xyz.X, pThis->xyz.Y, pThis->xyz.Z, 
      pThis->isXyzValid,
      pThis->type[0]);
  }
  else
  {
    // deal with isoleucine;
    if(strcmp(resiName, "ILE")==0 && strcmp(atomName, "CD")==0)
    {
      strcpy(atomName, "CD1");
    }
    fprintf(pFile, "%-6.6s%5d  %-4.4s%-4s%1.1s%4d     %7.3f %7.3f %7.3f  %d                    %c\n", 
      header, 
      atomIndex, 
      atomName, 
      resiName, 
      chainName, 
      resiIndex, 
      pThis->xyz.X, pThis->xyz.Y, pThis->xyz.Z, 
      pThis->isXyzValid,
      pThis->type[0]);
  }

  return Success;
}


int AtomArrayCreate(AtomArray* pThis)
{
  pThis->atoms = NULL;
  pThis->atomNum = 0;
  return Success;
}
void AtomArrayDestroy(AtomArray* pThis)
{
  int i;
  for(i=0;i<pThis->atomNum;i++)
  {
    AtomDestroy(&pThis->atoms[i]);
  }
  free(pThis->atoms);
  pThis->atoms = NULL;
  pThis->atomNum = 0;
}

int AtomArrayCopy(AtomArray* pThis, AtomArray* pOther)
{
  int i;
  AtomArrayDestroy(pThis);
  AtomArrayCreate(pThis);
  pThis->atomNum = pOther->atomNum;
  pThis->atoms = (Atom*)malloc(sizeof(Atom)*pThis->atomNum);
  for(i=0;i<pThis->atomNum;i++)
  {
    AtomCreate(&pThis->atoms[i]);
    AtomCopy(&pThis->atoms[i], &pOther->atoms[i]);
  }
  return Success;
}

int AtomArrayGetCount(AtomArray* pThis)
{
  return pThis->atomNum;
}

Atom* AtomArrayGet(AtomArray* pThis, int index)
{
  if(index<0 || index>=pThis->atomNum)
  {
        return NULL;
  }
  return pThis->atoms + index;
}

Atom* AtomArrayGetByName(AtomArray* pThis, char* atomName)
{
    int index = -1;
    int result;
    result = AtomArrayFind(pThis, atomName, &index);
    if(FAILED(result))
  {
        return NULL;
    }
    else
  {
        return AtomArrayGet(pThis, index);
    }
}

int AtomArrayFind(AtomArray* pThis, char* atomName, int* pIndex)
{
  int i;
  for(i=0;i<pThis->atomNum;i++)
  {
    if(strcmp(AtomGetName(&pThis->atoms[i]), atomName)==0)
    {
      *pIndex = i;
      return Success;
    }
  }
  return DataNotExistError;
}

int AtomArrayInsert(AtomArray* pThis, int index, Atom* pNewAtom)
{
  int newCount;
  int i;
  if(index<0 || index>pThis->atomNum)
  {
    int    result = IndexError;
    return result;
  }
  newCount = pThis->atomNum + 1;
  pThis->atoms = (Atom*)realloc(pThis->atoms, sizeof(Atom)*newCount);
  pThis->atomNum = newCount;

  AtomCreate(&pThis->atoms[newCount-1]);
  for(i=newCount-1;i>index;i--)
  {
    AtomCopy(&pThis->atoms[i], &pThis->atoms[i-1]);
  }
  return AtomCopy(&pThis->atoms[index], pNewAtom);
}

int AtomArrayRemove(AtomArray* pThis, int index)
{
  int i;
  if(index<0 || index>=pThis->atomNum)
  {
    int    result = IndexError;
    return result;
  }
  
  for(i=index;i<pThis->atomNum-1;i++)
  {
    AtomCopy(&pThis->atoms[i], &pThis->atoms[i+1]);
  }
  AtomDestroy(&pThis->atoms[pThis->atomNum-1]);
  (pThis->atomNum)--;
  return Success;
}

int AtomArrayRemoveByName(AtomArray* pThis, char* atomName)
{
    int index = -1;
    int    result;
    result = AtomArrayFind(pThis, atomName, &index);
    if(FAILED(result))
  {
        return result;
    }
    else
  {
        return AtomArrayRemove(pThis, index);
    }
}

int AtomArrayAppend(AtomArray* pThis, Atom* pNewAtom)
{
  return AtomArrayInsert(pThis, AtomArrayGetCount(pThis), pNewAtom);
}

double AtomArrayCalcTotalCharge(AtomArray* pThis)
{
  double totalCharge = 0.0;
  int i;
  for(i=0;i<AtomArrayGetCount(pThis);i++)
  {
    totalCharge += pThis->atoms[i].CHARMM_charge;
  }
  return totalCharge;
}

double AtomArrayCalcMinDistance(AtomArray* pThis, AtomArray* pOther)
{
  int i;
  int j;
  double minDist = 1e8;
  for(i=0; i<pThis->atomNum; i++){
    if(pThis->atoms[i].name[0] == 'H') continue;
    for(j=0; j<pOther->atomNum; j++){
      if(pOther->atoms[j].name[0] == 'H') continue;
      double dist = XYZDistance(&pThis->atoms[i].xyz, &pOther->atoms[j].xyz);
      if(dist < minDist){
        minDist = dist;
      }
    }
  }
  return minDist;
}

BOOL AtomArrayAllAtomXYZAreValid(AtomArray* pThis)
{
    int i;
    for(i=0;i<pThis->atomNum;i++)
  {
        if(pThis->atoms[i].isXyzValid==FALSE)
    {
            return FALSE;
        }
    }
    return TRUE;
}

int AtomArrayShowInPDBFormat(AtomArray* pThis, char* header, 
                char* resiName, char* chainName, 
                int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile)
{
  int i;
  for(i=0;i<AtomArrayGetCount(pThis);i++)
  {
    AtomShowInPDBFormat(&pThis->atoms[i], header, resiName, chainName, atomIndex+i, 
      resiIndex, showHydrogen, pFile);
  }
  return Success;
}

/**************************************************/
/* AtomParamsSet
/**************************************************/
int AtomParamsSetCreate(AtomParamsSet* pThis)
{
    StringArrayCreate(&pThis->residueNames);
    pThis->atomCount = NULL;
    pThis->atoms = NULL;
    return Success;
}

void AtomParamsSetDestroy(AtomParamsSet* pThis)
{
    int i, j;
    for(i=0;i<AtomParamsSetGetResidueCount(pThis);i++)
  {
        for(j=0;j<pThis->atomCount[i];j++)
    {
            AtomDestroy(&pThis->atoms[i][j]);
        }
        free(pThis->atoms[i]);
    }
    free(pThis->atoms);
    free(pThis->atomCount);
    StringArrayDestroy(&pThis->residueNames);
}

int AtomParamsSetAdd(AtomParamsSet* pThis, char* residueName, Atom* pNewAtom)
{
    int    result;
    int resiIndex;
    int atomIndex;

    if(strlen(residueName)>MAX_LENGTH_RESIDUE_NAME)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d, residue name is too long", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
    }

    result = StringArrayFind(&pThis->residueNames, residueName, &resiIndex);
    if(FAILED(result))
  {  // new residue
        int newResidueCount = AtomParamsSetGetResidueCount(pThis) + 1;
        pThis->atoms = (Atom**)realloc(pThis->atoms, sizeof(Atom*)*newResidueCount);
        pThis->atoms[newResidueCount-1] = NULL;

        StringArrayAppend(&pThis->residueNames, residueName);

        pThis->atomCount = (int*) realloc(pThis->atomCount, sizeof(int)*newResidueCount);
        pThis->atomCount[newResidueCount-1] = 0;

        resiIndex = AtomParamsSetGetResidueCount(pThis)-1;
    }

    for(atomIndex=0;atomIndex<pThis->atomCount[resiIndex];atomIndex++)
  {
        if(strcmp(pThis->atoms[resiIndex][atomIndex].name, pNewAtom->name) == 0)
    {
            break;
        }
    }
  // new atom
    if(atomIndex == pThis->atomCount[resiIndex])
  {
        int newAtomCount = pThis->atomCount[resiIndex] + 1;
        pThis->atoms[resiIndex] = (Atom*)realloc(pThis->atoms[resiIndex], sizeof(Atom)*newAtomCount);
        pThis->atomCount[resiIndex]++;
    }

  AtomCopy(&pThis->atoms[resiIndex][atomIndex], pNewAtom);
  //AtomShowParams(&pThis->atoms[resiIndex][atomIndex]);

  return Success;
}

int AtomParamsSetAddFromFile(AtomParamsSet* pThis, char* filePath)
{
    FileReader file;
    int    result;
    char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  char usrMsg[MAX_LENGTH_ERR_MSG+1];

    result = FileReaderCreate(&file, filePath);
  if(FAILED(result))
  {
    sprintf(usrMsg, "in file %s function %s() line %d, when opening: \n%s", 
      __FILE__, __FUNCTION__, __LINE__, filePath);
    TraceError(usrMsg, IOError);
    exit(IOError);
  }

    while(!FAILED(FileReaderGetNextLine(&file, line)))
  {
        Atom atom;
        StringArray params;
        char* resiName;
        char usrMsg[MAX_LENGTH_ERR_MSG+1];

        AtomCreate(&atom);
        StringArrayCreate(&params);
        StringArraySplitString(&params, line, ' ');

        resiName = StringArrayGet(&params, 0);
        if( resiName == NULL   ||
            FAILED(AtomSetParamsByStringArray(&atom, &params)) ||
            FAILED(AtomParamsSetAdd(pThis, resiName, &atom)) )
    {
      sprintf(usrMsg, "in file %s function %s() line %d, when reading:\n%s", 
          __FILE__, __FUNCTION__, __LINE__, line);
      TraceError(usrMsg, FormatError);
      exit(FormatError);
        }

        StringArrayDestroy(&params);
        AtomDestroy(&atom);
    }

    FileReaderDestroy(&file);
    return Success;
}

int AtomParamsSetGetResidueCount(AtomParamsSet* pThis)
{
    return StringArrayGetCount(&pThis->residueNames);
}

int AtomParamsSetGetResidueName(AtomParamsSet* pThis, int index, char* residueName)
{
    char* name = StringArrayGet(&pThis->residueNames, index);
    if(name==NULL)
  {
        return DataNotExistError;
    }
    else
  {
        strcpy(residueName, name);
        return Success;
    }
}

int AtomParamsSetGetAtomCount(AtomParamsSet* pThis, 
                                        char* residueName, int* pCount)
{
    int resiIndex;
    int    result;
    result = StringArrayFind(&pThis->residueNames, residueName, &resiIndex);
    if(FAILED(result))
  {
        return result;
    }
    *pCount = pThis->atomCount[resiIndex];
    return Success;
}

int AtomParamsSetGetAtomParam(AtomParamsSet* pThis, 
                                        char* residueName, int index, Atom* pDestAtom)
{
    int resiIndex;
    int    result;
    result = StringArrayFind(&pThis->residueNames, residueName, &resiIndex);
    if(FAILED(result))
  {
        return result;
    }

    if( index<0 || index>=pThis->atomCount[resiIndex] )
  {
        result = IndexError;
        return result;
    }

    AtomCopy(pDestAtom, &pThis->atoms[resiIndex][index]);
    return Success;  
}

int AtomParamsSetGetAtomParamByName(AtomParamsSet* pThis, 
                                              char* residueName, char* atomName, Atom* pDestAtom)
{
    int resiIndex, atomIndex;
    int    result;
    result = StringArrayFind(&pThis->residueNames, residueName, &resiIndex);
    if(FAILED(result))
  {
        return result;
    }
    for(atomIndex=0;atomIndex<pThis->atomCount[resiIndex];atomIndex++)
{
        if(strcmp(pThis->atoms[resiIndex][atomIndex].name, atomName)==0)
            break;
    }
    if(atomIndex==pThis->atomCount[resiIndex])
  { // atom not found
        result = DataNotExistError;
        return result;
    }
    AtomCopy(pDestAtom, &pThis->atoms[resiIndex][atomIndex]);
    return Success;
}


int AtomAssignEEF1Parameter(Atom *pThis, int eef1_atType, double refDG, double freeDG, double volume, double lamda)
{
  pThis->EEF1_atType = eef1_atType;
  pThis->EEF1_refDG_ = refDG;
  pThis->EEF1_freeDG = freeDG;
  pThis->EEF1_volume = volume;
  pThis->EEF1_lamda_ = lamda;

  return Success;
}
