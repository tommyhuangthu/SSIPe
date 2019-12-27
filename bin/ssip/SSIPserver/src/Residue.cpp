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

#include "Residue.h"
#include <string.h>
#include <ctype.h>

int BondCreate(Bond* pThis)
{
  strcpy(pThis->atomFromName, "");
  strcpy(pThis->atomToName, "");
  pThis->type = Type_Bond_None;
  return Success;
}

void BondDestroy(Bond* pThis)
{
  strcpy(pThis->atomFromName, "");
  strcpy(pThis->atomToName, "");
  pThis->type = Type_Bond_None;
}

int BondCopy(Bond* pThis, Bond* pOther)
{
  *pThis = *pOther;
  return Success;
}

char* BondGetFromName(Bond* pThis)
{
  return pThis->atomFromName;
}

int BondSetFromName(Bond* pThis, char* from)
{
  if(strlen(from)>MAX_LENGTH_ATOM_NAME)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->atomFromName, from);
  return Success;
}

char* BondGetToName(Bond* pThis)
{
  return pThis->atomToName;
}

int BondSetToName(Bond* pThis, char* to)
{
  if(strlen(to)>MAX_LENGTH_ATOM_NAME)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->atomToName, to);
  return Success;
}

Type_Bond BondGetType(Bond* pThis)
{
  return pThis->type;
}

int BondSetType(Bond* pThis, Type_Bond newType)
{
  if(newType<Type_Bond_Single || newType>Type_Bond_None)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  pThis->type = newType;
  return Success;
}

int BondShow(Bond* pThis)
{
  char bondType[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  switch(BondGetType(pThis))
  {
    case Type_Bond_None:
      strcpy(bondType, "NONE  ");break;
    case Type_Bond_Single:
      strcpy(bondType, "BOND  ");break;
    case Type_Bond_double:
      strcpy(bondType, "double");break;
    case Type_Bond_Triple:
      strcpy(bondType, "TRIPLE");break;
    default:
      strcpy(bondType, "");break;
  }
  printf("%s  %s  %s", 
    bondType, BondGetFromName(pThis), 
    BondGetToName(pThis)
    );
  return Success;
}

int BondSetCreate(BondSet* pThis)
{
  pThis->count = 0;
  pThis->bonds = NULL;
  return Success;
}

void BondSetDestroy(BondSet* pThis)
{
  int i;
  for(i=0;i<pThis->count;i++)
  {
    BondDestroy(&pThis->bonds[i]);
  }
  free(pThis->bonds);
  pThis->bonds = NULL;
  pThis->count = 0;
}

int BondSetCopy(BondSet* pThis, BondSet* pOther)
{
  int i;
  BondSetDestroy(pThis);
  BondSetCreate(pThis);
  pThis->count = pOther->count;
  pThis->bonds = (Bond*)malloc(sizeof(Bond)*pOther->count);
  for(i=0;i<pThis->count;i++)
  {
    BondCreate(&pThis->bonds[i]);
    BondCopy(&pThis->bonds[i], &pOther->bonds[i]);
  }
  return Success;
}

int BondSetGetCount(BondSet* pThis)
{
  return pThis->count;
}

int BondSetAdd(BondSet* pThis, char* atom1, char* atom2, Type_Bond bondType)
{
  if( BondSetFind(pThis, atom1, atom2) == Type_Bond_None )
  {
    Bond* pNewlyAddedBond;
    (pThis->count)++;
    pThis->bonds = (Bond*)realloc(pThis->bonds, sizeof(Bond)*pThis->count);
    pNewlyAddedBond = &pThis->bonds[pThis->count-1];
    BondCreate(pNewlyAddedBond);
    if( FAILED(BondSetFromName(pNewlyAddedBond, atom1)) ||
      FAILED(BondSetToName(pNewlyAddedBond, atom2))  ||
      FAILED(BondSetType(pNewlyAddedBond, bondType)) )
    {
      return ValueError;
    }
  }
  return Success;
}

int BondSetRemove(BondSet* pThis, char* atom1, char* atom2)
{
  int i;
  for(i=0;i<pThis->count;i++)
  {
    Bond* pCurBond = &pThis->bonds[i];
    if( (strcmp(BondGetFromName(pCurBond), atom1)==0 && strcmp(BondGetToName(pCurBond), atom2)==0) ||
      (strcmp(BondGetFromName(pCurBond), atom2)==0 && strcmp(BondGetToName(pCurBond), atom1)==0) )
    {
      /* Remove this bond, swap this bond with the last one in the set, and decrease the bond counter */
      /* Because these is at least one bond when reaching here, the operations are safe */
      BondCopy(pCurBond, &pThis->bonds[pThis->count-1]); 
      BondDestroy(&pThis->bonds[pThis->count-1]);
      (pThis->count)--;
      return Success;
    }
  }
  return DataNotExistError;
}

Type_Bond BondSetFind(BondSet* pThis, char* atom1, char* atom2)
{
  int i;
  for(i=0;i<pThis->count;i++)
  {
    Bond* pCurBond = &pThis->bonds[i];
    if( (strcmp(BondGetFromName(pCurBond), atom1)==0 && strcmp(BondGetToName(pCurBond), atom2)==0) ||
      (strcmp(BondGetFromName(pCurBond), atom2)==0 && strcmp(BondGetToName(pCurBond), atom1)==0) )
    {
      return BondGetType(pCurBond);
    }
  }
  return Type_Bond_None;
}
Bond* BondSetGet(BondSet* pThis, int index)
{
  if(index<0 || index>=pThis->count)
  {
    return NULL;
  }
  return &pThis->bonds[index];
}
int BondSetShow(BondSet* pThis)
{
  int i;
  for(i=0;i<pThis->count;i++)
  {
    BondShow(&pThis->bonds[i]);
    printf("\n");
  }
  return Success;
}
  

int CharmmICCreate(CharmmIC* pThis)
{
  int i;
  for(i=0;i<4;i++)
    strcpy(pThis->atomNames[i], "");
  return Success;
}
int CharmmICCreateFromStringArray(CharmmIC* pThis, StringArray* params)
{
  int i;
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];

  CharmmICCreate(pThis);
  if(StringArrayGetCount(params) != 10)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }

  for(i=0;i<4;i++)
  {
    if(strlen(StringArrayGet(params, i+1))>MAX_LENGTH_ATOM_NAME)
    {
      char usrMsg[MAX_LENGTH_ERR_MSG+1];
      int errorCode = ValueError;
      sprintf(usrMsg, "in file %s function %s() line %d, name is too long", 
        __FILE__, __FUNCTION__, __LINE__);
      TraceError(usrMsg, errorCode);
      return errorCode;
    }
    strcpy(pThis->atomNames[i], StringArrayGet(params, i+1));
  }

  for(i=0;i<5;i++)
  {
    pThis->icParam[i] = atof(StringArrayGet(params, i+5));
        if(i>=1 && i<=3)
    {
      pThis->icParam[i] = DegToRad(pThis->icParam[i]); /* Convert degree to arc */
        }  
  }

  strcpy(buffer, CharmmICGetAtomC(pThis));
  if(buffer[0]=='*')
  {
    pThis->torsionProperFlag = FALSE;   /* Improper torsion */
    strcpy(pThis->atomNames[2], buffer+1);  /* Discard the '*' in atomC's name */
  }
  else
  {
    pThis->torsionProperFlag = TRUE;   /* Proper torsion */
  }
  return Success;
}

void CharmmICDestroy(CharmmIC* pThis)
{
  int i;
  for(i=0;i<4;i++)
    strcpy(pThis->atomNames[i], "");
}

int CharmmICCopy(CharmmIC* pThis, CharmmIC* pOther)
{
  *pThis = *pOther;
  return Success;
}

char* CharmmICGetAtomA(CharmmIC* pThis)
{
  return pThis->atomNames[0];
}

char* CharmmICGetAtomB(CharmmIC* pThis)
{
  return pThis->atomNames[1];
}

char* CharmmICGetAtomC(CharmmIC* pThis)
{
  return pThis->atomNames[2];
}

char* CharmmICGetAtomD(CharmmIC* pThis)
{
  return pThis->atomNames[3];
}

double* CharmmICGetICParams(CharmmIC* pThis)
{
  return pThis->icParam;
}

int CharmmICSetTorsionAngle(CharmmIC* pThis, double newTorsionAngle)
{
    pThis->icParam[2] = newTorsionAngle;
    return Success;
}

BOOL CharmmICGetTorsionProperFlag(CharmmIC* pThis)
{
  return pThis->torsionProperFlag;
}

int CharmmICCalcXYZ(CharmmIC* pThis, AtomArray* atomArray, XYZ* pDestXYZ)
{
    int i;
    XYZ* pAtomsABC[3];

    for(i=0;i<3;i++)
  {
        Atom* pAtom = AtomArrayGetByName(atomArray, pThis->atomNames[i]);
        if(pAtom == NULL)
    {
            return DataNotExistError;
        }
        
        if(pAtom->isXyzValid==FALSE)
    {
            return DataNotExistError;
        }

        pAtomsABC[i] = &pAtom->xyz;
    }
    return GetFourthAtom(pAtomsABC[0], pAtomsABC[1], 
        pAtomsABC[2], pThis->icParam, pDestXYZ);
}

int CharmmICShow(CharmmIC* pThis)
{
  int i;
  printf("IC  %-7.7s %-7.7s %c%-7.7s %-7.7s  ", 
    CharmmICGetAtomA(pThis), CharmmICGetAtomB(pThis), 
    CharmmICGetTorsionProperFlag(pThis) ? ' ' : '*', 
    CharmmICGetAtomC(pThis), CharmmICGetAtomD(pThis)
    );
  printf("%8.3f ", pThis->icParam[0]);
  for(i=1;i<4;i++)
  {
    printf("%8.3f ", RadToDeg(pThis->icParam[i]));
  }
  printf("%8.3f", pThis->icParam[4]);
  return Success;
}


int ResidueTopologyCreate(ResidueTopology* pThis)
{
  StringArrayCreate(&pThis->atoms);
  StringArrayCreate(&pThis->deletes);
  BondSetCreate(&pThis->bonds);
  pThis->icCount = 0;
  pThis->ics = NULL;
  return Success;
}

int ResidueTopologyCreateFromFileReader(ResidueTopology* pThis, FileReader* pFileReader)
{
  BOOL residueNameAlreadySet = FALSE;

  ResidueTopologyCreate(pThis);

  while(!FileReaderEndOfFile(pFileReader))
  {
    char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char* keyword;
    StringArray wordsInLine;
    StringArrayCreate(&wordsInLine);

    FileReaderGetNextLine(pFileReader, line);
    StringArraySplitString(&wordsInLine, line, ' ');
    keyword = StringArrayGet(&wordsInLine, 0);

    if(strcmp(keyword, "RESI")==0 || strcmp(keyword, "PRES")==0)
    {
      if(residueNameAlreadySet)
      {
        FileReaderSetCurrentPos(pFileReader, FileReaderGetCurrentPos(pFileReader)-1);
        /* Don't forget to destroy 'wordsInLine' before quit: */
        StringArrayDestroy(&wordsInLine);
        return Success;
      }
      else
      {
        char* residueName = StringArrayGet(&wordsInLine, 1);
        if(  residueName == NULL  ||
           FAILED(ResidueTopologySetName(pThis, residueName))  )
        {
          char usrMsg[MAX_LENGTH_ERR_MSG+1];
          int errorCode = ValueError;
          sprintf(usrMsg, "in file %s function %s() line %d", 
            __FILE__, __FUNCTION__, __LINE__);
          TraceError(usrMsg, errorCode);
          return errorCode;
        }
        residueNameAlreadySet = TRUE;
      }
    }
    else if(strcmp(keyword, "ATOM")==0)
    {
      char* atomName = StringArrayGet(&wordsInLine, 1);
      if(  atomName == NULL   ||
         FAILED(StringArrayAppend(&pThis->atoms, atomName))  )
      {
        char usrMsg[MAX_LENGTH_ERR_MSG+1];
        int errorCode = ValueError;
        sprintf(usrMsg, "in file %s function %s() line %d", 
          __FILE__, __FUNCTION__, __LINE__);
        TraceError(usrMsg, errorCode);
        return errorCode;
      }
    }
    else if(strcmp(keyword, "DELETE")==0)
    {
      char* deleteAtomName = StringArrayGet(&wordsInLine, 2);
      if(  deleteAtomName == NULL ||
         FAILED(StringArrayAppend(&pThis->deletes, deleteAtomName))  )
      {
        char usrMsg[MAX_LENGTH_ERR_MSG+1];
        int errorCode = ValueError;
        sprintf(usrMsg, "in file %s function %s() line %d", 
          __FILE__, __FUNCTION__, __LINE__);
        TraceError(usrMsg, errorCode);
        return errorCode;
      }
    }
    else if( strcmp(keyword, "BOND")==0||strcmp(keyword, "DOUBLE")==0||strcmp(keyword, "TRIPLE")==0)
    {
      int i;
      for(i=1;i<StringArrayGetCount(&wordsInLine);i+=2)
      {
        char* atomFromName = StringArrayGet(&wordsInLine, i);
        char* atomToName   = StringArrayGet(&wordsInLine, i+1);
        Type_Bond bondType;

        if(  atomFromName == NULL ||
           atomToName == NULL   )
        {
          char usrMsg[MAX_LENGTH_ERR_MSG+1];
          int errorCode = ValueError;
          sprintf(usrMsg, "in file %s function %s() line %d", 
            __FILE__, __FUNCTION__, __LINE__);
          TraceError(usrMsg, errorCode);
          return errorCode;
        }

        if(strcmp(keyword, "BOND")==0)
        {
          bondType = Type_Bond_Single;
        }
        else if(strcmp(keyword, "DOUBLE")==0)
        {
          bondType = Type_Bond_double;
        }
        else if(strcmp(keyword, "TRIPLE")==0)
        { 
          bondType = Type_Bond_Triple;
        }
        else
        {
          char usrMsg[MAX_LENGTH_ERR_MSG+1];
          int errorCode = ValueError;
          sprintf(usrMsg, "in file %s function %s() line %d", 
            __FILE__, __FUNCTION__, __LINE__);
          TraceError(usrMsg, errorCode);
          return errorCode;
        }

        if(FAILED( BondSetAdd(&pThis->bonds, atomFromName, atomToName, bondType)) )
        {
          char usrMsg[MAX_LENGTH_ERR_MSG+1];
          int errorCode = ValueError;
          sprintf(usrMsg, "in file %s function %s() line %d", 
            __FILE__, __FUNCTION__, __LINE__);
          TraceError(usrMsg, errorCode);
          return errorCode;
        }
      }
    }
    else if(strcmp(keyword, "IC")==0)
    {
      CharmmIC newIC;
      if( FAILED(CharmmICCreateFromStringArray(&newIC, &wordsInLine)) ||
        FAILED(ResidueTopologyAddCharmmIC(pThis, &newIC)))
      {
        char usrMsg[MAX_LENGTH_ERR_MSG+1];
        int errorCode = ValueError;
        sprintf(usrMsg, "in file %s function %s() line %d", 
          __FILE__, __FUNCTION__, __LINE__);
        TraceError(usrMsg, errorCode);
        return errorCode;
      }
      CharmmICDestroy(&newIC);
      
    }

    else
    {
      /* Unrecognized keyword in ResidueTopology file */
    }

  StringArrayDestroy(&wordsInLine);
  }
  return Success;

}

void ResidueTopologyDestroy(ResidueTopology* pThis)
{
  int i;
  StringArrayDestroy(&pThis->atoms);
  StringArrayDestroy(&pThis->deletes);
  BondSetDestroy(&pThis->bonds);

  for(i=0;i<pThis->icCount;i++)
  {
    CharmmICDestroy(&pThis->ics[i]);
  }
  free(pThis->ics);
  pThis->ics = NULL;
  pThis->icCount = 0;

}

int ResidueTopologyCopy(ResidueTopology* pThis, ResidueTopology* pOther)
{
  int i;

  ResidueTopologyDestroy(pThis);
  ResidueTopologyCreate(pThis);

  ResidueTopologySetName(pThis, ResidueTopologyGetName(pOther));

  StringArrayCopy(&pThis->atoms, &pOther->atoms);
  StringArrayCopy(&pThis->deletes, &pOther->deletes);
  BondSetCopy(&pThis->bonds, &pOther->bonds);

  pThis->icCount = pOther->icCount;
  pThis->ics = (CharmmIC*)malloc(sizeof(CharmmIC)*pThis->icCount);
  for(i=0;i<pThis->icCount;i++)
  {
    CharmmICCreate(&pThis->ics[i]);
    CharmmICCopy(&pThis->ics[i], &pOther->ics[i]);
  }
  return Success;
}

char* ResidueTopologyGetName(ResidueTopology* pThis)
{
  return pThis->residueName;
}

int ResidueTopologySetName(ResidueTopology* pThis, char* newName)
{
  if(strlen(newName)>MAX_LENGTH_RESIDUE_NAME)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d, name is too long", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->residueName, newName);
  return Success;
}

StringArray* ResidueTopologyGetAtoms(ResidueTopology* pThis)
{
  return &pThis->atoms;
}

StringArray* ResidueTopologyGetDeletes(ResidueTopology* pThis)
{
  return &pThis->deletes;
}

BondSet* ResidueTopologyGetBonds(ResidueTopology* pThis)
{
  return &pThis->bonds;
}

int ResidueTopologyGetCharmmICCount(ResidueTopology* pThis)
{
    return pThis->icCount;
}

int ResidueTopologyGetCharmmIC(ResidueTopology* pThis, int index, CharmmIC* pDestIC)
{
    if(index<0 || index>=pThis->icCount)
        return IndexError;
    CharmmICCopy(pDestIC, &pThis->ics[index]);
    return Success;
}

int ResidueTopologyFindCharmmIC(ResidueTopology* pThis, char* atomDName, CharmmIC* pDestIC)
{
  int i;
  for(i=0;i<pThis->icCount;i++)
  {
    if(strcmp(CharmmICGetAtomD(&pThis->ics[i]), atomDName)==0)
    {
      CharmmICCopy(pDestIC, &pThis->ics[i]);
      return Success;
    }
  }
  return DataNotExistError;
}

int ResidueTopologyAddCharmmIC(ResidueTopology* pThis, CharmmIC* pNewIC)
{
  /* Will not check if it already exists*/
  int newICCount = pThis->icCount+1;
  pThis->ics = (CharmmIC*)realloc(pThis->ics, sizeof(CharmmIC)*newICCount);
  CharmmICCreate(&pThis->ics[newICCount-1]);
  CharmmICCopy(&pThis->ics[newICCount-1], pNewIC);
  pThis->icCount = newICCount;
  return Success;
}

int ResidueTopologyCalcAtomArrayXYZ(ResidueTopology* pThis, AtomArray* pAtomArray)
{
    int atomCount;
    int i;
    BOOL done;

    CharmmIC ic;
    CharmmICCreate(&ic);
    atomCount = AtomArrayGetCount(pAtomArray);

    done = FALSE;
    while(!done)
  {
        done = TRUE; /* If no new XYZ has been calculated in this 'while' loop, it will terminate */
        for(i=0;i<atomCount;i++)
    {
            Atom* pAtom = AtomArrayGet(pAtomArray, i);
            XYZ newXYZforAtom;

            if(pAtom->isXyzValid)
      {
                continue;
            }

            if( FAILED(ResidueTopologyFindCharmmIC(pThis, AtomGetName(pAtom), &ic)) ||
                FAILED(CharmmICCalcXYZ(&ic, pAtomArray, &newXYZforAtom)) )
      {
                 continue; /* continue 'for' */
            }
            pAtom->xyz = newXYZforAtom;
            pAtom->isXyzValid = TRUE;

            done = FALSE; /* New XYZ has been calculated, hope to find more in next 'while' loop*/
        }
    }

    CharmmICDestroy(&ic);
    return TRUE;
}

int ResidueTopologyShow(ResidueTopology* pThis)
{
  int i;

  printf("RESI  %s\n", ResidueTopologyGetName(pThis));
  for(i=0;i<StringArrayGetCount(&pThis->atoms);i++)
  {
    printf("ATOM   %s\n", StringArrayGet(&pThis->atoms, i));
  }
  BondSetShow(&pThis->bonds);
  for(i=0;i<pThis->icCount;i++)
  {
    CharmmICShow(&pThis->ics[i]);
    printf("\n");
  }
  return Success;
}


int ResiTopoCollectionCreate(ResiTopoSet* pThis)
{
  pThis->count = 0;
  pThis->topos = NULL;
  return Success;
}

void ResiTopoCollectionDestroy(ResiTopoSet* pThis)
{
  int i;
  for(i=0;i<pThis->count;i++)
  {
    ResidueTopologyDestroy(&pThis->topos[i]);
  }
  free(pThis->topos);
  pThis->topos = NULL;
  pThis->count = 0;
}

int ResiTopoCollectionCopy(ResiTopoSet* pThis, ResiTopoSet* pOther)
{
  int i;
  ResiTopoCollectionDestroy(pThis);
  pThis->count = pOther->count;
  pThis->topos = (ResidueTopology*)malloc(sizeof(ResidueTopology)*pThis->count);
  for(i=0;i<pThis->count;i++)
  {
    ResidueTopologyCreate(&pThis->topos[i]);
    ResidueTopologyCopy(&pThis->topos[i], &pOther->topos[i]);
  }
  return Success;
}

int ResiTopoCollectionGet(ResiTopoSet* pThis, char* resiName, ResidueTopology* pDestTopo)
{
  int i;
  for(i=0;i<pThis->count;i++)
  {
    if(strcmp(ResidueTopologyGetName(&pThis->topos[i]), resiName)==0)
    {
      ResidueTopologyCopy(pDestTopo, &pThis->topos[i]);
      return Success;
    }
  }
  return DataNotExistError;
}
int ResiTopoCollectionAdd(ResiTopoSet* pThis, ResidueTopology* pNewTopo)
{
  /* If already exist, replace the original */
  int i;
  int newCount;
  for(i=0;i<pThis->count;i++)
  {
    if(strcmp(ResidueTopologyGetName(&pThis->topos[i]), 
          ResidueTopologyGetName(pNewTopo))==0)
    {
      ResidueTopologyCopy(&pThis->topos[i], pNewTopo);
      return Success;
    }
  }
  
  newCount = pThis->count+1;
  pThis->topos = (ResidueTopology*)realloc(pThis->topos, sizeof(ResidueTopology)*newCount);
  ResidueTopologyCreate(&pThis->topos[newCount-1]);
  ResidueTopologyCopy(&pThis->topos[newCount-1], pNewTopo);
  pThis->count = newCount;
  return Success;
}

int ResiTopoCollectionAddFromFile(ResiTopoSet* pThis, char* filepath)
{
  FileReader file;
  int    result;
  
  result = FileReaderCreate(&file, filepath);
  if(FAILED(result))
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, result);
    return result;
  }

  while(!FileReaderEndOfFile(&file))
  {
    ResidueTopology newTopo;
    result = ResidueTopologyCreateFromFileReader(&newTopo, &file);
    if(FAILED(result))
    {
      char usrMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(usrMsg, "in file %s function %s() line %d", 
        __FILE__, __FUNCTION__, __LINE__);
      TraceError(usrMsg, result);
      return result;
    }
    ResiTopoCollectionAdd(pThis, &newTopo);
    ResidueTopologyDestroy(&newTopo);
  }
  FileReaderDestroy(&file);
  return Success;
}

int ResiTopoCollectionShow(ResiTopoSet* pThis)
{
  int i;
  for(i=0;i<pThis->count;i++)
  {
    printf("Residue No. %d\n", i);
    ResidueTopologyShow(&pThis->topos[i]);
  }
  return Success;
}

int ResidueCreate(Residue* pThis)
{
  strcpy(pThis->name, "");
  strcpy(pThis->chainName, "");
  pThis->posInChain = -1;
  pThis->designSiteFlag = FALSE;
  AtomArrayCreate(&pThis->atoms);
  StringArrayCreate(&pThis->patches);
  BondSetCreate(&pThis->bonds);
  pThis->residueTerminalFlag = Type_ResidueIsNotTerminal;
  pThis->numCBwithin8AtoCurResidue = 0;

  AtomCreate(&pThis->pseudoAtom);
  return Success;
}

void ResidueDestroy(Residue* pThis)
{
  BondSetDestroy(&pThis->bonds);
  AtomArrayDestroy(&pThis->atoms);
  StringArrayDestroy(&pThis->patches);

  AtomDestroy(&pThis->pseudoAtom);
}

int ResidueCopy(Residue* pThis, Residue* pOther)
{
  ResidueDestroy(pThis);
  strcpy(pThis->name, pOther->name);
  strcpy(pThis->chainName, pOther->chainName);
  pThis->posInChain = pOther->posInChain;
  pThis->designSiteFlag = pOther->designSiteFlag;
  pThis->numCBwithin8AtoCurResidue = pOther->numCBwithin8AtoCurResidue;
  pThis->DesignType = pOther->DesignType;

  AtomArrayCreate(&pThis->atoms);
  AtomArrayCopy(&pThis->atoms, &pOther->atoms);

  StringArrayCreate(&pThis->patches);
  StringArrayCopy(&pThis->patches, &pOther->patches);

  BondSetCopy(&pThis->bonds, &pOther->bonds);

  AtomCopy(&pThis->pseudoAtom, &pOther->pseudoAtom);
  pThis->residueTerminalFlag = pOther->residueTerminalFlag;

  return Success;
}

char* ResidueGetName(Residue* pThis)
{
  return pThis->name;
}

int ResidueSetName(Residue* pThis, char* newName)
{
  if(strlen(newName)>MAX_LENGTH_RESIDUE_NAME)
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

char* ResidueGetChainName(Residue* pThis)
{
  return pThis->chainName;
}

int ResidueSetSegName(Residue* pThis, char* newSegName)
{
  int i;
  if(strlen(newSegName)>MAX_LENGTH_CHAIN_NAME)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->chainName, newSegName);
  for(i=0;i<ResidueGetAtomCount(pThis);i++)
  {
    Atom* pCurAtom = ResidueGetAtom(pThis, i);
    AtomSetChainName(pCurAtom, newSegName);
  }
  return Success;
}

int ResidueGetPosInChain(Residue* pThis)
{
  return pThis->posInChain;
}

int ResidueSetPosInChain(Residue* pThis, int newSegPos)
{
  int i;
  if(newSegPos < 0)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  pThis->posInChain = newSegPos;
  for(i=0;i<ResidueGetAtomCount(pThis);i++)
  {
    Atom* pCurAtom = ResidueGetAtom(pThis, i);
    AtomSetPosInChain(pCurAtom, newSegPos);
  }

  return Success;
}

int ResidueSetDesignSiteFlag(Residue* pThis, BOOL newFlag)
{
  pThis->designSiteFlag = newFlag;
  return Success;
}

double ResidueGetCharge(Residue* pThis)
{
  return AtomArrayCalcTotalCharge(&pThis->atoms);
}

int ResidueGetPolarity(Residue* pThis, Type_ResiduePolarity* pPolarity)
{
  // if the total charge is not zero, the type is Type_ResiduePolarity_Charged;
  // if the total charge is zero, but any non-root atom has non-zero charge, the type is Type_ResiduePolarity_Polar;
  // if the total charge is zero and all non-root atoms have zero charge, the type is Type_ResiduePolarity_NonPolar;
  int i; 
  BOOL nonRootAtomCharged;
  BOOL residueCharged;
  
  nonRootAtomCharged = FALSE;
  for(i=0;i<ResidueGetAtomCount(pThis);i++)
  {
    Atom* pAtom;
    BOOL nonRoot;
    BOOL charged;

    pAtom = ResidueGetAtom(pThis, i);
    nonRoot = !(pAtom->isBBAtom);
    charged = fabs(pAtom->CHARMM_charge)>MIN_ZERO_TOLERANCE;

    if(nonRoot && charged)
      nonRootAtomCharged = TRUE;
  }

  residueCharged = fabs(ResidueGetCharge(pThis))>MIN_ZERO_TOLERANCE;

  if(residueCharged)
  {
    *pPolarity = Type_ResiduePolarity_Charged;
  }
  else if(nonRootAtomCharged)
  {
    *pPolarity = Type_ResiduePolarity_Polar;
  }
  else
  {
    *pPolarity = Type_ResiduePolarity_NonPolar;
  }

  return Success;
}

int ResidueGetAtomCount(Residue* pThis)
{
  return AtomArrayGetCount(&pThis->atoms);
}

Atom* ResidueGetAtom(Residue* pThis, int index)
{
    return AtomArrayGet(&pThis->atoms, index);
}

Atom* ResidueGetAtomByName(Residue* pThis, char* atomName)
{
  int index;
  int result;
  result = ResidueFindAtom(pThis, atomName, &index);
  if(FAILED(result))
  {
    return NULL;
  }
    else
  {
        return ResidueGetAtom(pThis, index);
    }
}

int ResidueFindAtom(Residue* pThis, char* atomName, int* pIndex)
{
  return AtomArrayFind(&pThis->atoms, atomName, pIndex);
}

int ResidueGetAtomXYZ(Residue* pThis, char* atomName, XYZ* pXYZ)
{
  Atom* pAtom = ResidueGetAtomByName(pThis, atomName);
    
  if( pAtom == NULL ||
        pAtom->isXyzValid==FALSE)
  {
    return DataNotExistError;
    }

  *pXYZ = pAtom->xyz;
  return Success;
}

AtomArray* ResidueGetAllAtoms(Residue* pThis)
{
  return &pThis->atoms;
}

int ResidueInsertAtom(Residue* pThis, int newIndex, Atom* pNewAtom)
{
  int index;

  if(newIndex < 0 || newIndex > AtomArrayGetCount(&pThis->atoms))
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = IndexError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  else if(FAILED(ResidueFindAtom(pThis, AtomGetName(pNewAtom), &index)))
  {
    return AtomArrayInsert(&pThis->atoms, newIndex, pNewAtom);
  }
  else
  {
    return AtomCopy(ResidueGetAtom(pThis, index), pNewAtom);
  }
}

int ResidueAddAtom(Residue* pThis, Atom* pNewAtom)
{
  int index;
  if(FAILED(ResidueFindAtom(pThis, AtomGetName(pNewAtom), &index)))
  {
    return AtomArrayAppend(&pThis->atoms, pNewAtom);
  }else
  {
    return AtomCopy(ResidueGetAtom(pThis, index), pNewAtom);
  }
}

int ResidueAddAtomsFromAtomParams(Residue* pThis, AtomParamsSet* pAtomParams)
{
  int    result;
  int count, i, j;
  int bbAtomCount = 0;

  result = AtomParamsSetGetAtomCount(pAtomParams, ResidueGetName(pThis), &count);
  if(FAILED(result))
  {
    return result;
  }

  // first add the backbone atoms
  for(i=0;i<count;i++)
  {
    BOOL atomExist = FALSE;
    Atom newAtom;
    AtomCreate(&newAtom);
    AtomParamsSetGetAtomParam(pAtomParams, ResidueGetName(pThis), i, &newAtom);
    for(j = 0; j < ResidueGetAtomCount(pThis); j++)
    {
      Atom *pAtom1 = ResidueGetAtom(pThis, j);
      if(strcmp(pAtom1->name, newAtom.name) == 0)
      {
        Atom tempAtom;
        AtomCreate(&tempAtom);
        AtomCopy(&tempAtom, pAtom1);
        AtomCopy(pAtom1, &newAtom);
        pAtom1->posInChain = tempAtom.posInChain;
        pAtom1->xyz = tempAtom.xyz;
        pAtom1->isXyzValid = tempAtom.isXyzValid;
        atomExist = TRUE;
        AtomDestroy(&tempAtom);
        break;
      }
    }
    if(atomExist == FALSE)
    {
      if(newAtom.isBBAtom == TRUE)
      {
        newAtom.posInChain = pThis->posInChain;
        result = ResidueInsertAtom(pThis, 0, &newAtom);
        bbAtomCount++;
      }
      else
      {
        newAtom.posInChain = pThis->posInChain;
        result = ResidueAddAtom(pThis, &newAtom);
      }
    }
    AtomDestroy(&newAtom);
  }
  pThis->bbAtomCount = bbAtomCount;

  return Success;

}
int ResidueDeleteAtom(Residue* pThis, char* atomName)
{
  int index;
  if(FAILED(ResidueFindAtom(pThis, atomName, &index)))
  {
    return DataNotExistError;
  }
  else
  {
    return AtomArrayRemove(&pThis->atoms, index);
  }
}


int ResidueReadXYZFromPDB(Residue* pThis, FileReader* pPDBFileReader, 
              AtomParamsSet* pAtomParam, 
              ResiTopoSet* pTopos)
{
  // the following three lines show the format of a pdb file;
  // type, serial, name, altLoc, resName, chainID, resSeq, ..., X,   Y,   Z
  // 0,    6,      12,   16,     17,      21,      22,     ..., 30,  38,  46
  //   6,    5,       4,    1,      4,       1,       5,   ...,    8,   8,  8
  BOOL firstLine = TRUE;
    char resSeqValue[MAX_LENGTH_ONE_LINE_IN_FILE+1] = "UNKNOWN";
  char line[MAX_LENGTH_ONE_LINE_IN_FILE+1];

  //strcpy(resSeqValue, "UNKNOWN");
  //resSeqValue[strlen(resSeqValue)]='\0';
  while(!FAILED(FileReaderGetNextLine(pPDBFileReader, line)))
  {
    char strType[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strAtomName[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strResName[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strResSeq[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strX[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strY[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    char strZ[MAX_LENGTH_ONE_LINE_IN_FILE+1];
    Atom* pAtom;
    XYZ  xyz;
    int  result;
    char usrMsg[MAX_LENGTH_ERR_MSG+1];

    ExtractTargetStringFromSourceString(strType, line, 0, 6);
    if(strcmp(strType, "ENDMDL")==0)
    {
      return Success;
    }
    if(strcmp(strType, "ATOM")!=0 && strcmp(strType, "HETATM")!=0)
    {
      continue;
    }

    ExtractTargetStringFromSourceString(strAtomName, line, 12, 4);
    ExtractTargetStringFromSourceString(strResName, line, 17, 4);
    ExtractTargetStringFromSourceString(strResSeq, line, 22, 5);
    ExtractTargetStringFromSourceString(strX, line, 30, 8);
    ExtractTargetStringFromSourceString(strY, line, 38, 8);
    ExtractTargetStringFromSourceString(strZ, line, 46, 8);

    // if the residue is HIS, use HSD or HSE
    if(strcmp(strResName, "HIS")==0  &&  strcmp(pThis->name, "HSD")==0)
    {
      strcpy(strResName, "HSD");
    }
    else if(strcmp(strResName, "HIS")==0  &&  strcmp(pThis->name, "HSE")==0)
    {
            strcpy(strResName, "HSE");
        }
    // the CD1 atom of ILE in pdb is altered into CD
    if(strcmp(strResName, "ILE")==0 && strcmp(strAtomName, "CD1")==0)
    {
      strcpy(strAtomName, "CD");
    }
    // skip the hydrogen atoms
    if(strAtomName[0]=='H')
      continue;
    // sometimes the hydrogen atom can have a name like '1HG1' in residue VAL
    else if(isdigit(strAtomName[0]) && strAtomName[1]=='H')
      continue;

    if(firstLine)
    {
      strcpy(resSeqValue, strResSeq);
      firstLine = FALSE;
    }
    else
    {
      // read the line of a new residue
      if(strcmp(resSeqValue, strResSeq) != 0)
      {
        FileReaderSetCurrentPos(pPDBFileReader, FileReaderGetCurrentPos(pPDBFileReader)-1);
        return Success;
      }
    }

    if(strcmp(strAtomName, "OXT") == 0)
    {
      //ResiduePatchNTERorCTER(pThis, "CTER", pAtomParam, pTopos);
      //continue;
      Atom atomOXT;
      AtomCreate(&atomOXT);
      AtomCopy(&atomOXT, ResidueGetAtomByName(pThis, "O"));
      atomOXT.xyz.X = atof(strX);
      atomOXT.xyz.Y = atof(strY);
      atomOXT.xyz.Z = atof(strZ);
      atomOXT.isXyzValid = TRUE;
      atomOXT.CHARMM_charge = -0.55;
      strcpy(atomOXT.name, "OXT");
      ResidueGetAtomByName(pThis, "O")->CHARMM_charge = -0.55;
      ResidueGetAtomByName(pThis, "C")->CHARMM_charge = 0.1;
      AtomArrayAppend(&pThis->atoms, &atomOXT);

      BondSetAdd(&pThis->bonds, "C", "OXT", Type_Bond_double);
      AtomDestroy(&atomOXT);
    }
    else
    {
      pAtom = ResidueGetAtomByName(pThis, strAtomName);
      if(pAtom == NULL)
      {
        result = FormatError;
        sprintf(usrMsg, "in file %s function %s() line %d, "
          "residue %s do not have atom %s when reading line\n%s", 
          __FILE__, __FUNCTION__, __LINE__, 
          strResName, strAtomName, line);
        TraceError(usrMsg, result);
        continue;
      }

      xyz.X = atof(strX); xyz.Y = atof(strY); xyz.Z = atof(strZ);
      pAtom->xyz = xyz;
      pAtom->isXyzValid = TRUE;
    }
  }
  return Success;
}

BondSet* ResidueGetBonds(Residue* pThis)
{
  return &pThis->bonds;
}
int ResidueAddBondsFromResiTopos(Residue* pThis, ResiTopoSet* pResiTopoCollection)
{
  ResidueTopology resiTopo;
  ResidueTopologyCreate(&resiTopo);
  if(FAILED(ResiTopoCollectionGet(pResiTopoCollection, pThis->name, &resiTopo)))
  {
    return DataNotExistError;
  }
  BondSetCopy(&pThis->bonds, ResidueTopologyGetBonds(&resiTopo));
  ResidueTopologyDestroy(&resiTopo);
  return Success;
}

int ResidueShowInPDBFormat(Residue* pThis, char* header, char* chainName, 
                                   int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile)
{
  AtomArrayShowInPDBFormat(&pThis->atoms, header, ResidueGetName(pThis), 
              chainName, atomIndex, resiIndex, showHydrogen, pFile);
  return Success;
}


int ResiduePatch(Residue* pThis, char* patchName, 
             AtomParamsSet* pAtomParam, 
             ResiTopoSet* pTopos)
{
  int    result;
  StringArray* deleteAtoms;
  char resiName[MAX_LENGTH_RESIDUE_NAME];
  int i;
  
  BondSet* pPatchBonds;
  BondSet newBonds;

  ResidueTopology patchResiTopo;
  ResidueTopologyCreate(&patchResiTopo);

  // delete old atoms
  result  = ResiTopoCollectionGet(pTopos, patchName, &patchResiTopo);
  if(FAILED(result))
    return result;
  deleteAtoms = ResidueTopologyGetDeletes(&patchResiTopo);
  for(i=0;i<StringArrayGetCount(deleteAtoms);i++)
  {
    ResidueDeleteAtom(pThis, StringArrayGet(deleteAtoms, i));
  }

  // Temporarily reset its name to that of the patch residue, to find topology in the topology collection
  strcpy(resiName, ResidueGetName(pThis));
  ResidueSetName(pThis, patchName);
  result = ResidueAddAtomsFromAtomParams(pThis, pAtomParam);
  ResidueSetName(pThis, resiName);

  if(FAILED(result))
    return result;

  // new patches are stored at the head
  StringArrayInsert(&pThis->patches, 0, patchName);

  // deal with the bonds
  pPatchBonds = ResidueTopologyGetBonds(&patchResiTopo);
  for(i=0;i<BondSetGetCount(pPatchBonds);i++)
  {
    Bond* pCurBond = BondSetGet(pPatchBonds, i);
    BondSetAdd(&pThis->bonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSetCreate(&newBonds);
  for(i=0;i<BondSetGetCount(&pThis->bonds);i++)
  {
    Bond* pCurBond = BondSetGet(&pThis->bonds, i);
    char* atomName;
    atomName = BondGetFromName(pCurBond);
    if( atomName[0]!='+' && atomName[0]!='-' && ResidueGetAtomByName(pThis, atomName)==NULL)
    {
      continue;
    }
    atomName = BondGetToName(pCurBond);
    if( atomName[0]!='+' && atomName[0]!='-' && ResidueGetAtomByName(pThis, atomName)==NULL)
    {
      continue;
    }
    BondSetAdd(&newBonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSetCopy(&pThis->bonds, &newBonds);
  BondSetDestroy(&newBonds);

  ResidueTopologyDestroy(&patchResiTopo);
  
  return Success;
}

int ResiduePatchCTER(Residue* pThis, char* patchName, 
         AtomParamsSet* pAtomParam, 
         ResiTopoSet* pTopos)
{
  int    result;
  char resiName[MAX_LENGTH_RESIDUE_NAME];
  int i;

  BondSet* pPatchBonds;
  BondSet newBonds;

  ResidueTopology patchResiTopo;
  ResidueTopologyCreate(&patchResiTopo);

  // delete old atoms
  result  = ResiTopoCollectionGet(pTopos, "CTER", &patchResiTopo);
  if(FAILED(result))
    return result;
  // do not delete atom O;
  /*deleteAtoms = ResidueTopologyGetDeletes(&patchResiTopo);
  for(i=0;i<StringArrayGetCount(deleteAtoms);i++)
  {
    ResidueDeleteAtom(pThis, StringArrayGet(deleteAtoms, i));
  }*/

  // Temporarily reset its name to that of the patch residue, to find topology in the topology collection
  strcpy(resiName, ResidueGetName(pThis));
  ResidueSetName(pThis, patchName);
  result = ResidueAddAtomsFromAtomParams(pThis, pAtomParam);
  ResidueSetName(pThis, resiName);

  if(FAILED(result))
    return result;

  // new patches are stored at the head
  StringArrayInsert(&pThis->patches, 0, patchName);

  // deal with the bonds
  pPatchBonds = ResidueTopologyGetBonds(&patchResiTopo);
  for(i=0;i<BondSetGetCount(pPatchBonds);i++)
  {
    Bond* pCurBond = BondSetGet(pPatchBonds, i);
    BondSetAdd(&pThis->bonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSetCreate(&newBonds);
  for(i=0;i<BondSetGetCount(&pThis->bonds);i++)
  {
    Bond* pCurBond = BondSetGet(&pThis->bonds, i);
    char* atomName;
    atomName = BondGetFromName(pCurBond);
    if( atomName[0]!='+' && atomName[0]!='-' && ResidueGetAtomByName(pThis, atomName)==NULL)
    {
      continue;
    }
    atomName = BondGetToName(pCurBond);
    if( atomName[0]!='+' && atomName[0]!='-' && ResidueGetAtomByName(pThis, atomName)==NULL)
    {
      continue;
    }
    BondSetAdd(&newBonds, BondGetFromName(pCurBond), BondGetToName(pCurBond), BondGetType(pCurBond));
  }
  BondSetCopy(&pThis->bonds, &newBonds);
  BondSetDestroy(&newBonds);

  ResidueTopologyDestroy(&patchResiTopo);

  return Success;
}

int ResiduePatchNTERorCTER(Residue* pThis, char* NTERorCTER, 
                   AtomParamsSet* pAtomParam, 
                   ResiTopoSet* pTopos)
{
  int result;
  char deleteBondWithPreviousOrNextResidue = ' ';

  if(strcmp(NTERorCTER, "NTER")==0)
  {
    if(strcmp(ResidueGetName(pThis), "GLY")==0)
    {
      result = ResiduePatch(pThis, "GLYP", pAtomParam, pTopos);
    }
    else if(strcmp(ResidueGetName(pThis), "PRO")==0)
    {
      result = ResiduePatch(pThis, "PROP", pAtomParam, pTopos);
    }
    else
    {
      result = ResiduePatch(pThis, "NTER", pAtomParam, pTopos);
    }
    deleteBondWithPreviousOrNextResidue = '-';
  }
  else if(strcmp(NTERorCTER, "CTER")==0)
  {
    result = ResiduePatch(pThis, "CTER", pAtomParam, pTopos);
    deleteBondWithPreviousOrNextResidue = '+';
  }
  else
  {
    result = ValueError;
  }

  if(FAILED(result))
  {
    return result;
  }
  else
  {
    Bond* pCurBond;
    char* fromAtomName;
    char* toAtomName;
    int i;
    for(i=0;i<BondSetGetCount(&pThis->bonds);i++)
    {
      pCurBond = BondSetGet(&pThis->bonds, i);
      fromAtomName = BondGetFromName(pCurBond);
      toAtomName = BondGetToName(pCurBond);
      if( fromAtomName[0] == deleteBondWithPreviousOrNextResidue || 
        toAtomName[0]   == deleteBondWithPreviousOrNextResidue)
      {
        BondSetRemove(&pThis->bonds, fromAtomName, toAtomName);
        break;
      }
    }
  }
  return Success;
}

StringArray* ResidueGetPatchingHistory(Residue* pThis)
{
  return &pThis->patches;
}


int ResidueCalcAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos, 
                 Residue* pPrevResi, Residue* pNextResi, 
                 char* atomName, XYZ* pDestXYZ)
{
  int i;
  char* namesOfAtomABC[3];
  XYZ xyzsOfAtomABC[3];
  BOOL icFound;

  CharmmIC ic;
  ResidueTopology topo;
  CharmmICCreate(&ic);
  ResidueTopologyCreate(&topo);

  icFound = FALSE;
  // find IC from patch topology orderly;
  for(i=0;i<StringArrayGetCount(&pThis->patches);i++)
  {

    ResiTopoCollectionGet(pResiTopos, StringArrayGet(&pThis->patches, i), &topo);
    
    if(!FAILED(ResidueTopologyFindCharmmIC(&topo, atomName, &ic)))
    {
      icFound = TRUE;
      break;
    }          
  }
  // find IC in the residue topology;
  if(icFound == FALSE)
  {
    ResiTopoCollectionGet(pResiTopos, ResidueGetName(pThis), &topo);
    if(FAILED(ResidueTopologyFindCharmmIC(&topo, atomName, &ic)))
    {
      // if the IC is still not found in the residue topology, return an error;
      int result = DataNotExistError;
      char usrMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(usrMsg, "in file %s function %s() line %d, cannot find charmm IC of "
        "atom '%s' in residue %s %d %s", __FILE__, __FUNCTION__, __LINE__,
        atomName, ResidueGetName(pThis), ResidueGetPosInChain(pThis), ResidueGetChainName(pThis));
      TraceError(usrMsg, result);
      return result;
    }
  }

  namesOfAtomABC[0] = CharmmICGetAtomA(&ic);
  namesOfAtomABC[1] = CharmmICGetAtomB(&ic);
  namesOfAtomABC[2] = CharmmICGetAtomC(&ic);

  // find atom coordinate from the current residue, preceding residue and next residue;
  for(i=0;i<3;i++)
  {
    if( namesOfAtomABC[i][0] == '-' && pPrevResi!=NULL && 
      !FAILED(ResidueGetAtomXYZ(pPrevResi, namesOfAtomABC[i]+1, &xyzsOfAtomABC[i])))
    {
        continue;
    }
    else if( namesOfAtomABC[i][0] == '+' && pNextResi!=NULL && 
      !FAILED(ResidueGetAtomXYZ(pNextResi, namesOfAtomABC[i]+1, &xyzsOfAtomABC[i])))
    {
        continue;
    }
    else if( namesOfAtomABC[i][0] != '+' && namesOfAtomABC[i][0] != '-' &&
      !FAILED(ResidueGetAtomXYZ(pThis, namesOfAtomABC[i], &xyzsOfAtomABC[i])))
    {
      continue;
    }
    else
    {
      ResidueTopologyDestroy(&topo);
      CharmmICDestroy(&ic);
      return DataNotExistError;
    }
  }

  GetFourthAtom(&xyzsOfAtomABC[0], &xyzsOfAtomABC[1], 
                   &xyzsOfAtomABC[2], CharmmICGetICParams(&ic), pDestXYZ);

  ResidueTopologyDestroy(&topo);
  CharmmICDestroy(&ic);
  return Success;
}

int ResidueCalcAllAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos, 
                  Residue* pPrevResi, Residue* pNextResi)
{
    BOOL allAtomsXYZAreValid = FALSE;
  int done = FALSE;
  while(!done)
  {
    int i;
    done = TRUE;  /* If nothing new has been calculated in this loop, It will terminate */
    
    allAtomsXYZAreValid = TRUE;
    for(i=0;i<ResidueGetAtomCount(pThis);i++)
    {
      XYZ newXYZ;
      int    result;
      Atom* pCurAtom = ResidueGetAtom(pThis, i);
            if(pCurAtom->isXyzValid)
      {   /* This atom already got Its XYZ */
        continue;
            }
  
      result = ResidueCalcAtomXYZ(pThis, pResiTopos, pPrevResi, pNextResi, 
        AtomGetName(pCurAtom), &newXYZ);
      if(FAILED(result))
      {
        allAtomsXYZAreValid = FALSE;
        continue;
      }
      else
      {
                pCurAtom->xyz = newXYZ;
                pCurAtom->isXyzValid = TRUE;
        done = FALSE;/* New atom XYZ has been calculated in this 'while' loop, go on and try to find more */
      }
    }
  }

  if(!allAtomsXYZAreValid)
  {
    int result = Warning;
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    sprintf(usrMsg, "in file %s function %s() line %d, "
      "not all atoms' XYZ can be calculated for residue %s in chain %s %d", 
      __FILE__, __FUNCTION__, __LINE__,
      ResidueGetName(pThis), ResidueGetChainName(pThis), ResidueGetPosInChain(pThis));
    TraceError(usrMsg, result);
    return result;
  }
  return Success;
}
