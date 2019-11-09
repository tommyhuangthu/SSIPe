#include "Chain.h"
#include <string.h>

int ChainTypeConvertFromString(char* typeName, Type_Chain* type)
{
  if(strcmp(typeName, "PROTEIN")==0)
  {
    *type = Type_Chain_Protein;
  }
  else if(strcmp(typeName, "SMALLMOL")==0)
  {
    *type = Type_Chain_SmallMol;
  }
  else if(strcmp(typeName, "METALION")==0)
  {
    *type = Type_Chain_MetalIon;
  }
  else if(strcmp(typeName, "NUCLEOTIDE")==0)
  {
    *type = Type_Chain_Nucleotide;
  }
  else if(strcmp(typeName, "WATER")==0)
  {
    *type = Type_Chain_Water;
  }
  else
  {
    return ValueError;
  }
  return Success;
}

Type_Chain ChainTypeIdentifiedFromResidueName(char *resiName)
{
  if(strcmp(resiName, "ALA") == 0
    || strcmp(resiName, "ARG") == 0
    || strcmp(resiName, "ASN") == 0
    || strcmp(resiName, "ASP") == 0
    || strcmp(resiName, "CYS") == 0
    || strcmp(resiName, "GLN") == 0
    || strcmp(resiName, "GLU") == 0
    || strcmp(resiName, "GLY") == 0
    || strcmp(resiName, "HIS") == 0
    || strcmp(resiName, "HSD") == 0
    || strcmp(resiName, "HSE") == 0
    || strcmp(resiName, "HSP") == 0
    || strcmp(resiName, "ILE") == 0
    || strcmp(resiName, "LEU") == 0
    || strcmp(resiName, "LYS") == 0
    || strcmp(resiName, "MET") == 0
    || strcmp(resiName, "PHE") == 0
    || strcmp(resiName, "PRO") == 0
    || strcmp(resiName, "SER") == 0
    || strcmp(resiName, "THR") == 0
    || strcmp(resiName, "TRP") == 0
    || strcmp(resiName, "TYR") == 0
    || strcmp(resiName, "VAL") == 0)
  {
    return Type_Chain_Protein;
  }
  else if(strcmp(resiName, "GUA") == 0
    || strcmp(resiName, "ADE") == 0
    || strcmp(resiName, "CYT") == 0
    || strcmp(resiName, "THY") == 0
    || strcmp(resiName, "URA") == 0)
  {
    return Type_Chain_Nucleotide;
  }
  else if(strcmp(resiName, "HOH") == 0
    || strcmp(resiName, "H2O") == 0
    || strcmp(resiName, "WAT") == 0)
  {
    return Type_Chain_Water;
  }
  else
  {
    return Type_Chain_SmallMol;
  }
}

int ChainCreate(Chain* pThis)
{
  ChainSetName(pThis, "");
  pThis->residueNum = 0;
  pThis->residues = NULL;
  return Success;
}

int ChainDestroy(Chain* pThis)
{
  int i;
  for(i=0;i<pThis->residueNum;i++)
  {
    ResidueDestroy(&pThis->residues[i]);
  }
  free(pThis->residues);
  pThis->residues = NULL;
  pThis->residueNum = 0;
  return Success;
}

int ChainCopy(Chain* pThis, Chain* pOther)
{
  int i;
  ChainDestroy(pThis);
  ChainSetName(pThis, pOther->name);
  ChainSetType(pThis, pOther->type);

  pThis->residueNum = pOther->residueNum;
  pThis->residues = (Residue*)malloc(sizeof(Residue)*pThis->residueNum);
  for(i=0;i<pThis->residueNum;i++)
  {
    int result;
    ResidueCreate(&pThis->residues[i]);
    result = ResidueCopy(&pThis->residues[i], &pOther->residues[i]);
    if(FAILED(result))
      return result;
  }
  return Success;
}

char*  ChainGetName(Chain* pThis)
{
  return pThis->name;
}

int ChainSetName(Chain* pThis, char* newName)
{
  if(strlen(newName)>MAX_LENGTH_CHAIN_NAME)
  {
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = ValueError;
    sprintf(usrMsg, "in file %s function %s() line %d, name is too long", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  strcpy(pThis->name, newName);
  return Success;
}

Type_Chain ChainGetType(Chain* pThis)
{
  return pThis->type;
}

int ChainSetType(Chain* pThis, Type_Chain newType)
{
  pThis->type = newType;
  return Success;
}

int ChainGetResidueCount(Chain* pThis)
{
  return pThis->residueNum;
}

Residue* ChainGetResidue(Chain* pThis, int index)
{
  if(index<0 || index>=ChainGetResidueCount(pThis))
  {
    return NULL;
  }
  else
  {
    return pThis->residues + index;
  }

}
int ChainInsertResidue(Chain* pThis, int index, Residue* pNewResi)
{
  int i, newCount;
  if(index<0 || index>pThis->residueNum)
    return IndexError;

  newCount = pThis->residueNum+1;
  pThis->residues = (Residue*)realloc(pThis->residues, sizeof(Residue)*newCount);

  ResidueCreate(&pThis->residues[newCount-1]);
  for(i=newCount-1;i>index;i--)
  {
    ResidueCopy(&pThis->residues[i], &pThis->residues[i-1]);
    ResidueSetPosInChain(&pThis->residues[i], pThis->residues[i-1].posInChain);
  }

  ResidueCopy(&pThis->residues[index], pNewResi);
  ResidueSetSegName(&pThis->residues[index], pThis->name);
  ResidueSetPosInChain(&pThis->residues[index], pNewResi->posInChain);
  pThis->residueNum = newCount;

  return Success;
}

int ChainRemoveResidue(Chain* pThis, int index)
{
  int i;
  if(index<0 || index>=pThis->residueNum)
    return IndexError;

  for(i=index;i<pThis->residueNum-1;i++)
  {
    ResidueCopy(&pThis->residues[i], &pThis->residues[i+1]);
    ResidueSetPosInChain(&pThis->residues[i], i);
  }
  ResidueDestroy(&pThis->residues[pThis->residueNum-1]);
  (pThis->residueNum)--;
  return Success;
}

int ChainAppendResidue(Chain* pThis, Residue* pNewResi)
{
  return ChainInsertResidue(pThis, pThis->residueNum, pNewResi);
}


int ChainCalcAllAtomXYZ(Chain* pThis, ResiTopoSet* topos)
{
  int i;
  Residue* prevResi, *nextResi, *curResi;

  for(i=0;i<ChainGetResidueCount(pThis);i++)
  {
    int    result;

    prevResi = ChainGetResidue(pThis, i-1);
    nextResi = ChainGetResidue(pThis, i+1);
    curResi = ChainGetResidue(pThis, i);

    result = ResidueCalcAllAtomXYZ(curResi, topos, prevResi, nextResi);
    if(FAILED(result))
    {
      char usrMsg[MAX_LENGTH_ERR_MSG+1];
      sprintf(usrMsg, "in file %s function %s() line %d", 
        __FILE__, __FUNCTION__, __LINE__);
      TraceError(usrMsg, result);
      return result;
    }
  }
  return Success;
}

int ChainShowInPDBFormat(Chain* pThis, int resiIndex, int atomIndex, BOOL showHydrogen, FILE* pFile)
{
  char header[7];
  char chainName[2];
  int i;
  int atomCounter;
  if(pThis->type==Type_Chain_Protein)
  {
    strcpy(header, "ATOM  ");
    strcpy(chainName, " ");
    chainName[0] = pThis->name[strlen(pThis->name)-1];
  }
  else
  {
    strcpy(header, "HETATM");
    strcpy(chainName, " ");
  }
  atomCounter = atomIndex;
  for(i=0;i<pThis->residueNum;i++)
  {
    ResidueShowInPDBFormat(&pThis->residues[i], header, chainName, atomCounter, 
      ResidueGetPosInChain(&pThis->residues[i])+resiIndex, showHydrogen, pFile);
    atomCounter += ResidueGetAtomCount(&pThis->residues[i]);
  }
  return Success;
}

