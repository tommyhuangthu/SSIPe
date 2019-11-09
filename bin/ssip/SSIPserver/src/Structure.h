#ifndef STRUCTURE_H
#define STRUCTURE_H

#include "Chain.h"


typedef struct _Structure
{
  char name[MAX_LENGTH_STRUCTURE_NAME+1];
  int chainNum;
  Chain* chains;
} Structure;

int StructureCreate(Structure* pThis);
int StructureDestroy(Structure* pThis);
int StructureCopy(Structure* pThis, Structure* pOther);
char* StructureGetName(Structure* pThis);
int StructureSetName(Structure* pThis, char* newName);
int StructureGetChainCount(Structure* pThis);
Chain* StructureGetChain(Structure* pThis, int index);
Chain* StructureGetChainByName(Structure* pThis, char* segName);
int StructureFindChain(Structure* pThis, char* segName, int* index);
int StructureFindSmallMol(Structure* pThis, Residue** ppSmallMol);
int StructureAddChain(Structure* pThis, Chain* newSeg);
int StructureDeleteChain(Structure* pThis, char* segName);
int StructureShowInPDBFormat(Structure* pThis, BOOL showHydrogen, FILE* pFile);
int StructureGetDesignSiteCount(Structure* pThis);


#endif // STRUCTURE_H
