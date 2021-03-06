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

#ifndef RESIDUE_H
#define RESIDUE_H

#include "Atom.h"

typedef enum _Type_Bond
{
  Type_Bond_Single, 
  Type_Bond_double, 
  Type_Bond_Triple, 
  Type_Bond_None
}Type_Bond ;

typedef struct _Bond
{
  Type_Bond type;
  char atomFromName[MAX_LENGTH_ATOM_NAME+1];
  char atomToName[MAX_LENGTH_ATOM_NAME+1];
} Bond;
int         BondCreate(Bond* pThis);
void        BondDestroy(Bond* pThis);
int         BondCopy(Bond* pThis, Bond* pOther);
char* BondGetFromName(Bond* pThis);
int         BondSetFromName(Bond* pThis, char* from);
char* BondGetToName(Bond* pThis);
int         BondSetToName(Bond* pThis, char* to);
Type_Bond   BondGetType(Bond* pThis);
int         BondSetType(Bond* pThis, Type_Bond newType);
int         BondShow(Bond* pThis);

typedef struct _BondSet{
  int count;
  Bond* bonds;
} BondSet;

int         BondSetCreate(BondSet* pThis);
void        BondSetDestroy(BondSet* pThis);
int         BondSetCopy(BondSet* pThis, BondSet* pOther);
int         BondSetGetCount(BondSet* pThis);
int         BondSetAdd(BondSet* pThis, char* atom1, char* atom2, Type_Bond bondType);
int         BondSetRemove(BondSet* pThis, char* atom1, char* atom2);
Type_Bond   BondSetFind(BondSet* pThis, char* atom1, char* atom2);
Bond*       BondSetGet(BondSet* pThis, int index);
int         BondSetShow(BondSet* pThis);

typedef struct _CharmmIC
{
  /* All Members are private */
  char atomNames[4][MAX_LENGTH_ATOM_NAME+1];
  BOOL torsionProperFlag;
  double icParam[5]; 
  /*    Proper dihedral angle (A, B, C, D): five parameters are R(AB), Theta(ABC), Phi(ABCD), Theta(BCD), R(CD) */
  /*    Improper dihedral angle (A, B, *C, D); five parameters are R(AC), Theta(ACB), Phi(ABCD), Theta(BCD), R(CD) */
}CharmmIC;
int           CharmmICCreate(CharmmIC* pThis);
int           CharmmICCreateFromStringArray(CharmmIC* pThis, StringArray* params);
void          CharmmICDestroy(CharmmIC* pThis);
int           CharmmICCopy(CharmmIC* pThis, CharmmIC* pOther);
char*   CharmmICGetAtomA(CharmmIC* pThis);
char*   CharmmICGetAtomB(CharmmIC* pThis);
char*   CharmmICGetAtomC(CharmmIC* pThis);
char*   CharmmICGetAtomD(CharmmIC* pThis);
double* CharmmICGetICParams(CharmmIC* pThis);
int           CharmmICSetTorsionAngle(CharmmIC* pThis, double newTorsionAngle);
BOOL          CharmmICGetTorsionProperFlag(CharmmIC* pThis);
int           CharmmICCalcXYZ(CharmmIC* pThis, AtomArray* atomArray, XYZ* pDestXYZ);
int           CharmmICShow(CharmmIC* pThis);



typedef struct _ResidueTopology
{
  char residueName[MAX_LENGTH_RESIDUE_NAME+1];
  
  StringArray atoms;

  StringArray deletes;

  BondSet bonds;

  int icCount;
  CharmmIC* ics;

} ResidueTopology;

int          ResidueTopologyCreate(ResidueTopology* pThis);
int          ResidueTopologyCreateFromFileReader(ResidueTopology* pThis, FileReader* pFileReader);
void         ResidueTopologyDestroy(ResidueTopology* pThis);
int          ResidueTopologyCopy(ResidueTopology* pThis, ResidueTopology* pOther);

char*  ResidueTopologyGetName(ResidueTopology* pThis);
int          ResidueTopologySetName(ResidueTopology* pThis, char* newName);

StringArray* ResidueTopologyGetAtoms(ResidueTopology* pThis);
StringArray* ResidueTopologyGetDeletes(ResidueTopology* pThis);
BondSet*     ResidueTopologyGetBonds(ResidueTopology* pThis);

int          ResidueTopologyGetCharmmICCount(ResidueTopology* pThis);
int          ResidueTopologyGetCharmmIC(ResidueTopology* pThis, int index, CharmmIC* pDestIC);
int          ResidueTopologyFindCharmmIC(ResidueTopology* pThis, char* atomDName, CharmmIC* pDestIC);
int          ResidueTopologyAddCharmmIC(ResidueTopology* pThis, CharmmIC* pNewIC);
int          ResidueTopologyCalcAtomArrayXYZ(ResidueTopology* pThis, AtomArray* pAtomArray);

int          ResidueTopologyShow(ResidueTopology* pThis);



typedef struct _ResiTopoCollection{
  int count;
  ResidueTopology* topos;
} ResiTopoSet;

int     ResiTopoCollectionCreate(ResiTopoSet* pThis);
void    ResiTopoCollectionDestroy(ResiTopoSet* pThis);
int     ResiTopoCollectionCopy(ResiTopoSet* pThis, ResiTopoSet* pOther);
int     ResiTopoCollectionGet(ResiTopoSet* pThis, char* resiName, ResidueTopology* pDestTopo);
int     ResiTopoCollectionAdd(ResiTopoSet* pThis, ResidueTopology* pNewTopo);
int     ResiTopoCollectionAddFromFile(ResiTopoSet* pThis, char* filepath);
int     ResiTopoCollectionShow(ResiTopoSet* pThis);


/* enum types associated with Residue */
typedef enum _Type_ResiduePolarity
{
  Type_ResiduePolarity_Charged, 
  Type_ResiduePolarity_Polar, 
  Type_ResiduePolarity_NonPolar
}Type_ResiduePolarity;

typedef enum _Type_ResiduePosition
{
  Type_ResiduePosition_Buried, 
  Type_ResiduePosition_Intermediate, 
  Type_ResiduePosition_Exposed
}Type_ResiduePosition;


typedef enum _Type_ResidueDesign
{
  Type_ResidueDesign_Catalytic, 
  Type_ResidueDesign_PrimaryContact, 
  Type_ResidueDesign_SecondaryContact
}Type_ResidueDesign;


typedef enum _Type_ResidueChainType
{
  Type_ResidueChainType_Protein, 
  Type_ResidueChainType_Ligand, 
  Type_ResidueChainType_DNA, 
  Type_ResidueChainType_RNA, 
  Type_ResidueChainType_Water, 
  Type_ResidueChainType_Metal, 
} Type_ResidueChainType;

typedef enum _Type_ResidueIsTerminal
{
  Type_ResidueIsNter,
  Type_ResidueIsCter,
  Type_ResidueIsNotTerminal
} Type_ResidueIsTerminal;

typedef struct _Residue
{
  char name[MAX_LENGTH_RESIDUE_NAME+1];
  char chainName[MAX_LENGTH_CHAIN_NAME+1];
  int posInChain;
  AtomArray atoms;
  BondSet bonds;
  // residue patches; new patches are stored at the head, and old patches at the tail;
  StringArray patches;
  BOOL designSiteFlag;
  Atom pseudoAtom;
  int bbAtomCount;

  int numCBwithin8AtoCurResidue;

  Type_ResidueIsTerminal residueTerminalFlag;
  Type_ResidueDesign DesignType;
  Type_ResidueChainType resiChainType;

} Residue;

int ResidueCreate(Residue* pThis);
void ResidueDestroy(Residue* pThis);
int ResidueCopy(Residue* pThis, Residue* pOther);

char* ResidueGetName(Residue* pThis);
int ResidueSetName(Residue* pThis, char* newName);

char* ResidueGetChainName(Residue* pThis);
int ResidueSetSegName(Residue* pThis, char* newSegName);

int ResidueGetPosInChain(Residue* pThis);
int ResidueSetPosInChain(Residue* pThis, int newSegPos);

int ResidueSetDesignSiteFlag(Residue* pThis, BOOL newFlag);

double ResidueGetCharge(Residue* pThis);
int ResidueGetPolarity(Residue* pThis, Type_ResiduePolarity* pPolarity);

int ResidueGetAtomCount(Residue* pThis);
Atom* ResidueGetAtom(Residue* pThis, int index);
Atom* ResidueGetAtomByName(Residue* pThis, char* atomName);
int ResidueFindAtom(Residue* pThis, char* atomName, int* pIndex);
int ResidueGetAtomXYZ(Residue* pThis, char* atomName, XYZ* pXYZ);
AtomArray* ResidueGetAllAtoms(Residue* pThis);
int ResidueInsertAtom(Residue* pThis, int newIndex, Atom* pNewAtom);
int ResidueAddAtom(Residue* pThis, Atom* pNewAtom);
int ResidueDeleteAtom(Residue* pThis, char* atomName);
int ResidueReadXYZFromPDB(Residue* pThis, FileReader* pPDBFileReader,
                  AtomParamsSet* pAtomParam, 
                  ResiTopoSet* pTopos);
int ResidueAddAtomsFromAtomParams(Residue* pThis, AtomParamsSet* pAtomParams);


BondSet* ResidueGetBonds(Residue* pThis);
int ResidueAddBondsFromResiTopos(Residue* pThis, ResiTopoSet* pResiTopoCollection);


int ResidueShowInPDBFormat(Residue* pThis, char* header, char* chainName, 
                                   int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile);

int ResiduePatch(Residue* pThis, char* patchName, 
             AtomParamsSet* pAtomParam, 
             ResiTopoSet* pTopos);
int ResiduePatchCTER(Residue* pThis, char* patchName, 
           AtomParamsSet* pAtomParam, 
           ResiTopoSet* pTopos);
int ResiduePatchNTERorCTER(Residue* pThis, char* NTERorCTER, 
             AtomParamsSet* pAtomParam, 
             ResiTopoSet* pTopos);
StringArray* ResidueGetPatchingHistory(Residue* pThis);

int ResidueCalcAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos, 
                 Residue* pPrevResi, Residue* pNextResi, 
                 char* atomName, XYZ* pDestXYZ);
int ResidueCalcAllAtomXYZ(Residue* pThis, ResiTopoSet* pResiTopos, 
                  Residue* pPrevResi, Residue* pNextResi);

#endif //RESIDUE_H
