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

#ifndef ATOM_H
#define ATOM_H

#include "ErrorHandling.h"
#include "Utility.h"
#include "GeometryCalc.h"

#define MAX_LENGTH_ATOM_NAME          5
#define MAX_LENGTH_ATOM_TYPE          5
#define MAX_LENGTH_ATOM_HYBRIDTYPE    5
#define MAX_LENGTH_ATOM_DONOR         5
#define MAX_LENGTH_ATOM_ACCEPTOR      5
#define MAX_LENGTH_RESIDUE_NAME       5
#define MAX_LENGTH_CHAIN_NAME         5
#define MAX_LENGTH_STRUCTURE_NAME    10

typedef enum _Type_AtomPolarity
{
  Type_AtomPolarity_P, 
  Type_AtomPolarity_C, 
  Type_AtomPolarity_NPAliphatic, 
  Type_AtomPolarity_NPAromatic
} Type_AtomPolarity;

typedef enum _Type_AtomHydrogen
{
  Type_AtomHydrogen_PolarH, 
  Type_AtomHydrogen_NPolarH, 
  Type_AtomHydrogen_Heavy,
  Type_AtomHydrogen_United
} Type_AtomHydrogen;

typedef enum _Type_Backbone_AtomType
{
  Type_BB_Atom_HN,
  Type_BB_Atom_HA,
  Type_BB_Atom_N,
  Type_BB_Atom_CA,
  Type_BB_Atom_C,
  Type_BB_Atom_O,
  Type_NotBBAtom,
} Type_Backbone_AtomType;

typedef enum _Type_AtomHybridType
{
  Type_AtomHybridType_SP3,
  Type_AtomHybridType_SP2,
  Type_AtomHybridType_SP,
  Type_AtomHybridType_None,
} Type_AtomHybridType;


typedef struct _Atom
{
  // attributes from parameter file;
  char name[MAX_LENGTH_ATOM_NAME+1];
  char type[MAX_LENGTH_ATOM_TYPE+1];
  char hbHorA[MAX_LENGTH_ATOM_DONOR+1];
  char hbDorB[MAX_LENGTH_ATOM_ACCEPTOR+1];
  char hbB2[MAX_LENGTH_ATOM_ACCEPTOR+1];
  Type_AtomPolarity polarity;
  Type_AtomHybridType hybridType;
  Type_Backbone_AtomType bbAtomType;

  // attributes from PDB file;
  char chainName[MAX_LENGTH_CHAIN_NAME+1];
  int  posInChain;
  XYZ  xyz;
  
  BOOL isXyzValid;
  BOOL isBBAtom;
  BOOL isInHBond;
  BOOL isHBatomH;
  BOOL isHBatomA;

  // CHARMM parameters
  double CHARMM_epsilon;
  double CHARMM_radius;
  double CHARMM_charge;

  // parameters used in LK model
  int    EEF1_atType;
  double EEF1_volume;
  double EEF1_lamda_;
  double EEF1_radius;
  double EEF1_refDG_;
  double EEF1_freeDG;

  // parameters used in FoldEF model
  double FOLDEF_charge;
  double FOLDEF_volume;
  double FOLDEF_radius;
  double FOLDEF_OccCal;
  double FOLDEF_Occmin;
  double FOLDEF_Occmax;
  double FOLDEF_VDWene;
  double FOLDEF_SolEne;
} Atom;

int AtomCreate(Atom* pThis);
void AtomDestroy(Atom* pThis);
int AtomCopy(Atom* pThis, Atom* pOther);
char* AtomGetName(Atom* pThis);
int AtomSetName(Atom* pThis, char* newName);
char* AtomGetType(Atom* pThis);
Type_AtomHybridType AtomGetHybridType(Atom* pThis);
char* AtomGetHbHorA(Atom* pThis);
char* AtomGetHbDorB(Atom* pThis);
char* AtomGetHbB2(Atom* pThis);
BOOL AtomIsHydrogen(Atom* pThis);
int AtomSetParamsByStringArray(Atom* pThis, StringArray* pParams);
int AtomShowParams(Atom* pThis);

char* AtomGetChainName(Atom* pThis);
int AtomSetChainName(Atom* pThis, char* newChainName);
int AtomGetPosInChain(Atom* pThis);
int AtomSetPosInChain(Atom* pThis, int newChainPos);
int AtomShowInPDBFormat(Atom* pThis, char* header, char* resiName, char* chainName, int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile);


typedef struct _AtomArray
{
    Atom* atoms;
  int atomNum;
} AtomArray;

int AtomArrayCreate(AtomArray* pThis);
void AtomArrayDestroy(AtomArray* pThis);
int AtomArrayCopy(AtomArray* pThis, AtomArray* pOther);
int AtomArrayGetCount(AtomArray* pThis);
Atom* AtomArrayGet(AtomArray* pThis, int index);
Atom* AtomArrayGetByName(AtomArray* pThis, char* atomName);
int AtomArrayFind(AtomArray* pThis, char* atomName, int* pIndex);
int AtomArrayInsert(AtomArray* pThis, int index, Atom* pNewAtom);
int AtomArrayRemove(AtomArray* pThis, int index);
int AtomArrayRemoveByName(AtomArray* pThis, char* atomName);
int AtomArrayAppend(AtomArray* pThis, Atom* pNewAtom);
double AtomArrayCalcTotalCharge(AtomArray* pThis);
double AtomArrayCalcMinDistance(AtomArray* pThis, AtomArray* pOther);
BOOL AtomArrayAllAtomXYZAreValid(AtomArray* pThis);
int AtomArrayShowInPDBFormat(AtomArray* pThis, char* header, char* resiName, char* chainName, int atomIndex, int resiIndex, BOOL showHydrogen, FILE* pFile);

typedef struct _AtomParamsSet
{
    StringArray residueNames;
    int*   atomCount;
    Atom** atoms;
} AtomParamsSet;

int AtomParamsSetCreate(AtomParamsSet* pThis);
void AtomParamsSetDestroy(AtomParamsSet* pThis);
int AtomParamsSetAdd(AtomParamsSet* pThis, char* residueName, Atom* pNewAtom);
int AtomParamsSetAddFromFile(AtomParamsSet* pThis, char* filepath);
int AtomParamsSetGetResidueCount(AtomParamsSet* pThis);
int AtomParamsSetGetResidueName(AtomParamsSet* pThis, int index, char* residueName);
int AtomParamsSetGetAtomCount(AtomParamsSet* pThis, char* residueName, int* pCount);
int AtomParamsSetGetAtomParam(AtomParamsSet* pThis, char* residueName, int index, Atom* pDestAtom);
int AtomParamsSetGetAtomParamByName(AtomParamsSet* pThis, char* residueName, char* atomName, Atom* pDestAtom);

int AtomAssignEEF1Parameter(Atom *pThis, int eef1_atType, double refDG, double freeDG, double volume, double lamda);
int AtomparamsSetAssignEEF1Parameters(AtomParamsSet* pThis);

#endif // ATOM_H
