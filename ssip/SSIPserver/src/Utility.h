#ifndef UTILITY_H
#define UTILITY_H

#include <stdio.h>
#include <math.h>
#include <time.h>

typedef struct _StringArray
{
    int stringCount;
    int capacity;
    char** strings;
} StringArray;

int StringArrayCreate(StringArray* pThis);
void StringArrayDestroy(StringArray* pThis);
int StringArrayCopy(StringArray* pThis, StringArray* pOther);
int StringArrayGetCount(StringArray* pThis);
int StringArraySet(StringArray* pThis, int index, char* srcString);
char* StringArrayGet(StringArray* pThis, int index);
int StringArrayFind(StringArray* pThis, char* srcString, int* pos);
int StringArrayInsert(StringArray* pThis, int pos, char* srcString);
int StringArrayRemove(StringArray* pThis, int pos);
int StringArrayAppend(StringArray* pThis, char* srcString);
int StringArraySplitString(StringArray* pThis, char* srcStr, char splitter);
int StringArrayShow(StringArray* pThis);
int StringArrayTester();


#define MAX_LENGTH_ONE_LINE_IN_FILE 1024

typedef enum _Type_CoordinateFile
{
  Type_CoordinateFile_PDB, 
  Type_CoordinateFile_MOL2, 
  Type_CoordinateFile_Unrecognized
} Type_CoordinateFile ;

typedef struct _FileReader
{
  StringArray lines;
  int position;
} FileReader;

int FileReaderCreate(FileReader* pThis, char* path);
void FileReaderDestroy(FileReader* pThis);
int FileReaderGetLineCount(FileReader* pThis);
int FileReaderGetLine(FileReader* pThis, int index, char* dest);
int FileReaderGetCurrentPos(FileReader* pThis);
int FileReaderSetCurrentPos(FileReader* pThis, int index);
int FileReaderGetNextLine(FileReader* pThis, char* dest);
int FileReaderEndOfFile(FileReader* pThis);
Type_CoordinateFile FileReaderRecognizeCoordinateFileType(FileReader* pThis);
int FileReaderTester(char* path);

typedef struct _IntArray
{
  int length;
  int capacity;
  int* content;
} IntArray;

int IntArrayCreate(IntArray* pThis, int length);
void IntArrayDestroy(IntArray* pThis);
int IntArrayCopy(IntArray* pThis, IntArray* pOther);
int IntArrayResize(IntArray* pThis, int newLength);
int IntArrayGetLength(IntArray* pThis);
int IntArrayGet(IntArray* pThis, int index);
int IntArraySet(IntArray* pThis, int index, int newValue);
int* IntArrayGetAll(IntArray* pThis);
int IntArraySetAll(IntArray* pThis, int* pNewContent);
int IntArrayInsert(IntArray* pThis, int index, int newValue);
int IntArrayRemove(IntArray* pThis, int index);
int IntArrayAppend(IntArray* pThis, int newValue);
int IntArrayShow(IntArray* pThis);
int IntArrayTester();
int IntArrayFind(IntArray *pThis, int num);

typedef struct _doubleArray
{
  int length;
  int capacity;
  double* content;
} DoubleArray;

int DoubleArrayCreate(DoubleArray* pThis, int length);
void DoubleArrayDestroy(DoubleArray* pThis);
int DoubleArrayCopy(DoubleArray* pThis, DoubleArray* pOther);
int DoubleArrayResize(DoubleArray* pThis, int newLength);
int DoubleArrayGetLength(DoubleArray* pThis);
double DoubleArrayGet(DoubleArray* pThis, int index);
int DoubleArraySet(DoubleArray* pThis, int index, double newValue);
double* DoubleArrayGetAll(DoubleArray* pThis);
int DoubleArraySetAll(DoubleArray* pThis, double* pNewContent);
int DoubleArrayInsert(DoubleArray* pThis, int index, double newValue);
int DoubleArrayRemove(DoubleArray* pThis, int index);
int DoubleArrayAppend(DoubleArray* pThis, double newValue);
int DoubleArrayShow(DoubleArray* pThis);
int DoubleArrayTester();

double DoubleArrayInnerProduct(DoubleArray *pThis, DoubleArray *pOther);
int DoubleArrayScale(DoubleArray *pThis, double scale);
double DoubleArrayNorm(DoubleArray *pThis);
int  DoubleArrayMinus(DoubleArray *pThis, DoubleArray *pOther);

int ExtractTargetStringFromSourceString(char* dest, char* src, int start, int length);
int ExtractFirstStringFromSourceString(char* dest, char* src);
int ExtractFirstStringFromSourceStringNew(char* dest, char* src, char ch);

void Model(int i, FILE* pFile);
void EndModel(FILE* pFile);

void ShowProgress(int width, double percentage);
void SpentTimeShow(time_t ts, time_t te);

#endif // UTILITY_H
