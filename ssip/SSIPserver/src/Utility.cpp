#pragma warning(disable:4996)

#include "Utility.h"
#include "ErrorHandling.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>

#include <errno.h>


int StringArrayCreate(StringArray* pThis){
  pThis->stringCount = 0;
  pThis->capacity = 1;
  pThis->strings = (char**)malloc(sizeof(char*)*(pThis->capacity));
  return Success;
}

void StringArrayDestroy(StringArray* pThis){
  int i;
  for(i=0;i<pThis->stringCount;i++){
    free(pThis->strings[i]);
  }
  free(pThis->strings);
  pThis->stringCount = pThis->capacity = 0;
  pThis->strings = NULL;
}

int StringArrayCopy(StringArray* pThis, StringArray* pOther){
  int i;
  StringArrayDestroy(pThis);
  pThis->capacity = pOther->capacity;
  pThis->stringCount = pOther->stringCount;
  pThis->strings = (char**)malloc(sizeof(char*)*pThis->capacity);
  for(i=0;i<pThis->stringCount;i++){
    pThis->strings[i] = (char*)malloc( sizeof(char) * (strlen(pOther->strings[i])+1) );
    strcpy(pThis->strings[i], pOther->strings[i]);
  }
  return Success;
}

int StringArrayGetCount(StringArray* pThis){
  return pThis->stringCount;
}

int StringArraySet(StringArray* pThis, int index, char* srcString){
  if(index<0 || index>=pThis->stringCount){
    return IndexError;
  }
  if(strlen(srcString) > strlen(pThis->strings[index]) ){
    pThis->strings[index] = (char*)realloc(pThis->strings[index], 
      sizeof(char)*(strlen(srcString)+1));
  }
  strcpy(pThis->strings[index], srcString);
  return Success;
}

char* StringArrayGet(StringArray* pThis, int index){
  if(index<0 || index>=pThis->stringCount){
    return NULL;
  }
  return pThis->strings[index];
}

int StringArrayFind(StringArray* pThis, char* srcString, int* pos){
  int index;
  for(index=0;index<pThis->stringCount;index++){
    if(strcmp(pThis->strings[index], srcString)==0)
      break;
  }
  if(index==pThis->stringCount){
    return DataNotExistError;
  }
  *pos = index;
  return Success;
}

int StringArrayInsert(StringArray* pThis, int pos, char* srcString){
  int i;
  char* newString;
  if(pos<0 || pos>pThis->stringCount){
    return IndexError;
  }
  
  if(pThis->stringCount == pThis->capacity){
    int newCap = pThis->capacity*2;
    pThis->strings = (char**)realloc(pThis->strings, sizeof(char*) * newCap);
    pThis->capacity = newCap;
  }
  
  newString = (char*)malloc(
    sizeof(char)*( strlen(srcString)+1 ));
  strcpy(newString, srcString);

  for(i=pThis->stringCount; i>pos; i--){
    pThis->strings[i] = pThis->strings[i-1];
  }
  pThis->strings[pos] = newString;
  pThis->stringCount++;
  return Success;
}

int StringArrayRemove(StringArray* pThis, int pos){
  int i;
  if(pos<0 || pos>=pThis->stringCount){
    return IndexError;
  }
  free(pThis->strings[pos]);
  for(i=pos;i < pThis->stringCount-1;i++){
    pThis->strings[i] = pThis->strings[i+1];
  }
  pThis->stringCount--;
  return Success;  
}

int StringArrayAppend(StringArray* pThis, char* srcString){
  return StringArrayInsert(pThis, StringArrayGetCount(pThis), srcString);
}

int StringArraySplitString(StringArray* pThis, char* srcStr, char splitter){
  char* buffer;
  int beg, end, length;
  StringArrayDestroy(pThis);
  StringArrayCreate(pThis);
  length = (int)strlen(srcStr);
  buffer = (char*)malloc(sizeof(char)*(length+1));

  beg = 0;
  while(beg<length){
    end = beg;
    while(end<length){
      char endChar = srcStr[end];
      if(  splitter==endChar || (isspace(splitter)&&isspace(endChar))  ){
        break;
      }
      else{
        end++;
      }
    }
    if(end>beg){
      strcpy(buffer, srcStr+beg);
      buffer[end-beg] = '\0';
      StringArrayAppend(pThis, buffer);
    }
    beg = end+1;
  }
  free(buffer);
  return Success;
}


int FileReaderCreate(FileReader* pThis, char* path)
{
  char buffer[MAX_LENGTH_ONE_LINE_IN_FILE+1];
  char usrMsg[MAX_LENGTH_ERR_MSG+1];
  FILE* pFile;

  StringArrayCreate(&pThis->lines);
  pFile = fopen(path, "r");
  if(pFile==NULL){
    int    result = IOError;
    sprintf(usrMsg, "in file %s function %s() line %d, when opening:\n%s", 
      __FILE__, __FUNCTION__, __LINE__, path);
    TraceError(usrMsg, result);
    return result;
  }
  while(fgets(buffer, MAX_LENGTH_ONE_LINE_IN_FILE, pFile)){
    int length, i;
    length = 0;
    while(length< (int)(strlen(buffer)) && buffer[length] != '!')
      length++;
    while(length>0 && buffer[length-1]>0 && isspace(buffer[length-1]))
      length--;
    if(length==0)
      continue;
    buffer[length]='\0';
    for(i=0;i<length;i++){
      if(buffer[i] < 0){
        int    result = FormatError;
        sprintf(usrMsg, "in file %s function %s() line %d, only ASCII characters "
          "are allowed in non-comment Area : \n%s", 
          __FILE__, __FUNCTION__, __LINE__, buffer);
        TraceError(usrMsg, result);
        return result;
      }
    }
    StringArrayAppend(&pThis->lines, buffer);
  }
  pThis->position = 0;

  fclose(pFile);
  return Success;
}

void FileReaderDestroy(FileReader* pThis){
  StringArrayDestroy(&pThis->lines);
}

int FileReaderGetLineCount(FileReader* pThis){
  return StringArrayGetCount(&pThis->lines);
}

int FileReaderGetLine(FileReader* pThis, int index, char* dest){
  char* line = StringArrayGet(&pThis->lines, index);
  if(line==NULL){
    return IndexError;
  }
  else{
    strcpy(dest, line);
    return Success;
  }
}

int FileReaderGetCurrentPos(FileReader* pThis){
  return pThis->position;
}

int FileReaderSetCurrentPos(FileReader* pThis, int index){
  if(index<0 || index>=FileReaderGetLineCount(pThis)){
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = IndexError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  pThis->position = index;
  return Success;
}

int FileReaderGetNextLine(FileReader* pThis, char* dest){
  int    result;
  result = FileReaderGetLine(pThis, pThis->position, dest);
  if(FAILED(result)){
    return result;
  }
  else{
    (pThis->position)++;
    return Success;
  }
}

int FileReaderEndOfFile(FileReader* pThis){
  if(pThis->position == FileReaderGetLineCount(pThis))
    return 1;
  else
    return 0;
}


int IntArrayCreate(IntArray* pThis, int length){
  if(length<0)
    return IndexError;

  pThis->content = (int*)calloc(length, sizeof(int));
  pThis->length = length;
  pThis->capacity = length;
  return Success;
}

void IntArrayDestroy(IntArray* pThis){
  free(pThis->content);
  pThis->content = NULL;
  pThis->length = 0;
  pThis->capacity = 0;
}

int IntArrayCopy(IntArray* pThis, IntArray* pOther){
  IntArrayDestroy(pThis);
  pThis->content = (int*)calloc(pOther->length, sizeof(int));
  pThis->length = pOther->length;
  pThis->capacity = pThis->length;
  IntArraySetAll(pThis, IntArrayGetAll(pOther));
  return Success;
}

int IntArrayResize(IntArray* pThis, int newLength){
  int i;
  if(newLength<0)
    return IndexError;

  pThis->content = (int*)realloc(pThis->content, sizeof(int)*newLength);
  for(i=pThis->length;i<newLength;i++){
    pThis->content[i] = 0;
  }
  pThis->length = newLength;
  pThis->capacity = pThis->length;
  return Success;
}

int IntArrayGetLength(IntArray* pThis){
  return pThis->length;
}

int IntArrayGet(IntArray* pThis, int index){
  if(index<0 || index>=pThis->length){
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = IndexError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  return pThis->content[index];
}

int IntArraySet(IntArray* pThis, int index, int newValue){
  if(index<0 || index>=pThis->length)
    return IndexError;
  pThis->content[index] = newValue;
  return Success;
}

int* IntArrayGetAll(IntArray* pThis){
  return pThis->content;
}

int IntArraySetAll(IntArray* pThis, int* pNewContent){
  int i;
  if(pNewContent==NULL)
    return ValueError;
  for(i=0;i<pThis->length;i++){
    pThis->content[i] = pNewContent[i];
  }
  return Success;
}

int IntArrayInsert(IntArray* pThis, int index, int newValue){
  int i;
  if(index<0 || index>pThis->length)
    return IndexError;
  if(pThis->capacity == pThis->length){
    pThis->capacity = pThis->capacity*2 + 1;
    pThis->content = (int*)realloc(pThis->content, sizeof(int)*pThis->capacity);
  }
  for(i=pThis->length;i>index;i--){
    pThis->content[i] = pThis->content[i-1];
  }
  pThis->content[index] = newValue;

  (pThis->length)++;
  return Success;
}
  
int IntArrayRemove(IntArray* pThis, int index){
  int i;
  if(index<0 || index>pThis->length){
    return IndexError;
  }

  for(i=index;i<pThis->length-1;i++){
    pThis->content[i] = pThis->content[i+1];
  }
  (pThis->length)--;
  return Success;
}

int IntArrayAppend(IntArray* pThis, int newValue){
  return IntArrayInsert(pThis, IntArrayGetLength(pThis), newValue);
}

int IntArrayFind(IntArray *pThis, int num){
  int i;
  for(i = 0;i < pThis->length; i++){
    if(pThis->content[i] == num){
      return i;
    }
  }
  return -1;
}


int DoubleArrayCreate(DoubleArray* pThis, int length){
  int i;
  if(length<0)
    return IndexError;

  pThis->content = (double*)calloc(length, sizeof(double));
  pThis->length = length;
  pThis->capacity = length;
  for(i=0;i<pThis->length;i++){
    pThis->content[i] = 0.0;
  }
  return Success;
}

void DoubleArrayDestroy(DoubleArray* pThis){
  free(pThis->content);
  pThis->content = NULL;
  pThis->length = 0;
  pThis->capacity = 0;
}

int DoubleArrayCopy(DoubleArray* pThis, DoubleArray* pOther){
  if(pThis->capacity >= pOther->length){
    DoubleArraySetAll(pThis, DoubleArrayGetAll(pOther));
    return Success;
  }
  DoubleArrayDestroy(pThis);
  pThis->content = (double*)calloc(pOther->length, sizeof(double));
  pThis->length = pOther->length;
  pThis->capacity = pThis->length;
  DoubleArraySetAll(pThis, DoubleArrayGetAll(pOther));
  return Success;
}

int DoubleArrayResize(DoubleArray* pThis, int newLength){
  int i;
  if(newLength<0)
    return IndexError;

  pThis->content = (double*)realloc(pThis->content, sizeof(double)*newLength);
  for(i=pThis->length;i<newLength;i++){
    pThis->content[i] = 0.0;
  }
  pThis->length = newLength;
  pThis->capacity = pThis->length;
  return Success;
}

int DoubleArrayGetLength(DoubleArray* pThis){
  return pThis->length;
}

double DoubleArrayGet(DoubleArray* pThis, int index){
  if(index<0 || index>=pThis->length){
    char usrMsg[MAX_LENGTH_ERR_MSG+1];
    int errorCode = IndexError;
    sprintf(usrMsg, "in file %s function %s() line %d", 
      __FILE__, __FUNCTION__, __LINE__);
    TraceError(usrMsg, errorCode);
    return errorCode;
  }
  return pThis->content[index];
}

int DoubleArraySet(DoubleArray* pThis, int index, double newValue){
  if(index<0 || index>=pThis->length)
    return IndexError;
  pThis->content[index] = newValue;
  return Success;
}

double* DoubleArrayGetAll(DoubleArray* pThis){
  return pThis->content;
}

int DoubleArraySetAll(DoubleArray* pThis, double* pNewContent){
  int i;
  if(pNewContent==NULL)
    return ValueError;
  for(i=0;i<pThis->length;i++){
    pThis->content[i] = pNewContent[i];
  }
  return Success;
}

int DoubleArrayInsert(DoubleArray* pThis, int index, double newValue){
  int i;
  if(index<0 || index>pThis->length)
    return IndexError;
  if(pThis->capacity == pThis->length){
    pThis->capacity = pThis->capacity*2 + 1;
    pThis->content = (double*)realloc(pThis->content, sizeof(double)*pThis->capacity);
  }
  for(i=pThis->length;i>index;i--){
    pThis->content[i] = pThis->content[i-1];
  }
  pThis->content[index] = newValue;

  (pThis->length)++;
  return Success;
}

int DoubleArrayRemove(DoubleArray* pThis, int index){
  int i;
  if(index<0 || index>pThis->length)
    return IndexError;

  for(i=index;i<pThis->length-1;i++){
    pThis->content[i] = pThis->content[i+1];
  }
  (pThis->length)--;
  return Success;
}

int DoubleArrayAppend(DoubleArray* pThis, double newValue){
  return DoubleArrayInsert(pThis, DoubleArrayGetLength(pThis), newValue);
}

double DoubleArrayInnerProduct(DoubleArray *pThis, DoubleArray *pOther){
  int i = 0;
  double result = 0.0;
  if(pThis->length == pOther->length){
    for(i = 0; i < pThis->length; i++){
      result += pThis->content[i]*pOther->content[i];
    }
    return result;
  }
  printf("In doubleArrayInnerProduct, Dim error\n");
  return 0;
}

int DoubleArrayScale(DoubleArray *pThis, double scale){
  int i = 0;
  for(i = 0; i < pThis->length; i++)
  {
    pThis->content[i] *= scale;
  }
  return Success;
}

double DoubleArrayNorm(DoubleArray *pThis)
{
  int i = 0;
  double result = 0.0;
  for(i = 0; i < pThis->length; i++)
  {
    result += pThis->content[i]*pThis->content[i];
  }
  return sqrt(result);
}

int  DoubleArrayMinus(DoubleArray *pThis, DoubleArray *pOther)
{
  int i = 0;
  if(pThis->length != pOther->length)
  {
    printf("In doubleArrayMinus, length error\n");
    return 1;
  }
  for(i = 0; i < pThis->length; i++)
  {
    pThis->content[i] -= pOther->content[i];
  }
  return Success;
}

int ExtractTargetStringFromSourceString(char* dest, char* src, int start, int length)
{
  // the function strncpy() from <string.h> has the following format;
  // char *strncpy(char *strDest, char *strSource, size_t count);
  // the function copies characters with a maximum number of 'count' from string 'strSource' to string 'strDest', 
  // if 'count' <= strlen(strSource), '\0' will not be added to the end of 'strDest'; else added.
  // given the parameters 'beg' and 'length';
  // the following operation:
  // strncpy(strDest, strSource+beg, length);
  // strDest[length] = '\0';
  // equals to
  // copying a string with a length of 'length' from the 'beg' character of 'strSource' to 'strDest'.

  char* temp = (char*)malloc(sizeof(char)*(length+1));
  strncpy(temp, src+start, length);
  temp[length]='\0';
  // discard redundant space characters;
  sscanf(temp, "%s", dest);
  free(temp);
  return Success;
}

int ExtractFirstStringFromSourceString(char* dest, char* src)
{

    int from, to, i;
    if(sscanf(src, "%s", dest)==EOF)
        return ValueError;
  // the number of space characters before the first string;
    from = 0;
    while(isspace(src[from]))
        from++;
  // add the length of the first string;
    from += (int)strlen(dest);
  // add the number of space characters after the first string;
    while(isspace(src[from]))
        from++;
    to = strlen(src);
  // the left string is still stored in string 'src';
    for(i=from;i<=to;i++)
  {
        src[i-from] = src[i];
    }
    return Success;
}

int ExtractFirstStringFromSourceStringNew(char* dest, char* src, char ch)
{

  int i;
  strcpy(dest, src);
  for(i = 0; i < (int)strlen(dest); i++){
    if(dest[i] == ch){
      dest[i] = '\0';
      break;
    }
  }
  for(int from = i; from <= (int)strlen(src); from++){
    src[from - i] = src[from];
  }
  if(src[0] == ch){
    for(int i = 0; i < (int)strlen(src); i++){
      src[i] = src[i+1];
    }
  }
  return Success;
}

void Model(int i, FILE* pFile)
{
    if(pFile==NULL)
  {
        pFile = stdout;
    }
    fprintf(pFile, "MODEL     %d\n", i);
}

void EndModel(FILE* pFile)
{
    if(pFile==NULL)
  {
        pFile = stdout;
    }
    fprintf(pFile, "ENDMDL\n");
}

void ShowProgress(int width, double percentage)
{
    int i;
    int complete;
    if(width<0)
  {
        width=0;
    }
    complete = (int)(width*percentage/100.0);
    printf("|");
    for(i=0;i<complete;i++)
  {
        printf(">");
    }
    for(i=complete;i<width;i++)
  {
        printf("=");
    }
    printf("|");
    printf("%5.2f%% ", percentage);
}

void SpentTimeShow(time_t ts, time_t te)
{
  printf("Total time spent: %f\n", (double)(te-ts)/CLOCKS_PER_SEC);
}

