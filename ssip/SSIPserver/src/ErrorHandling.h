#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#define _CRT_SECURE_NO_WARNINGS

#define MAX_LENGTH_ERR_MSG 1024
#define BOOL               char
#define TRUE               1
#define FALSE              0

typedef enum _ErrorCode
{
  Success           =   0, 
  Warning           =   1, 
  NameError         =  -1, 
  IOError           =  -2, 
  FormatError       =  -3, 
  IndexError        =  -4, 
  ValueError        =  -5, 
  ZeroDivisonError  =  -6, 
  DataNotExistError =  -7, 
  Exception         = -10, 
} ErrorCode;

BOOL FAILED(int errorCode);
void TraceError(char* userMsg, int errorCode);

#endif // ERROR_HANDLING_H
