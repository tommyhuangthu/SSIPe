#include "ErrorHandling.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

BOOL FAILED(int errorCode)
{
  if(errorCode == Success || errorCode == Warning)
  {
    return FALSE;
  }
  return TRUE;
}

void TraceError(char* userMsg, int errorCode)
{
  char errMsg[MAX_LENGTH_ERR_MSG+1];

  if(errorCode == Success)
  {
    return;
  }
  else if (errorCode == Warning)
  {
    strcpy(errMsg, "Warning message,");
    printf("--------------------------------------------\n");
    printf("Warning %d : %s %s.\n", errorCode, errMsg, userMsg);
    return;
  }

  switch(errorCode)
  {
    case IOError:
      strcpy(errMsg, "File cannot be opened,"); break;
    case FormatError:
      strcpy(errMsg, "File format is wrong,"); break;
    case IndexError:
      strcpy(errMsg, "Index not in range,"); break;
    case ValueError:
      strcpy(errMsg, "Invalid value used,"); break;
    case ZeroDivisonError:
      strcpy(errMsg, "Denominator is near zero,"); break;
    case DataNotExistError:
      strcpy(errMsg, "Data not exist in records,"); break;
    case NameError:
      strcpy(errMsg, "Parameter name is wrong,"); break;
    default:
      strcpy(errMsg, "Undefined error type,");
  }

  printf("--------------------------------------------\n");
  printf("Error %d : %s %s.\n", errorCode, errMsg, userMsg);
}
