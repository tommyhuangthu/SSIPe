#include "AminoacidOrder.h"
#include "ErrorHandling.h"
#include <stdio.h>
#include <stdlib.h>

SIP::AminoacidOrder::AminoacidOrder(){
  ;
}

SIP::AminoacidOrder::~AminoacidOrder(){
  ;
}

int SIP::AminoacidOrder::get_aminoacid_order(char letter){
  switch(letter){
    case 'A': aaOrder = 0; break;
    case 'C': aaOrder = 1; break;
    case 'D': aaOrder = 2; break;
    case 'E': aaOrder = 3; break;
    case 'F': aaOrder = 4; break;
    case 'G': aaOrder = 5; break;
    case 'H': aaOrder = 6; break;
    case 'I': aaOrder = 7; break;
    case 'K': aaOrder = 8; break;
    case 'L': aaOrder = 9; break;
    case 'M': aaOrder = 10; break;
    case 'N': aaOrder = 11; break;
    case 'P': aaOrder = 12; break;
    case 'Q': aaOrder = 13; break;
    case 'R': aaOrder = 14; break;
    case 'S': aaOrder = 15; break;
    case 'T': aaOrder = 16; break;
    case 'V': aaOrder = 17; break;
    case 'W': aaOrder = 18; break;
    case 'Y': aaOrder = 19; break;
    case '-': aaOrder = 20; break;
    default:
      printf("Amino-acid type %c is not identified, please check!\n", letter);
      exit(FormatError);
  }
  return Success;
}