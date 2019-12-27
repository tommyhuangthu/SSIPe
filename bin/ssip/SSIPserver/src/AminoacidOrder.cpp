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