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

#ifndef MUTANT_H
#define MUTANT_H

#include <stdio.h>
#include "Utility.h"

namespace SIP{
  class SingleMutant{
    public:
    char native_aa;
    char chain_name;
    int  pos;
    char mutant_aa;
    int index_in_alignment;
    double smddg;
   
    SingleMutant();
    ~SingleMutant();
    int single_mutant_initialize();
    int single_mutant_create(char *mutStr);
    int single_mutant_destroy();
    int single_mutant_print();
  };

  class MultipleMutant{
    public:
    int nmut;
    SingleMutant *mutants;
    double mmddg;
   
    MultipleMutant();
    ~MultipleMutant();
    int multiple_mutant_initialize();
    int multiple_mutant_create(StringArray *mutStrs);
    int multiple_mutant_destroy();
    int multiple_mutant_print();
  };

  class MutantSet{
    public:
    int nset;
    MultipleMutant *mutset;
   
    MutantSet();
    ~MutantSet();
    int mutant_set_initialize();
    int mutant_set_create(char* mutantfile);
    int mutant_set_append(char *multipleMutStr);
    int mutant_set_destroy();
    int mutant_set_get_count();
    int mutant_set_print();
  };

}






#endif
