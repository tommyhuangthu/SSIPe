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
