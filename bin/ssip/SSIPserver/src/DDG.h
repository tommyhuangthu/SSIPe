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

#ifndef DDG_H
#define DDG_H

#include <stdio.h>
#include "InterfaceAlignment.h"
#include "Mutant.h"
#include "Parameter.h"

namespace SIP{
  class DDG{
    public:
    SIP::MutantSet ms;
    SIP::InterfaceAlignment ia_structure;
    SIP::InterfaceAlignment ia_sequence;

    DDG();
    ~DDG();
    int ddg_initialize();
    int find_index_for_single_mutation(SIP::SingleMutant *sm);

    // calculate DDG with oberved and pseudo amino-acid counts
    double calc_aat_pseudo_count_for_single_mutation(char aatype, char chain_name, int aln_index, SIP::Parameter *par);
    double calc_gap_pseudo_count_for_single_mutation(char aatype, char chain_name, int aln_index, SIP::Parameter *par);
    double calc_evo_pseudo_count_for_single_mutation(char aatype, char chain_name, int aln_index, SIP::Parameter *par);
    int calc_ddg_for_single_mutation(SIP::SingleMutant *sm, SIP::Parameter *par);
    int calc_ddg_for_multiple_mutation(SIP::MultipleMutant *mm, SIP::Parameter *par);
    int calc_ddg_for_mutation_set(SIP::Parameter *par);


    int calc_ddg_for_multiple_mutation_new(SIP::MultipleMutant *mm, SIP::Parameter *par);
    int calc_ddg_for_mutation_set_new(SIP::Parameter *par);
    int write_ddg_for_mutation_set(char* filepath);

    // calculate DDG with henikoff weighted amino-acid counts
    double calc_gap_pseudo_count_for_single_mutation_henikoff(char aatype, char chain_name, int alnindex, SIP::Parameter *par);
    double calc_evo_pseudo_count_for_single_mutation_henikoff(char aatype, char chainname, int alnindex, SIP::Parameter *par);
    int calc_ddg_for_single_mutation_henikoff(SIP::SingleMutant *sm, SIP::Parameter *par);
    int calc_ddg_for_multiple_mutation_henikoff(SIP::MultipleMutant *mm, SIP::Parameter *par);
    int calc_ddg_for_mutation_set_henikoff(SIP::Parameter *par);

    // calculate DDG using combined interface alignment
    double calc_gap_pseudo_count_for_single_mutation12(char aatype, char chain_name, int alnindex, SIP::Parameter *par);
    double calc_evo_pseudo_count_for_single_mutation12(char aatype, char chainname, int alnindex, SIP::Parameter *par);
    int calc_ddg_for_single_mutation12(SIP::SingleMutant *sm, SIP::Parameter *par);
    int calc_ddg_for_multiple_mutation12(SIP::MultipleMutant *mm, SIP::Parameter *par);
    int calc_ddg_for_mutation_set12(SIP::Parameter *par);


    // combine all profiles: ialign profile, string interface alignment profile & tmalign profile
    double calc_gap_pseudo_count_for_single_mutation1234(char aatype, char chain_name, int alnindex, SIP::Parameter *par);
    double calc_evo_pseudo_count_for_single_mutation1234(char aatype, char chainname, int alnindex, SIP::Parameter *par);
    int calc_ddg_for_single_mutation1234(SIP::SingleMutant *sm, SIP::Parameter *par);
    int calc_ddg_for_multiple_mutation1234(SIP::MultipleMutant *mm, SIP::Parameter *par);
    int calc_ddg_for_mutation_set1234(SIP::Parameter *par);
  };

}


#endif
