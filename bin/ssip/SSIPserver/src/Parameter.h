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

#ifndef PARAMETER_H
#define PARAMETER_H

#include "PseudoCount.h"
#include "Utility.h"

namespace SIP{
  class Parameter{
  public:
    //parameters for pseudo count and lambda
    int fixed_count;
    SIP::PseudoCount pc;
    double lamda0;
    //parameters used to control the profile
    double ial_weight;
    double isscore_cutoff;
    int ial_count_cutoff;
    //psiblast
    double psi_weight;
    double psi_linkscore_cutoff;
    double seqid_high_cutoff;
    double seqid_low_cutoff;
    int psi_count_cutoff;

    Parameter();
    ~Parameter();
    int parameter_initialize();
    int parameter_initialize_random();
    int parameter_copy(Parameter *newpar);
    int parameter_randomize_pseudo_count(int nmut);
    int parameter_randomize_lamda(int nmut);
    int parameter_print();
    int parameter_add(Parameter *newpar);
    int parameter_average(int count);
    int parameter_read(char* parameterfile);
    int parameter_randomize_4par_model(int nmut);
    int parameter_randomize_21par_model_aat(int nmut);
    int parameter_randomize_41par_model_gap(int nmut);
    int parameter_randomize_41par_model_evo(int nmut);
    int parameter_randomize_61par_model(int nmut);
  };
}



#endif // end parameter.h