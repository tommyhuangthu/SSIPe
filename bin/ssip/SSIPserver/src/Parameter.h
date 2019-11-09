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