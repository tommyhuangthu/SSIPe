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

#ifndef PSEUDO_COUNT_H
#define PSEUDO_COUNT_H

#define AAT_COUNT_ALA    1
#define AAT_COUNT_CYS    1
#define AAT_COUNT_ASP    1
#define AAT_COUNT_GLU    1
#define AAT_COUNT_PHE    1
#define AAT_COUNT_GLY    1
#define AAT_COUNT_HIS    1
#define AAT_COUNT_ILE    1
#define AAT_COUNT_LYS    1
#define AAT_COUNT_LEU    1
#define AAT_COUNT_MET    1
#define AAT_COUNT_ASN    1
#define AAT_COUNT_PRO    1
#define AAT_COUNT_GLN    1
#define AAT_COUNT_ARG    1
#define AAT_COUNT_SER    1
#define AAT_COUNT_THR    1
#define AAT_COUNT_VAL    1
#define AAT_COUNT_TRP    1
#define AAT_COUNT_TYR    1

#define GAP_COUNT_ALA    1
#define GAP_COUNT_CYS    1
#define GAP_COUNT_ASP    1
#define GAP_COUNT_GLU    1
#define GAP_COUNT_PHE    1
#define GAP_COUNT_GLY    1
#define GAP_COUNT_HIS    1
#define GAP_COUNT_ILE    1
#define GAP_COUNT_LYS    1
#define GAP_COUNT_LEU    1
#define GAP_COUNT_MET    1
#define GAP_COUNT_ASN    1
#define GAP_COUNT_PRO    1
#define GAP_COUNT_GLN    1
#define GAP_COUNT_ARG    1
#define GAP_COUNT_SER    1
#define GAP_COUNT_THR    1
#define GAP_COUNT_VAL    1
#define GAP_COUNT_TRP    1
#define GAP_COUNT_TYR    1

#define EVO_COUNT_ALA    1
#define EVO_COUNT_CYS    1
#define EVO_COUNT_ASP    1
#define EVO_COUNT_GLU    1
#define EVO_COUNT_PHE    1
#define EVO_COUNT_GLY    1
#define EVO_COUNT_HIS    1
#define EVO_COUNT_ILE    1
#define EVO_COUNT_LYS    1
#define EVO_COUNT_LEU    1
#define EVO_COUNT_MET    1
#define EVO_COUNT_ASN    1
#define EVO_COUNT_PRO    1
#define EVO_COUNT_GLN    1
#define EVO_COUNT_ARG    1
#define EVO_COUNT_SER    1
#define EVO_COUNT_THR    1
#define EVO_COUNT_VAL    1
#define EVO_COUNT_TRP    1
#define EVO_COUNT_TYR    1


namespace SIP{
  //const double mat[20][20] = {
  //  {0.297, 0.061, 0.033, 0.039, 0.036, 0.115, 0.044, 0.039, 0.040, 0.044, 0.057, 0.041, 0.065, 0.043, 0.036, 0.113, 0.066, 0.059, 0.026, 0.029},
  //  {0.013, 0.510, 0.003, 0.003, 0.004, 0.004, 0.004, 0.008, 0.004, 0.007, 0.009, 0.009, 0.017, 0.004, 0.004, 0.014, 0.014, 0.010, 0.004, 0.004},
  //  {0.017, 0.007, 0.419, 0.057, 0.008, 0.010, 0.025, 0.012, 0.029, 0.012, 0.015, 0.042, 0.016, 0.032, 0.025, 0.028, 0.024, 0.013, 0.006, 0.010},
  //  {0.047, 0.016, 0.128, 0.318, 0.024, 0.014, 0.046, 0.046, 0.128, 0.042, 0.052, 0.058, 0.025, 0.145, 0.094, 0.060, 0.065, 0.044, 0.018, 0.030},
  //  {0.018, 0.009, 0.008, 0.010, 0.392, 0.009, 0.032, 0.025, 0.012, 0.024, 0.040, 0.013, 0.014, 0.015, 0.011, 0.011, 0.014, 0.020, 0.107, 0.109},
  //  {0.083, 0.013, 0.014, 0.009, 0.014, 0.654, 0.023, 0.008, 0.010, 0.008, 0.013, 0.023, 0.022, 0.014, 0.012, 0.037, 0.015, 0.011, 0.015, 0.012},
  //  {0.013, 0.006, 0.014, 0.011, 0.019, 0.009, 0.414, 0.009, 0.011, 0.008, 0.009, 0.020, 0.011, 0.015, 0.016, 0.011, 0.012, 0.008, 0.014, 0.025},
  //  {0.038, 0.035, 0.022, 0.038, 0.050, 0.011, 0.029, 0.213, 0.054, 0.083, 0.094, 0.072, 0.017, 0.046, 0.037, 0.031, 0.050, 0.130, 0.025, 0.032},
  //  {0.037, 0.019, 0.051, 0.099, 0.022, 0.012, 0.034, 0.050, 0.207, 0.047, 0.043, 0.056, 0.020, 0.095, 0.124, 0.043, 0.046, 0.045, 0.015, 0.022},
  //  {0.115, 0.084, 0.058, 0.093, 0.125, 0.027, 0.071, 0.220, 0.135, 0.447, 0.236, 0.133, 0.028, 0.116, 0.088, 0.084, 0.118, 0.180, 0.077, 0.111},
  //  {0.020, 0.014, 0.010, 0.015, 0.028, 0.006, 0.011, 0.033, 0.016, 0.031, 0.161, 0.015, 0.005, 0.022, 0.015, 0.013, 0.018, 0.022, 0.012, 0.018},
  //  {0.028, 0.029, 0.054, 0.033, 0.018, 0.021, 0.047, 0.050, 0.041, 0.035, 0.030, 0.247, 0.014, 0.036, 0.035, 0.038, 0.043, 0.054, 0.010, 0.026},
  //  {0.017, 0.022, 0.008, 0.006, 0.008, 0.008, 0.010, 0.005, 0.006, 0.003, 0.004, 0.005, 0.597, 0.005, 0.008, 0.012, 0.014, 0.008, 0.007, 0.006},
  //  {0.030, 0.014, 0.042, 0.086, 0.021, 0.014, 0.036, 0.033, 0.073, 0.031, 0.046, 0.037, 0.014, 0.221, 0.061, 0.032, 0.038, 0.028, 0.015, 0.019},
  //  {0.026, 0.012, 0.035, 0.058, 0.017, 0.012, 0.041, 0.028, 0.100, 0.025, 0.033, 0.038, 0.023, 0.064, 0.318, 0.031, 0.037, 0.025, 0.016, 0.022},
  //  {0.072, 0.041, 0.034, 0.032, 0.014, 0.033, 0.025, 0.021, 0.030, 0.021, 0.024, 0.036, 0.029, 0.029, 0.027, 0.303, 0.065, 0.026, 0.011, 0.026},
  //  {0.046, 0.047, 0.032, 0.039, 0.019, 0.014, 0.029, 0.036, 0.035, 0.032, 0.037, 0.045, 0.037, 0.037, 0.035, 0.071, 0.274, 0.044, 0.015, 0.021},
  //  {0.065, 0.052, 0.026, 0.041, 0.044, 0.016, 0.032, 0.146, 0.054, 0.076, 0.070, 0.088, 0.032, 0.044, 0.037, 0.045, 0.069, 0.257, 0.016, 0.034},
  //  {0.005, 0.003, 0.002, 0.003, 0.040, 0.004, 0.009, 0.005, 0.003, 0.005, 0.007, 0.003, 0.005, 0.004, 0.004, 0.003, 0.004, 0.003, 0.535, 0.023},
  //  {0.013, 0.009, 0.009, 0.011, 0.097, 0.007, 0.038, 0.014, 0.010, 0.019, 0.023, 0.017, 0.009, 0.012, 0.013, 0.018, 0.013, 0.014, 0.054, 0.420}
  //};

  class PseudoCount{
    public:
    int aat_pseudo_count[20]; // amino acid type dependent pseudo counts, ACDEFGHIKLMNPQRSTVWY
    int gap_pseudo_count[20];
    int evo_pseudo_count[20];
    static double iPTM[20][20];

    PseudoCount();
    ~PseudoCount();
    int set_aat_pseudo_count();
    int set_aat_pseudo_count(int index, int count);
    int set_gap_pseudo_count();
    int set_gap_pseudo_count(int index, int count);
    int set_evo_pseudo_count();
    int set_evo_pseudo_count(int index, int count);
    int set_random_pseudo_count();

    int get_aat_pseudo_count(int aaindex);
    int get_gap_pseudo_count(int aaindex);
    int get_evo_pseudo_count(int aaindex);
    int read_interface_probability_transition_matrix();

    int pseudo_counts_copy(PseudoCount *newdata);
  };

}



#endif
