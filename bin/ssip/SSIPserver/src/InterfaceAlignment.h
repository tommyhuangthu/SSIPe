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

#ifndef INTERFACE_ALIGNMENT_H
#define INTERFACE_ALIGNMENT_H

#include <stdio.h>
#include "Utility.h"

#define MAX_LENGTH_FILE_NAME     100
#define MAX_NUM_AA_TYPE          21 // ACDEFGHIKLMNPQRSTVWY-
#define MAX_NUM_INTERFACE_RESI   100
#define MAX_NUM_INTERFACE_ALN    500
#define MAX_NUM_EFF_INF_ALN      500


namespace SIP{
  class InterfaceAlignment{
    public:
    char   alignment_filename[100];
    int    lig_alignment_num;
    int    rec_alignment_num;
    char   lig_chn_name[2];
    char   rec_chn_name[2];
    int    lig_res_count;
    int    rec_res_count;
    int    lig_res_index[MAX_NUM_INTERFACE_RESI];
    int    rec_res_index[MAX_NUM_INTERFACE_RESI];
    char   native_lig_aas[MAX_NUM_INTERFACE_RESI];
    char   native_rec_aas[MAX_NUM_INTERFACE_RESI];
    double lig_aligntable[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE];
    double rec_aligntable[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE];
    double eff_ligalign[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE];
    double eff_recalign[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE];
    
    InterfaceAlignment();
    ~InterfaceAlignment();
    int interface_alignment_count_amino_acid_number_in_string(double table[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE], char* string, int length);
    int interface_alignment_count_amino_acid_number_in_string_with_weight(double table[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE], char* string, int length, double seq_weight);
    double interface_alignment_get_gap_count(char chnid, int alnindex);
    double interface_alignment_get_gap_ratio(char chnid, int alnindex);
    int interface_alignment_initialize();
    int interface_alignment_print_statistics();

    // read Ialign interface alignment file
    int interface_alignment_read_bindprofx(char* ialignfile, double seq_weight, double cutoff_isscore);
    int interface_alignment_read_bindprofx_with_cutoff(char* ialignfile, double seq_weight, double cutoff_isscore, int seq_count_cutoff);
    int interface_alignment_read_bindprofx_with_henikoff_weight(char* ialignfile, double cutoff_isscore);

    // read String sequence-based interface alignment file
    int interface_alignment_read_sip_with_cutoff(char* ialignfile, double seq_weight, double cutoff_linkscore, double seqid_high, double seqid_low, int seq_count_cutoff);
    int interface_alignment_read_sip_with_henikoff_weight(char* ialignfile, double cutoff_linkscore);

    // read ligand and receptor monomer alignment separately
  };

}



#endif
