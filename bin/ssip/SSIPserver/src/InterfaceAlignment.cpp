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

#pragma warning(disable:4996)
#include "InterfaceAlignment.h"
#include "Utility.h"
#include "ErrorHandling.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "AminoacidOrder.h"

int SIP::InterfaceAlignment::interface_alignment_count_amino_acid_number_in_string(double table[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE], char* string, int length){
  for(int i = 0; i < length; i++){
    switch(string[i]){
      case 'A': table[i][0] += 1.0;  break;
      case 'C': table[i][1] += 1.0;  break;
      case 'D': table[i][2] += 1.0;  break;
      case 'E': table[i][3] += 1.0;  break;
      case 'F': table[i][4] += 1.0;  break;
      case 'G': table[i][5] += 1.0;  break;
      case 'H': table[i][6] += 1.0;  break;
      case 'I': table[i][7] += 1.0;  break;
      case 'K': table[i][8] += 1.0;  break;
      case 'L': table[i][9] += 1.0;  break;
      case 'M': table[i][10] += 1.0; break;
      case 'N': table[i][11] += 1.0; break;
      case 'P': table[i][12] += 1.0; break;
      case 'Q': table[i][13] += 1.0; break;
      case 'R': table[i][14] += 1.0; break;
      case 'S': table[i][15] += 1.0; break;
      case 'T': table[i][16] += 1.0; break;
      case 'V': table[i][17] += 1.0; break;
      case 'W': table[i][18] += 1.0; break;
      case 'Y': table[i][19] += 1.0; break;
      case '-': table[i][20] += 1.0; break;
      default: break;
    }
  }
  return 0;
}

int SIP::InterfaceAlignment::interface_alignment_count_amino_acid_number_in_string_with_weight(double table[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE], char* string, int length, double seq_weight){
  for(int i = 0; i < length; i++){
    switch(string[i]){
      case 'A': table[i][0] += seq_weight;  break;
      case 'C': table[i][1] += seq_weight;  break;
      case 'D': table[i][2] += seq_weight;  break;
      case 'E': table[i][3] += seq_weight;  break;
      case 'F': table[i][4] += seq_weight;  break;
      case 'G': table[i][5] += seq_weight;  break;
      case 'H': table[i][6] += seq_weight;  break;
      case 'I': table[i][7] += seq_weight;  break;
      case 'K': table[i][8] += seq_weight;  break;
      case 'L': table[i][9] += seq_weight;  break;
      case 'M': table[i][10] += seq_weight; break;
      case 'N': table[i][11] += seq_weight; break;
      case 'P': table[i][12] += seq_weight; break;
      case 'Q': table[i][13] += seq_weight; break;
      case 'R': table[i][14] += seq_weight; break;
      case 'S': table[i][15] += seq_weight; break;
      case 'T': table[i][16] += seq_weight; break;
      case 'V': table[i][17] += seq_weight; break;
      case 'W': table[i][18] += seq_weight; break;
      case 'Y': table[i][19] += seq_weight; break;
      case '-': table[i][20] += seq_weight; break;
      default: break;
    }
  }
  return 0;
}

SIP::InterfaceAlignment::InterfaceAlignment(){
  lig_alignment_num=0;
  rec_alignment_num=0;
}


SIP::InterfaceAlignment::~InterfaceAlignment(){
  ;
}

int SIP::InterfaceAlignment::interface_alignment_initialize(){
  lig_alignment_num = 0;
  rec_alignment_num = 0;
  for(int i=0; i<MAX_NUM_INTERFACE_RESI; i++){
    for(int j=0; j<MAX_NUM_AA_TYPE; j++){
      lig_aligntable[i][j]=0;
      rec_aligntable[i][j]=0;
      eff_ligalign[i][j]=0;
      eff_recalign[i][j]=0;
    }
  }
  return Success;
}



int SIP::InterfaceAlignment::interface_alignment_print_statistics(){
  // test and output the information read from ialign.txt
  printf("The ligand interface residues are:\n");
  for(int i = 0; i < lig_res_count; i++){
    printf("%s%d%c\n", lig_chn_name, lig_res_index[i], native_lig_aas[i]);
  }
  printf("The receptor interface residues are:\n");
  for(int i = 0; i < rec_res_count; i++){
    printf("%s%d%c\n", rec_chn_name, rec_res_index[i], native_rec_aas[i]);
  }
  printf("Alignment table:\n");
  printf("---------Ligand protein-------------\n");
  printf("       A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y      -\n");
  for(int j = 0; j < lig_res_count; j++){
    printf("%c ", native_lig_aas[j]);
    for(int k = 0; k < 21; k++){
      printf("%6.3f ", lig_aligntable[j][k]);
    }
    printf("\n");
  }
  printf("---------Receptor protein-------------\n");
  printf("       A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y      -\n");
  for(int j = 0; j < rec_res_count; j++){
    printf("%c ", native_rec_aas[j]);
    for(int k = 0; k < 21; k++){
      printf("%6.3f ", rec_aligntable[j][k]);
    }
    printf("\n");
  }

  //printf("Effective amino acid number:\n");
  //printf("---------Ligand protein-------------\n");
  //printf("       A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y      -\n");
  //for(int j = 0; j < lig_res_count; j++){
  //  printf("%c ", native_lig_aas[j]);
  //  for(int k = 0; k < 21; k++){
  //    printf("%6.3f ", eff_ligalign[j][k]);
  //  }
  //  printf("\n");
  //}
  //printf("---------Receptor protein-------------\n");
  //printf("       A      C      D      E      F      G      H      I      K      L      M      N      P      Q      R      S      T      V      W      Y      -\n");
  //for(int j = 0; j < rec_res_count; j++){
  //  printf("%c ", native_rec_aas[j]);
  //  for(int k = 0; k < 21; k++){
  //    printf("%6.3f ", eff_recalign[j][k]);
  //  }
  //  printf("\n");
  //}

  return 0;
}

double SIP::InterfaceAlignment::interface_alignment_get_gap_count(char chnid, int alnindex){
  if(alnindex < 0){
    return 0;
  }
  if(chnid == lig_chn_name[0]){
    return lig_aligntable[alnindex][20];
  }
  else if(chnid == rec_chn_name[0]){
    return rec_aligntable[alnindex][20];
  }

  return 0;
}

double SIP::InterfaceAlignment::interface_alignment_get_gap_ratio(char chnid, int alnindex){
  if(alnindex < 0){
    return 0;
  }
  if(chnid == lig_chn_name[0]){
    double tot = 0;
    for(int i = 0; i < 21; i++){
      tot += lig_aligntable[alnindex][i];
    }
    if(tot > 0){
      return (double)lig_aligntable[alnindex][20]/tot;
    }
    else{
      return 0.0;
    }
  }
  else if(chnid == rec_chn_name[0]){
    double tot = 0;
    for(int i = 0; i < 21; i++){
      tot += rec_aligntable[alnindex][i];
    }
    if(tot > 0){
      return (double)rec_aligntable[alnindex][20]/tot;
    }
    else{
      return 0.0;
    }
  }

  return 0;
}




// read the interface alignment file with given linkscore cutoff from 0 to 1;
// and sequence identity cutoff from 0 to 1
// the link score cutoff should not be too low, link score > 0.5?
// the sequence identity cutoff should not be too high to remove redundant sequence alignment, identity cutoff <= 0.5?
int SIP::InterfaceAlignment::interface_alignment_read_sip_with_cutoff(char* ialignfile, double seq_weight, double cutoff_linkscore, double seqid_high, double seqid_low, int seq_count_cutoff){
  FileReader fr;
  if(FAILED(FileReaderCreate(&fr, ialignfile))){
    printf("failed open file %s\n", ialignfile);
    exit(IOError);
  }
  strcpy(alignment_filename, ialignfile);

  char line[1000];
  char eff_ligstrings[MAX_NUM_EFF_INF_ALN][MAX_NUM_INTERFACE_RESI];
  char eff_recstrings[MAX_NUM_EFF_INF_ALN][MAX_NUM_INTERFACE_RESI];
  int eff_aln_num = 0;
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, line, ' ');
    if(strcmp(StringArrayGet(&strings, 0), "chain_ID_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_A:") == 0){
      strcpy(lig_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "chain_ID_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_B:") == 0){
      strcpy(rec_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_A:") == 0){
      lig_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        lig_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < lig_res_count; j++){
        for(int k = 0; k < 21; k++){
          lig_aligntable[j][k] = 0;
          eff_ligalign[j][k] = 0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_B:") == 0){
      rec_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        rec_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < rec_res_count; j++){
        for(int k = 0; k < 21; k++){
          rec_aligntable[j][k] = 0;
          eff_recalign[j][k] = 0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "query") == 0 || strcmp(StringArrayGet(&strings, 0), "target") == 0){
      char ligstring[200];
      char recstring[200];
      strcpy(ligstring, StringArrayGet(&strings, 6));
      strcpy(recstring, StringArrayGet(&strings, 7));
      interface_alignment_count_amino_acid_number_in_string_with_weight(lig_aligntable, ligstring, strlen(ligstring), seq_weight);
      interface_alignment_count_amino_acid_number_in_string_with_weight(rec_aligntable, recstring, strlen(recstring), seq_weight);
      strcpy(native_lig_aas, ligstring);
      strcpy(native_rec_aas, recstring);
      strcpy(eff_ligstrings[eff_aln_num], ligstring);
      strcpy(eff_recstrings[eff_aln_num], recstring);
      lig_alignment_num++;
      eff_aln_num++;
    }
    else if(strstr(StringArrayGet(&strings, 0), "aln") != 0 || isdigit(StringArrayGet(&strings, 0)[0])){
      double linkscore = atof(StringArrayGet(&strings, 5));
      int inf_num = atoi(StringArrayGet(&strings, 3));
      int aln_num = atoi(StringArrayGet(&strings, 4));
      //if((double)aln_num/inf_num < 0.90) continue;
      // only consider the alignment whose linkscore is larger than cutoff
      if(linkscore >= cutoff_linkscore){
      //if(linkscore*aln_num/inf_num >= cutoff_linkscore){
        char ligstring[MAX_NUM_INTERFACE_RESI];
        char recstring[MAX_NUM_INTERFACE_RESI];
        strcpy(ligstring, StringArrayGet(&strings, 6));
        strcpy(recstring, StringArrayGet(&strings, 7));
        // compare current alignment with previous stored alignments
        int alignment_is_effective = 1;
        for(int x = 0; x < eff_aln_num; x++){
          char* str1 = eff_ligstrings[x];
          char* str2 = eff_recstrings[x];
          int num_sameaa = 0;
          for(int j = 0; j < (int)strlen(str1); j++){
            if(ligstring[j] == str1[j]) num_sameaa++;
          }
          double ratio = (double)num_sameaa/(strlen(str1));
          if(ratio > seqid_high || ratio < seqid_low){
            alignment_is_effective = 0;
            break;
          }
          num_sameaa = 0;
          for(int j = 0; j < (int)strlen(str2); j++){
            if(recstring[j] == str2[j]) num_sameaa++;
          }
          ratio = (double)num_sameaa/(strlen(str2));
          if(ratio > seqid_high || ratio < seqid_low){
            alignment_is_effective = 0;
            break;
          }
        }
        if(alignment_is_effective == 1){
          interface_alignment_count_amino_acid_number_in_string_with_weight(lig_aligntable, ligstring, strlen(ligstring), seq_weight);
          interface_alignment_count_amino_acid_number_in_string_with_weight(rec_aligntable, recstring, strlen(recstring), seq_weight);
          lig_alignment_num++;
          // save the alignment as an effective alignment
          strcpy(eff_ligstrings[eff_aln_num], ligstring);
          strcpy(eff_recstrings[eff_aln_num], recstring);
          eff_aln_num++;
        }
      }
    }
    StringArrayDestroy(&strings);
    if(eff_aln_num >= seq_count_cutoff) break;
  }
  FileReaderDestroy(&fr);

  rec_alignment_num = lig_alignment_num;

  return Success;
}

// read the interface alignment file with given linkscore cutoff from 0 to 1;
// the link score cutoff should not be too low, link score > 0.5?
int SIP::InterfaceAlignment::interface_alignment_read_sip_with_henikoff_weight(char* ialignfile, double cutoff_linkscore){
  FileReader fr;
  if(FAILED(FileReaderCreate(&fr, ialignfile))){
    printf("failed open file %s\n", ialignfile);
    exit(IOError);
  }
  strcpy(alignment_filename, ialignfile);

  StringArray ligstrings;
  StringArray recstrings;
  StringArrayCreate(&ligstrings);
  StringArrayCreate(&recstrings);
  char line[2000];
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, line, ' ');
    if(strcmp(StringArrayGet(&strings, 0), "chain_ID_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_A:") == 0){
      strcpy(lig_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "chain_ID_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_B:") == 0){
      strcpy(rec_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_A:") == 0){
      lig_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        lig_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < lig_res_count; j++){
        for(int k = 0; k < 21; k++){
          lig_aligntable[j][k] = 0;
          eff_ligalign[j][k] = 0.0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_B:") == 0){
      rec_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        rec_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < rec_res_count; j++){
        for(int k = 0; k < 21; k++){
          rec_aligntable[j][k] = 0;
          eff_recalign[j][k] = 0.0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "query") == 0 || strcmp(StringArrayGet(&strings, 0), "target") == 0){
      char ligstring[200];
      char recstring[200];
      strcpy(ligstring, StringArrayGet(&strings, 6));
      strcpy(recstring, StringArrayGet(&strings, 7));
      interface_alignment_count_amino_acid_number_in_string(lig_aligntable, ligstring, lig_res_count);
      interface_alignment_count_amino_acid_number_in_string(rec_aligntable, recstring, rec_res_count);
      strcpy(native_lig_aas, ligstring);
      strcpy(native_rec_aas, recstring);
      StringArrayAppend(&ligstrings, ligstring);
      StringArrayAppend(&recstrings, recstring);
      lig_alignment_num++;
    }
    else if(strstr(StringArrayGet(&strings, 0), "aln") != 0 || isdigit(StringArrayGet(&strings, 0)[0])){
      double linkscore = atof(StringArrayGet(&strings, 5));
      // only consider the alignment whose linkscore is larger than cutoff
      if(linkscore >= cutoff_linkscore){
        char ligstring[MAX_NUM_INTERFACE_RESI];
        char recstring[MAX_NUM_INTERFACE_RESI];
        strcpy(ligstring, StringArrayGet(&strings, 6));
        strcpy(recstring, StringArrayGet(&strings, 7));
        interface_alignment_count_amino_acid_number_in_string(lig_aligntable, ligstring, lig_res_count);
        interface_alignment_count_amino_acid_number_in_string(rec_aligntable, recstring, rec_res_count);
        StringArrayAppend(&ligstrings, ligstring);
        StringArrayAppend(&recstrings, recstring);
        lig_alignment_num++;
      }
    }
    StringArrayDestroy(&strings);
  }
  FileReaderDestroy(&fr);
  rec_alignment_num = lig_alignment_num;

  // calculate w_ij, i is the index of interface residue position, j is the index of sequence
  double weights_aatype_lig[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE];
  double weights_aatype_rec[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE];
  for(int i = 0; i < lig_res_count; i++){
    int num_aatype = 0;
    for(int j = 0; j < MAX_NUM_AA_TYPE; j++){
      if(lig_aligntable[i][j] > 0){
        num_aatype++;
      }
    }
    for(int j = 0; j < MAX_NUM_AA_TYPE; j++){
      if(lig_aligntable[i][j] > 0){
        weights_aatype_lig[i][j] = 1.0/(num_aatype*lig_aligntable[i][j]);
      }
      else{
        weights_aatype_lig[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < rec_res_count; i++){
    int num_aatype = 0;
    for(int j = 0; j < MAX_NUM_AA_TYPE; j++){
      if(rec_aligntable[i][j] > 0){
        num_aatype++;
      }
    }
    for(int j = 0; j < MAX_NUM_AA_TYPE; j++){
      if(rec_aligntable[i][j] > 0){
        weights_aatype_rec[i][j] = 1.0/(num_aatype*rec_aligntable[i][j]);
      }
      else{
        weights_aatype_rec[i][j] = 0.0;
      }
    }
  }
  // calculate the weights for each alignment sequence
  // here, '1' stands for the query sequence
  for(int j = 0; j < lig_alignment_num; j++){
    double weight_seq = 0.0;
    char* ligstr = StringArrayGet(&ligstrings, j);
    for(int i = 0; i < lig_res_count; i++){
      char aa = ligstr[i];
      AminoacidOrder aao;
      aao.get_aminoacid_order(aa);
      weight_seq += weights_aatype_lig[i][aao.aaOrder];
    }
    char* recstr = StringArrayGet(&recstrings, j);
    for(int i = 0; i < rec_res_count; i++){
      char aa = recstr[i];
      AminoacidOrder aao;
      aao.get_aminoacid_order(aa);
      weight_seq += weights_aatype_rec[i][aao.aaOrder];
    }
    weight_seq = weight_seq*lig_alignment_num/(lig_res_count + rec_res_count);
    interface_alignment_count_amino_acid_number_in_string_with_weight(eff_ligalign, ligstr, lig_res_count, weight_seq);
    interface_alignment_count_amino_acid_number_in_string_with_weight(eff_recalign, recstr, rec_res_count, weight_seq);
  }

  return Success;
}


// read the interface alignment file with given isscore cutoff from 0 to 1;
// the optimal isscore cutoff is 0.55 according to the original bindprofx paper
int SIP::InterfaceAlignment::interface_alignment_read_bindprofx(char* ialignfile, double seq_weight, double cutoff_isscore){
  FILE* pf=fopen(ialignfile,"r");
  if(pf==NULL) return IOError;
  fclose(pf);
  FileReader fr;
  if(FAILED(FileReaderCreate(&fr, ialignfile))){
    printf("failed open file %s\n", ialignfile);
    exit(IOError);
  }
  strcpy(alignment_filename, ialignfile);

  char line[1000];
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, line, ' ');
    if(strcmp(StringArrayGet(&strings, 0), "chain_ID_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_A:") == 0){
      strcpy(lig_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "chain_ID_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_B:") == 0){
      strcpy(rec_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_A:") == 0){
      lig_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        lig_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < lig_res_count; j++){
        for(int k = 0; k < 21; k++){
          lig_aligntable[j][k] = 0;
          eff_ligalign[j][k] = 0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_B:") == 0){
      rec_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        rec_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < rec_res_count; j++){
        for(int k = 0; k < 21; k++){
          rec_aligntable[j][k] = 0;
          eff_recalign[j][k] = 0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "query") == 0 || strcmp(StringArrayGet(&strings, 0), "target") == 0){
      char ligstring[MAX_NUM_INTERFACE_RESI];
      char recstring[MAX_NUM_INTERFACE_RESI];
      strcpy(ligstring, StringArrayGet(&strings, 6));
      strcpy(recstring, StringArrayGet(&strings, 7));
      interface_alignment_count_amino_acid_number_in_string_with_weight(lig_aligntable, ligstring, strlen(ligstring), seq_weight);
      interface_alignment_count_amino_acid_number_in_string_with_weight(rec_aligntable, recstring, strlen(recstring), seq_weight);
      strcpy(native_lig_aas, ligstring);
      strcpy(native_rec_aas, recstring);
      lig_alignment_num++;
    }
    else if(strstr(StringArrayGet(&strings, 0), "aln") != 0 || isdigit(StringArrayGet(&strings, 0)[0])){
      double isscore = atof(StringArrayGet(&strings, 5));
      // only consider the alignment whose linkscore is larger than cutoff
      if(isscore > cutoff_isscore){
        if(fabs(isscore-1.0)<1e-4) continue;
        char ligstring[MAX_NUM_INTERFACE_RESI];
        char recstring[MAX_NUM_INTERFACE_RESI];
        strcpy(ligstring, StringArrayGet(&strings, 6));
        strcpy(recstring, StringArrayGet(&strings, 7));
        interface_alignment_count_amino_acid_number_in_string_with_weight(lig_aligntable, ligstring, strlen(ligstring), seq_weight);
        interface_alignment_count_amino_acid_number_in_string_with_weight(rec_aligntable, recstring, strlen(recstring), seq_weight);
        lig_alignment_num++;
      }
    }
    StringArrayDestroy(&strings);
  }
  FileReaderDestroy(&fr);
  rec_alignment_num = lig_alignment_num;

  return Success;
}

int SIP::InterfaceAlignment::interface_alignment_read_bindprofx_with_cutoff(char* ialignfile, double seq_weight, double cutoff_isscore, int seq_count_cutoff){
  FileReader fr;
  if(FAILED(FileReaderCreate(&fr, ialignfile))){
    printf("failed open file %s\n", ialignfile);
    exit(IOError);
  }
  strcpy(alignment_filename, ialignfile);

  char line[1000];
  int eff_aln_num = 0;
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, line, ' ');
    if(strcmp(StringArrayGet(&strings, 0), "chain_ID_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_A:") == 0){
      strcpy(lig_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "chain_ID_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_B:") == 0){
      strcpy(rec_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_A:") == 0){
      lig_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        lig_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < lig_res_count; j++){
        for(int k = 0; k < 21; k++){
          lig_aligntable[j][k] = 0;
          eff_ligalign[j][k] = 0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_B:") == 0){
      rec_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        rec_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < rec_res_count; j++){
        for(int k = 0; k < 21; k++){
          rec_aligntable[j][k] = 0;
          eff_recalign[j][k] = 0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "query") == 0 || strcmp(StringArrayGet(&strings, 0), "target") == 0){
      char ligstring[MAX_NUM_INTERFACE_RESI];
      char recstring[MAX_NUM_INTERFACE_RESI];
      strcpy(ligstring, StringArrayGet(&strings, 6));
      strcpy(recstring, StringArrayGet(&strings, 7));
      interface_alignment_count_amino_acid_number_in_string_with_weight(lig_aligntable, ligstring, strlen(ligstring), seq_weight);
      interface_alignment_count_amino_acid_number_in_string_with_weight(rec_aligntable, recstring, strlen(recstring), seq_weight);
      strcpy(native_lig_aas, ligstring);
      strcpy(native_rec_aas, recstring);
      lig_alignment_num++;
      eff_aln_num++;
    }
    else if(strstr(StringArrayGet(&strings, 0), "aln") != 0 || isdigit(StringArrayGet(&strings, 0)[0])){
      double isscore = atof(StringArrayGet(&strings, 5));
      // only consider the alignment whose linkscore is larger than cutoff
      if(isscore > cutoff_isscore){
        if(fabs(isscore-1.0)<1e-4) continue;
        char ligstring[MAX_NUM_INTERFACE_RESI];
        char recstring[MAX_NUM_INTERFACE_RESI];
        strcpy(ligstring, StringArrayGet(&strings, 6));
        strcpy(recstring, StringArrayGet(&strings, 7));
        interface_alignment_count_amino_acid_number_in_string_with_weight(lig_aligntable, ligstring, strlen(ligstring), seq_weight);
        interface_alignment_count_amino_acid_number_in_string_with_weight(rec_aligntable, recstring, strlen(recstring), seq_weight);
        lig_alignment_num++;
        eff_aln_num++;
      }
    }
    StringArrayDestroy(&strings);
    if(eff_aln_num >= seq_count_cutoff) break;
  }
  FileReaderDestroy(&fr);
  rec_alignment_num = lig_alignment_num;

  return Success;
}


// read the interface alignment file with given isscore cutoff from 0 to 1;
// the optimal isscore cutoff is 0.55 according to the original bindprofx paper
int SIP::InterfaceAlignment::interface_alignment_read_bindprofx_with_henikoff_weight(char* ialignfile, double cutoff_isscore){
  FileReader fr;
  if(FAILED(FileReaderCreate(&fr, ialignfile))){
    printf("failed open file %s\n", ialignfile);
    exit(IOError);
  }
  strcpy(alignment_filename, ialignfile);

  StringArray ligstrings;
  StringArray recstrings;
  StringArrayCreate(&ligstrings);
  StringArrayCreate(&recstrings);

  char line[1000];
  while(!FAILED(FileReaderGetNextLine(&fr, line))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, line, ' ');
    if(strcmp(StringArrayGet(&strings, 0), "chain_ID_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_A:") == 0){
      strcpy(lig_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "chain_ID_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "CHAIN_ID_B:") == 0){
      strcpy(rec_chn_name, StringArrayGet(&strings, 1));
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_A:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_A:") == 0){
      lig_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        lig_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < lig_res_count; j++){
        for(int k = 0; k < 21; k++){
          lig_aligntable[j][k] = 0;
          eff_ligalign[j][k] = 0.0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "residue_B:") == 0 || strcmp(StringArrayGet(&strings, 0), "RESIDUE_B:") == 0){
      rec_res_count = StringArrayGetCount(&strings) - 1;
      for(int j = 1; j < StringArrayGetCount(&strings); j++){
        rec_res_index[j-1] = atoi(StringArrayGet(&strings, j));
      }
      for(int j = 0; j < rec_res_count; j++){
        for(int k = 0; k < 21; k++){
          rec_aligntable[j][k] = 0;
          eff_recalign[j][k] = 0.0;
        }
      }
    }
    else if(strcmp(StringArrayGet(&strings, 0), "query") == 0 || strcmp(StringArrayGet(&strings, 0), "target") == 0){
      char ligstring[200];
      char recstring[200];
      strcpy(ligstring, StringArrayGet(&strings, 6));
      strcpy(recstring, StringArrayGet(&strings, 7));
      strcpy(native_lig_aas, ligstring);
      strcpy(native_rec_aas, recstring);
      interface_alignment_count_amino_acid_number_in_string(lig_aligntable, ligstring, lig_res_count);
      interface_alignment_count_amino_acid_number_in_string(rec_aligntable, recstring, rec_res_count);
      StringArrayAppend(&ligstrings, ligstring);
      StringArrayAppend(&recstrings, recstring);
      lig_alignment_num++;
    }
    else if(strstr(StringArrayGet(&strings, 0), "aln") != 0 || isdigit(StringArrayGet(&strings, 0)[0])){
      double isscore = atof(StringArrayGet(&strings, 5));
      // only consider the alignment whose linkscore is larger than cutoff
      if(isscore > cutoff_isscore){
        //if(fabs(isscore-1.0)<1e-4) continue;
        char ligstring[MAX_NUM_INTERFACE_RESI];
        char recstring[MAX_NUM_INTERFACE_RESI];
        strcpy(ligstring, StringArrayGet(&strings, 6));
        strcpy(recstring, StringArrayGet(&strings, 7));
        interface_alignment_count_amino_acid_number_in_string(lig_aligntable, ligstring, lig_res_count);
        interface_alignment_count_amino_acid_number_in_string(rec_aligntable, recstring, rec_res_count);
        StringArrayAppend(&ligstrings, ligstring);
        StringArrayAppend(&recstrings, recstring);
        lig_alignment_num++;
      }
    }
    StringArrayDestroy(&strings);
  }
  FileReaderDestroy(&fr);
  rec_alignment_num = lig_alignment_num;

  // calculate w_ij, i is the index of interface residue position, j is the index of sequence
  double weights_aatype_lig[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE];
  double weights_aatype_rec[MAX_NUM_INTERFACE_RESI][MAX_NUM_AA_TYPE];
  for(int i = 0; i < lig_res_count; i++){
    int num_aatype = 0;
    for(int j = 0; j < MAX_NUM_AA_TYPE; j++){
      if(lig_aligntable[i][j] > 0){
        num_aatype++;
      }
    }
    for(int j = 0; j < MAX_NUM_AA_TYPE; j++){
      if(lig_aligntable[i][j] > 0){
        weights_aatype_lig[i][j] = 1.0/(num_aatype*lig_aligntable[i][j]);
      }
      else{
        weights_aatype_lig[i][j] = 0.0;
      }
    }
  }
  for(int i = 0; i < rec_res_count; i++){
    int num_aatype = 0;
    for(int j = 0; j < MAX_NUM_AA_TYPE; j++){
      if(rec_aligntable[i][j] > 0){
        num_aatype++;
      }
    }
    for(int j = 0; j < MAX_NUM_AA_TYPE; j++){
      if(rec_aligntable[i][j] > 0){
        weights_aatype_rec[i][j] = 1.0/(num_aatype*rec_aligntable[i][j]);
      }
      else{
        weights_aatype_rec[i][j] = 0.0;
      }
    }
  }
  // calculate the weights for each alignment sequence
  // here, '1' stands for the query sequence
  for(int j = 0; j < lig_alignment_num; j++){
    double weight_seq = 0.0;
    char* ligstr = StringArrayGet(&ligstrings, j);
    for(int i = 0; i < lig_res_count; i++){
      char aa = ligstr[i];
      AminoacidOrder aao;
      aao.get_aminoacid_order(aa);
      weight_seq += weights_aatype_lig[i][aao.aaOrder];
    }
    char* recstr = StringArrayGet(&recstrings, j);
    for(int i = 0; i < rec_res_count; i++){
      char aa = recstr[i];
      AminoacidOrder aao;
      aao.get_aminoacid_order(aa);
      weight_seq += weights_aatype_rec[i][aao.aaOrder];
    }
    weight_seq = weight_seq*lig_alignment_num/(lig_res_count + rec_res_count);
    interface_alignment_count_amino_acid_number_in_string_with_weight(eff_ligalign, ligstr, lig_res_count, weight_seq);
    interface_alignment_count_amino_acid_number_in_string_with_weight(eff_recalign, recstr, rec_res_count, weight_seq);
  }

  return Success;
}

