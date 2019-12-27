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
#include "CrossValidation.h"
#include "ParameterOptimization.h"
#include "Utility.h"
#include "ErrorHandling.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>


SIP::CrossValidation::CrossValidation(){
  ;
}

SIP::CrossValidation::~CrossValidation(){
  ;
}

int SIP::CrossValidation::cross_validation_divide_train_and_test_set(char* alldatafile, double testsetratio){
  FileReader fr;
  FileReaderCreate(&fr, alldatafile);
  int tot_data_num = StringArrayGetCount(&fr.lines);
  int test_num = (int)(tot_data_num * testsetratio);
  int *flag = (int*)malloc(sizeof(int)*tot_data_num);
  for(int i = 0; i < tot_data_num; i++){
    flag[i] = 0;
  }
  // randomly select split dataset into training and test
  // 1 stands for the test set entry
  int counter = 0;
  while(counter < test_num){
    int index = rand()%tot_data_num;
    if(flag[index] == 0){
      flag[index] = 1;
      counter++;
      if(counter == test_num) break;
    }
  }
  //
  char newfilename[1024];
  sprintf(newfilename, "%s.train", alldatafile);
  FILE* trainfile = fopen(newfilename, "w");
  sprintf(newfilename, "%s.test", alldatafile);
  FILE* testfile = fopen(newfilename, "w");
  for(int i = 0; i < tot_data_num; i++){
    if(flag[i] == 1){
      fprintf(testfile, "%s\n", StringArrayGet(&fr.lines, i));
    }
    else{
      fprintf(trainfile, "%s\n", StringArrayGet(&fr.lines, i));
    }
  }
  fclose(testfile);
  fclose(trainfile);

  FileReaderDestroy(&fr);
  free(flag);
  flag = NULL;

  return Success;
}

int SIP::CrossValidation::cross_validation_divide_train_set_into_kfold(char* traindatafile, int kfold){
  if(kfold <= 1){
    printf("your kfold = %d,for cross validation, kfold should be a number greater than one\n", kfold);
    printf("it has been changed into 10 instead\n");
    kfold = 10;
  }
  FileReader fr;
  FileReaderCreate(&fr, traindatafile);
  int tot_data_num = StringArrayGetCount(&fr.lines);
  int data_num_each_fold;
  if(tot_data_num % kfold == 0){
    data_num_each_fold = (int)(tot_data_num/kfold);
  }
  else{
    data_num_each_fold = (int)(tot_data_num/kfold)+1;
  }
  int *flag = (int*)malloc(sizeof(int)*tot_data_num);
  for(int i = 0; i < tot_data_num; i++){
    flag[i] = 0;
  }
  // randomly select split dataset into training and test
  // 1 stands for the test set entry
  // the trainset is split into kfold small sets, marked from 0 to k-1
  for(int i = 1; i < kfold; i++){
    int counter = 0;
    while(counter < data_num_each_fold){
      int index = rand()%tot_data_num;
      if(flag[index] == 0){
        flag[index] = i;
        counter++;
        if(counter == data_num_each_fold) break;
      }
    }
  }

  for(int i = 0; i < kfold; i++){
    char newfilename[1024];
    sprintf(newfilename, "%s.%d", traindatafile, i);
    FILE* newfile = fopen(newfilename, "w");
    for(int j = 0; j < tot_data_num; j++){
      if(flag[j] == i){
        fprintf(newfile, "%s\n", StringArrayGet(&fr.lines, j));
      }
    }
    fclose(newfile);
  }

  FileReaderDestroy(&fr);
  free(flag);
  flag = NULL;

  return Success;
}


int SIP::CrossValidation::cross_validation_run1234(char* traindatafile, int kfold, char* working_path, char* ialign_bindprofx, char* ialign_sip, char* lig_tmalign, char* rec_tmalign, char* parameterfile){
  double cvpreddg[3000];
  double cvexpddg[3000];
  int    cvndata=0;
  StringArray cvmutants[2];
  StringArrayCreate(&cvmutants[0]);
  StringArrayCreate(&cvmutants[1]);

  
  double avg_pearson_training = 0.0;
  double avg_rmse_training = 0.0;
  for(int i = 0; i < kfold; i++){
    SIP::MCoptimization mco;
    SIP::Parameter par;
    char tempfile[1024];
    sprintf(tempfile, "%s.temp", traindatafile);
    FILE *crossfile = fopen(tempfile, "w");
    for(int j = 0; j < kfold; j++){
      // read the training data and fit model
      if(j != i){
        char tempfile1[1024];
        sprintf(tempfile1, "%s.%d", traindatafile, j);
        FileReader fr;
        FileReaderCreate(&fr, tempfile1);
        for(int i = 0; i < StringArrayGetCount(&fr.lines); i++){
          fprintf(crossfile, "%s\n", StringArrayGet(&fr.lines, i));
        }
        FileReaderDestroy(&fr);
      }
    }
    fclose(crossfile);

    //optimize parameter on the combined training dataset
    mco.optimization_by_sa1234(&par, tempfile, working_path, ialign_bindprofx, ialign_sip, lig_tmalign, rec_tmalign, parameterfile);
    double pearson, rmse;
    mco.optimization_get_pearson_and_rmse1234(&par, &pearson, &rmse);
    //par.parameter_print();
    avg_pearson_training+=pearson;
    avg_rmse_training+=rmse;
    mco.optimization_release_memory();

    //calculate ddg on the sub cv dataset and add collect them together
    sprintf(tempfile, "%s.%d", traindatafile, i);
    cross_validation_calculate_cvddg1234(&par, tempfile, working_path, ialign_bindprofx, ialign_sip, lig_tmalign, rec_tmalign, cvexpddg, cvpreddg, &cvndata, &cvmutants[0], &cvmutants[1]);
  }

  avg_pearson_training/=kfold;
  avg_rmse_training/=kfold;

  double cvpearson=0.0, cvrmse=1000.0;
  SIP::calculate_pearson_and_rmse(cvndata,cvexpddg,cvpreddg,&cvpearson,&cvrmse);
  printf("pearson & rmse of cross validation training dataset (average): %f %f\n", avg_pearson_training, avg_rmse_training);
  printf("pearson & rmse of cross validation testing dataset: %f %f\n", cvpearson, cvrmse);
  for(int k=0; k<cvndata; k++){
    printf("%-10s %7.3f %7.3f %s\n", StringArrayGet(&cvmutants[0], k), cvexpddg[k], cvpreddg[k], StringArrayGet(&cvmutants[1], k));
  }
  StringArrayDestroy(&cvmutants[0]);
  StringArrayDestroy(&cvmutants[1]);


  return Success;
}

int SIP::CrossValidation::cross_validation_calculate_pearson_and_rmse(SIP::Parameter *par, char* cvdatafile, char*workingpath, char*alignmentfilename, double *pearson, double *rmse){
  SIP::MCoptimization mco;
  mco.optimization_read_data_from_file1or2(par,cvdatafile, workingpath, alignmentfilename);
  mco.optimization_get_pearson_and_rmse(par, pearson, rmse);
  mco.optimization_release_memory();
  return Success;
}

int SIP::CrossValidation::cross_validation_calculate_cvddg1234(SIP::Parameter *par, char* cvdatafile, char*workingpath, char*ialign_bindprofx, char* ialign_sip, char* lig_tmalign, char* rec_tmalign, double* expddg, double* preddg, int *ndata, StringArray* pdbs, StringArray* mutants){
  SIP::MCoptimization mco;
  mco.optimization_read_data_from_file1234(par,cvdatafile, workingpath, ialign_bindprofx, ialign_sip, lig_tmalign, rec_tmalign);
  //mco.optimization_get_pearson_and_rmse(par, pearson, rmse);
  mco.optimization_get_cvddg1234(par,expddg,preddg,ndata,pdbs,mutants);
  mco.optimization_release_memory();
  return Success;
}

int SIP::CrossValidation::cross_validation_check_test_performance1234(SIP::Parameter *par, char* testdatafile, char*workingpath, char*ialignfile_bindprofx, char* ialignfile_sip, char* tmalignfile_lig, char* tmalignfile_rec, double *pearson, double *rmse){
  SIP::MCoptimization mco;
  mco.optimization_read_data_from_file1234(par, testdatafile, workingpath, ialignfile_bindprofx, ialignfile_sip, tmalignfile_lig, tmalignfile_rec);
  mco.optimization_get_pearson_and_rmse1234(par, pearson, rmse);
  mco.optimization_print_ddgs();
  mco.optimization_release_memory();
  return Success;
}