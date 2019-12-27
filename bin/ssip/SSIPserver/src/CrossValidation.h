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

#ifndef CROSS_VALIDATION_H
#define CROSS_VALIDATION_H

#include "Parameter.h"

namespace SIP{
  class CrossValidation{
  public:
    SIP::Parameter par0;

    CrossValidation();
    ~CrossValidation();
    int cross_validation_divide_train_and_test_set(char* alldatafile, double testsetratio);
    int cross_validation_divide_train_set_into_kfold(char* traindatafile, int kfold);
    int cross_validation_calculate_pearson_and_rmse(SIP::Parameter *par, char* cvdatafile, char*workingpath, char*alignmentfilename, double *pearson, double *rmse);
    int cross_validation_check_test_performance1234(SIP::Parameter *par, char* testdatafile, char*workingpath, char*ialignfile_bindprofx, char* ialignfile_sip, char* tmalignfile_lig, char* tmalignfile_rec, double *pearson, double *rmse);
    int cross_validation_run1234(char* traindatafile, int kfold, char* working_path, char* ialign_bindprofx, char* ialign_sip, char* lig_tmalign, char* rec_tmalign, char* parameterfile);

    int cross_validation_calculate_cvddg1234(SIP::Parameter *par, char* cvdatafile, char*workingpath, char*ialign_bindprofx, char* ialign_sip, char* lig_tmalign, char* rec_tmalign, double* expddg, double* preddg, int *ndata, StringArray* pdbs, StringArray* mutants);
  };
}
#endif