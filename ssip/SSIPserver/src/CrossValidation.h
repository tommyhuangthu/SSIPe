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