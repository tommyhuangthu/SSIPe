#ifndef PARAMETER_OPTIMIZATION_H
#define PARAMETER_OPTIMIZATION_H

#include <stdio.h>
#include "PseudoCount.h"
#include "Utility.h"
#include "Mutant.h"
#include "Parameter.h"
#include "DDG.h"

// scale of dataset;
#define MAX_NUM_PDB_ENTRY       500
#define MAX_NUM_DDG_VALUE       5000


#define MONTE_CARLO_MAX_STEP    200000
#define MONTE_CARLO_TEMPERATURE 1e-4

#define SIMUL_ANNEA_TMAX        1e-3
#define SIMUL_ANNEA_TMIN        1e-4
#define SIMUL_ANNEA_TDECREASE   0.5
#define SIMUL_ANNEA_STEP        200000
#define SIMUL_ANNEA_CYCLE       1000


#define REPLICA_EXCH_REPLICAS   20
#define REPLICA_EXCH_NSWAP      100
#define REPLICA_EXCH_NMETROPLIS 500
#define REPLICA_EXCH_TMIN       1e-5
#define REPLICA_EXCH_TMAX       1e-4

namespace SIP{
  int generate_mutant_file(char* mutantfile, char* working_path);
  int calculate_pearson_and_rmse(int nentry, double *exp_ddgs, double *cal_ddgs, double *pearson, double *rmse);
  int metropolis_criteria(double old_val, double new_val, double temperature);
  int swap_criteria(double val1, double val2, double temp1, double temp2);
  double define_objective_function(double pearson, double rmse, SIP::Parameter *par);
  int print_pearson_and_rmse(double pearson, double rmse);

  // read data from file for only once
  class Data{
    public:
      int npdb;
      char pdblist[MAX_NUM_PDB_ENTRY][10];
      int nmutants;
      double expddgs[MAX_NUM_DDG_VALUE];
      StringArray pdbs;
      StringArray mutants;

      Data();
      ~Data();
  };

  class MCoptimization{
    public:
      Parameter best_par;
      Data expdata;
      SIP::DDG *ddgs;
      int nentry;
      double calddgs[MAX_NUM_DDG_VALUE];
      double best_pearson;
      double best_rmse;
      double bestobj;

      MCoptimization();
      ~MCoptimization();
      int optimization_read_data_from_file1or2(SIP::Parameter *par,char* mutantfile, char* working_path, char* ialignfile_bindprofx);
      int optimization_append_calcddg(SIP::DDG *ddg);
      int optimization_get_pearson_and_rmse(SIP::Parameter *par, double* pearson, double *rmse);
      int optimization_print_pearson_and_rmse();
      int optimization_print_ddgs();
      int optimization_release_memory();
      int remc_local_movement(SIP::Parameter *par, double *pearson, double *rmse, double temperature, int step);

      int optimization_by_mc1or2(SIP::Parameter *par, char *mutantfile, char *working_path, char* alignmentfilename);
      int optimization_by_sa1or2(SIP::Parameter *par, char *mutantfile, char *working_path, char* alignmentfilename);
      int optimization_by_remc1or2(SIP::Parameter *par, char *mutantfile, char *working_path, char* alignmentfilename);


      // other functions
      int optimization_show_interface_alignment_gap_count();

      // optimizing parameters with combined interface alignment
      int optimization_read_data_from_file12(SIP::Parameter *par, char* mutantfile, char* working_path, char* alignfile_bindprofx, char *alignfile_sip);
      int optimization_by_mc12(SIP::Parameter *par, char *mutantfile, char *working_path, char* ialign_bindprofx, char* ialign_sip);
      int optimization_get_pearson_and_rmse12(SIP::Parameter *par, double* pearson, double *rmse);
      int optimization_by_sa12(SIP::Parameter *par, char *mutantfile, char *working_path, char* alignfile_bindprofx, char* alignfile_sip);

      //
      int optimization_read_data_from_file1234(SIP::Parameter* par,char* mutantfile, char* working_path, char* alignfile_bindprofx, char *alignfile_sip, char* lig_alignment, char* rec_alignment);
      int optimization_by_sa1234(SIP::Parameter *par, char *mutantfile, char *working_path, char* ialignfile_bindprofx, char* ialignfile_sip, char* tmalignfile_lig, char* tmalignfile_rec, char* parameterfile);
      int optimization_get_pearson_and_rmse1234(SIP::Parameter *par, double* pearson, double *rmse);
      
      //for cross validation
      int optimization_get_cvddg1234(SIP::Parameter *par, double* exp, double* pre, int* ndata, StringArray* pdbs, StringArray* mutants);
  };

}

#endif
