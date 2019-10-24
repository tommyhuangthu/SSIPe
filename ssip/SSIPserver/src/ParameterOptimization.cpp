#pragma warning(disable:4996)
#include "ParameterOptimization.h"
#include "Utility.h"
#include "ErrorHandling.h"
#include "PseudoCount.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdio.h>

// general functions 
int SIP::generate_mutant_file(char* mutantfile, char* working_path){
  FileReader fr;
  FileReaderCreate(&fr, mutantfile);
  char buffer[1024];
  StringArray pdbs;
  StringArrayCreate(&pdbs);
  while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, buffer, ' ');
    int index = -1;
    StringArrayFind(&pdbs, StringArrayGet(&strings, 0), &index);
    if(index == -1){
      StringArrayAppend(&pdbs, StringArrayGet(&strings, 0));
    }
    StringArrayDestroy(&strings);
  }

  // read file for the second time to generate mutant file
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    FileReaderSetCurrentPos(&fr, 0);
    char *pdbid = StringArrayGet(&pdbs, i);
    FILE *fout = NULL, *fout2 = NULL;
    char filename[1024], filename2[1024];
    sprintf(filename, "%s/%s/mutation.txt", working_path, pdbid, pdbid);
    sprintf(filename2, "%s/%s/ddg.txt", working_path, pdbid, pdbid);
    fout = fopen(filename, "w");
    fout2 = fopen(filename2, "w");
    while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
      StringArray strings;
      StringArrayCreate(&strings);
      StringArraySplitString(&strings, buffer, ' ');
      if(strcmp(pdbid, StringArrayGet(&strings, 0)) == 0){
        // type 0: show mutations in the form of "A C100A"
        // type 1: show mutations in the form of "CA100A;DA101A;"
        // type 2: show mutations in the form of "CA100A_DA101A"
        // else  : show mutations in the form of "CA100A,DA101A;"
        int type = 2;
        if(type == 0){
          char *mutstr = StringArrayGet(&strings, 5);
          char chain = mutstr[1];
          mutstr[1] = mutstr[0];
          mutstr[0] = ' ';
          fprintf(fout, "%c%s\n", chain, mutstr);
          fprintf(fout2, "%s\n", StringArrayGet(&strings, 4));
        }
        else if(type == 1){
          char *mutstr = StringArrayGet(&strings, 5);
          char temp[1000];
          strcpy(temp, mutstr);
          for(int j = 0; j < (int)strlen(temp); j++){
            if(temp[j] == ',') temp[j] = ';';
          }
          fprintf(fout, "%s;\n", temp);
          fprintf(fout2, "%s\n", StringArrayGet(&strings, 4));
        }
        else if(type == 2){
          char *mutstr = StringArrayGet(&strings, 5);
          char temp[1000];
          strcpy(temp, mutstr);
          for(int j = 0; j < (int)strlen(temp); j++){
            if(temp[j] == ',') temp[j] = '_';
          }
          fprintf(fout, "%s\n", temp);
          fprintf(fout2, "%s\n", StringArrayGet(&strings, 4));
        }
        else{
          char *mutstr = StringArrayGet(&strings, 5);
          fprintf(fout, "%s;\n", mutstr);
          fprintf(fout2, "%s\n", StringArrayGet(&strings, 4));
        }
        
        
      }
      StringArrayDestroy(&strings);
    }
    fclose(fout);
    fclose(fout2);
  }
  StringArrayDestroy(&pdbs);
  FileReaderDestroy(&fr);
  return Success;
}


// general functions
int SIP::calculate_pearson_and_rmse(int nentry, double *exp_ddgs, double *cal_ddgs, double *pearson, double *rmse){
  if(nentry <= 1){
    printf("cannot calculate the pearson and rmse, because nentry <= 1\n");
    exit(ZeroDivisonError);
  }
  double xsum = 0.0, xavg = 0.0, ysum = 0.0, yavg = 0.0;
  for(int i = 0; i < nentry; i++){
    xsum += exp_ddgs[i];
    ysum += cal_ddgs[i];
  }
  xavg = xsum/nentry;
  yavg = ysum/nentry;
  double nominator = 0.0;
  for(int i = 0; i < nentry; i++){
    nominator += (exp_ddgs[i] - xavg)*(cal_ddgs[i] - yavg);
  }
  double denominator = 0.0;
  double sumxx = 0.0, sumyy = 0.0;
  for(int i = 0; i < nentry; i++){
    sumxx += (exp_ddgs[i] - xavg)*(exp_ddgs[i] - xavg);
    sumyy += (cal_ddgs[i] - yavg)*(cal_ddgs[i] - yavg);
  }
  denominator = sqrt(sumxx*sumyy);
  if(denominator < 1e-6) denominator = 1e-6;
  *pearson = nominator/denominator;

  double ddg2 = 0;
  for(int i = 0; i < nentry; i++){
    ddg2 += (exp_ddgs[i] - cal_ddgs[i])*(exp_ddgs[i] - cal_ddgs[i]);
  }
  *rmse = sqrt(ddg2/nentry);

  return Success;
}

int SIP::metropolis_criteria(double old_val, double new_val, double temperature){
  if(new_val - old_val <= 0){
    return 1;
  }
  else{
    double val1 = exp((old_val - new_val)/temperature);
    double val2 = rand()/(RAND_MAX+1.0);
    if(val1 > val2){
      return 1;
    }
    else{
      return 0;
    }
  }

  return 1;
}

int SIP::swap_criteria(double val1, double val2, double temp1, double temp2){
  double db = 1/temp2 - 1/temp1;
  double dv = val2 - val1;
  if(db*dv >= 0){
    return 1;
  }
  else{
    double rdm = rand()/(RAND_MAX+1.0);
    if(exp(db*dv) > rdm){
      return 1;
    }
    else{
      return 0;
    }
  }

  return 1;
}

double SIP::define_objective_function(double pearson, double rmse, SIP::Parameter *par){
  double obj = rmse;
  return obj;
}

int SIP::print_pearson_and_rmse(double pearson, double rmse){
  printf("pearson & rmse: %7.3f %7.3f\n\n", pearson, rmse);
  return Success;
}



SIP::Data::Data(){
  npdb = 0;
  nmutants = 0;
  StringArrayCreate(&mutants);
  StringArrayCreate(&pdbs);
}

SIP::Data::~Data(){
  npdb = 0;
  nmutants = 0;
  StringArrayDestroy(&mutants);
  StringArrayDestroy(&pdbs);
}

SIP::MCoptimization::MCoptimization(){
  ddgs = NULL;
  bestobj = 1e8;
}

SIP::MCoptimization::~MCoptimization(){
  ;
}

// alignfilename can be ialign.txt, ialign_bindprofx.txt
int SIP::MCoptimization::optimization_read_data_from_file1or2(SIP::Parameter *par,char* mutantfile, char* working_path, char* alignfilename){
  FileReader fr;
  FileReaderCreate(&fr, mutantfile);
  char buffer[1024];
  StringArray pdbs;
  StringArrayCreate(&pdbs);
  while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, buffer, ' ');
    int index = -1;
    StringArrayFind(&pdbs, StringArrayGet(&strings, 0), &index);
    if(index == -1){
      StringArrayAppend(&pdbs, StringArrayGet(&strings, 0));
    }
    StringArrayDestroy(&strings);
  }
  FileReaderDestroy(&fr);

  // check if interface alignment file exists, remove the pdbid if not existing
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    char ialignfile[MAX_LENGTH_ONE_LINE_IN_FILE];
    sprintf(ialignfile, "%s/%s/%s_%s", working_path, StringArrayGet(&pdbs, i), StringArrayGet(&pdbs, i), alignfilename);
    FileReader fr;
    if(FAILED(FileReaderCreate(&fr, ialignfile))){
      printf("neglect interface alignment file %s and the corresponding PDB entry\n", ialignfile);
      StringArrayRemove(&pdbs, i);
      i--;
    }
    FileReaderDestroy(&fr);
  }

  // read the interface alignment file
  expdata.npdb = StringArrayGetCount(&pdbs);
  ddgs = (SIP::DDG*)malloc(sizeof(SIP::DDG)*expdata.npdb);
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    // initialize the ddg data structure
    ddgs[i].ddg_initialize();
    strcpy(expdata.pdblist[i], StringArrayGet(&pdbs, i));
  }

  // read the interface alignment file
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    char ialignfile[MAX_LENGTH_ONE_LINE_IN_FILE];
    sprintf(ialignfile, "%s/%s/%s_%s", working_path, StringArrayGet(&pdbs, i), StringArrayGet(&pdbs, i), alignfilename);
    ddgs[i].ia_structure.interface_alignment_read_bindprofx_with_cutoff(ialignfile, par->ial_weight, par->isscore_cutoff, par->ial_count_cutoff);
  }

  // read the mutant file
  FileReaderCreate(&fr, mutantfile);
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    FileReaderSetCurrentPos(&fr, 0);
    char *pdbid = StringArrayGet(&pdbs, i);
    while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
      StringArray strings;
      StringArrayCreate(&strings);
      StringArraySplitString(&strings, buffer, ' ');
      if(strcmp(pdbid, StringArrayGet(&strings, 0)) == 0){
        // remove the data that has <= 10 alignments with the wild-type sequence included
        //if(ddgs[i].ia_from_bindprofx.alignment_num > 16){
        //  StringArrayDestroy(&strings);
        //  continue;
        //}

        // extract the experimental ddg data
        expdata.expddgs[expdata.nmutants] = atof(StringArrayGet(&strings, 4));
        ddgs[i].ms.mutant_set_append(StringArrayGet(&strings, 5));
        expdata.nmutants++;
      }
      StringArrayDestroy(&strings);
    }
  }
  StringArrayDestroy(&pdbs);
  FileReaderDestroy(&fr);
  return Success;
}


// alignmentfileaname can be ialign.txt, ialign_bindprofx.txt
int SIP::MCoptimization::optimization_read_data_from_file12(SIP::Parameter *par, char* mutantfile, char* working_path, char* structure_alignfile, char *sequence_ialignfile){
  FileReader fr;
  FileReaderCreate(&fr, mutantfile);
  char buffer[1024];
  StringArray pdbs;
  StringArrayCreate(&pdbs);
  while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, buffer, ' ');
    int index = -1;
    StringArrayFind(&pdbs, StringArrayGet(&strings, 0), &index);
    if(index == -1){
      StringArrayAppend(&pdbs, StringArrayGet(&strings, 0));
    }
    StringArrayDestroy(&strings);
  }
  FileReaderDestroy(&fr);

  // check if interface alignment file exists, remove the pdbid if not existing
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    char ialignfile[MAX_LENGTH_ONE_LINE_IN_FILE];
    sprintf(ialignfile, "%s/%s/%s_%s", working_path, StringArrayGet(&pdbs, i), StringArrayGet(&pdbs, i), structure_alignfile);
    FileReader fr;
    if(FAILED(FileReaderCreate(&fr, ialignfile))){
      printf("neglect interface alignment file %s and the corresponding PDB entry\n", ialignfile);
      StringArrayRemove(&pdbs, i);
      i--;
    }
    FileReaderDestroy(&fr);
  }

  // read the interface alignment file
  expdata.npdb = StringArrayGetCount(&pdbs);
  ddgs = (SIP::DDG*)malloc(sizeof(SIP::DDG)*expdata.npdb);
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    // initialize the ddg data structure
    ddgs[i].ddg_initialize();
    strcpy(expdata.pdblist[i], StringArrayGet(&pdbs, i));
  }

  // read the interface alignment file
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    char ialignfile[MAX_LENGTH_ONE_LINE_IN_FILE];
    sprintf(ialignfile, "%s/%s/%s_%s", working_path, StringArrayGet(&pdbs, i), StringArrayGet(&pdbs, i), structure_alignfile);
    ddgs[i].ia_structure.interface_alignment_read_bindprofx_with_cutoff(ialignfile, par->ial_weight,par->isscore_cutoff, par->ial_count_cutoff);
    sprintf(ialignfile, "%s/%s/%s_%s", working_path, StringArrayGet(&pdbs, i), StringArrayGet(&pdbs, i), sequence_ialignfile);
    ddgs[i].ia_sequence.interface_alignment_read_sip_with_cutoff(ialignfile, par->psi_weight,par->psi_linkscore_cutoff, par->seqid_high_cutoff, par->seqid_low_cutoff, par->psi_count_cutoff);
  }

  // read the mutant file
  FileReaderCreate(&fr, mutantfile);
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    FileReaderSetCurrentPos(&fr, 0);
    char *pdbid = StringArrayGet(&pdbs, i);
    while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
      StringArray strings;
      StringArrayCreate(&strings);
      StringArraySplitString(&strings, buffer, ' ');
      if(strcmp(pdbid, StringArrayGet(&strings, 0)) == 0){
        // remove the data that has <= 10 alignments with the wild-type sequence included
        //if(ddgs[i].ia.alignment_num > 0){
        //  StringArrayDestroy(&strings);
        //  continue;
        //}

        // extract the experimental ddg data
        expdata.expddgs[expdata.nmutants] = atof(StringArrayGet(&strings, 4));
        ddgs[i].ms.mutant_set_append(StringArrayGet(&strings, 5));
        expdata.nmutants++;
      }
      StringArrayDestroy(&strings);
    }
  }
  StringArrayDestroy(&pdbs);
  FileReaderDestroy(&fr);
  return Success;
}


int SIP::MCoptimization::optimization_read_data_from_file1234(SIP::Parameter* par, char* mutantfile, char* working_path, char* structure_alignfile, char *seq_alignfile, char* tmalignfile_lig, char* tmalignfile_rec){
  FileReader fr;
  FileReaderCreate(&fr, mutantfile);
  char buffer[1024];
  StringArray pdbs;
  StringArrayCreate(&pdbs);
  while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
    StringArray strings;
    StringArrayCreate(&strings);
    StringArraySplitString(&strings, buffer, ' ');
    int index = -1;
    StringArrayFind(&pdbs, StringArrayGet(&strings, 0), &index);
    if(index == -1){
      StringArrayAppend(&pdbs, StringArrayGet(&strings, 0));
    }
    StringArrayDestroy(&strings);
  }
  FileReaderDestroy(&fr);

  // check if interface alignment file exists, remove the pdbid if not existing
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    char ialignfile[MAX_LENGTH_ONE_LINE_IN_FILE];
    sprintf(ialignfile, "%s/%s/%s_%s", working_path, StringArrayGet(&pdbs, i), StringArrayGet(&pdbs, i), structure_alignfile);
    if(structure_alignfile!=NULL && strcmp(structure_alignfile,"")!=0){
      FileReader fr;
      if(FAILED(FileReaderCreate(&fr, ialignfile))){
        printf("neglect interface alignment file %s and the corresponding PDB entry\n", ialignfile);
        StringArrayRemove(&pdbs, i);
        i--;
      }
      FileReaderDestroy(&fr);
    }
  }

  // read the interface alignment file
  expdata.npdb = StringArrayGetCount(&pdbs);
  ddgs = (SIP::DDG*)malloc(sizeof(SIP::DDG)*expdata.npdb);
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    // initialize the ddg data structure
    ddgs[i].ddg_initialize();
    strcpy(expdata.pdblist[i], StringArrayGet(&pdbs, i));
  }

  // read the interface alignment file
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    int result1=Success, result2=Success;
    char ialignfile[MAX_LENGTH_ONE_LINE_IN_FILE];
    if(structure_alignfile!=NULL && strcmp(structure_alignfile,"")!=0){
      sprintf(ialignfile, "%s/%s/%s_%s", working_path, StringArrayGet(&pdbs, i), StringArrayGet(&pdbs, i), structure_alignfile);
      result1 = ddgs[i].ia_structure.interface_alignment_read_bindprofx_with_cutoff(ialignfile, par->ial_weight, par->isscore_cutoff, par->ial_count_cutoff);
      if(FAILED(result1)){
        printf("cannot find bindprofx interface alignment file\n");
        exit(IOError);
      }
    }

    if(seq_alignfile!=NULL && strcmp(seq_alignfile,"")!=0){
      sprintf(ialignfile, "%s/%s/%s_%s", working_path, StringArrayGet(&pdbs, i), StringArrayGet(&pdbs, i), seq_alignfile);
      result2 = ddgs[i].ia_sequence.interface_alignment_read_sip_with_cutoff(ialignfile,par->psi_weight, par->psi_linkscore_cutoff, par->seqid_high_cutoff, par->seqid_low_cutoff, par->psi_count_cutoff);
    }
    
  }

  // read the mutant file
  FileReaderCreate(&fr, mutantfile);
  for(int i = 0; i < StringArrayGetCount(&pdbs); i++){
    FileReaderSetCurrentPos(&fr, 0);
    char *pdbid = StringArrayGet(&pdbs, i);
    while(!FAILED(FileReaderGetNextLine(&fr, buffer))){
      StringArray strings;
      StringArrayCreate(&strings);
      StringArraySplitString(&strings, buffer, ' ');
      if(strcmp(pdbid, StringArrayGet(&strings, 0)) == 0){
        // remove the data that has <= 10 alignments with the wild-type sequence included
        //if(ddgs[i].ia.alignment_num > 0){
        //  StringArrayDestroy(&strings);
        //  continue;
        //}

        // extract the experimental ddg data
        expdata.expddgs[expdata.nmutants] = atof(StringArrayGet(&strings, 4));
        ddgs[i].ms.mutant_set_append(StringArrayGet(&strings, 5));
        expdata.nmutants++;
        StringArrayAppend(&expdata.pdbs, StringArrayGet(&strings,0));
        StringArrayAppend(&expdata.mutants, StringArrayGet(&strings,5));
      }
      StringArrayDestroy(&strings);
    }
  }
  StringArrayDestroy(&pdbs);
  FileReaderDestroy(&fr);
  return Success;
}



int SIP::MCoptimization::optimization_by_mc1or2(SIP::Parameter *par, char *mutantfile, char *working_path, char* alignfile){
  SIP::Parameter old_par, new_par;
  double old_pearson, old_rmse;
  double new_pearson, new_rmse;
  double bestobj = 1e8;

  srand((int)time(NULL));

  // read mutant and alignment data file file
  optimization_read_data_from_file1or2(par,mutantfile, working_path, alignfile);
  old_par.parameter_initialize();
  //old_par.parameter_initialize_random();
  optimization_get_pearson_and_rmse(&old_par, &old_pearson, &old_rmse);
  best_par.parameter_copy(&old_par);
  best_pearson = old_pearson;
  best_rmse = old_rmse;
  bestobj = SIP::define_objective_function(best_pearson, best_rmse, &best_par);
  // show parameters for best pearson correlation and rmse value
  printf("Initialization:\n");
  old_par.parameter_print();
  SIP::print_pearson_and_rmse(old_pearson, old_rmse);

  printf("\nIteration:\n");
  double t = MONTE_CARLO_TEMPERATURE;
  int max_iter = MONTE_CARLO_MAX_STEP;
  int naccept = 0;
  int ndecrease = 0;
  //printf("optimize lamda values, cycle %d ...\n", j);
  for(int i = 0; i < max_iter; i++){
    new_par.parameter_copy(&old_par);
    new_par.parameter_randomize_21par_model_aat(1);
    optimization_get_pearson_and_rmse(&new_par, &new_pearson, &new_rmse);
    double oldobj = SIP::define_objective_function(old_pearson, old_rmse, &old_par);
    double newobj = SIP::define_objective_function(new_pearson, new_rmse, &new_par);
    if(SIP::metropolis_criteria(oldobj, newobj, t)){
      old_par.parameter_copy(&new_par);
      old_pearson = new_pearson;
      old_rmse = new_rmse;
      naccept++;
      if(newobj-oldobj < 0) ndecrease++;
      if(newobj - bestobj < 0){
        best_par.parameter_copy(&new_par);
        best_pearson = new_pearson;
        best_rmse = new_rmse;
        bestobj = newobj;
        printf("Iter: %d\n", i);
        best_par.parameter_print();
        SIP::print_pearson_and_rmse(best_pearson, best_rmse);
        printf("objective function: %f\n", newobj);
      }
    }
  }
  printf("acceptance rate: %.3f\n", (double)naccept/max_iter);
  printf("decreasing rate: %.3f\n", (double)ndecrease/naccept);

  printf("***********************************************\n");
  printf("Final best parameters:\n");
  best_par.parameter_print();
  optimization_get_pearson_and_rmse(&best_par, &new_pearson, &new_rmse);
  SIP::print_pearson_and_rmse(best_pearson, best_rmse);
  optimization_print_ddgs();

  par->parameter_copy(&best_par);

  return Success;
}

int SIP::MCoptimization::optimization_by_mc12(SIP::Parameter *par, char *mutantfile, char *working_path, char* stru_alignfile, char* seq_alignfile){
  SIP::Parameter old_par, new_par;
  double old_pearson, old_rmse;
  double new_pearson, new_rmse;
  double bestobj = 1e8;

  srand((int)time(NULL));

  // read mutant and alignment data file file
  optimization_read_data_from_file12(par, mutantfile, working_path, stru_alignfile, seq_alignfile);

  // initialize parameters
  old_par.parameter_initialize();
  optimization_get_pearson_and_rmse12(&old_par, &old_pearson, &old_rmse);
  best_par.parameter_copy(&old_par);
  best_pearson = old_pearson;
  best_rmse = old_rmse;
  bestobj = SIP::define_objective_function(best_pearson, best_rmse, &best_par);
  // show parameters for best pearson correlation and rmse value
  printf("Initialization:\n");
  old_par.parameter_print();
  SIP::print_pearson_and_rmse(old_pearson, old_rmse);

  printf("\nIteration:\n");
  double t = MONTE_CARLO_TEMPERATURE;
  int max_iter = MONTE_CARLO_MAX_STEP;
  int naccept = 0;
  int ndecrease = 0;
  for(int i = 0; i < max_iter; i++){
    new_par.parameter_copy(&old_par);
    new_par.parameter_randomize_21par_model_aat(1);
    optimization_get_pearson_and_rmse12(&new_par, &new_pearson, &new_rmse);
    double oldobj = SIP::define_objective_function(old_pearson, old_rmse, &old_par);
    double newobj = SIP::define_objective_function(new_pearson, new_rmse, &new_par);
    if(SIP::metropolis_criteria(oldobj, newobj, t)){
      old_par.parameter_copy(&new_par);
      old_pearson = new_pearson;
      old_rmse = new_rmse;
      naccept++;
      if(newobj-oldobj < 0) ndecrease++;
      if(newobj - bestobj < 0){
        best_par.parameter_copy(&new_par);
        best_pearson = new_pearson;
        best_rmse = new_rmse;
        bestobj = newobj;
        printf("Iter: %d\n", i);
        best_par.parameter_print();
        SIP::print_pearson_and_rmse(best_pearson, best_rmse);
        printf("objective function: %f\n", newobj);
      }
    }
  }
  printf("acceptance rate: %.3f\n", (double)naccept/max_iter);
  printf("decreasing rate: %.3f\n", (double)ndecrease/naccept);

  printf("***********************************************\n");
  printf("Final best parameters:\n");
  best_par.parameter_print();
  optimization_get_pearson_and_rmse(&best_par, &new_pearson, &new_rmse);
  SIP::print_pearson_and_rmse(best_pearson, best_rmse);
  //optimization_print_ddgs();

  par->parameter_copy(&best_par);

  return Success;
}




int SIP::MCoptimization::optimization_by_sa1or2(SIP::Parameter *par, char *mutantfile, char *working_path, char* alignfile){
  SIP::Parameter old_par, new_par;
  double old_pearson, old_rmse;
  double new_pearson, new_rmse;
  double bestobj = 1e8;

  srand((int)time(NULL));

  // read mutant and alignment data file file
  optimization_read_data_from_file1or2(par,mutantfile, working_path, alignfile);

  old_par.parameter_initialize();
  //old_par.parameter_initialize_random();
  optimization_get_pearson_and_rmse(&old_par, &old_pearson, &old_rmse);
  best_par.parameter_copy(&old_par);
  best_pearson = old_pearson;
  best_rmse = old_rmse;
  bestobj = SIP::define_objective_function(best_pearson, best_rmse, &best_par);
  // show parameters for best pearson correlation and rmse value
  printf("Initialization:\n");
  old_par.parameter_print();
  SIP::print_pearson_and_rmse(old_pearson, old_rmse);

  int iter = 1;
  printf("\nIteration:\n");
  int cycle_max = SIMUL_ANNEA_CYCLE;
  int naccept = 0;
  for(int i = 0; i < SIMUL_ANNEA_STEP; i++){
    double x = (double)i/SIMUL_ANNEA_STEP;
    double n = 1.0+(double)i*SIMUL_ANNEA_CYCLE/SIMUL_ANNEA_STEP;
    double m = 1.0;
    double miu = 1.0;
    double pai = 3.1415926;
    double t = SIMUL_ANNEA_TMIN + (SIMUL_ANNEA_TMAX - SIMUL_ANNEA_TMIN)*(1-pow(x, miu))*pow(cos(n*pai*x), 2.0*m);
    new_par.parameter_copy(&old_par);
    new_par.parameter_randomize_21par_model_aat(1);
    optimization_get_pearson_and_rmse(&new_par, &new_pearson, &new_rmse);
    double oldobj = SIP::define_objective_function(old_pearson, old_rmse, &old_par);
    double newobj = SIP::define_objective_function(new_pearson, new_rmse, &new_par);
    if(SIP::metropolis_criteria(oldobj, newobj, t)){
      // accept the move and copy data back
      old_par.parameter_copy(&new_par);
      old_pearson = new_pearson;
      old_rmse = new_rmse;
      naccept++;
      if(newobj - bestobj < 0){
        best_par.parameter_copy(&new_par);
        best_pearson = new_pearson;
        best_rmse = new_rmse;
        bestobj = newobj;
        printf("Iter: %d\n", iter);
        best_par.parameter_print();
        SIP::print_pearson_and_rmse(best_pearson, best_rmse);
        printf("objective function: %f, temperature: %f\n", newobj, t);
      }
    }
    iter++;
  }

  printf("***********************************************\n");
  printf("Final best parameters:\n");
  best_par.parameter_print();
  SIP::print_pearson_and_rmse(best_pearson, best_rmse);
  //optimization_print_ddgs();

  par->parameter_copy(&best_par);

  return Success;
}



int SIP::MCoptimization::optimization_by_sa12(SIP::Parameter *par, char *mutantfile, char *working_path, char* stru_alignfile, char* seq_alignfile){
  SIP::Parameter old_par, new_par;
  double old_pearson, old_rmse;
  double new_pearson, new_rmse;
  double bestobj = 1e8;

  srand((int)time(NULL));

  // read mutant and alignment data file file
  optimization_read_data_from_file12(par, mutantfile, working_path, stru_alignfile, seq_alignfile);

  old_par.parameter_initialize();
  //old_par.parameter_initialize_random();
  optimization_get_pearson_and_rmse12(&old_par, &old_pearson, &old_rmse);
  best_par.parameter_copy(&old_par);
  best_pearson = old_pearson;
  best_rmse = old_rmse;
  bestobj = SIP::define_objective_function(best_pearson, best_rmse, &best_par);
  // show parameters for best pearson correlation and rmse value
  printf("Initialization:\n");
  old_par.parameter_print();
  SIP::print_pearson_and_rmse(old_pearson, old_rmse);

  int iter = 1;
  printf("\nIteration:\n");
  int cycle_max = SIMUL_ANNEA_CYCLE;
  int naccept = 0;
  for(int i = 0; i < SIMUL_ANNEA_STEP; i++){
    double x = (double)i/SIMUL_ANNEA_STEP;
    double n = 1.0+(double)i*SIMUL_ANNEA_CYCLE/SIMUL_ANNEA_STEP;
    double m = 1.0;
    double miu = 1.0;
    double pai = 3.1415926;
    double t = SIMUL_ANNEA_TMIN + (SIMUL_ANNEA_TMAX - SIMUL_ANNEA_TMIN)*(1-pow(x, miu))*pow(cos(n*pai*x), 2.0*m);
    //double t = SIMUL_ANNEA_TMIN + (SIMUL_ANNEA_TMAX - SIMUL_ANNEA_TMIN)*(1-pow(x, miu))*pow(sin(n*pai*x), 2.0*m);
    new_par.parameter_copy(&old_par);
    new_par.parameter_randomize_21par_model_aat(1);
    optimization_get_pearson_and_rmse12(&new_par, &new_pearson, &new_rmse);
    double oldobj = SIP::define_objective_function(old_pearson, old_rmse, &old_par);
    double newobj = SIP::define_objective_function(new_pearson, new_rmse, &new_par);
    if(SIP::metropolis_criteria(oldobj, newobj, t)){
      // accept the move and copy data back
      old_par.parameter_copy(&new_par);
      old_pearson = new_pearson;
      old_rmse = new_rmse;
      naccept++;
      if(newobj - bestobj < 0){
        best_par.parameter_copy(&new_par);
        best_pearson = new_pearson;
        best_rmse = new_rmse;
        bestobj = newobj;
        printf("Iter: %d\n", iter);
        best_par.parameter_print();
        SIP::print_pearson_and_rmse(best_pearson, best_rmse);
        printf("objective function: %f, temperature: %f\n", newobj, t);
      }
    }
    iter++;
  }

  printf("***********************************************\n");
  printf("Final best parameters:\n");
  best_par.parameter_print();
  SIP::print_pearson_and_rmse(best_pearson, best_rmse);
  //optimization_print_ddgs();

  par->parameter_copy(&best_par);

  return Success;
}

int SIP::MCoptimization::optimization_by_sa1234(SIP::Parameter *par, char *mutantfile, char *working_path, char* stru_alignfile, char* seq_alignfile, char* tmalignfile_lig, char* tmalignfile_rec, char* parameterfile){
  SIP::Parameter new_par;
  double old_pearson = 1e8, old_rmse = 1e8;
  double new_pearson = 1e8, new_rmse = 1e8;
  double bestobj = 1e8;
  BOOL DEBUG_SHOW_PARAMETERS = TRUE;

  srand((int)time(NULL));

  // read mutant and alignment data file file
  if(FAILED(par->parameter_read(parameterfile))){
    par->parameter_initialize();
  }
  optimization_read_data_from_file1234(par, mutantfile, working_path, stru_alignfile, seq_alignfile, tmalignfile_lig, tmalignfile_rec);

  optimization_get_pearson_and_rmse1234(par, &old_pearson, &old_rmse);
  best_par.parameter_copy(par);
  best_pearson = old_pearson;
  best_rmse = old_rmse;
  bestobj = SIP::define_objective_function(best_pearson, best_rmse, &best_par);
  //show parameters for best pearson correlation and rmse value
  if(DEBUG_SHOW_PARAMETERS){
    printf("Initialization:\n");
    par->parameter_print();
    SIP::print_pearson_and_rmse(old_pearson, old_rmse);
  }

  int iter = 1;
  if(DEBUG_SHOW_PARAMETERS){
    printf("Iteration:\n");
  }
  int cycle_max = SIMUL_ANNEA_CYCLE;
  int naccept = 0;
  for(int i = 0; i < SIMUL_ANNEA_STEP; i++){
    double x = (double)i/SIMUL_ANNEA_STEP;
    double n = 1.0+(double)i*SIMUL_ANNEA_CYCLE/SIMUL_ANNEA_STEP;
    double m = 1.0;
    double miu = 1.0;
    double pai = 3.1415926;
    double t = SIMUL_ANNEA_TMIN + (SIMUL_ANNEA_TMAX - SIMUL_ANNEA_TMIN)*(1-pow(x, miu))*pow(cos(n*pai*x), 2.0*m);
    //double t = SIMUL_ANNEA_TMIN + (SIMUL_ANNEA_TMAX - SIMUL_ANNEA_TMIN)*(1-pow(x, miu))*pow(sin(n*pai*x), 2.0*m);
    new_par.parameter_copy(par);
    new_par.parameter_randomize_21par_model_aat(1);
    optimization_get_pearson_and_rmse1234(&new_par, &new_pearson, &new_rmse);
    double oldobj = SIP::define_objective_function(old_pearson, old_rmse, par);
    double newobj = SIP::define_objective_function(new_pearson, new_rmse, &new_par);
    if(SIP::metropolis_criteria(oldobj, newobj, t)){
      // accept the move and copy data back
      par->parameter_copy(&new_par);
      old_pearson = new_pearson;
      old_rmse = new_rmse;
      naccept++;
      if(newobj - bestobj < 0){
        best_par.parameter_copy(&new_par);
        best_pearson = new_pearson;
        best_rmse = new_rmse;
        bestobj = newobj;
        if(DEBUG_SHOW_PARAMETERS){
          printf("Iter: %d\n", iter);
          best_par.parameter_print();
          SIP::print_pearson_and_rmse(best_pearson, best_rmse);
        }
      }
    }
    iter++;
  }

  if(DEBUG_SHOW_PARAMETERS){
    printf("Final best parameters:\n");
    best_par.parameter_print();
    printf("Final best pearson & rmse: %7.3f %7.3f\n\n", best_pearson, best_rmse);
    //optimization_print_ddgs();
  }

  par->parameter_copy(&best_par);
  return Success;
}


int SIP::MCoptimization::remc_local_movement(SIP::Parameter *par, double *pearson, double *rmse, double temperature, int step){
  SIP::Parameter new_par;
  double new_pearson, new_rmse;
  int nacc = 0;
  for(int i = 0; i < step; i++){
    new_par.parameter_copy(par);
    // randomly change pseudo counts
    new_par.parameter_randomize_61par_model(1);
    optimization_get_pearson_and_rmse(&new_par, &new_pearson, &new_rmse);
    double oldobj = SIP::define_objective_function(*pearson, *rmse, par);
    double newobj = SIP::define_objective_function(new_pearson, new_rmse, &new_par);
    if(SIP::metropolis_criteria(oldobj, newobj, temperature)){
      // accept the move and copy data back
      par->parameter_copy(&new_par);
      *pearson = new_pearson;
      *rmse = new_rmse;
      nacc++;
      //printf("objective function: %f, delta_e: %f, temperature: %f\n", newobj, newobj - oldobj, temperature);
      if(newobj - bestobj < 0){
        best_par.parameter_copy(&new_par);
        best_pearson = new_pearson;
        best_rmse = new_rmse;
        bestobj = newobj;
        best_par.parameter_print();
        SIP::print_pearson_and_rmse(best_pearson, best_rmse);
        printf("objective function: %f, temperature: %f\n", newobj, temperature);
      }
    }
  }
  printf("acceptance rate: %f @ temperature: %f\n", (double)nacc/step, temperature);

  return Success;
}



int SIP::MCoptimization::optimization_by_remc1or2(SIP::Parameter *par, char *mutantfile, char *working_path, char* alignfile){
  SIP::Parameter old_par, new_par;

  srand((int)time(NULL));

  // read mutant and alignment data file file
  optimization_read_data_from_file1or2(par,mutantfile, working_path, alignfile);

  int NREPLICA = REPLICA_EXCH_REPLICAS;
  double tmin = REPLICA_EXCH_TMIN, tmax = REPLICA_EXCH_TMAX;
  SIP::Parameter *parameters = (Parameter*)malloc(sizeof(Parameter)*NREPLICA);
  double *pearsons = (double*)malloc(sizeof(double)*NREPLICA);
  double *rmses = (double*)malloc(sizeof(double)*NREPLICA);
  double *temperatures = (double*)malloc(sizeof(double)*NREPLICA);
  for(int i = 0; i < NREPLICA; i++){
    parameters[i].parameter_initialize();
    optimization_get_pearson_and_rmse(&parameters[i], &pearsons[i], &rmses[i]);
    // lower index means lower temperature
    temperatures[i] = tmin*pow((tmax/tmin),(double)i/(NREPLICA-1));
  }

  int nglobal = REPLICA_EXCH_NSWAP;
  int nlocal = REPLICA_EXCH_NMETROPLIS;
  for(int i = 1; i < nglobal; i++){
    // do local movement for replicas (metropolis)
    for(int j = 0; j < NREPLICA; j++){
      remc_local_movement(&parameters[j], &pearsons[j], &rmses[j], temperatures[j], nlocal);
    }

    // global swap between replicas
    // odd index, try to swap 0<->1, 2<->3, ... 38<->39
    if(i%2 == 1){
      for(int j = 0; j < NREPLICA; j +=2){
        double val1 = SIP::define_objective_function(pearsons[j], rmses[j], &parameters[j]);
        double val2 = SIP::define_objective_function(pearsons[j+1], rmses[j+1], &parameters[j+1]);
        double dval = val2 - val1;
        double dtem = temperatures[j+1] - temperatures[j];
        if(SIP::swap_criteria(val1, val2, temperatures[j], temperatures[j+1])){
          printf("successful swap between temperature %f and %f, pearson, rmse: %f %f <=> %f %f\n",
            temperatures[j], temperatures[j+1], pearsons[j], rmses[j], pearsons[j+1], rmses[j+1]);
          SIP::Parameter temp;
          temp.parameter_copy(&parameters[j]);
          parameters[j].parameter_copy(&parameters[j+1]);
          parameters[j+1].parameter_copy(&temp);
          double p = pearsons[j];
          pearsons[j] = pearsons[j+1];
          pearsons[j+1] = p;
          double r = rmses[j];
          rmses[j] = rmses[j+1];
          rmses[j+1] = rmses[j];
        }
      }
    }
    // even index, try to swap 2<->3, 4<->5, ... 37<->38
    else{
      for(int j = 1; j < NREPLICA-1; j +=2){
        double val1 = SIP::define_objective_function(pearsons[j], rmses[j], &parameters[j]);
        double val2 = SIP::define_objective_function(pearsons[j+1], rmses[j+1], &parameters[j+1]);
        double dval = val2 - val1;
        double dtem = temperatures[j+1] - temperatures[j];
        if(SIP::swap_criteria(val1, val2, temperatures[j], temperatures[j+1])){
          printf("successful swap between temperature %f and %f, pearson, rmse: %f %f <=>  %f %f\n",
            temperatures[j], temperatures[j+1], pearsons[j], rmses[j], pearsons[j+1], rmses[j+1]);
          Parameter temp;
          temp.parameter_copy(&parameters[j]);
          parameters[j].parameter_copy(&parameters[j+1]);
          parameters[j+1].parameter_copy(&temp);
          double p = pearsons[j];
          pearsons[j] = pearsons[j+1];
          pearsons[j+1] = p;
          double r = rmses[j];
          rmses[j] = rmses[j+1];
          rmses[j+1] = rmses[j];
        }
      }
    }
  }

  printf("***********************************************\n");
  printf("Final best parameters:\n");
  best_par.parameter_print();
  SIP::print_pearson_and_rmse(best_pearson, best_rmse);
  //optimization_print_ddgs();

  par->parameter_copy(&best_par);

  return Success;
}





int SIP::MCoptimization::optimization_print_ddgs(){
  int expddgindex = 0;
  for(int i = 0; i < expdata.npdb; i++){
    SIP::DDG *ddg = &ddgs[i];
    for(int j = 0; j < ddg->ms.nset; j++){
      SIP::MultipleMutant *mm = &ddg->ms.mutset[j];
      printf("%s %7.3f %7.3f ", expdata.pdblist[i], expdata.expddgs[expddgindex], mm->mmddg);
      expddgindex++;
      for(int k = 0; k < mm->nmut; k++){
        SIP::SingleMutant *sm = &mm->mutants[k];
        sm->single_mutant_print();
        if(k != mm->nmut - 1){
          printf(",");
        }
        else{
          printf("\n");
        }
      }
    }
  }
  return Success;
}



int SIP::MCoptimization::optimization_append_calcddg(SIP::DDG *ddg){
  for(int i = 0; i < ddg->ms.nset; i++){
    calddgs[nentry+i] = ddg->ms.mutset[i].mmddg;
  }
  nentry += ddg->ms.nset;
  return Success;
}



int SIP::MCoptimization::optimization_get_pearson_and_rmse(SIP::Parameter *par, double* pearson, double *rmse){
  nentry = 0;
  for(int i = 0; i < expdata.npdb; i++){
    SIP::DDG *ddg = &ddgs[i];
    ddg->calc_ddg_for_mutation_set(par);
    optimization_append_calcddg(ddg);
  }
  SIP::calculate_pearson_and_rmse(expdata.nmutants, expdata.expddgs, calddgs, pearson, rmse);
  return Success;
}

int SIP::MCoptimization::optimization_get_pearson_and_rmse12(SIP::Parameter *par, double* pearson, double *rmse){
  nentry = 0;
  for(int i = 0; i < expdata.npdb; i++){
    SIP::DDG *ddg = &ddgs[i];
    ddg->calc_ddg_for_mutation_set12(par);
    optimization_append_calcddg(ddg);
  }
  SIP::calculate_pearson_and_rmse(expdata.nmutants, expdata.expddgs, calddgs, pearson, rmse);
  return Success;
}

int SIP::MCoptimization::optimization_get_pearson_and_rmse1234(SIP::Parameter *par, double* pearson, double *rmse){
  nentry = 0;
  for(int i = 0; i < expdata.npdb; i++){
    SIP::DDG *ddg = &ddgs[i];
    ddg->calc_ddg_for_mutation_set1234(par);
    optimization_append_calcddg(ddg);
  }
  SIP::calculate_pearson_and_rmse(expdata.nmutants, expdata.expddgs, calddgs, pearson, rmse);
  return Success;
}


int SIP::MCoptimization::optimization_get_cvddg1234(SIP::Parameter *par, double* exp, double* pre, int* ndata, StringArray* pdbs, StringArray* mutants){
  nentry = 0;
  for(int i = 0; i < expdata.npdb; i++){
    SIP::DDG *ddg = &ddgs[i];
    ddg->calc_ddg_for_mutation_set1234(par);
    optimization_append_calcddg(ddg);
  }
  //record the ddg of sub cv dataset into bigger array
  for(int i=0; i<expdata.nmutants; i++){
    exp[*ndata+i]=expdata.expddgs[i];
    pre[*ndata+i]=calddgs[i];
    StringArrayAppend(pdbs, StringArrayGet(&expdata.pdbs, i));
    StringArrayAppend(mutants, StringArrayGet(&expdata.mutants, i));
  }
  *ndata+=expdata.nmutants;

  return Success;
}


int SIP::MCoptimization::optimization_release_memory(){
  if(ddgs != NULL){
    for(int i = 0; i < expdata.npdb; i++){
      ddgs[i].ms.mutant_set_destroy();
    }
    free(ddgs);
    ddgs = NULL;
  }
  nentry = 0;
  expdata.nmutants = 0;
  expdata.npdb = 0;

  return Success;
}



int SIP::MCoptimization::optimization_show_interface_alignment_gap_count(){
  for(int i = 0; i < expdata.npdb; i++){
    SIP::DDG *ddg = &ddgs[i];
    SIP::MutantSet *ms = &ddgs->ms;
    for(int j = 0; j < ms->nset; j++){
      SIP::MultipleMutant *mm = &ms->mutset[j];
      for(int k = 0; k < mm->nmut; k++){
        SIP::SingleMutant *sm = &mm->mutants[k];
        int aln_index = ddg->find_index_for_single_mutation(sm);
        double ratio = ddg->ia_structure.interface_alignment_get_gap_ratio(sm->chain_name, aln_index);
        if(ratio > 0.0){
          printf("single mutant ");
            sm->single_mutant_print();
            printf(" gap/naln = %.3f in the alignment file %s\n", ratio, ddg->ia_structure.alignment_filename);
        }
      }
    }
    
  }
  return Success;
}

