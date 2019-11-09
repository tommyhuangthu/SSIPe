//###################################################################################################
//# build_complex_alignment
//#
//# This is a perl script used to build the protein complex MSA from two monomer MSAs via the STRING 
//# protein links database.
//# Input:  <-- protein.links.txt
//#         <-- pdb_lig.msa
//#         <-- pdb_rec.msa
//# Output: --> protein_cpx_alignment.txt
//#
//# Xiaoqiang Huang
//# xiaoqiah@umich.edu
//###################################################################################################

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_SPE_PAIR  1000
#define MAX_PRO_PAIR  1000

#define MAX_DEF_LENG  100
#define MAX_LIN_LENG  10000
#define MAX_SEQ_LENG  10000
#define MAX_NAM_LENG  100


char *linkfile = NULL;
char *msa1file = NULL;
char *msa2file = NULL;
char *cpxoutfile = NULL;


typedef struct _Comm{
  int spe;
  int count;
  char content[MAX_PRO_PAIR][MAX_DEF_LENG];
}Comm;


int make_new_files(){
  char cmd[1000];
  sprintf(cmd, "cp -rf %s msa1\n", msa1file);
  system(cmd);
  sprintf(cmd, "cp -rf %s msa2\n", msa2file);
  system(cmd);
  return 0;
}


int sort_msa(){
  printf("Sorting ligand and receptor MSA files ... \n");
  system("sort -k2 msa1 > sorted_msa1\n");
  system("sort -k2 msa2 > sorted_msa2\n");
  printf("Done.\n");
  return 0;
}

// find common organisms of the two MSAs and build a pseudolink 
// between two proteins for each organism
int find_common_species_and_build_links(){
  Comm* comms=(Comm*)malloc(sizeof(Comm)*MAX_SPE_PAIR);
  printf("Find common species of two MSAs and build links between proteins ... \n");
  // set a maximum common specie number: 1000
  for(int i = 0; i < MAX_SPE_PAIR; i++){
    comms[i].count = 0;
    // for each common specie, set the maximum protein pair: 1000
    for(int j = 0; j < MAX_PRO_PAIR; j++){
      strcpy(comms[i].content[j], "");
    }
  }
  int comm_spe_count = 0;
  int prespe = -10000;
  FILE *f1, *f2;
  f1 = fopen("sorted_msa1", "r");
  f2 = fopen("sorted_msa2", "r");
  char buffer1[MAX_LIN_LENG],buffer2[MAX_LIN_LENG];
  while(fgets(buffer1, MAX_LIN_LENG, f1)){
    char msa1[MAX_SEQ_LENG], def1[MAX_DEF_LENG];
    int spe1;
    sscanf(buffer1, "%s  %d.%s", msa1, &spe1, def1);
    //printf("from ligmsa: %s %d.%s\n", msa1, spe1, def1);
    while(fgets(buffer2, MAX_LIN_LENG, f2)){
      char msa2[MAX_SEQ_LENG], def2[MAX_DEF_LENG];
      int spe2;
      sscanf(buffer2, "%s  %d.%s", msa2, &spe2, def2);
      //printf("from recmsa: %s %d.%s\n", msa2, spe2, def2);
      // a pair with new species
      if(spe1 == spe2 && spe1 != prespe){
        comms[comm_spe_count].count = 0;
        prespe = spe1;
        comms[comm_spe_count].spe = spe1;
	char temp[MAX_DEF_LENG];
	sprintf(temp, "%d.%s %d.%s", spe1, def1, spe1, def2);
	//printf("%s\n", temp);
	strcpy(comms[comm_spe_count].content[comms[comm_spe_count].count], temp);
	comms[comm_spe_count].count++;
	comm_spe_count++;
	if(comm_spe_count == MAX_SPE_PAIR){
	  break;
	}
      }
      // a pair with existed species
      else if(spe1 == spe2 && spe1 == prespe){
	char temp[MAX_DEF_LENG];
	sprintf(temp, "%d.%s %d.%s", spe1, def1, spe1, def2);
        //printf("%s\n", temp);
	strcpy(comms[comm_spe_count-1].content[comms[comm_spe_count-1].count], temp);
        comms[comm_spe_count-1].count++;
      }
    }
    fseek(f2, 0, SEEK_SET);
    if(comm_spe_count == MAX_SPE_PAIR){
      break;
    }
  }
  fclose(f1); fclose(f2);
  // show the comm species
  for(int index = 0; index < comm_spe_count; index++){
    printf("The following are common protein pairs for species: %d\n", comms[index].spe);
    for(int index2 = 0; index2 < comms[index].count; index2++){
      printf("%s\n", comms[index].content[index2]);
    }
  }
  printf("Done.\n");

  // find the real link and the link score if possible via the species based link files
  // this is the most time-consuming part of this program
  printf ("Build complex links from the common species ... \n");
  FILE *fout = fopen("cpx_link.txt", "w");
  for(int index = 0; index < comm_spe_count; index++){
    char newlinkfile[MAX_NAM_LENG];
    sprintf(newlinkfile, "%s.%d", linkfile, comms[index].spe);
    printf("searching specie pair %d from string links %s ...\n", comms[index].spe, newlinkfile);
    FILE *fin = fopen(newlinkfile, "r");
    char buffer[MAX_LIN_LENG];
    int startindex = 0;
    while(fgets(buffer, MAX_LIN_LENG, fin)){
      for(int index2 = startindex; index2 < comms[index].count; index2++){
        if(strstr(buffer, comms[index].content[index2]) != 0){
          fprintf(fout, "%s", buffer);
          // next time will start searching from the next protein pair, not from the beginning
          startindex = index2+1;
          break;
        }
      }
    }
    fclose(fin);
  }
  fclose(fout);
  printf("Done.\n");

  return 0;
}

// step 4: sort the real link result with sore ranked from high to low
// select the link with the highest score when there're several
// links between one protein with several different proteins
int sort_links(){
  printf("Sort complex links and select only the top 1 complex link for each protein ... \n");
  system("sort -k1,1 -k3nr,3  cpx_link.txt > cpx_link_sort13.txt");
  char initDef[MAX_DEF_LENG] = "XXX";
  char line[MAX_LIN_LENG];
  FILE *f1 = fopen("cpx_link_sort13.txt", "r");
  FILE *f2 = fopen("cpx_link_sort13_top1.txt", "w");
  while(fgets(line, MAX_LIN_LENG, f1)){
    char def1[MAX_DEF_LENG];
    sscanf(line, "%s ", def1);
    if(strcmp(def1, initDef) != 0){
      fprintf(f2, "%s", line);
      strcpy(initDef, def1);
    }
  }
  fclose(f1); fclose(f2);

  system("sort -k2,2 -k3nr,3 cpx_link_sort13_top1.txt > cpx_link_sort23.txt");
  strcpy(initDef, "XXX");
  f1 = fopen("cpx_link_sort23.txt", "r");
  f2 = fopen("top_cpx_link.txt", "w");
  while(fgets(line, MAX_LIN_LENG, f1)){
    char def1[MAX_DEF_LENG], def2[MAX_DEF_LENG];
    sscanf(line, "%s %s ", def1, def2);
    if(strcmp(def2, initDef) != 0){
      fprintf(f2, "%s", line);
      strcpy(initDef, def2);
    }
  }
  fclose(f1); fclose(f2);
  printf("Done.\n");
  return 0;
}


// build the MSA file for the whole protein complex 
// the default medium-high cutoff link score is 400 [400, 700]
// the high link score is [700, 900]
// the ultra high link score cutoff is [900, 999]
int build_cpx_msa(){
  printf("Build complex MSA for protein complex ... \n");
  int CUTOFF = 1;
  FILE *f1 = fopen("top_cpx_link.txt", "r");
  FILE *f2 = fopen(cpxoutfile, "w");
  char buffer[MAX_LIN_LENG];
  while(fgets(buffer, MAX_LIN_LENG, f1)){
    char def1[MAX_DEF_LENG], def2[MAX_DEF_LENG];
    int score;
    sscanf(buffer, "%s %s %d", def1, def2, &score);
    char msa1[MAX_LIN_LENG], msa2[MAX_LIN_LENG];

    FILE *f3 = fopen("msa1", "r");
    char line[MAX_LIN_LENG];
    while(fgets(line, MAX_LIN_LENG, f3)){
      char msa[MAX_LIN_LENG], def[MAX_DEF_LENG];
      sscanf(line, "%s  %s", msa, def);
      if(strcmp(def, def1) == 0){
        strcpy(msa1, msa);
      }
    }
    fclose(f3);

    f3 = fopen("msa2", "r");
    while(fgets(line, MAX_LIN_LENG, f3)){
      char msa[MAX_LIN_LENG], def[MAX_DEF_LENG];
      sscanf(line, "%s  %s", msa, def);
      if(strcmp(def, def2) == 0){
        strcpy(msa2, msa);
      }
    }
    fclose(f3);

    if(score > CUTOFF){
      fprintf(f2, "%s %s %.3f\n", msa1, msa2, (double)score/1000.0);
    }
  }
  fclose(f1); fclose(f2);
  printf("Done.\n");
  return 0;
}

int delete_files(){
  system("rm msa1 msa2 sorted_msa1 sorted_msa2");
  system("rm -rf cpx_link_sort13.txt");
  system("rm -rf cpx_link_sort13_top1.txt");
  system("rm -rf cpx_link_sort23.txt");
  system("rm -rf cpx_link.txt");
  system("rm -rf top_cpx_link.txt");
  return 0;
}

int main(int argc, char **argv){
  if(argc != 9){
    printf("------------------------------------ERROR-----------------------------------------------------------\n");
    printf("This program is used to build sequence alignment for a dimeric protein, the values in [] are required, ");
    printf("and the ligand and receptor msa file should be put in the current directory\n");
    printf("Usage: build_cpx_alignment -link [STRING_protein_links] -ligmsa [ligand.msa] -recmsa [receptor.msa] -out [protein_cpx_alignment.txt]\n");
    exit(-1);
  }
  
  linkfile = argv[2];
  msa1file = argv[4];
  msa2file = argv[6];
  cpxoutfile = argv[8];


  make_new_files();
  sort_msa();
  find_common_species_and_build_links();
  sort_links();
  build_cpx_msa();
  delete_files();

  return 0;
}

