
/*
Copyright (C) 2008
Mu Gao
Georgia Institute of Technology
*/

#include "structures.h"
#include "sidechain.h"

#define  MAX_NUM_CHAIN        100
#define  MAX_NUM_CONTACT      10000000
#define  MAX_NUM_RES_CONT     500

/******  contact data structure  ******/
struct Contact{
  int    r_resid;
  int    l_resid;
  int    r_atomid;
  int    l_atomid;
  float  distance;
};




/*********************************************************************************/
/*--------------------------->    Subroutine declarations -----------------------*/
/*********************************************************************************/
long calc_contact( struct Structure *structure, char *re_chain, char *li_chain, float cutoff2, struct Contact *contact, int flag );
void print_usage();


void calc_sec_struct( struct Structure *structure );
float dist( double r1[], double r2[] );
void atomcoor_to_vector ( struct Structure myStructure, char *atom_name, double **vector );


struct Structure mk_cg_struct( struct Structure aa_structure, int rebuilt_flag );
void get_scom( int* res_type, int length, float* cax, float* cay, float* caz, float* scmx, float* scmy, float* scmz );

char** read_residue_file( char *residue_file_name, int *num_sel_residue );
struct Structure sel_struct( struct Structure aa_structure, char** sel_residue, int );







/*********************************************************************************/
/*--------------------------->    Main Procedure   <-----------------------------*/
/*********************************************************************************/

int main( int argc , char *argv[] ) {


  int i,j,ii,jj;
  long k;

  /* File stuff */
  char	*complex_file_name;
  char  *contact_file_name;
  char  *interface_file_name;
  char  *residue_file_name;

  FILE  *output_file;

  /* Structures */
  struct Structure   complex_structure,  interface;
  char   re_chain[MAX_NUM_CHAIN], li_chain[MAX_NUM_CHAIN];
  int    re_chain_len[MAX_NUM_CHAIN], li_chain_len[MAX_NUM_CHAIN];
  int    num_re_chain = 0, num_li_chain = 0;
  int    r_res, l_res, r_atom, l_atom, residue, atom;
  int    *is_int, num_int_res = 0, num_sel_residue;


  /* contact stuff */
  struct Contact *contact;
  long   num_contact = 0;
  int    *contres_num, **contres_lst, tot_num_res_cont;
  float  cont_dist2_cutoff = 20.25;    /* default contact distance 4.5 A squared */

  char **selected_residue;

  /* option flags */
  int sec_flag     = 0;   /* calculate secondary structure if flag = 1 */
  int ca_flag      = 0;   /* output only ca coordiates */
  int scom_flag    = 0;   /* calculate side-chain center of mass contacts instead of all-atom contacts */
  int rebuilt_flag = 1;   /* rebuilt side-chain center of mass from Calpha only structures */
  int dist_flag    = 1;   /* use residue specific cutoffs instead of a constant cutoff, only effective when scom_flag is 1 */

  int write_allres_flag = 0;  /* write both interfacial and non-interfacial residues */

  /*------------------>  End varaiable declaration <-------------------*/



  /*-----------------> Memory allocations and initialization <--------------------------*/

  if( ( ( complex_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL )  ||
      ( ( contact_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL )  ||
      ( ( residue_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL )  ||
      ( ( interface_file_name  = ( char * ) malloc ( 500 * sizeof( char ) ) ) == NULL) ) {
    GENERAL_MEMORY_PROBLEM
  }

  /* initialize the contact list */
  if( ( contact = malloc( MAX_NUM_CONTACT * sizeof( struct Contact ) ) ) == NULL ) { GENERAL_MEMORY_PROBLEM }



  /*------------------> defaults for command line options <---------------------------*/
  strcpy( complex_file_name,   ""   );
  strcpy( contact_file_name,   ""   );
  strcpy( interface_file_name, ""   );
  strcpy( residue_file_name,   " "  );


  /*------------------> read command line options <----------------------------------*/
  for( i = 1 ; i < argc ; i ++ ) {

    if( strcmp( argv[i] , "-s" ) == 0 ) {
      i++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
        print_usage();
        exit( EXIT_FAILURE ) ;
      }
      strcpy( complex_file_name , argv[i] ) ;
    }
    else if( strcmp( argv[i] , "-l" ) == 0 ) {
      i++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
        print_usage();
	exit( EXIT_FAILURE ) ;
      }
      strcpy( li_chain, argv[i] );
      num_li_chain += strlen( argv[i] );
      while( ( (i+1) < argc ) && ( strncmp( argv[i+1] , "-" , 1 ) != 0 ) ) {
	if( num_li_chain == ( MAX_NUM_CHAIN - 1 ) ) {
	  printf( "Error: Too many chains\nDying\n" ) ;
	  exit( EXIT_FAILURE ) ;
	}
	i++;
	num_li_chain += strlen( argv[i] );
	strcat( li_chain, argv[i] ) ;
      }
    }
    else if( strcmp( argv[i] , "-r" ) == 0 ) {
      i++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
	print_usage();
	exit( EXIT_FAILURE ) ;
      }
      strcpy( re_chain, argv[i] );
      num_re_chain += strlen( argv[i] );
      while( ( (i+1) < argc ) && ( strncmp( argv[i+1] , "-" , 1 ) != 0 ) ) {
	if( num_li_chain == ( MAX_NUM_CHAIN - 1 ) ) {
	  printf( "Error: Too many chains\nDying\n" ) ;
	  exit( EXIT_FAILURE ) ;
	}
	i++;
	num_re_chain += strlen( argv[i] );
	strcat( re_chain, argv[i] ) ;
      }
    }
    else if( strcmp( argv[i] , "-c" ) == 0 ) {
      i++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
	print_usage();
	exit( EXIT_FAILURE ) ;
      }
      strcpy( contact_file_name , argv[i] );
    }
    else if( strcmp( argv[i] , "-res" ) == 0 ) {
      i++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
	print_usage();
	exit( EXIT_FAILURE ) ;
      }
      strcpy( residue_file_name , argv[i] );
    }
    else if( strcmp( argv[i] , "-i" ) == 0 ) {
      i++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
	print_usage();
	exit( EXIT_FAILURE ) ;
      }
      strcpy( interface_file_name , argv[i] );
    }
    else if( strcmp( argv[i] , "-e" ) == 0 ) {
      sec_flag = 1;
    }
    else if( strcmp( argv[i] , "-ca" ) == 0 ) {
      ca_flag = 1;
    }
    else if( strcmp( argv[i] , "-scom" ) == 0 ) {
      scom_flag = 1;
    }
    else if( strcmp( argv[i] , "-norebuilt" ) == 0 ) {
      rebuilt_flag = 0;
    }
    else if( strcmp( argv[i] , "-wa" ) == 0 ) {
      write_allres_flag = 1;
    }
    else if( strcmp( argv[i] , "-d" ) == 0 ) {
      i++ ;
      if( ( i == argc ) || ( strncmp( argv[i] , "-" , 1 ) == 0 ) ) {
	print_usage();
	exit( EXIT_FAILURE ) ;
      }
      cont_dist2_cutoff = atof( argv[i] );
      cont_dist2_cutoff = cont_dist2_cutoff * cont_dist2_cutoff;
      dist_flag = 0;
      if ( cont_dist2_cutoff < 0.01 ) {
	printf( "Error: contact distance cutoff too small\n");
	exit( EXIT_FAILURE );
      }
    }
    else if( strcmp( argv[i] , "-h" ) == 0 ) {
      print_usage();
      exit( EXIT_FAILURE ) ;
    }
    else {
      print_usage();
      exit( EXIT_FAILURE ) ;
    }


  }

  if( scom_flag == 0 ) { dist_flag = 0; }

  if( num_li_chain == 0 || num_re_chain == 0 ) {
    print_usage();
    exit( EXIT_FAILURE );
  }

  /*---------------> end of reading command line options <----------------------*/



  /********  Read structures from pdb files **********/
  complex_structure = read_pdb_to_structure( complex_file_name );
  printf( "Complex size: %d residues\n", complex_structure.length );


  /******* Read a list of selected residues, only output residues on this list *******/
  selected_residue = NULL;
  if( strcmp( residue_file_name, " " ) != 0 ) {
    selected_residue = read_residue_file( residue_file_name, &num_sel_residue );
  }




  /*---------------> secondary structure <------------------*/
  if(sec_flag) {
    calc_sec_struct( &complex_structure );
  }


  /*---------------> coarse-grained structure with only Calpha and sidechain center of mass <------------------*/
  if( scom_flag == 1 ) {
    complex_structure = mk_cg_struct( complex_structure, rebuilt_flag );
  }

  if( selected_residue != NULL ) {
    complex_structure = sel_struct( complex_structure, selected_residue, num_sel_residue ); 

    for(i=0; i<num_sel_residue; i++) {
      free( selected_residue[i] );
    }
    free( selected_residue );
  }


  /*---------------> get chain length  <------------------*/
  for( i = 0; i < num_re_chain; i++ ) {
    re_chain_len[i] = 0;
    for( j = 1; j <= complex_structure.length; j++ ) {
      if( strchr( complex_structure.Residue[j].chainID, re_chain[i] ) != NULL ) {
	re_chain_len[i]++;
      }
    }
    printf( "Receptor chain %c length %d\n", re_chain[i], re_chain_len[i] );
  }
  for( i = 0; i < num_li_chain; i++ ) {
    li_chain_len[i] = 0;
    for( j = 1; j <= complex_structure.length; j++ ) {
      if( strchr( complex_structure.Residue[j].chainID, li_chain[i] ) != NULL ) {
	li_chain_len[i]++;
      }
    }
    printf( "Ligand chain %c length %d\n", li_chain[i], li_chain_len[i] );
  }





  /* Line buffering only - so user can see whats going on a bit better */
  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;
  printf( "Calculating contacts ...\n" ) ;

  /******** Calculate atomic contacts between receptor and ligand  *******/
  num_contact = calc_contact ( &complex_structure, re_chain, li_chain, cont_dist2_cutoff, contact, dist_flag );






  /*---------------> extract interfact residues from the contacts <------------------*/
  if ( (     is_int = malloc ( (complex_structure.length + 1) * sizeof(int) )) == NULL ) { GENERAL_MEMORY_PROBLEM }
  if ( (contres_num = malloc ( (complex_structure.length + 1) * sizeof(int) )) == NULL ) { GENERAL_MEMORY_PROBLEM }
  if ( (contres_lst = malloc ( (complex_structure.length + 1) * sizeof(int *) )) == NULL ) { GENERAL_MEMORY_PROBLEM }
  for( i = 0; i <= complex_structure.length; i++ ) {
    if ( (contres_lst[i] = malloc ( MAX_NUM_RES_CONT * sizeof(int) )) == NULL ) { GENERAL_MEMORY_PROBLEM; }
  }


  /* interfacial residues: is_int = 1, non-interfacial residues: is_int = 0 */
  for( i = 0; i <= complex_structure.length; i++ ) { is_int[i] = 0; }
  for( k = 0; k < num_contact; k++ ) {
    r_res  = contact[k].r_resid;
    l_res  = contact[k].l_resid;
    is_int[r_res] = 1;
    is_int[l_res] = 1;
  }

  /* get the total number of interfacial residues */
  for( i = 0; i <= complex_structure.length; i++ ) { 
    if ( is_int[i] == 1 ) num_int_res++;
  }

  /***** get the list of residue-residue contact for each interface residue ****/
  for( i = 0; i <= complex_structure.length; i++ )
    for( j = 0; j < MAX_NUM_RES_CONT; j++ ) contres_lst[i][j] = -1;

  for( k = 0; k < num_contact; k++ ) {
    r_res  = contact[k].r_resid;
    l_res  = contact[k].l_resid;
    for( j = 0; j < MAX_NUM_RES_CONT; j++ ) {
      if( contres_lst[r_res][j] == l_res ) break;
      if( contres_lst[r_res][j] == -1 ) {
	contres_lst[r_res][j] = l_res;
	break;
      }
    }
    for( j = 0; j < MAX_NUM_RES_CONT; j++ ) {
      if( contres_lst[l_res][j] == r_res ) break;
      if( contres_lst[l_res][j] == -1 ) {
	contres_lst[l_res][j] = r_res;
	break;
      }
    }
  }

  tot_num_res_cont = 0;
  for( i = 0; i <= complex_structure.length; i++ ) { 
    contres_num[i] = 0;
    for( j = 0; j < MAX_NUM_RES_CONT ; j++ ) {
      if( contres_lst[i][j] > -1 ) {
	contres_num[i]++;
	tot_num_res_cont++;
      }
    }
    if( contres_num[i] == MAX_NUM_RES_CONT ) {
      printf( "Warning: Residue %d has resdiue contacts more than the maximum (%d) allowed\n", i, MAX_NUM_RES_CONT ); 
    }
  }
  tot_num_res_cont = tot_num_res_cont / 2;   /* each residue-residue contact counted exactly twice */
  printf( "Found %ld atomic contacts, %d residue-residue contacts\n", num_contact, tot_num_res_cont );


  /**** create a structure contains only interfacial residues ****/
  strcpy( interface.ident, "interface" );

  if( ( interface.Residue = ( struct Amino_Acid * ) malloc ( ( complex_structure.length  + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

  i = 0;
  for( residue = 1; residue <= complex_structure.length; residue++ ) {
    if( is_int[residue] == 0 && write_allres_flag == 0 ) continue;

    if( strchr( re_chain, complex_structure.Residue[residue].chainID[0] ) == NULL &&
	strchr( li_chain, complex_structure.Residue[residue].chainID[0] ) == NULL ) continue;  /* must be selected chains */


    i++;
    is_int[residue]=i;  /* index of the residue in the pdb file of the interface */
    /*printf("chain %s residue %s\n", complex_structure.Residue[residue].chainID, complex_structure.Residue[residue].res_seq_plus_iCode ); */
    interface.Residue[i] = complex_structure.Residue[residue] ;
    if( ( interface.Residue[i].Atom = ( struct Atom * ) malloc ( ( complex_structure.Residue[residue].size + 1 ) * sizeof_Atom ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }

    for( atom = 1 ; atom <= complex_structure.Residue[residue].size ; atom++ ) {
      interface.Residue[i].Atom[atom] = complex_structure.Residue[residue].Atom[atom] ;
    }
  }
  interface.length = i;

  if (interface_file_name[0]){
    if(ca_flag) {
      write_structure_to_pdbca( interface, interface_file_name );
    }
    else {
      write_structure_to_pdb( interface, interface_file_name );
    }
  }

  /*  write_structure_to_pdb( complex_structure, "test.pdb" ); */
  printf("%s %d interface residues extracted\n", complex_structure.ident, num_int_res); 

  /*---------------> end of extract interfact residues <------------------*/







  /*------------------------->  Open output data file <-----------------------------*/
  if( contact_file_name[0] ){
    if(( output_file = fopen( contact_file_name , "w" ) ) == NULL ) {
      printf( "Could not open %s for writing.\nDying\n\n" , contact_file_name ) ;
      exit( EXIT_FAILURE );
    }
  }
  else{
      output_file=stdout;
  }
  /* Line buffering only - so user can see whats going on a bit better */
  setvbuf( output_file , (char *)NULL , _IOLBF , 0 ) ;


  /* Write out command line options */
  fprintf( output_file, "Atomic and residue contact lists\n\n" );
  fprintf( output_file, "Input PDB file                :: %s\n", complex_file_name );
  fprintf( output_file, "Receptor Chain(s)             ::");
  for( i = 0; i < num_re_chain; i++ ) fprintf( output_file, " %c", re_chain[i] );
  fprintf( output_file, "\nLigand Chain(s)               ::");
  for( i = 0; i < num_li_chain; i++ ) fprintf( output_file, " %c", li_chain[i] );
  fprintf( output_file, "\nContact output                :: %s\n", contact_file_name );

  if( dist_flag == 0 ) {
    fprintf( output_file, "Contact distance cutoff       :: %-6.2f A\n",   sqrt(cont_dist2_cutoff) );
  }
  else {
    fprintf( output_file, "Contact distance cutoff       :: Residue specific\n" );
  }

  fprintf( output_file, "Secondary structure flag      :: %-6d\n",  sec_flag );

  if( scom_flag == 1 ) {
    fprintf( output_file, "Contacts using CA and SCOM    :: %-6d\n",  scom_flag );
    fprintf( output_file, "Rebuilding sidechain COM      :: %-6d\n",  rebuilt_flag );
  }

  fprintf( output_file, "Total atomic contacts found   :: %-ld\n",   num_contact );
  fprintf( output_file, "Total interacting residues    :: %-d\n",   num_int_res );

  /* atomic contacts */
  if( ca_flag == 0 ) {
    fprintf( output_file , "\nR_ch R_resid R_name R_atom   L_ch L_resid L_name L_atom  Distance\n" ) ;
  
    for ( k = 0; k < num_contact; k++ ) {
      r_res  = contact[k].r_resid;
      l_res  = contact[k].l_resid;
      r_atom = contact[k].r_atomid;
      l_atom = contact[k].l_atomid;
      fprintf( output_file, "%4s %7s %-6s %6s   %4s %7s %-6s %6s %8.3f\n",
	       complex_structure.Residue[r_res].chainID,  complex_structure.Residue[r_res].res_seq_plus_iCode,
	       complex_structure.Residue[r_res].res_name, complex_structure.Residue[r_res].Atom[r_atom].atom_name,
	       complex_structure.Residue[l_res].chainID,  complex_structure.Residue[l_res].res_seq_plus_iCode,
	       complex_structure.Residue[l_res].res_name, complex_structure.Residue[l_res].Atom[l_atom].atom_name, contact[k].distance );
    }
  }

  /* residue-residue contacts */
  fprintf( output_file, "\nTotal residue-residue contacts found  :: %-d\n",  tot_num_res_cont );
  fprintf( output_file, " ResIndex  NumCont Contact_Residues (starts 1 from the first interface residue)\n" );
  for( i = 1; i <= complex_structure.length; i++ ) { 
    if( contres_num[i] == 0 && write_allres_flag == 0 ) continue;
    ii = is_int[i];
    fprintf( output_file, "RES  %4d%4d", ii, contres_num[i] );
    for( j = 0; j < contres_num[i] ; j++ ) {
      if( contres_lst[i][j] > -1 ) {
	jj = is_int[ contres_lst[i][j] ];
	fprintf( output_file, "%5d", jj );
      }
    }
    fprintf( output_file, "\n" );
  }


  fclose( output_file ) ;



  printf( "Finished\n\n" ) ;

  return( 0 ) ;

}

/*------------------------> end main <-----------------------------------*/




/*********************************************************************************/
/*--------------------------->    Subroutines      <-----------------------------*/
/*********************************************************************************/


/*** calculate contacts between ligand and receptor   ******************/
long  calc_contact ( struct Structure *structure, char *re_chain, char *li_chain,
		     float cutoff2, struct Contact *contact, int dist_flag ) {

  int l_residue, l_atom;
  int r_residue, r_atom;

  int l_resind, r_resind;

  long num_contact = 0;

  float  distance2, x, y, z;

  /*------------------>  End varaiable declaration <-------------------*/


  /*------------------>  contact calucations at atomic level <-------------------*/
  for( r_residue = 1; r_residue <= structure->length; r_residue++ ) {
    if( strchr( re_chain, structure->Residue[r_residue].chainID[0] ) == NULL ) continue;

    r_resind = structure->Residue[r_residue].nc;
    if( r_resind < 0 || r_resind > 19 ) { r_resind = 1; }

    for( r_atom = 1; r_atom <= structure->Residue[r_residue].size; r_atom++ ) {

      for( l_residue = 1; l_residue <= structure->length; l_residue++ ) {
	if( strchr( li_chain, structure->Residue[l_residue].chainID[0] ) == NULL ) continue;

	l_resind = structure->Residue[l_residue].nc;
	if( l_resind < 0 || l_resind > 19 ) { l_resind = 1; }

	for( l_atom = 1; l_atom <= structure->Residue[l_residue].size; l_atom++ ) {
	  x = structure->Residue[r_residue].Atom[r_atom].coord[1] - structure->Residue[l_residue].Atom[l_atom].coord[1];
	  y = structure->Residue[r_residue].Atom[r_atom].coord[2] - structure->Residue[l_residue].Atom[l_atom].coord[2];
	  z = structure->Residue[r_residue].Atom[r_atom].coord[3] - structure->Residue[l_residue].Atom[l_atom].coord[3];
	  distance2 = x*x + y*y + z*z;

	  if( dist_flag == 0 ) {
	    if( distance2 < cutoff2 ) {
	      contact[num_contact].r_resid  = r_residue;
	      contact[num_contact].l_resid  = l_residue;
	      contact[num_contact].r_atomid = r_atom;
	      contact[num_contact].l_atomid = l_atom;
	      contact[num_contact].distance = sqrt( distance2 );
	      num_contact++;
	      if( num_contact == MAX_NUM_CONTACT - 1 ) {
		printf( "Error: too many contacts\n");
		exit( EXIT_FAILURE );;
	      }
	    }
	  }
	  else {
	    if( distance2 < CG_CF2[r_resind][l_resind] ) {
	      contact[num_contact].r_resid  = r_residue;
	      contact[num_contact].l_resid  = l_residue;
	      contact[num_contact].r_atomid = r_atom;
	      contact[num_contact].l_atomid = l_atom;
	      contact[num_contact].distance = sqrt( distance2 );
	      num_contact++;
	      if( num_contact == MAX_NUM_CONTACT - 1 ) {
		printf( "Error: too many contacts\n");
		exit( EXIT_FAILURE );;
	      }
	    }
	  }
	}
      }

    }
  }

  return num_contact;
}


/***** print program usage *****/
void print_usage() {

  printf("\nUsage: extint <options>\n");
  printf("Options:\n");
  printf("\t-s  <pdbfile of the complex input> default is read from stdin\n");
  printf("\t-r  <receptor chain IDs> required\n");
  printf("\t-l  <ligand chain IDs> required\n");
  printf("\t-c  <contact output file> default is write to stdout\n");
  printf("\t-ca output only coordinates of ca atoms\n");
  printf("\t-i  <interface pdb file> default is do not write anything\n");
  printf("\t-e  calculate secondary structure default off\n");
  printf("\t    saved to the occupancy field: 1->coil, 2->helix, 3->turn, 4->strand\n");
  printf("\t-d  <contact distance cutoff> default 4.5 Angstrom\n");
  printf("\t-res  read a list files of residues, only consider contacts of these residues.\n");
  printf("\t     (file format: Index ChaindId Residue_id (5 chars including insertion code)\n");
  printf("\t-scom  contacts of sidechain centers of mass (default all-atom contacts)\n");
  printf("\t-norebuilt  do not rebuild sidechain COM, using original sidechain atom coordinates\n");
  printf("\t-wa write all residues including non-interfacial ones\n");
  printf("\t-h this help message\n\n");
  printf("Example: ./extint -s mycomplex.pdb -r A B -l X Y Z -d 6\n");
  printf("This calculate contacts within 6A between chains A B and chains X Y Z of mycomplex.pdb.\n\n");
}



/******************************************************************/
/*** calculate secondary structure according to Ca coordinates ****/
/*** 1->coil, 2->helix, 3->turn, 4->strand                     ****/
/******************************************************************/
void calc_sec_struct( struct Structure *structure ) {

  int *isec, atom;
  int i,j, j1,j2,j3,j4,j5;
  float dis13,dis14,dis15,dis24,dis25,dis35;
  float delta;

  double **ca_coor;

  /***** initialization ****/
  if( ( ca_coor = malloc( structure->length * sizeof(double *) ) ) == NULL ) { GENERAL_MEMORY_PROBLEM }
  for( i = 0; i < structure->length; i++ ) {
    if( ( ca_coor[i] = malloc( 3*sizeof(double) ) ) == NULL ) { GENERAL_MEMORY_PROBLEM }
  }

  if( ( isec = malloc( structure->length * sizeof(int *) ) ) == NULL ) { GENERAL_MEMORY_PROBLEM }

  /*************************/

  atomcoor_to_vector( *structure, " CA ", ca_coor );

  for( i = 0; i < structure->length; i++ ) {
    j1 = i-2;
    j2 = i-1;
    j3 = i;
    j4 = i+1;
    j5 = i+2;

    isec[i] = 1;
    if( j1 >= 0 && j5 < structure->length ) {
      if ( strcmp (structure->Residue[j1+1].chainID, structure->Residue[j3+1].chainID) == 0 && \
	   strcmp (structure->Residue[j5+1].chainID, structure->Residue[j3+1].chainID) == 0 ) {
	dis13 = dist( ca_coor[j1], ca_coor[j3] );
	dis14 = dist( ca_coor[j1], ca_coor[j4] );
	dis15 = dist( ca_coor[j1], ca_coor[j5] );
	dis24 = dist( ca_coor[j2], ca_coor[j4] );
	dis25 = dist( ca_coor[j2], ca_coor[j5] );
	dis35 = dist( ca_coor[j3], ca_coor[j5] );

	if(dis15 < 8) isec[i] = 3;

	delta=2.1;
	if(fabs(dis15-6.37) < delta)
	  if(fabs(dis14-5.18) < delta)
            if(fabs(dis25-5.18) < delta)
	      if(fabs(dis13-5.45) < delta)
		if(fabs(dis24-5.45) < delta)
		  if(fabs(dis35-5.45) < delta) {
		    isec[i] = 2;
		  } /* helix */

	delta=1.42;
	if(fabs(dis15-13) < delta && fabs(dis14-10.4) < delta && fabs(dis25-10.4) < delta && \
	   fabs(dis13-6.1)< delta && fabs(dis24-6.1 ) < delta && fabs(dis35-6.1)  < delta) {
	  isec[i] = 4; } /* strand */
      }
    }

    /*    printf("CA %4d %8.3f %8.3f %8.3f %2d\n", i+1, ca_coor[i][0], ca_coor[i][1], ca_coor[i][2], sec); */
  }


  /***** smooth secondary structure assignment ******/
  for( i = 0; i < structure->length; i++ ) {
    /*** boundary check, prevent cross chains ****/
    if( i-2 < 0 || i+2 >= structure->length ) continue;
    if ( strcmp (structure->Residue[i].chainID, structure->Residue[i-2].chainID) != 0 || \
	 strcmp (structure->Residue[i].chainID, structure->Residue[i+2].chainID) != 0 ) continue;

    if(isec[i] == 2 || isec[i] == 4) {
      j = isec[i];
      if( (isec[i-2] != j) && (isec[i-1] != j) && (isec[i+1] != j) && (isec[i+2] != j) ) isec[i]=1;
    }
  }

  for( i = 0; i < structure->length; i++ ) {
    if( i+5 >= structure->length ) continue;
    if ( strcmp (structure->Residue[i].chainID, structure->Residue[i+5].chainID) != 0 ) continue;

    if( (isec[i] != 2) && (isec[i+1] != 2) && (isec[i+2] == 2) && (isec[i+3] == 2) && \
	(isec[i+4] != 2) && (isec[i+5] != 2) ) {
           isec[i+2] = 1;
	   isec[i+3] = 1;
    }

    if( (isec[i] != 4) && (isec[i+1] != 4) && (isec[i+2] == 4) && (isec[i+3] == 4) && \
	(isec[i+4] != 4) && (isec[i+5] != 4) ) {
         isec[i+2]=1;
         isec[i+3]=1;
    }
  }

  for( i = 0; i < structure->length; i++ ) {
    if( i+2 >= structure->length ) continue;
    if ( strcmp (structure->Residue[i].chainID, structure->Residue[i+2].chainID) != 0 ) continue;

    if( (isec[i] == 2) && (isec[i+1] != 2) && (isec[i+2] == 2) ) isec[i+1] = 2;
    if( (isec[i] == 4) && (isec[i+1] != 4) && (isec[i+2] == 4) ) isec[i+1] = 4;
  }

  /********************************************************************/


  /***** tag secondary structure to the occupancy factor field ******/
  for( i = 0; i < structure->length; i++ ) {
    /*    printf("CA %4d %8.3f %8.3f %8.3f %2d\n", i+1, ca_coor[i][0], ca_coor[i][1], ca_coor[i][2], isec[i]); */
    for( atom = 1; atom <= structure->Residue[i+1].size; atom++ )
      structure->Residue[i+1].Atom[atom].occupancy = (float)isec[i];
  }

}



float dist ( double r1[3], double r2[3] ) {

  double dist;

  dist  = (r1[0] - r2[0])*(r1[0] - r2[0]) + (r1[1] - r2[1])*(r1[1] - r2[1]) + (r1[2] - r2[2])*(r1[2] - r2[2]);
  return (float) (sqrt( dist ));
}



void atomcoor_to_vector ( struct Structure myStructure, char *atom_name, double **vector )
{
  int i,j;

  for( i = 1 ; i <= myStructure.length ; i ++ ) {

    for( j = 1; j <= myStructure.Residue[i].size; j++ ) {
      if ( strcmp (myStructure.Residue[i].Atom[j].atom_name, atom_name) != 0 ) { continue; }
      vector[i-1][0] = myStructure.Residue[i].Atom[j].coord[1];
      vector[i-1][1] = myStructure.Residue[i].Atom[j].coord[2];
      vector[i-1][2] = myStructure.Residue[i].Atom[j].coord[3];
      break;
    }
  }

}



/******************************************************************/
/*** make a coarse-grained structure with Calpha and sidechain ****/
/*** center of mass, rebuilt if necessary                      ****/
/******************************************************************/
struct Structure mk_cg_struct( struct Structure aa_structure, int rebuilt_flag ) {
  struct Structure cg_struct;

  char atom_name[5];
  float scomx, scomy, scomz;

  int residue, atom;
  int tot_mass = 0, mass = 0;
  int i = 0, j = 0, counter;

  int *res_type, *id2resid;

  
  float *cax, *cay, *caz, *scmx, *scmy, *scmz;
  

 
  if( ( cg_struct.Residue = ( struct Amino_Acid * ) malloc ( ( aa_structure.length  + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

  if( ( res_type = malloc ( (aa_structure.length + 1) * sizeof(int)  )) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
  if( ( id2resid = malloc ( (aa_structure.length + 1) * sizeof(int)  )) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

  if( ( cax = malloc ( (aa_structure.length + 2) * sizeof(float)  )) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
  if( ( cay = malloc ( (aa_structure.length + 2) * sizeof(float)  )) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
  if( ( caz = malloc ( (aa_structure.length + 2) * sizeof(float)  )) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
  if( ( scmx = malloc ( (aa_structure.length + 2) * sizeof(float)  )) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
  if( ( scmy = malloc ( (aa_structure.length + 2) * sizeof(float)  )) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }
  if( ( scmz = malloc ( (aa_structure.length + 2) * sizeof(float)  )) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }



/*** find residues with Calpha atom, get their coordinates and numeric residue names ***/
  counter = 0;
  for( residue = 1; residue <= aa_structure.length; residue++ ) {
    for( atom = 1 ; atom <= aa_structure.Residue[residue].size ; atom++ ) {
      strcpy( atom_name, aa_structure.Residue[residue].Atom[atom].atom_name );
      if( strcmp( atom_name, " CA " ) == 0 ) {
	counter++;
	cax[counter] = aa_structure.Residue[residue].Atom[atom].coord[1];
	cay[counter] = aa_structure.Residue[residue].Atom[atom].coord[2];
	caz[counter] = aa_structure.Residue[residue].Atom[atom].coord[3];

	res_type[counter] = aa_structure.Residue[residue].nc;
	if( res_type[counter] < 0 || res_type[counter] > 19 ) { res_type[counter] = 1; } /* alanine for unknowns */ 

        id2resid[counter] = residue;

	break;
      }
    }
  }

  strcpy( cg_struct.ident, "coarse_grained" );
  cg_struct.length = counter;   /* ignore residues without Calpha */


/*** create a corase-graind structure contains only Calphas and COM of sidechains ***/

  for( j = 1; j <= cg_struct.length; j++ ) {
    residue = id2resid[j];
    cg_struct.Residue[j] = aa_structure.Residue[residue];
    /* initialize cg structure, note +2 (as compared to +1) is important to allocate a memory space for ca only models */
    if( ( cg_struct.Residue[j].Atom =
      ( struct Atom * ) malloc ( ( aa_structure.Residue[residue].size + 2 ) * sizeof_Atom ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }

    scomx = 0;
    scomy = 0;
    scomz = 0;
    tot_mass = 0;

    i = 1;

    for( atom = 1 ; atom <= aa_structure.Residue[residue].size ; atom++ ) {
      strcpy( atom_name, aa_structure.Residue[residue].Atom[atom].atom_name );
      if( strcmp( atom_name, " CA " ) == 0 ) {
	cg_struct.Residue[j].Atom[i] = aa_structure.Residue[residue].Atom[atom];
	/* printf("resid=%d, counter=%d, k=%d\n",residue, j, res_type[j]); */

	/* deal with GLY */
	if (res_type[j] == 0) {
	  scomx = aa_structure.Residue[residue].Atom[atom].coord[1];
	  scomy = aa_structure.Residue[residue].Atom[atom].coord[2];
	  scomz = aa_structure.Residue[residue].Atom[atom].coord[3];
	}

	continue;
      }
      else if( strcmp( atom_name, " N  " ) == 0 ) {
	continue;
      }
      else if( strcmp( atom_name, " C  " ) == 0 ) {
	continue;
      }
      else if( strcmp( atom_name, " O  " ) == 0 ) {
	continue;
      }
      else {
	mass = 0;
	if      ( strchr( atom_name, 'C' ) != NULL ) { mass = 12; }
	else if ( strchr( atom_name, 'N' ) != NULL ) { mass = 14; }
	else if ( strchr( atom_name, 'O' ) != NULL ) { mass = 16; }
	else if ( strchr( atom_name, 'S' ) != NULL ) { mass = 32; }

	scomx += mass * aa_structure.Residue[residue].Atom[atom].coord[1];
	scomy += mass * aa_structure.Residue[residue].Atom[atom].coord[2];
	scomz += mass * aa_structure.Residue[residue].Atom[atom].coord[3];
	tot_mass += mass;
      }
    }

    if( tot_mass > 0 ) {
      scomx = scomx / tot_mass;
      scomy = scomy / tot_mass;
      scomz = scomz / tot_mass;
    }

    i++;
    cg_struct.Residue[j].Atom[i] = aa_structure.Residue[residue].Atom[1];
    strcpy( cg_struct.Residue[j].Atom[i].atom_name, "SCOM" );
    cg_struct.Residue[j].Atom[i].coord[1]  = scomx;
    cg_struct.Residue[j].Atom[i].coord[2]  = scomy;
    cg_struct.Residue[j].Atom[i].coord[3]  = scomz;

    cg_struct.Residue[j].size = i;
  }




/*** rebuilt sidechain center of mass ***/
  if( rebuilt_flag == 1 ) {
    printf("Rebuilting side chain center of mass using only Calpha atoms %d\n", counter);
    get_scom( res_type, cg_struct.length, cax, cay, caz, scmx, scmy, scmz );

    for( i = 1; i <= counter; i++ ) {
      residue = id2resid[i];
      for( atom = 1 ; atom <= cg_struct.Residue[residue].size ; atom++ ) {
	strcpy( atom_name, cg_struct.Residue[residue].Atom[atom].atom_name );
	if( strcmp( atom_name, "SCOM" ) == 0 ) {
	  cg_struct.Residue[residue].Atom[atom].coord[1] = scmx[residue];
	  cg_struct.Residue[residue].Atom[atom].coord[2] = scmy[residue];
	  cg_struct.Residue[residue].Atom[atom].coord[3] = scmz[residue];
	}
      }
   }
  }


  free(cax ); free(cay ); free(caz );
  free(scmx); free(scmy); free(scmz);

  return cg_struct;
}




/******************************************************************/
/*** rebuilt side-chain center of mass from Calpha atoms       ****/
/******************************************************************/
void get_scom( int* res_type, int length, float* cax, float* cay, float* caz, float* scmx, float* scmy, float* scmz ) {

  int i,j,jm,jp,k;
  float amx,amy,amz,aaa,apx,apy,apz;
  float aaax,aaay,aaaz,ax,ay,az;
  float ccx,ccy,ccz,cx,cy,cz;
  float bbx,bby,bbz,bx,by,bz;
  float dx,dy,dz;
  float ang;


  cax[0]=cax[1]+(cax[2]-cax[3]);
  cay[0]=cay[1]+(cay[2]-cay[3]);
  caz[0]=caz[1]+(caz[2]-caz[3]);


  cax[length+1]=cax[length]+(cax[length-1]-cax[length-2]);         
  cax[length+1]=cay[length]+(cay[length-1]-cay[length-2]);
  cax[length+1]=caz[length]+(caz[length-1]-caz[length-2]);


  for(i=1; i<=length; i++) {
    j=i;
    jm=i-1;
    jp=i+1;
    amx=cax[j]-cax[jm];
    amy=cay[j]-cay[jm];
    amz=caz[j]-caz[jm];

    aaa=sqrt(amx*amx+amy*amy+amz*amz);
    amx=amx/aaa;
    amy=amy/aaa;
    amz=amz/aaa;
    apx=cax[jp]-cax[j];
    apy=cay[jp]-cay[j];
    apz=caz[jp]-caz[j];
    aaa=sqrt(apx*apx+apy*apy+apz*apz);
    apx=apx/aaa;
    apy=apy/aaa;
    apz=apz/aaa;

    ang=acos(-(amx*apx+amy*apy+amz*apz))*180/3.1415926;


    aaax=amx+apx;
    aaay=amy+apy;
    aaaz=amz+apz;
    aaa=sqrt(aaax*aaax+aaay*aaay+aaaz*aaaz);
    ax=aaax/aaa;
    ay=aaay/aaa;
    az=aaaz/aaa;


    ccx=amx-apx;
    ccy=amy-apy;
    ccz=amz-apz;
    aaa=sqrt(ccx*ccx+ccy*ccy+ccz*ccz);
    cx=ccx/aaa;
    cy=ccy/aaa;
    cz=ccz/aaa;


    bbx=amy*apz-amz*apy;
    bby=amz*apx-amx*apz;
    bbz=amx*apy-amy*apx;
    aaa=sqrt(bbx*bbx+bby*bby+bbz*bbz);
    bx=bbx/aaa;
    by=bby/aaa;
    bz=bbz/aaa;

    k=res_type[i];

    if(ang < 105) { /* alpha-helix or turn like*/
      dx = SC_PARA_A[k][0]*ax + SC_PARA_A[k][1]*bx + SC_PARA_A[k][2]*cx;
      dy = SC_PARA_A[k][0]*ay + SC_PARA_A[k][1]*by + SC_PARA_A[k][2]*cy;
      dz = SC_PARA_A[k][0]*az + SC_PARA_A[k][1]*bz + SC_PARA_A[k][2]*cz;
    }
    else {              /* beta-sheet */
      dx = SC_PARA_B[k][0]*ax + SC_PARA_B[k][1]*bx + SC_PARA_B[k][2]*cx;
      dy = SC_PARA_B[k][0]*ay + SC_PARA_B[k][1]*by + SC_PARA_B[k][2]*cy;
      dz = SC_PARA_B[k][0]*az + SC_PARA_B[k][1]*bz + SC_PARA_B[k][2]*cz;
    }


    scmx[i]=cax[j]+dx; 
    scmy[i]=cay[j]+dy;
    scmz[i]=caz[j]+dz; 

  }


}



/***********************************************************/
/******   read a list of selected residue               ****/ 
/***********************************************************/

char** read_residue_file( char *residue_file_name, int *num_sel_residue ) {

  char** sel_residue;
  int ind, num = 0;

  char  line_buffer[200];
  char  chainId[2], aa_name[4], resid[6], chain_resid[8], string[20];

  FILE  *res_file;


  sel_residue = NULL;

  /* read the list file */
  if( ( res_file = fopen( residue_file_name, "r" ) ) == NULL ) {
    printf( "Could not open %s for reading.\nDying\n\n" , residue_file_name ) ;
    exit( EXIT_FAILURE ) ;
  }


  while( fgets( line_buffer , 199 , res_file ) != NULL ) {
    sscanf( line_buffer, "%3s", string );
    if( strcmp( string, "###" ) == 0 ) continue;
    if( strcmp( string, "Ind" ) == 0 ) continue;
    else {
      sscanf( line_buffer, "%d %c%*c%5c %3c", &ind, chainId, resid, aa_name );

      strcpy( chain_resid, chainId );
      strcat( chain_resid, "." );
      strcat( chain_resid, resid );

      if( (sel_residue = realloc(sel_residue, sizeof(char *) * (num + 1) )) == NULL ) {
	GENERAL_MEMORY_PROBLEM
      }
      if( (sel_residue[num] = malloc(sizeof(char) * 8 )) == NULL ) {
	GENERAL_MEMORY_PROBLEM
      }

      strcpy( sel_residue[num], chain_resid );

      /* printf("num = %d, chain_resid = %s, sel_residue = %s\n",num, chain_resid, sel_residue[num]);   */

      num++;
    }
  }

  *num_sel_residue = num;

  return sel_residue;
}



/*******************************************************************/
/*** select a part of complex structure according to pre-selected **/
/*** residues                                                     **/
/*******************************************************************/
struct Structure sel_struct( struct Structure aa_structure, char** selected_residue, int num_sel_residue ) {
  struct Structure new_struct;

  char chain_resid[8];
  int j, residue, atom;
  int counter = 0, hit_flag;

  if( ( new_struct.Residue = ( struct Amino_Acid * ) malloc ( ( aa_structure.length  + 1 ) * sizeof_Amino_Acid ) ) == NULL ) {
    GENERAL_MEMORY_PROBLEM
  }

  for( residue = 1; residue <= aa_structure.length; residue++ ) {

    strcpy( chain_resid, aa_structure.Residue[residue].chainID );
    strcat( chain_resid, "." );
    strcat( chain_resid, aa_structure.Residue[residue].res_seq_plus_iCode );

    hit_flag = 0;
    for( j = 0; j < num_sel_residue; j++ ) {
      if( strcmp( selected_residue[j], chain_resid ) == 0 ) {
	hit_flag = 1;
	break;  /* hit a selected residue */
      }
    }

    if( hit_flag == 0 ) continue;  /* did not hit any selected residue */

    counter++;
    new_struct.Residue[counter] = aa_structure.Residue[residue];
    if( ( new_struct.Residue[counter].Atom =
      ( struct Atom * ) malloc ( ( aa_structure.Residue[residue].size + 1 ) * sizeof_Atom ) ) == NULL ) {
      GENERAL_MEMORY_PROBLEM
    }

    for( atom = 1 ; atom <= aa_structure.Residue[residue].size ; atom++ ) {
      new_struct.Residue[counter].Atom[atom] = aa_structure.Residue[residue].Atom[atom];
    }
  }

  strcpy( new_struct.ident, aa_structure.ident );
  new_struct.length = counter;   /* only considered selected residues */

  return new_struct;
}
