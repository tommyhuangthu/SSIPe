/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund

Copyright (C) 2008 Mu Gao
Georgia Institute of Technology

*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


/************/

/* These values directly below may be altered, and the programs rebuilt */

#define MAX_ROTATIONS  100000
#define NUMBER_TO_KEEP 10000
#define NUMBER_OF_CONSTRAINTS 50
#define SAVED_HEADER_LINES 1000

/* I do not advise messing with anything below here */

/************/


#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

#define GENERAL_MEMORY_PROBLEM printf( "You do not have enough memory ([m|re]alloc failure)\nDying\n\n" ) ; exit( EXIT_FAILURE ) ;

/************/

/* The structures comprising a Structure (representation of an organic molecule in 3D) */

struct Atom{
	int		serial ;
	char		atom_name[5] ;
	float		coord[4] ;
	float		occupancy ;
	float		temp_factor ;
	float		charge ;
} ;

struct Amino_Acid{	
	char		res_name[4] ;
	char		chainID[2] ;
	char		res_seq_plus_iCode[6] ;
	char		olc[2] ;
	int		nc ;
	int		size ;
	struct Atom	*Atom ;
} ;

struct Structure{
	char			ident[256] ;
	int			length ;
	struct Amino_Acid	*Residue ;	
} ;








/************/

/* Memory allocation sizes */

#define sizeof_Atom		sizeof( struct Atom )
#define sizeof_Amino_Acid	sizeof( struct Amino_Acid )
#define sizeof_Structure	sizeof( struct Structure )

/************/

extern struct Structure read_pdb_to_structure( char *pdb_file_name ) ;
extern void write_structure_to_pdb( struct Structure This_Structure , char *pdb_file_name ) ;
extern void write_structure_to_pdbca( struct Structure This_Structure , char *pdb_file_name ) ;
extern void write_structure_to_trajectory_frame ( struct Structure This_Structure , FILE *traj_file, int frame );
extern struct Structure duplicate_structure( struct Structure This_Structure ) ;
extern struct Structure translate_structure( struct Structure This_Structure , float x_shift , float y_shift , float z_shift ) ;
extern struct Structure translate_structure_onto_origin( struct Structure This_Structure ) ;
extern struct Structure rotate_structure( struct Structure This_Structure , int z_twist , int theta , int phi ) ;
extern struct Structure merge_structures( struct Structure Structure_One , struct Structure Structure_Two ) ;
extern float radius_of_structure( struct Structure This_Structure ) ;
extern float total_span_of_structures( struct Structure Structure_1 , struct Structure Structure_2 ) ;

extern void  free_structure( struct Structure *pStr );



extern void atom_coor_to_vector ( struct Structure myStructure, char *atom_name, double vector[][3] );
extern void transform_vector  ( double vector[][3], int dim, double trans[3], double rotation[3][3] );



