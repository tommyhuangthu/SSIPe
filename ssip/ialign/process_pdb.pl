#!/usr/bin/perl -w

# Parse PDB molecule
# Copyright (C) 1997-2000 Gidon Moont
# 
# Biomolecular Modelling Laboratory
# Imperial Cancer Research Fund
#
# Copyright (C) 2008 Mu Gao
# Georgia Institute of Technology
#

use strict ;
use File::Basename;

my $sys_com ;
my $record ;
my $i ;
my $j ;
my $info ;

#############
#
# ensure messages appear asap

select STDOUT ;
$| = 1 ;

#############
#
# Other programs and files used by this one

BEGIN{

  my $full_program_name = $0 ;
  my $relative_path ;
  ( $relative_path ) = ( $full_program_name =~ /^(.+)\/.+$/ ) ;

  unshift( @INC , ( "$relative_path" ) ) ;

}

use PDB_Types ;
use PDB_Parse ;

print STDOUT "$ARGV[0] : $ARGV[1]\n" ;

#############
#
# Command line options

my $pdb_file = 'undefined' ;
my $warn = 1 ;
my $multidock_flag = 0 ;
my $out_file;

while( @ARGV ) {

  my $bit = shift @ARGV ;

  if( $bit eq '-nowarn' ) {
    $warn = 0 ;
  }

  if( $bit eq '-multidock' ) {
    $multidock_flag = 1 ;
  }

  if( $bit eq '-pdb' ) {
    $pdb_file = shift @ARGV ;
  }

  if( $bit eq '-o' ) {
    $out_file = shift @ARGV ;
  }
}

#############
#
# Screen notices

my $full_program_name = $0 ;
my $program_name ;
( $program_name ) = ( $full_program_name =~ /^.+\/(.+)$/ ) ;

print STDOUT "\nRunning $program_name ...\n" ;

print STDOUT "pdb file     :: ".$pdb_file."\n" ;
print STDOUT "warn         :: ".$warn."\n" ;
#print STDOUT "multidock    :: ".$multidock_flag."\n" ;

#print STDOUT "\n" ;

#############
#
# Get ATOM records

my @atomrecords = () ;
open( PDB_File , $pdb_file ) || die "Could not open molecule file\n" ;
while( defined( $record = <PDB_File> ) ) {

  last if ( $record =~ /^ENDMDL/ );  ### Mu: only consider the first model, e.g., from multiple NMR models

  if( $record =~ /^ATOM  / ) {
    push( @atomrecords , $record ) ;
  }
}
close( PDB_File ) ;

#############
#
# Parse

my @correct_records = PDB_Parse::parse( $warn , $multidock_flag , @atomrecords ) ;

my $correct_sequence = pop @correct_records ;

if( $multidock_flag == 0 ) {
  @correct_records = PDB_Parse::add_types( @correct_records ) ;
}

#############
#
# Fasta style file dump

#my $name ;
#( $name ) = ( $pdb_file =~ /^(.+)\..+/ ) ;

my ($name,$dir,$ext) = fileparse($pdb_file, qr/\..*/);

#open( Dump , "> ".$name.".fasta" ) || die "Could not open for writing\n" ;
#print Dump "> ".$name."\n" ;
#print Dump $correct_sequence ;
#close( Dump ) ;

#############
#
# PDB style file dump

unless (defined $out_file) { $out_file = "$name.parsed"; }

open( Dump , "> ", $out_file ) || die "Could not open for writing\n" ;
print Dump @correct_records ;
close( Dump ) ;

#############
#
# Finished

exit ;
