#! /usr/bin/perl -w
## A perl script to transform a pdb structure given a transformation matrix.
## Mu Gao

use strict;
use File::Basename;
use Getopt::Long;

###################### Input & options ###################
my $pdb_file;
my $out_file;
my $ialn_file;
my $help;
my @trans_vec;
my @rot_vec;
my @rot_vec1;

if ( @ARGV > 0 ) {
  GetOptions('i|ialign=s'   => \$ialn_file,   ### output from iAlign
             'p|pdbfile=s'  => \$pdb_file,    ### pdb file to be transformed
             'o|output=s'   => \$out_file,    ### file name of transformed pdb
             't|trans=f{3}' => \@trans_vec,   ### translation vector, 3 elements
             'r|rot=f{9}'   => \@rot_vec,     ### rotation matrix, 9 elements
             'r1|rot1=f{3}' => \@rot_vec1,    ### rotation vector, roll, pitch, yaw representation
             'h|help'       => \$help)        ### print the help message
  or printHelp();
}

printHelp() if( $help );

my @rot_mat = ();


if( defined $ialn_file ) {
  my ($ialn_trans, $ialn_rot) = readiAlignFile( $ialn_file );
  @trans_vec = @$ialn_trans;
  @rot_mat = @$ialn_rot;
}
elsif( scalar @rot_vec == 9 ) {
  my $k = 0;
  for(my $i=0; $i<3; $i++) {
    my @a = ();
    for(my $j=0; $j<3; $j++) {
      push( @a, $rot_vec[$k] );
      $k++;
    }
    push( @rot_mat, \@a );
  }
}
elsif( scalar @rot_vec1 == 3 ) {
  my $a = $rot_vec1[2];   ## yaw,   counterclockwise rotation of a about the z-axis
  my $b = $rot_vec1[1];   ## pitch, counterclockwise rotation of b about the y-axis
  my $c = $rot_vec1[0];   ## roll,  counterclockwise rotation of c about the x-axis

  my $r11 = cos($a)*cos($b);
  my $r12 = cos($a)*sin($b)*sin($c) - sin($a)*cos($c);
  my $r13 = cos($a)*sin($b)*cos($c) + sin($a)*sin($c);

  my $r21 = sin($a)*cos($b);
  my $r22 = sin($a)*sin($b)*sin($c) + cos($a)*cos($c);
  my $r23 = sin($a)*sin($b)*cos($c) - cos($a)*sin($c);

  my $r31 = -sin($b);
  my $r32 =  cos($b)*sin($c);
  my $r33 =  cos($b)*cos($c);

  @rot_mat = ( [$r11, $r12, $r13], [$r21, $r22, $r23], [$r31, $r32, $r33] );

  #printf "rotation matrix: %.3f %.3f %.3f %.3f\n",$r11, $r12, $r13, $rot_mat[2][0];
}


my $num1 = scalar @trans_vec;
my $num2 = scalar @rot_mat;
unless ( $num1 == 3 && $num2 == 3 && defined $pdb_file ) {
  printHelp();
}

unless( -s $pdb_file ) {
  die "Error: could not locate the pdb file: $pdb_file\n";
}


unless( defined $out_file ) {
  my ($file, $dir, $ext) = fileparse( $pdb_file, qr/\..*/ );
  $out_file = $dir . "$file\_rot.pdb";
}


####################################################################

if( transform_pdb( $pdb_file, $out_file, \@rot_mat, \@trans_vec ) ) {
  print "The transformed structure of $pdb_file was saved to $out_file.\n";
}









#####################################################################################
############################         Subroutine        ##############################
#####################################################################################


########### transform a pdb structure given a transformation ##############
sub transform_pdb {
  my ($inp_file, $out_file, $rot_mat, $trans_vec) = @_;


  open PDB, "<$inp_file" or die "Error: could not open $inp_file\n";
  open OUT, ">$out_file" or die "Error: could not write $out_file\n";
  while (my $line = <PDB>) {
    next if ( $line =~ /^REMARK/ );
    last if ( $line =~ /^ENDMDL/ );  #### only consider the first model

    #### PDB format: http://www.wwpdb.org/documentation/format23/sect9.html ####
    if ( $line =~ /^(ATOM  |HETATM)/ ) {
      (my $atomnum = substr($line, 6, 5)) =~ s/ //g;  # have to remove spaces
      (my $atomname = substr($line, 12, 4))=~ s/ //g;
      (my $altloc = substr($line, 16, 1));
      (my $resname = substr($line,17,3)) =~ s/ //g;
      (my $chain = substr($line, 21, 1));
      (my $resid = substr($line, 22, 4)) =~ s/ //g;
      (my $icode = substr($line, 26, 1)) =~ s/ //g;
      (my $x = substr($line,30,8)) =~ s/ //g;
      (my $y = substr($line,38,8)) =~ s/ //g;
      (my $z = substr($line,46,8)) =~ s/ //g;

      unless( defined $x and defined $y and defined $z ) {
	die "Error: cuold not find coordinates\n";
      }

      my $tx = $$trans_vec[0] + $x*$$rot_mat[0][0] + $y*$$rot_mat[0][1] + $z*$$rot_mat[0][2];
      my $ty = $$trans_vec[1] + $x*$$rot_mat[1][0] + $y*$$rot_mat[1][1] + $z*$$rot_mat[1][2];
      my $tz = $$trans_vec[2] + $x*$$rot_mat[2][0] + $y*$$rot_mat[2][1] + $z*$$rot_mat[2][2];

      substr($line,30,8) = sprintf("%8.3f", $tx);
      substr($line,38,8) = sprintf("%8.3f", $ty);
      substr($line,46,8) = sprintf("%8.3f", $tz);

      print OUT $line;
    }
    elsif( $line =~ /^TER/ ) {
      print OUT $line;
    }
  }
  close PDB;
  close OUT;

  return 1;
}


#### read transformation matrix from iAlign output ####
sub readiAlignFile {
  my $file = shift @_;

  my @trans_vec = ();
  my @rot_mat   = ();

  open INP, "<$file" or die "Error: could not open $file\n";
  while(<INP>){
    if( /Transformation matrix/ or /rotation matrix/ ) {
      my $line0 = <INP>;
      for(my $i=0; $i<3; $i++) {
	my $line   = <INP>;
	my @fields = split(' ',$line);
	$trans_vec[$i] = $fields[1];
	push( @rot_mat, [ ($fields[2], $fields[3], $fields[4]) ] );
      }
    }
  }
  close INP;


  return( \@trans_vec, \@rot_mat );
}



sub printHelp {

  print "Usage:\n";
  print "\tMatrix Representation:\n";
  print "\t$0 -p <pdb file>",
    " -i <transformation input file> -o <output file>\n\n";

  print "\t$0 -p <pdb file> -t <translation vector> -r <rotation matrix 9 elements> -o <output file>\n\n";

  print "\tRoll-Pitch-Yaw Representation:\n";
  print "\t$0 -p <pdb file>",
    " -t <translation vector> -r1 <rotation vector, roll (x), pitch (y), yaw (z)> -o <output file>\n\n";

  exit;
}
