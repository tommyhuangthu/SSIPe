#!/usr/bin/perl -w
######################################################
#this script is used to calculate the EvoEF ddG scores
# Input:  complex.pdb
#         mutList.txt
# Output: mutant models
#         evoef_score.txt
######################################################

use strict;
use warnings;

use File::Basename;
use Cwd 'abs_path';
my $bin_path=dirname(abs_path(__FILE__));
my $evoef = "$bin_path/EvoEF/src/EvoEF";
my $par_dir = "$bin_path/EvoEF/src/data";

if(@ARGV != 3){
  printf "Usage: run_evoef.pl  complex.pdb  mutList.txt  evoef_score.txt\n";
  exit;
}


my $complex   = $ARGV[0];
my $mutList   = $ARGV[1];
my $scorefile = $ARGV[2];

#read mutList file
open MUTANT, "<$mutList";
chomp(my @mutants = <MUTANT>);
close MUTANT;

#copy evoef parameter folder to the current folder
if(not -d "data"){
  `cp -r $par_dir .`;
}

my $individual_list = "individual_list.txt";
open EVOSCORE, ">$scorefile";
for(my $i=0; $i<@mutants; $i++){
  #check mutation
  my $mutation=$mutants[$i];
  $mutation =~ s/;//g;
  my @new=();
  my @temp=split(/,/,$mutation);
  for(my $j=0;$j<@temp;$j++){
    my $first=substr $temp[$j],0,1;
    my $last=substr $temp[$j],-1;
    if($first ne $last){
     push @new, $temp[$j];
    }
  }

  #if the mutation is that amino acids are mutated into themselves
  if(@new==0){
    my $preddg=0;
    printf EVOSCORE "%.3f\n",$preddg;
    next;
  }
  #compile the effective mutations
  my $newmutation='';
  for(my $j=0;$j<@new;$j++){
    if($j==0){
      $newmutation=$newmutation.$new[$j];
    }
    else{
      $newmutation=$newmutation.",".$new[$j];
    }
  }

  #generate a temp $individual_list file
  open INDI,">$individual_list";
  #print INDI "$mutants[$i]\n";
  print INDI "$newmutation;\n";
  close INDI;
  #repair and minimize the initial structure
  my $cpx1 = $complex;
  $cpx1 =~ s/.pdb//g;
  if(not -e "$cpx1\_Repair.pdb"){
    `$evoef --command=RepairStructure --pdb=$complex`;
  }
  #generate the models and calculate the binding energy
  my $string = $mutants[$i];
  $string =~ s/;//g;
  $string =~ s/,/_/g;
  if(not -e "$cpx1\_$string.pdb"){
    `$evoef --command=BuildMutant --pdb=$cpx1\_Repair.pdb --mutant-file=$individual_list`;
    `mv $cpx1\_Repair_Model_1.pdb $cpx1\_$string.pdb`;
  }
  my $preddg = 0;
  if(-e "$cpx1\_Repair.pdb" and -e "$cpx1\_$string.pdb"){
    my $dgwt = evalCpx("$cpx1\_Repair.pdb");
    my $dgmt = evalCpx("$cpx1\_$string.pdb");
    $preddg = $dgmt - $dgwt;
  }
  else{
    #print "EvoEF failed generating wild-type or mutant model for mutant $mutants[$i]\n";
    print "EvoEF failed generating wild-type or mutant model for mutant $newmutation\n";
    $preddg = 0;
  }
  printf EVOSCORE "%.3f\n",$preddg;
  if(-e $individual_list){
    `rm -rf $individual_list`;
  }
}
close EVOSCORE;

#clean the directory
`rm -rf data`;
exit;

sub evalCpx{
  my $pdb=shift;
  my $temp = `$evoef --command=ComputeBinding --pdb=$pdb | grep \"Total                 =\" | tail -n 1`;
  my @buff=split(/\s+/, $temp);
  return $buff[$#buff];
}
