#!/usr/bin/perl

###############################################################################
#
# This script is used to calculate SSIPe ddG by calling the other programs 
# located in bin/.
# 
# Author: Xiaoqiang Huang;
# contact: xiaoqiah@umich.edu
################################################################################

use strict;
use warnings;

# PROGRAM PATH
use File::Basename;
use Cwd 'abs_path';
my $bin_path = dirname(abs_path(__FILE__));

my $run_ialign="$bin_path/bin/ssip/run_ialign.pl";
my $run_seqalign="$bin_path/bin/ssip/run_seqalign.pl";
my $get_ssip_score="$bin_path/bin/ssip/get_ssip_score.pl";
my $run_evoef="$bin_path/bin/evoef/run_evoef.pl";

#################################################
use Benchmark qw(:hireswallclock);
my $starttime = Benchmark->new;
my $finishtime;
my $timespent;

######################################
my $timeStart = localtime();

######################################
# INPUT-OPTIONS-DEFAULTS
my $pdb="complex.pdb";
my $mutList="mutList.txt";
#three forcefields for choice:
#SSIPE/SSIP/EVOEF
my $forcefield="SSIPE";
my $isscore=0.5;

if(@ARGV != 8){
  printf "######################################################################################################\n";
  printf "#                                                                                                    #\n";
  printf "#                                              SSIPe                                                 #\n";
  printf "#                                                                                                    #\n";
  printf "#   SSIPe is a method for accurate binding affinity change (ddGbind) prediction upon mutation at     #\n";
  printf "# protein-protein interface, using structural and sequence evolutionary profiles in combination with #\n";
  printf "# an optimized physical energy function                                                              #\n";
  printf "#                                                                                                    #\n";
  printf "#   Usage: run_SSIPe.pl -pdb complex.pdb -mutlist mutList.txt -forcefield FF -isscore SCORE          #\n";
  printf "#   FF can be chosen from SSIPE, SSIP or EVOEF                                                       #\n";
  printf "#   SCORE can be selected from [0.45, 0.55], default = 0.50                                          #\n";
  printf "######################################################################################################\n";
  exit;
}
else{
  $pdb=$ARGV[1];
  $mutList=$ARGV[3];
  $forcefield=$ARGV[5];
  $isscore=$ARGV[7];
}

if($forcefield eq 'SSIPE'){
  `$run_ialign $pdb $isscore structure_align.out`;
  `$run_seqalign $pdb structure_align.out sequence_align.out`;
  `$get_ssip_score $mutList structure_align.out sequence_align.out $isscore ssip_score.txt`;
  `$run_evoef $pdb $mutList evoef_score.txt`;
  my $wssip=0.734;
  my $wevoef=0.341;
  my $offset=0.205;
  open FILE,"<$mutList";
  chomp(my @mutants=<FILE>);
  close FILE;
  open FILE,"<ssip_score.txt";
  chomp(my @ssipscores=<FILE>);
  close FILE;
  open FILE,"<evoef_score.txt";
  chomp(my @evoefscores=<FILE>);
  close FILE;
  open FILE,">result.txt";
  printf FILE "DDG     SSIPscore EvoEFscore mutations\n";
  for(my $i=0;$i<@mutants;$i++){
    if($evoefscores[$i]==0.0 and $ssipscores[$i]==0.0){
      printf FILE "%-7.3f %-9.3f %-10.3f %s\n",0.0,0.0,0.0,$mutants[$i];
    }
    else{
      my $ddg=$wssip*$ssipscores[$i]+$wevoef*$evoefscores[$i]+$offset;
      printf FILE "%-7.3f %-9.3f %-10.3f %s\n",$ddg,$ssipscores[$i],$evoefscores[$i],$mutants[$i];
    }
  }
  close FILE;
}
elsif($forcefield eq 'SSIP'){
  `$run_ialign $pdb $isscore structure_align.out`;
  `$run_seqalign $pdb structure_align.out sequence_align.out`;
  `$get_ssip_score $mutList structure_align.out sequence_align.out $isscore ssip_score.txt`;
  open FILE,"<$mutList";
  chomp(my @mutants=<FILE>);
  close FILE;
  open FILE,"<ssip_score.txt";
  chomp(my @ssipscores=<FILE>);
  close FILE;
  open FILE,">result.txt";
  printf FILE "DDG     SSIPscore mutations\n";
  for(my $i=0;$i<@mutants;$i++){
    printf FILE "%-7.3f %-9.3f %s\n",$ssipscores[$i],$ssipscores[$i],$mutants[$i];
  }
  close FILE;
}
else{
  `$run_evoef $pdb $mutList evoef_score.txt`;
  open FILE,"<$mutList";
  chomp(my @mutants=<FILE>);
  close FILE;
  open FILE,"<evoef_score.txt";
  chomp(my @evoefscores=<FILE>);
  close FILE;
  open FILE,">result.txt";
  printf FILE "DDG     EvoEFscore mutations\n";
  for(my $i=0;$i<@mutants;$i++){
    printf FILE "%-7.3f %-10.3f %s\n",$evoefscores[$i],$evoefscores[$i],$mutants[$i];
  }
  close FILE;
}

# timestamping to end program
$finishtime = Benchmark->new;
$timespent = timediff($finishtime,$starttime);
print "time elapsed for SSIPe prediction: ". timestr($timespent),"\n";
exit;
