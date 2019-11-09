#!/usr/bin/perl -w
###################################################################################################
# build_interface_alignment.pl
#
# This is a perl script used to build the protein complex MSA from two monomer MSAs via the STRING 
# protein links database.
# Input:
#         <-- protein.links.txt
#         <-- pdb_lig.msa
#         <-- pdb_rec.msa
# Output:
#         --> cpx_alignment.txt
#
# Xiaoqiang Huang
# xiaoqiah@umich.edu
# 08/01/2018
###################################################################################################

use 5.010;
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

if(@ARGV != 8){
  print "This script build_cpx_alignment.pl is used to build sequence alignment for a dimeric protein, the values in [] are required\n";
  print "Usage: build_cpx_alignment.pl -link [STRING_protein_links] -ligmsa [ligand.msa] -recmsa [receptor.msa] -out [protein_cpx_alignment.txt]\n";
  exit;
}

my $filepath = dirname(abs_path(__FILE__));

my $link = $ARGV[1];
my $msa1 = $ARGV[3];
my $msa2 = $ARGV[5];
my $cpxout = $ARGV[7];

# step 1: sort and read the multiple sequence alignment data
print "Sorting ligand and receptor MSA files ... ";
system "sort -k2 $msa1 > sorted_$msa1";
system "sort -k2 $msa2 > sorted_$msa2";
open MSA1, "<sorted_$msa1";
open MSA2, "<sorted_$msa2";
my @lines1 = <MSA1>;
my @lines2 = <MSA2>;
chomp(@lines1);
chomp(@lines2);
close MSA1;
close MSA2;
print "Done.\n";

# step 2: find common organisms of the two MSAs and build a pseudolink 
# between two proteins for each organism
print "Find common species of two MSAs and build links between proteins ... ";
my @comm = ();
for(my $k = 0; $k < @lines1; $k++){
  my $def1 = $lines1[$k];
  my @buff1 = split(/\s+/, $def1);
  $def1 = $buff1[1];
  my @buff2 = split(/\./,$def1);
  my $spe1 = $buff2[0];
  for(my $j = 0; $j < @lines2; $j++){
    my $def2 = $lines2[$j];
    my @buff3 = split(/\s+/, $def2);
    $def2 = $buff3[1];
    my @buff4 = split(/\./,$def2);
    my $spe2 = $buff4[0];
    if($spe1 == $spe2){
      push @comm, "$def1 $def2";
    }
  }
}
print "Done.\n";

# step 3: find the real link and the link score if possible via the species based link files
# this is the most time-consuming part of this program
print "Build complex links from the common species ... ";
my $firstLine = 1;
for(my $k = 0; $k < @comm; $k++){
  my $temp = $comm[$k];
  my @buff = split(/\./, $temp);
  my $speindex = $buff[0];
  my @buff2 = split(/\s+/, $temp);
  my $def1 = $buff2[0];
  if($firstLine == 1){
    `grep -r "$temp" $link.$speindex > cpx_link.out`;
    $firstLine = 0;
  }
  else{
    `grep -r "$temp" $link.$speindex >> cpx_link.out`;
  }
}
print "Done.\n";

# step 4: sort the real link result with sore ranked from high to low
# select the link with the highest score when there're several
# links between one protein with several different proteins
print "Sort complex links and select only the top 1 complex link for each protein ... ";
system "sort -k1,1 -k3nr,3  cpx_link.out > cpx_link_sort13.out";
open FILE, "<cpx_link_sort13.out";
my $initDef = "XXX";
$firstLine = 1;
while(my $line = <FILE>){
  chomp($line);
  my @buff = split(/\s+/, $line);
  my $def1 = $buff[0];
  if($def1 ne $initDef){
    if($firstLine == 1){
      `echo $line > cpx_link_sort13_top1.out`;
      $initDef = $def1;
      $firstLine = 0;
    }
    else{
      `echo $line >> cpx_link_sort13_top1.out`;
      $initDef = $def1;
    }
  }
}
close FILE;
system "sort -k2,2 -k3nr,3 cpx_link_sort13_top1.out > cpx_link_sort23.out";
open FILE, "<cpx_link_sort23.out";
$initDef = "XXX";
$firstLine = 1;
while(my $line = <FILE>){
  chomp($line);
  my @buff = split(/\s+/, $line);
  my $def1 = $buff[1];
  if($def1 ne $initDef){
    if($firstLine == 1){
      `echo $line > top_cpx_link.out`;
      $initDef = $def1;
      $firstLine = 0;
    }
    else{
      `echo $line >> top_cpx_link.out`;
      $initDef = $def1;
    }
  }
}
close FILE;
print "Done.\n";

# step 5: build the MSA file for the whole protein complex 
# the default medium-high cutoff link score is 0.4 [0.4, 0.7]
# the high link score is (0.7, 0.9]
# the ultra high link score cutoff is (0.9, 1.0)
print "Build complex MSA for protein complex ... ";
my $CUTOFF = 0.7;
open FILE, "<top_cpx_link.out";
open FILE2, ">$cpxout";
while(my $line = <FILE>){
  chomp($line);
  my @buff = split(/\s+/, $line);
  my $def1 = $buff[0];
  my $def2 = $buff[1];
  my $line1 = `grep -r $def1 $msa1`;
  my @temp1 = split(/\s+/, $line1);
  my $seq1 = $temp1[0];
  my $line2 = `grep -r $def2 $msa2`;
  my @temp2 = split(/\s+/, $line2);
  my $seq2 = $temp2[0];
  if($buff[2] > $CUTOFF){
    print FILE2 "$seq1 $seq2\n";
  }
}
close FILE;
close FILE2;
print "Done.\n";

system "rm sorted_$msa1 sorted_$msa2";
system "rm -rf cpx_link_sort13.out";
system "rm -rf cpx_link_sort13_top1.out";
system "rm -rf cpx_link_sort23.out";
system "rm -rf cpx_link.out";
system "rm -rf top_cpx_link.out";
