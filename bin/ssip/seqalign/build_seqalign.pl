#!/usr/bin/env perl

###########################################################################################################
# this script is used to extract sequence-based interface alignment from structure alignment file
# Input: 
#        cpx_alignment.txt     # whole sequence alignment
#        structure_align.out   # structure alignment file which contains interface residues in two chains
#        lig.pdb               # pdb file of the first chain
#        rec.pdb               # pdb file of the second chain
# Output:
#        sequence_align.out    # sequence alignment file with similar format to structure_align.out
###########################################################################################################

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

if(@ARGV != 5){
  printf "Usage: build_seqalign.pl cpx_alignment.txt structure_align.out lig.pdb rec.pdb sequence_align.out\n";
  exit;
}

my $cpx_alignment   = $ARGV[0];
my $structure_align = $ARGV[1];
my $ligpdb          = $ARGV[2];
my $recpdb          = $ARGV[3];
my $sequence_align  = $ARGV[4];

if(not -e $structure_align){
  printf("iAlign interface alignment file $structure_align doesn't exist, cannot build sequence interface alignment\n");
  exit;
}

#sort the cpx_alignment file
#`sort -k3,3nr -k1,1r -k2,2r $cpx_alignment > $cpx_alignment.sort`;
#`mv $cpx_alignment.sort $cpx_alignment`;

my %AminoList3to1 = (
  "GLY" => "G", "PRO" => "P", "ASP" => "D", "GLU" => "E", "LYS" => "K",
  "ARG" => "R", "HIS" => "H", "SER" => "S", "THR" => "T", "ASN" => "N",
  "GLN" => "Q", "ALA" => "A", "MET" => "M", "TYR" => "Y", "TRP" => "W",
  "VAL" => "V", "ILE" => "I", "LEU" => "L", "PHE" => "F", "CYS" => "C",
);

#my $ligpdb = "lig.pdb";
#my $recpdb = "rec.pdb";

#read in all the amino acid index and name for ligand protein from PDB
my (@ligindex, @ligaas, @ligchainid);
open FILE1, "<$ligpdb";
my $iniIndex = -100;
while(my $line = <FILE1>){
  chomp($line);
  my $key = substr $line, 0, 4;
  if($key ne "ATOM"){next;}
  my $res = substr $line, 17, 3;
  my $index = substr $line, 22, 4;
  my $chainid = substr $line, 21, 1;
  $index =~ s/\s+//g;
  if($index != $iniIndex){
    push @ligaas, $AminoList3to1{$res};
    push @ligindex, $index;
    push @ligchainid, $chainid;
    $iniIndex = $index;
  }
}
close FILE1;
#read in all the amino acid index and name for receptor protein from PDB
my (@recindex, @recaas, @recchainid);
open FILE1, "<$recpdb";
$iniIndex = -100;
while(my $line = <FILE1>){
  chomp($line);
  my $key = substr $line, 0, 4;
  if($key ne "ATOM"){next;}
  my $res = substr $line, 17, 3;
  my $index = substr $line, 22, 4;
  my $chainid = substr $line, 21, 1;
  $index =~ s/\s+//g;
  if($index != $iniIndex){
    push @recaas, $AminoList3to1{$res};
    push @recindex, $index;
    push @recchainid, $chainid;
    $iniIndex = $index;
  }
}
close FILE1;

#read in the INTERFACE amino acid index and name from the bindprofx interface alignment file
#note that:
#(1) chain_ID_A in bindprofx can be different from that in sip
#(2) residue_A in bindprofx can be different from that in sip
#(3) chain_ID_B in bindprofx can be different from that in sip
#(4) residue_B in bindprofx can be different from that in sip
my $bindprofx_align = $structure_align;
open FILE, "<$bindprofx_align";
chomp(my @lines = <FILE>);
close FILE;
my (@lig_inf_index, @lig_inf_index_from_zero);
my (@rec_inf_index, @rec_inf_index_from_zero);
for(my $j = 0; $j < @lines; $j++){
  my $line = $lines[$j];
  my @temp = split(/\s+/, $line);
  # find the interface residues for ligand
  if($line =~ /chain_ID_A:/){
    my $line1 = $lines[$j+1];
    my @temp1 = split(/\s+/, $line1);
    if($temp[1] eq $ligchainid[0]){
      my $from = 0;
      for(my $k = 1; $k < @temp1; $k++){
        my $index = $temp1[$k];
        for(my $s = $from; $s < @ligindex; $s++){
          if($index == $ligindex[$s]){
            push @lig_inf_index_from_zero, $s;
            push @lig_inf_index, $index;
            $from = $s+1;
          }
        }
      }
    }
    else{
      my $from = 0;
      for(my $k = 1; $k < @temp1; $k++){
        my $index = $temp1[$k];
        for(my $s = $from; $s < @recindex; $s++){
	  if($index == $recindex[$s]){
	    push @rec_inf_index_from_zero, $s;
	    push @rec_inf_index, $index;
	    $from = $s+1;
	  }
        }
      }
    }
  }
  elsif($line =~ /chain_ID_B:/){
    my $line1 = $lines[$j+1];
    my @temp1 = split(/\s+/, $line1);
    if($temp[1] eq $ligchainid[0]){
      my $from = 0;
      for(my $k = 1; $k < @temp1; $k++){
        my $index = $temp1[$k];
        for(my $s = $from; $s < @ligindex; $s++){
          if($index == $ligindex[$s]){
            push @lig_inf_index_from_zero, $s;
            push @lig_inf_index, $index;
            $from = $s+1;
          }
        }
      }
    }
    else{
      my $from = 0;
      for(my $k = 1; $k < @temp1; $k++){
        my $index = $temp1[$k];
        for(my $s = $from; $s < @recindex; $s++){
          if($index == $recindex[$s]){
            push @rec_inf_index_from_zero, $s;
            push @rec_inf_index, $index;
            $from = $s+1;
          }
        }
      }
    }
  }
}

#build the BindProfX-like interface alignment file
open ALIGN, ">$sequence_align";
printf ALIGN "Interface residues:\n";
printf ALIGN "chain_ID_A: $ligchainid[0]\n";
printf ALIGN "residue_A:";
for(my $i = 0; $i < @lig_inf_index; $i++){
  printf( ALIGN "%5d ", $lig_inf_index[$i]);
}
printf ALIGN "\n";
printf ALIGN "chain_ID_B: $recchainid[0]\n";
printf ALIGN "residue_B:";
for(my $i = 0; $i < @rec_inf_index; $i++){
  printf( ALIGN "%5d ", $rec_inf_index[$i]);
}
printf ALIGN "\n";

printf ALIGN "Alignments:\n";
printf ALIGN "Name   chA chB ifResNum alnResNum linkscore ifSeqA ifSeqB\n";
# wildtype sequence
my $interfacenum = scalar(@lig_inf_index) + scalar(@rec_inf_index);
my $alignnum = $interfacenum;
printf ALIGN "%-7s%3s %3s %4d %4d %6.3f  ", "target", $ligchainid[0], $recchainid[0], $interfacenum, $alignnum, $alignnum/$interfacenum;
for(my $k = 0; $k < @lig_inf_index_from_zero; $k++){
  printf ALIGN "%s", $ligaas[$lig_inf_index_from_zero[$k]];
}
printf ALIGN " ";
for(my $k = 0; $k < @rec_inf_index_from_zero; $k++){
  printf ALIGN "%s", $recaas[$rec_inf_index_from_zero[$k]];
}
printf ALIGN "\n";

#other sequences from STRING database
open STRING, "<$cpx_alignment";
@lines = <STRING>;
chomp(@lines);
for(my $k = 0; $k < @lines; $k++){
  my $line = $lines[$k];
  my @strings = split(/\s+/, $line);
  my $lig = $strings[0];
  my $rec = $strings[1];
  my $linkscore = $strings[2];
  my @ligseq = split(//,$lig);
  my @recseq = split(//,$rec);
  my $name = "aln$k";
  #my $name = sprintf "aln%03d",$k+1;
  my $num = 0;
  for(my $j = 0; $j < @lig_inf_index_from_zero; $j++){
    if($ligseq[$lig_inf_index_from_zero[$j]] ne '-'){
      $num++;
    }
  }
  for(my $j = 0; $j < @rec_inf_index_from_zero; $j++){
    if($recseq[$rec_inf_index_from_zero[$j]] ne '-'){
      $num++;
    }
  }
  printf ALIGN "%-7s%3s %3s %4d %4d %6.3f  ", $name, $ligchainid[0], $recchainid[0], $interfacenum, $num, $linkscore;
  for(my $j = 0; $j < @lig_inf_index_from_zero; $j++){
    printf ALIGN "%s", $ligseq[$lig_inf_index_from_zero[$j]];
  }
  printf ALIGN " ";
  for(my $j = 0; $j < @rec_inf_index_from_zero; $j++){
    printf ALIGN "%s", $recseq[$rec_inf_index_from_zero[$j]];
  }
  printf ALIGN "\n";
}
close STRING;
close ALIGN;

