#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use Cwd 'abs_path';
use File::Basename;

my $complex = $ARGV[0];
my $ligpdb  = $ARGV[1];
my $recpdb  = $ARGV[2];

open FILE, "<$complex";
open FILE1, ">$ligpdb";
open FILE2, ">$recpdb";
my $firstLine = 1;
my $lig_chnid = 'A';
while(<FILE>){
  chomp(my $line = $_);
  my $key = substr $line, 0, 4;
  if($key ne "ATOM"){next;}
  my $chnid = substr $line, 21, 1;
  my $conf = substr $line, 16, 1;
  my $resi = substr $line, 17, 3;
  if($firstLine == 1){
    $lig_chnid = $chnid;
    $firstLine = 0;
  }
  if($chnid eq $lig_chnid){
    if($conf eq 'A' or $conf eq '1'){
      $line =~ s/$conf$resi/ $resi/g;
      printf FILE1 "$line\n";
    }
    elsif($conf eq ' '){
      printf FILE1 "$line\n";
    }
  }
  else{
    if($conf eq 'A' or $conf eq '1'){
      $line =~ s/$conf$resi/ $resi/g;
      printf FILE2 "$line\n";
    }
    elsif($conf eq ' '){
      printf FILE2 "$line\n";
    }
  }
}
close FILE1;
close FILE2;
close FILE;
exit;

