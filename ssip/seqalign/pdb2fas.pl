#!/usr/bin/env perl
#############################################################
# this script is used to convert pdb to fasta
#############################################################
use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

#programs
my $binpath = dirname(abs_path(__FILE__));
my $pdb2fas = "$binpath/PDB2FAS";

my $pdb = $ARGV[0];
my $fas = $ARGV[1];
my $fasseq = `$pdb2fas  $pdb`;
my $firstline = `head -1 $pdb`;
my $chnid  = substr $firstline, 21, 1;
`echo ">$chnid" > $fas`;
`echo $fasseq >> $fas`;
exit;

