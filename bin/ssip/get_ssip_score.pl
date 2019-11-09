#!/usr/bin/perl -w
#################################################################################
# this script is used to call the SSIP program to calculate the profile score
# Input:  
#         mutList.txt          # contains the mutants that
#         structure_align.out  # interface structure alignment file by iAlign 
#         sequence_align.out   # interface sequence alignment file by PSI-BLAST
#         isscore              # a float value (0, 1] for scoring
# Output: 
#         ssip_score.txt       #predicted scores, one ddG each line
#################################################################################

use strict;
use warnings;

use File::Basename;
use Cwd 'abs_path';

if(@ARGV != 5){
  printf "Usage: get_ssip_score.pl mutList.txt structure_align.out sequence_align.out isscore ssip_score.txt\n";
  exit(-1);
}


my $bin_path=dirname(abs_path(__FILE__));

my $mut_list        = $ARGV[0];
my $structure_align = $ARGV[1];
my $sequence_align  = $ARGV[2];
my $isscore         = $ARGV[3];
my $ssip_score      = $ARGV[4];



my $ssippath      = "$bin_path/SSIPserver";
my $sort_align = "$bin_path/sort_align.pl"; 

#sort the structure alignment
`$sort_align $structure_align`;
`$sort_align $sequence_align`;
`cp -r $ssippath/parameter .`;
`$ssippath/SSIP --command=DDG  --isscore=$isscore  --mutant_file=$mut_list --structure_align=$structure_align --sequence_align=$sequence_align --ssipscore_file=$ssip_score`;
if(-d "./parameter"){
  `rm -rf parameter`;
}
exit;
