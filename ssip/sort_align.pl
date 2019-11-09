#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use File::Basename;

if(@ARGV != 1){
  print "Usage: ./sort_ialign.pl unsorted_ialign.txt\n";
  exit(-1);
}

my $unsorted=$ARGV[0];
`head -8 $unsorted > X_tempfile0`;
`sed -n '9,\$p' $unsorted | sort -k6,6nr > X_tempfile1`;
`cat X_tempfile0 X_tempfile1 > $unsorted`;
`rm -rf X_tempfile0 X_tempfile1`;
exit;
