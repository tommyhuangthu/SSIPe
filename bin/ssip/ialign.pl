#!/usr/bin/perl -w

### The master perl script of iAlign.
### Mu Gao
### E-mail: <mu.gao@gatech.edu>
###
### changed by Chengxin to avoid writing to library directory

use strict;
use File::Basename;
use File::Path;
use Getopt::Long;
use Pod::Usage;
use Time::Local;
use Cwd 'abs_path';


############ Global Parameters ############
#### Do not modify these parameters unless you know what you are doing
my $CONT_CUTOFF  = 4.5;    #### 4.5 Angstrom distance cutoff for contact calculations
my $MIN_NUM_RES  = 25;     #### minimum number of protein residues, ignore short peptide in a PDB file
my $MIN_INT_RES  = 20;     #### minimum number of protein residues of a interface

my $VERSION = '1.0b7';  ### version number of this package

###################### Input & options ###################
my $pdb_file1;
my $pdb_file2;
my $pchains1 = '';
my $pchains2 = '';

my $work_path = '.';
my $pdb_path;

my $out_file;
my $set_file1;
my $set_file2;

my $help;
my $manual;
my $version;

my $seq_flag;
my $tra_flag;
my $aln_flag = 2;
my $lib_flag;
my $sel_flag;
my $measure  = 'is';
my $norm;
my $nrd_cf;
my $speed_mode;
my $vmd;

  GetOptions(
	     #### input files
	     'l|liblst1|l1=s'=> \$set_file1,    ### a list of entries
	     'l2|liblst2=s'  => \$set_file2,    ### a list of entries
	     'p1|pdb1=s'     => \$pdb_file1,    ### pdb file1
	     'p2|pdb2=s'     => \$pdb_file2,    ### pdb file2
	     'c1|chains1=s'  => \$pchains1,     ### list of protein chains of pdb1
	     'c2|chains2=s'  => \$pchains2,     ### list of protein chains of pdb2
	     #### paths
	     'p|pdbpath=s'   => \$pdb_path,     ### path to original pdb files
	     'w|workpath=s'  => \$work_path,    ### path to input files of interfaces
	     #### output file
             'o|output=s'    => \$out_file,     ### save results to an output file
	     #### options
             'g|mklib'       => \$lib_flag,     ### make a library of interfaces w/o alignment
             'nr=f'          => \$nrd_cf,       ### remove redundant interfaces within a pdb record
             't|trans'       => \$tra_flag,     ### the first structure aligned to the second
             's|nonseq'      => \$seq_flag,     ### non-sequential mode (default off)
             'a|aln=i'       => \$aln_flag,     ### choose a format of alignment printout
             'do_self'       => \$sel_flag,     ### allow alignment against itself
             'e|measure=s'   => \$measure,      ### similarity measure: TM-score or IS-score
             'n|norm=s'      => \$norm,         ### score normalization methods
             'q|quick=i'     => \$speed_mode,   ### speed mode
             'vmd=s'         => \$vmd,          ### generate tcl script for VMD
	     #### advanced options
	     'dc=f'          => \$CONT_CUTOFF,  ### distance cutoff for extracting interface
	     'minp=i'        => \$MIN_NUM_RES,  ### minimum number of protein residues
	     'mini=i'        => \$MIN_INT_RES,  ### minimum number of interface residues
	     #### help
             'h|help'        => \$help,         ### print the help message
             'm|manual'      => \$manual,       ### print the manual page
             'v|version'     => \$version       ### print the version number;
            ) or pod2usage(0);


  pod2usage(-verbose => 0) if $help;
  pod2usage(-verbose => 2) if $manual;

  if ($version) {
    print "iAlign version $VERSION\n" ;
    exit;
  }



  my $num_arg = $#ARGV + 1;
  if( $num_arg == 4 ) {
    $pdb_file1 = $ARGV[0];
    $pchains1  = $ARGV[1];
    $pdb_file2 = $ARGV[2];
    $pchains2  = $ARGV[3];
  }
  elsif( $num_arg == 2 and not defined $set_file1 ) {
    $pdb_file1 = $ARGV[0];
    $pdb_file2 = $ARGV[1];
  }
  elsif( $num_arg == 2 and defined $set_file1 ) {
    $pdb_file2 = $ARGV[0];   ### template is the first, target/query is the second
    $pchains2  = $ARGV[1];
  }
  elsif( $num_arg == 1 and not defined $set_file1 ) {
    $set_file1 = $ARGV[0];
  }
  elsif( $num_arg == 1 and defined $set_file1 ) {
    $pdb_file2 = $ARGV[0];
  }
  elsif( $num_arg > 0 ) {
    pod2usage(-verbose =>0, -msg=>'Error: invalid argument.');
  }


  unless( (defined $pdb_file1 and defined $pdb_file2) or
	(defined $pdb_file1 and defined $set_file1) or
	(defined $set_file1) ) {
    pod2usage(-verbose =>0, -msg=>'Error: invalid argument.');
  }


  if( defined $vmd ) {
    $aln_flag = 2;
    if( $vmd ne 'cartoon' and $vmd ne 'ribbon' ) {
      $vmd = 'cartoon';
    }
  }




#### creating two sets of structures for IS-align runs  ####
my @pdb_set1 = ();
my @pdb_set2 = ();

if( defined $pdb_file1 ) { $pdb_set1[0] = {'pdb_file'=>$pdb_file1,'pchains'=>$pchains1}; }
if( defined $pdb_file2 ) { $pdb_set2[0] = {'pdb_file'=>$pdb_file2,'pchains'=>$pchains2}; }
if( defined $set_file1 ) { readSetFile( $set_file1, \@pdb_set1 ); }
if( defined $set_file2 ) { readSetFile( $set_file1, \@pdb_set2 ); }

### get the path of executables
#my ($bin_base, $bin_path, $bin_ext) = fileparse( $0, qr/\..*/ );
my $bin_path=dirname(abs_path(__FILE__))."/ialign";

### create the work_path if it does not exist
unless( -d $work_path ) {
  mkpath( $work_path ) or die "Error: could not create the directory $work_path\n";
}


### options passing to IS-align
my $switch = '';
if( defined $seq_flag ) { $switch .= ' -s '; }
if( $measure eq 'tm' )  { $switch .= ' -t '; }
elsif( $measure ne 'is' ) {
  die "Error: the measure option can only be \"tm\" or \"is\".\n";
}



if( $aln_flag == 0 or $aln_flag == 1 or $aln_flag == 2 ) {
  $switch .= " -v $aln_flag ";
}
else { die "Error: invalid option for alignment printout format\n"; }

if( defined $norm ) {
  if    ( $norm eq 'short' ) { $switch .= ' -b '; }
  elsif ( $norm eq 'long'  ) { $switch .= ' -c '; }
  elsif ( $norm eq 'ave'   ) { $switch .= ' -a '; }
  else { die "Error: the norm option can only be: \"short\",\"long\", or \"ave\"\n"; }
}
if( defined $speed_mode ) {
  if ( $speed_mode >= 0 && $speed_mode <= 2 ) { $switch .= " -q $speed_mode"; }
  else { die "Error: the speed mode can only be: 0, 1 (default), 2\n"; }
}

##-------------------------- end of input & options   -------------------------##


#### determine where to ouput results ####
my $out;
if( defined $out_file ) {
  open $out, ">$out_file" or die "Error: could not open $out_file\n";
}
else {
  open $out, ">&STDOUT";
}
print_banner( $out );


print $out "iAlign starts at ", timeStamp(), "\n";

######################################################################
######### Stage 1: processing raw PDB files, sanity check   ##########
######################################################################
print $out "\n+++++ Stage 1: Processing PDB files ...\n";
print "+++++ Stage 1: Processing PDB files ...\n" if(defined $out_file);

my %parsed = ();   ### save information of processed PDB files
for (@pdb_set1) {
  my $name = $_->{pdb_file} . '_' . $_->{pchains};
  $parsed{$name} = $_;
}
for (@pdb_set2) {
  my $name = $_->{pdb_file} . '_' . $_->{pchains};
  $parsed{$name} = $_;
}

for my $name (sort keys %parsed) {
  my ($parsed_pdb_file, $pchains) = parsePDBFile( $parsed{$name}, $pdb_path, $work_path, $bin_path, $out );
  $parsed{$name}->{'parsed_file'} = $parsed_pdb_file;
  $parsed{$name}->{'chains'}      = $pchains;
}




######################################################################
######### Stage 2:  extracting protein-protein interfaces  ###########
######################################################################
print $out "\n+++++ Stage 2: Extracting protein-protein interfaces ...\n";
print $out "Mininum number of interface residue required: $MIN_INT_RES AAs\n";
print "+++++ Stage 2: Extracting protein-protein interfaces ...\n" if (defined $out_file);

if( defined $nrd_cf ) {
  if( $nrd_cf < 0.01 or $nrd_cf > 1 ) { $nrd_cf = 0.9; }  ### default cutoff is TM-score 0.9
  print $out "Redundant interfaces within a PDB will be removed at a TM-score cutoff of $nrd_cf\n";
}

for my $name (sort keys %parsed) {
  my $parsed_pdb_file = $parsed{$name}->{parsed_file};
  my $pchains         = $parsed{$name}->{chains};

  $parsed{$name}->{ppi} =
    extractProtProtInt( $bin_path, $parsed_pdb_file, $pchains, $out );

  remove_redundant( $parsed{$name}->{ppi}, $nrd_cf, $out ) if( defined $nrd_cf );
}

if( defined $lib_flag ) {
  print $out "\n*** List of PPI(s) found ***\n";
  for my $name (sort keys %parsed) {
    my $ppi = $parsed{$name}->{ppi};
    foreach (@$ppi) {
      print $out "$_->{name} $_->{size}\n";
    }
  }
  close $out;
  exit;
}






#####################################################################
############ Stage 3:  performing interface alignment  ##############
#####################################################################
print $out "\n+++++ Stage 3: Performing protein-protein interface alignment ...\n";
print "+++++ Stage 3: Performing protein-protein interface alignment ...\n" if(defined $out_file);

my $set1 = \@pdb_set1;
my $set2 = \@pdb_set2;
if (scalar @pdb_set2 == 0 ) { $set2 = \@pdb_set1; }  ### all-against-all

my $num1 = scalar @$set1;
my $num2 = scalar @$set2;
my %complete = ();
my $counter  = 1;

for(my $m=0; $m<$num1; $m++) {
  my $name1    = $$set1[$m]->{pdb_file} . '_' . $$set1[$m]->{pchains};
  my $ppi1     = $parsed{$name1}->{ppi};
  my $num_ppi1 = scalar @$ppi1;

  for(my $i=0; $i<$num_ppi1; $i++) {

    for(my $k=0; $k<$num2; $k++) {
      my $name2     = $$set2[$k]->{pdb_file} . '_' . $$set2[$k]->{pchains};
      my $ppi2      = $parsed{$name2}->{ppi};
      my $num_ppi2  = scalar @$ppi2;

      for(my $j=0; $j<$num_ppi2; $j++) {

	my $pair     = $$ppi1[$i]->{name} . '_' . $$ppi2[$j]->{name};
	my $pair_rev = $$ppi2[$j]->{name} . '_' . $$ppi1[$i]->{name};

	next if( exists $complete{$pair_rev} || exists $complete{$pair} );  ### avoid repeat
	next if( $$ppi1[$i]->{name} eq $$ppi2[$j]->{name} and not defined $sel_flag );  ### against itself

	my ($result, $trans_vec, $rot_mat, $parsed_out) = callISalign( $$ppi1[$i], $$ppi2[$j], $switch, $out, $tra_flag );
	if( defined $out_file ) { print ">>>$counter $result\n"; }
	if( defined $vmd ) {
	  my $script_file = "$pair.vmd";
	  my $pdbname1    = $$set1[$m]->{pdb_file};
	  my $pdbname2    = $$set2[$k]->{pdb_file};
	  mk_vmdscript( $script_file, $pdbname1, $pdbname2, $trans_vec, $rot_mat, $parsed_out, $vmd );
	}

	$counter++;
	$complete{$pair} = 1;
	$complete{$pair_rev} = 1;
      } ### j interface 2
    } ### k protein 2

  } ### i interface 1
} ### m protein 1



print $out "iAlign ends at ", timeStamp(), ".\n\n";

close $out;

####------------------------ End of main routine --------------------------------####





#####################################################################################
############################         Subroutine        ##############################
#####################################################################################

#### remove redundant entries
sub remove_redundant {
  my ($ppi, $cutoff, $out) = @_;

  for(my $i=0; $i<scalar @$ppi; $i++) {
    my $name1     = $$ppi[$i]->{name};
    my $int_file1 = $$ppi[$i]->{intfile};
    my $con_file1 = $$ppi[$i]->{confile};

    for(my $j=$i+1; $j<scalar @$ppi;) {
      my $name2     = $$ppi[$j]->{name};
      my $int_file2 = $$ppi[$j]->{intfile};
      my $con_file2 = $$ppi[$j]->{confile};

      my $alnout = `$bin_path/IS-align -t -q 2 $int_file2 $int_file1 $con_file2 $con_file1`;

      ### parsing output of IS-align
      my ($trans_vec, $rot_mat, $score_ln, $parsed_out) = parseISalignOut( $alnout );

      my $tms = 0;
      if( $score_ln =~ /^\w{2}-score =\s*(\d+\.\d+)/ ) { $tms = $1; }

      print $out "$name2 vs $name1 -- TM-score = $tms";
      if( $tms > $cutoff ) {
	splice( @$ppi, $j, 1 );
	print $out " redundant entry $name2 removed\n";
      }
      else {
	print $out "\n";
	$j++;
      }
    }
  }
}



#### perform interface alignment with IS-align
sub callISalign {
  my ($ppi1, $ppi2, $switch, $out, $tra_flag) = @_;

  my $name1     = $ppi1->{name};
  my $int_file1 = $ppi1->{intfile};
  my $con_file1 = $ppi1->{confile};

  my $name2     = $ppi2->{name};
  my $int_file2 = $ppi2->{intfile};
  my $con_file2 = $ppi2->{confile};

  my $pair = "$name1 vs $name2";
  print $out ">>>$name1 vs $name2\n";

  ### call IS-align
  my $alnout = `$bin_path/IS-align $switch $int_file1 $int_file2 $con_file1 $con_file2`;

  ### parsing output of IS-align
  my ($trans_vec, $rot_mat, $score_ln, $parsed_out) = parseISalignOut( $alnout );

  ### outputing results
  print $out "\n$parsed_out\n\n";

  ### transform original pdbfile1 to give the best inteface alignment to pdbfile2
  if( defined $tra_flag ) {
    my $org_pdbfile = $ppi1->{pdbfile};

    my $len1 = length( $name1 );
    my $org_name = substr( $name1, 0, $len1-2 );
    my $tra_pdbfile = "$name1\_to\_$name2.pdb";

    if( transform_pdb( $org_pdbfile, $tra_pdbfile, $rot_mat, $trans_vec ) ) {
      print $out "Transformed $org_name was saved to $tra_pdbfile\n";
    }
  }

  my $result = "$pair: $score_ln";
  return ($result, $trans_vec, $rot_mat, $parsed_out);
}



#### read a list of pdb files ####
sub readSetFile {
  my ($set_file, $pdb_set) = @_;

  open SET, "<$set_file" or die "Error: could not open $set_file\n";
  while(<SET>){
    next if /^#/;
    chomp;
    my @fields   = split(' ',$_);
    my $pdb_file = $fields[0];
    my $pchains  = '';
    if( defined $fields[1] and not ($fields[1] =~ /\W/) ) { $pchains = $fields[1]; }

    push( @$pdb_set, { 'pdb_file'=>$pdb_file, 'pchains'=>$pchains } );
  }
  close SET;
}



######### parse a raw PDB file  ######
sub parsePDBFile {
  my ($parsed_rec, $pdb_path, $work_path, $bin_path, $outfile) = @_;

  my $pdb_file   = $parsed_rec->{pdb_file};
  my $pchain_lst = $parsed_rec->{pchains};

  my ($base, $path, $ext) = fileparse( $pdb_file, qr/\..*/ );
  my $parsed_pdb_file = "$work_path/$base.parsed";

  my @pchains = ();
  unless ( -s $parsed_pdb_file ) {
    if( not (-e $pdb_file) and defined $pdb_path ) {
      $pdb_file = "$pdb_path/$pdb_file";
    }

    unless( -s $pdb_file ) {
      if( -s "$pdb_file.pdb" )    { $pdb_file = "$pdb_file.pdb"; }
      elsif( -s "$pdb_file.ent" ) { $pdb_file = "$pdb_file.ent"; }
      else {
	die "Error: could not find the pdb file $pdb_file\n";
      }
    }

    my $out = `perl $bin_path/process_pdb.pl -pdb $pdb_file -o $parsed_pdb_file -nowarn`;
  }

  @pchains = split(//,$pchain_lst);
  if( scalar @pchains == 0 ) { @pchains = getProtChains( $parsed_pdb_file ); }

  print $outfile "$pdb_file:  found ", scalar @pchains, " protein chains ", join(' ',@pchains), "\n";

  return( $parsed_pdb_file, \@pchains );
}






########### get protein chains from  Calpha atoms ##############
sub getProtChains {
  my $pdbfile = shift @_;

  open PDB, "<$pdbfile" or die "Error: could open file $pdbfile\n";
  my %pchain = ();
  while (my $line = <PDB>) {
    next if ( $line =~ /^REMARK/ );
    last if ( $line =~ /^ENDMDL/ );  #### only consider the first model

    #### PDB format: http://www.wwpdb.org/documentation/format23/sect9.html ####
    if ( $line =~ /^(ATOM  |HETATM)/ ) {
      (my $atomname = substr($line, 12, 4))=~ s/ //g;
      (my $altloc = substr($line, 16, 1));
      (my $resname = substr($line,17,3)) =~ s/ //g;
      (my $chain = substr($line, 21, 1));
      (my $resid = substr($line, 22, 4)) =~ s/ //g;
      (my $icode = substr($line, 26, 1)) =~ s/ //g;

      if( $atomname eq 'CA' && convAA( $resname ) ne 'X' ) {

	if ( $chain eq '' ) { $chain = '_'; }

	unless( exists $pchain{$chain} ) {
	  $pchain{$chain} = 0;
	}

	$pchain{$chain}++;
      }
    }

  }

  close PDB;

  my @protein_chains = ();
  foreach my $chain (sort keys %pchain) {
    if( $pchain{$chain} > $MIN_NUM_RES ) {
      #push( @protein_chains, { 'id'=>$chain, 'len'=>$pchain{$chain} } );
      push( @protein_chains, $chain );
    }
  }

  return @protein_chains;
}



######## convert residues from three-letter to one-letter and vice visa ############
sub convAA {
  my $res = shift @_;

  my %ts=(
     'GLY'=>'G',     'ALA'=>'A',     'VAL'=>'V',     'LEU'=>'L',
     'ILE'=>'I',     'SER'=>'S',     'THR'=>'T',     'CYS'=>'C',
     'MET'=>'M',     'PRO'=>'P',     'ASP'=>'D',     'ASN'=>'N',
     'GLU'=>'E',     'GLN'=>'Q',     'LYS'=>'K',     'ARG'=>'R',
     'HIS'=>'H',     'PHE'=>'F',     'TYR'=>'Y',     'TRP'=>'W',
  );

  if (exists $ts{$res}) {
    return $ts{$res};
  } else {
    return "X";
  }
}



################

sub extractProtProtInt {
  my ($bin_path, $parsed_pdb_file, $pchains, $outfile) = @_;

  my ($base, $path, $ext) = fileparse( $parsed_pdb_file, qr/\..*/ );

  my @ppi = ();
  my $num_pchains = scalar @$pchains;
  for(my $i=0; $i<$num_pchains; $i++) {
    my $pchain_a = $$pchains[$i];
    for(my $j=$i+1; $j<$num_pchains; $j++) {
      my $pchain_b = $$pchains[$j];
      my $pair     = $pchain_a . $pchain_b;
      my $name     = $base . $pair;
      my $int_file = "$path$name\_int.pdb";
      my $con_file = "$path$name\_con.lst";

      my $num_int_res = 0;
      if( -s $int_file and -s $con_file ) {
	my $output = `grep "^Total interacting residues" $con_file`;
	if( $output =~ /Total interacting residues.* (\d+)/ ) {
	  $num_int_res = $1;
	}
	else { die "Error: file $con_file seems corrupted\n"; }
      }
      else {
	my $output = `$bin_path/extint -s $parsed_pdb_file -r $pchain_a -l $pchain_b -d $CONT_CUTOFF -i $int_file -c $con_file -e -ca`;
	if( $output =~ / (\d+) interface residues extracted/ ) {
	  $num_int_res = $1;
	}
	else { die "Error: extraction of interfaces failed\n"; }
      }

      if( $num_int_res > $MIN_INT_RES ) {
        push( @ppi, { 'name'=>$name, 'size'=>$num_int_res, 'pair'=>$pair,
		      'intfile'=>$int_file, 'confile'=>$con_file, 'pdbfile'=>$parsed_pdb_file } );
      }
    }
  }

  printf $outfile "$parsed_pdb_file: found %d valid PPI(s)", scalar @ppi;
  for(my $i=0; $i<scalar @ppi; $i++) {
    print $outfile " $ppi[$i]->{pair} $ppi[$i]->{size}";
  }
  print $outfile "\n";

  return \@ppi;
}


##### get translation vector and rotation matrix from IS-align output
sub parseISalignOut {
  my $alnout = shift @_;

  my @output = split(/\n/, $alnout);

  ### remove useless information for end user
  my $num = scalar @output;
  for(my $i=0; $i < $num; $i++) {
    last if ( $output[0] =~ /^Structure 1/ );
    shift @output;
  }

  ### find measurement of the alignment
  my ($alnlen, $rmsd, $sid, $col, $pvalue, $score);
  my $score_ln;
  my @trans_vec = ();
  my @rot_mat = ();
  for (my $i=0; $i < scalar @output; $i++)  {
    my $line = $output[$i];
    if ( $line =~ /^\w{2}-score =\s*(\d+\.\d+), P-value = \s*(\S+)/ ) {
      $score    = $1;
      $pvalue   = $2;
      $score_ln = $line;
    }
    elsif ( $line =~ /^\w{2}-score =\s*(\d+\.\d+)/ ) {
      $score    = $1;
      $score_ln = $line;
    }
    elsif ( $line =~ /^Number of aligned residues  =\s*(\d+)/ ) {
      $alnlen = $1;
    }
    elsif ( $line =~ /^Number of aligned contacts  =\s*(\d+)/ ) {
      $col = $1;
    }
    elsif ( $line =~ /^RMSD =\s*(\d+\.\d+), Seq identity  =\s*(\d+\.\d+)/ ) {
      $rmsd = $1;
      $sid  = $2;
    }
    elsif( $line =~ /Transformation matrix/ ) {
      for(my $j=0; $j<3; $j++) {
        my @fields = split(' ',$output[$i+$j+2]);
        push( @trans_vec, $fields[1] );
        push( @rot_mat, [ ($fields[2], $fields[3], $fields[4]) ] );
      }
      last;
    }
  }

  my $out = join("\n", @output);
  return( \@trans_vec, \@rot_mat, $score_ln, $out );
}


########### transform a pdb structure given a transformation ##############
sub transform_pdb {
  my ($inp_file, $out_file, $rot_mat, $trans_vec) = @_;


  open PDB, "<$inp_file" or die "Error: could not open $inp_file\n";
  open OUT, ">$out_file" or die "Error: could not write $out_file\n";
  print OUT "REMARK  Generated by iAlign.\n";
  while (my $line = <PDB>) {
    next if ( $line =~ /^REMARK/ );
    last if ( $line =~ /^ENDMDL/ );  #### only consider the first model

    #### PDB format: http://www.wwpdb.org/documentation/format23/sect9.html ####
    if ( $line =~ /^(ATOM  |HETATM)/ ) {
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


sub print_banner {
  my $outfile = shift @_;
  print $outfile "iAlign version $VERSION\n";
}


##### Time stamp
sub timeStamp {
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  my $curr_time = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
    $year+1900,$mon+1,$mday,$hour,$min,$sec;

  return $curr_time;
}


sub mk_vmdscript {
  my ($script_file, $pdb_file1, $pdb_file2, $trans_vec, $rot_mat, $result, $rep ) = @_;

  my @lines = split(/\n/, $result);

  ######## decide the location of pdb files  #######
  if( not (-e $pdb_file1) and defined $pdb_path ) {
    $pdb_file1 = "$pdb_path/$pdb_file1";
  }

  unless( -s $pdb_file1 ) {
    if( -s "$pdb_file1.pdb" )    { $pdb_file1 = "$pdb_file1.pdb"; }
    elsif( -s "$pdb_file1.ent" ) { $pdb_file1 = "$pdb_file1.ent"; }
    else {
      print "Warning: could not find the pdb file $pdb_file1\n";
    }
  }


  if( not (-e $pdb_file2) and defined $pdb_path ) {
    $pdb_file2 = "$pdb_path/$pdb_file2";
  }

  unless( -s $pdb_file2 ) {
    if( -s "$pdb_file2.pdb" )    { $pdb_file2 = "$pdb_file2.pdb"; }
    elsif( -s "$pdb_file2.ent" ) { $pdb_file2 = "$pdb_file2.ent"; }
    else {
      print "Warning: could not find the pdb file $pdb_file1\n";
    }
  }
  ###--------------------------------------------###


  my @int1 = ( {'chain'=>'-'} );
  my @int2 = ( {'chain'=>'-'} );
  my $flag = 0;

  ####### find aligned interfacial residues
  foreach (@lines) {
    if( /^ Index Ch1 Resid1/ ) { $flag = 1; }
    next unless ( $flag );
    next unless( /^\s*\d+ / );

    my @fields = split(' ',$_);
    my $ch1   = $fields[1];
    my $res1  = $fields[2];
    my $ch2   = $fields[4];
    my $res2  = $fields[5];

    if ( $ch1 ne $int1[-1]->{chain} ) {
      push( @int1, { 'chain'=>$ch1, 'lst'=>[], } );
    }
    if ( $ch2 ne $int2[-1]->{chain} ) {
      push( @int2, { 'chain'=>$ch2, 'lst'=>[], } );
    }

    push( @{$int1[-1]->{lst}}, $res1 );
    push( @{$int2[-1]->{lst}}, $res2 );
  }

  my $ch1a = $int1[1]->{chain};
  my $ch1b = $int1[2]->{chain};
  my $ch2a = $int2[1]->{chain};
  my $ch2b = $int2[2]->{chain};

  my $aln1a = join(' ',@{$int1[1]->{lst}});
  my $aln1b = join(' ',@{$int1[2]->{lst}});
  my $aln2a = join(' ',@{$int2[1]->{lst}});
  my $aln2b = join(' ',@{$int2[2]->{lst}});

  my $vmd_rep;
  if( $rep eq 'cartoon' ) {
    $vmd_rep = 'NewCartoon 0.3 6.0 4.5 0';
  }
  else {
    $vmd_rep = 'NewRibbons 0.3 6.0 3.0 0';
  }

  ###############################################################
  open TCL, ">$script_file" or die "Error: could not write $script_file\n";
  print TCL "#!/usr/local/bin/vmd\n";
  print TCL "### VMD script generated by iAlign\n\n";

  ##### script for loading the first molecule
  print TCL "mol new $pdb_file1 type pdb\n";
  print TCL "mol delrep 0 top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 10\n";
  print TCL "mol selection {chain $ch1a and same residue as protein within 4.5 of chain $ch1b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 10\n";
  print TCL "mol selection {chain $ch1a and not same residue as protein within 4.5 of chain $ch1b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 3\n";
  print TCL "mol selection {chain $ch1b and same residue as protein within 4.5 of chain $ch1a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 3\n";
  print TCL "mol selection {chain $ch1b and not same residue as protein within 4.5 of chain $ch1a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  #### show aligned Ca in VDW rep
  print TCL "mol representation VDW 0.6 8.0\n";
  print TCL "mol color ColorID 10\n";
  print TCL "mol selection {chain $ch1a and resid $aln1a and name CA}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation VDW 0.6 8.0\n";
  print TCL "mol color ColorID 3\n";
  print TCL "mol selection {chain $ch1b and resid $aln1b and name CA}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";


  ##### script for transformation
  print TCL "\nset transMat {}\n";
  print TCL "lappend transMat {$$rot_mat[0][0] $$rot_mat[0][1] $$rot_mat[0][2] $$trans_vec[0]}\n";
  print TCL "lappend transMat {$$rot_mat[1][0] $$rot_mat[1][1] $$rot_mat[1][2] $$trans_vec[1]}\n";
  print TCL "lappend transMat {$$rot_mat[2][0] $$rot_mat[2][1] $$rot_mat[2][2] $$trans_vec[2]}\n";
  print TCL "lappend transMat {0.0 0.0 0.0 1.0}\n";
  print TCL "set myMol [atomselect top all]\n";
  print TCL "\$myMol move \$transMat\n\n";


  ##### script for loading the second molecule
  print TCL "mol new $pdb_file2 type pdb\n";
  print TCL "mol delrep 0 top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 0\n";
  print TCL "mol selection {chain $ch2a and same residue as protein within 4.5 of chain $ch2b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 0\n";
  print TCL "mol selection {chain $ch2a and not same residue as protein within 4.5 of chain $ch2b}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 1\n";
  print TCL "mol selection {chain $ch2b and same residue as protein within 4.5 of chain $ch2a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation $vmd_rep\n";
  print TCL "mol color ColorID 1\n";
  print TCL "mol selection {chain $ch2b and not same residue as protein within 4.5 of chain $ch2a}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  #### show aligned Ca in VDW rep
  print TCL "mol representation VDW 0.6 8.0\n";
  print TCL "mol color ColorID 0\n";
  print TCL "mol selection {chain $ch2a and resid $aln2a and name CA}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  print TCL "mol representation VDW 0.6 8.0\n";
  print TCL "mol color ColorID 1\n";
  print TCL "mol selection {chain $ch2b and resid $aln2b and name CA}\n";
  print TCL "mol material BrushedMetal\n";
  print TCL "mol addrep top\n";

  close TCL;
}






















 __END__


###################### PODIATRISTS ##############################

=head1 NAME

iAlign - Structural alignment of protein-protein interfaces.

=head1 DESCRIPTION

iAlign is a tool for structural comparison of protein-protein interfaces (PPIs).
One can align two PPIs, one PPI against a set of PPIs, or one set of PPIs against another set.
Two scoring functions, Interface Similarity score (IS-score) and Template Modeling score (TM-score),
are implemented to evaluate the similarity of two PPIs. By default, the alignment is sequence
order dependent (sequential), but one has an option to allow non-sequential alignment.


=head1 SYNOPSIS

=over 4

=item
ialign.pl [options] F<pdb_file1> chain_list1 F<pdb_file2> chain_list2

=item
ialign.pl [options] F<pdb_file1> F<pdb_file2>

=item
ialign.pl [options] F<pdb_list_file>

=back

Use the -m option to see the full documentation with examples.

Basic Options:

=over 4

=item
-a  <0/1/2>  0 - no alignment printout, 1 - concise, 2 - detailed

=item
-c1 <string> chain list for PDB file1

=item
-c2 <string> chain list for PDB file2

=item
-e  <tm/is>  score selection: tm - TM-score, is - IS-score (default)

=item
-g           extract PPIs from PDBs

=item
-h           print this help message

=item
-l  <file>   PDB list 1

=item
-l2 <file>   PDB list 2

=item
-m           print the manual with examples

=item
-nr <value>  TM-score cutoff for non-redundant entries within a PDB

=item
-n  <string> nomalization method: average, short, long

=item
-o  <file>   output file for saving the result

=item
-p1 <file>   PDB file 1

=item
-p2 <file>   PDB file 2

=item
-p  <path>   path to PDB files

=item
-t           transform PDB1 onto PDB2 according to the optimal alignment

=item
-w  <path>   workpath or path to parsed PDB files

=item
-vmd <style> generate VMD script in ribbon (default) or cartoon styles

=item
-v           print the version number

=back



Advanced Options:

=over 4

=item
-dc   <value> distance cutoff for an interfacial contact, default 4.5 A

=item
-do_self      allow self alignment.

=item
-minp <value> minimum number of residues for a protein chain

=item
-mini <value> minimum number of residues for an interface

=item
-q    <1/2>   speed mode: 1 - normal (default), 2 - fast

=item
-s            allow non-sequential alignment


=back



=head1 OPTIONS

=over 4

=item B<-a, -aln>

Choose a printout format for alignment. 0 - no printout, 1 - concise (sequential mode only), 2 - detailed (default).

=item B<-c1, -chains1>

Specify a list of protein chains in pdb file1, e.g., "AB", "ABC", without any space between chain IDs.
If more than two chains are specified, the program will attempt to extract all pairwise
combinations of these chains. If a list of chains is not provided, all chains
in pdb file1 will be examined for extracting protein-protein interfaces.

=item B<-c2, -chains2>

Specify a list of protein chains in pdb file2. Usage is the same as described in c1.


=item B<-dc>

Specify the distance cutoff (in Angstrom) for interfacial contact. If a pair of heavy atoms of two
residues from two separate proteins has a distance less than the cutoff, the two residues form
an interfacial contact. The default is 4.5 Angstroms. Note that specify a value other than 4.5 may
lead to an inaccurate estimation of statistical significance for TM/IS-scores.

=item B<-do_self>

Allow align an interface against itself.


=item B<-e, -measure>

Select a scoring function for measuring interface similarity: "is" IS-score (default) and "tm" TM-score.

=item B<-g, -mklib>

Parse PDB files and extract protein-protein interfaces, but skip the step of interface
alignment. Use this option to generate a library of protein-protein interfaces.

=item B<-h, -help>

show the brief help information.

=item B<-l, -liblst1>

Specify a file that contains a list of PDB entries. The format of the list is:

=over 12

=item pdbfile1 <chainIDs>

=item pdbfile2 <chainIDs>

=item ...

=back

ChainIDs are optional. If not specified, all protein chains will be examined. One can give full path
in "pdbfile". If not, ialign will search specified pdb files under the directory specified by the
option -pdbpath, or the current directory otherwise.

=item B<-l2, -liblst2>

Specify a list of query PDB files. The format is the same as the option -l. When -l2
is specified, the set file specified by -l1 serves as the list of templates. One can
compare two sets of interfaces by using both -l1 and -l2 options.

=item B<-m, -manual>

Show the manual with examples.

=item B<-minp>

Specify the minimum number of amino acids for a peptide chain to be considered as a valid protein.
Peptide chains with AAs lower than this number will be ignored. The default is 25 AAs.

=item B<-mini>

Specify the minimum number of amino acids for a protein-protein interface to be considered as a valid.
Interfaces with AAs lower than this number will be ignored. The default is 20 AAs.

=item B<-nr <cutoff>>

Remove redundant interfaces within a PDB record. The criterion for redundancy is defined by the
cutoff value for the TM-score.

=item B<-n, -norm>

Select a normalization method for the IS-score/TM-score. By default, iAlign uses the length of the
query interface (the second structure) to normalize its TM/IS-score. Setting an option with -n changes
the length for normalization: "average" - the average of two interfaces, "short" - the shorter interface,
and "long" - the longer interface.

=item B<-o, -output>

Specify the output file for saving alignment results.

=item B<-p1, -pdb1>

Specify the PDB file of structure 1.

=item B<-p2, -pdb2>

Specify the PDB file of structure 2.

=item B<-p, -pdbpath>

Specify the Path to unprocessed PDB files.

=item B<-q, -quick <num>>

Specify the speed mode: 1 - normal (default), 2 - fast, which is about two times faster but less accurate.

=item B<-s, -sequential>

Allow non-sequential alignment. By default, the alignment is sequential (sequence-order-dependent).

=item B<-t, -trans>

Superpose pdb1 onto pdb2 using the transformation matrix that generates the optimal interface
alignment. The output file is named as "pdb1chain1achain1b_to_pdb2chain2achain2b.pdb".

=item B<-w, -workpath>

Specify the path to processed PDB files and input interface files.

=item B<-vmd <style>>

Generate TCL script for visualizations with VMD. Two styles are available: ribbon and cartoon.

=item B<-v, -version>

Show the version number of this release.

=back




=head1 EXAMPLES

=over 4

=item Example 1

Align interfaces from two protein complexes: chains A and B of "foo1.pdb",
and chain C and D of "foo2.pdb". The second command also transforms foo1.pdb according to the optimal alignment.

ialign.pl foo1.pdb AB foo2.pdb CD

ialign.pl -t -p1 foo1.pdb -c1 AB -p2 foo2.pdb -c2 CD

=item Example 2

Scan all interfaces of "foo.pdb" against interfaces listed in goo.lst.
Parsed pdb and extracted interface files are saved to the working directory "/scratch/ialignlib/".

ialign.pl -w /scratch/ialignlib -p /database/pdb/20090101/  -l goo.lst foo.pdb

=item Example 3

All-against-all of a list of interfaces, and save results to the file "results.dat".

ialign.pl -w /scratch/ialignlib -o results.dat goo.lst

If these interfaces has been parsed by ialign and the parsed files are located in /scratch/ialighlib/,
ialign will retrieve the parsed files and avoid re-parsing the pdb files. Otherwise ialign will
look for the PDB files under the current directory, process the files and save extracted interfaces
to "/scratch/ialignlib".

=item Example 4

All-against-all of a list of interfaces against a list of interfaces.

ialign.pl -w /tmp/ialign -o results.dat -l foo1.lst -l2 foo2.lst

=item Example 5

Allow non-sequential alignment between two interfaces, and output a tcl script for VMD cartoon visualization.

ialign.pl -s -vmd cartoon foo1.pdb AB foo2.pdb CD

=item Example 6

Generate a library of protein-protein interfaces by scanning a list specificed in "goo.lst".
Remove redundant interfaces within a PDB record at a TM-score of 0.8. The interfaces files are
saved under "/scratch/ialignlib".

ialign.pl -mklib -nr 0.8 -w /scratch/ialignlib -l goo.lst

=back


=head1 AUTHOR

Mu Gao  <mu.gao@gatech.edu>.

=head1 DATE

10-May-2010

=cut
