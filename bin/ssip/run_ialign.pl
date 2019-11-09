#!/usr/bin/perl -w
use strict;
use File::Basename;
use Cwd 'abs_path';

#####HELP MESSAGE############################
if (@ARGV<3)
{
  die "run_align.pl DIMER_PDB_FILE IS-SCORE_CUTOFF ALIGNMENT_OUTPUT_FILE\n";
}
#############################################

my $maxLenth=16;
my @a = (0..9,'a'..'z','A'..'Z');

my $path=dirname(abs_path(__FILE__));

my $randStr = join '', map { $a[int rand @a] } 0..($maxLenth-1);
my $randName = "X_".$randStr.".pdb";
my $tmpResult = "X_".$randStr.".ial";


my $queryDimer=$ARGV[0];
my $is_cutoff=$ARGV[1];
my $outputFile=$ARGV[2];

my $ialign_type="is";
my $libList="$path/NIL/dimer.list";
my $libPath="$path/NIL";
my $IALIGN="$path/ialign/ialign.pl";
my $parseFile="X_".$randStr.".parsed";

#####PARAMETERS FROM LAUNCH INTERFACE########
my $mini='7';
my $minp='7';
my $dist_cutoff=4.5;
#############################################

`cp $queryDimer $randName`;

print "$IALIGN -a 2 -n ave -s -minp $minp -mini $mini -dc $dist_cutoff -e $ialign_type -w $libPath -l $libList $randName >$tmpResult\n";
`$IALIGN -a 2 -n ave -s -minp $minp -mini $mini -dc $dist_cutoff -e $ialign_type -w $libPath -l $libList $randName >$tmpResult`;


my $chainAB="XX";

open(IALTMP, "<$tmpResult");
while(my $line = <IALTMP>)
{
  if ( $line =~ /$parseFile/ )
  {
    chomp(my $ifNum=(split /[\s+,]/, $line)[2]);
    if($ifNum == 0)
    {
      die "no interface for ".$queryDimer;
    }
    elsif($ifNum > 1)
    {
      die "too many interface for ".$queryDimer;
    }
    chomp($chainAB=(split /[\s+,]/, $line)[5]);
    last;
  }
}

#print "$chainAB\n";
close IALTMP;

my $chainA=substr($chainAB,0,1);
my $chainB=substr($chainAB,1,1);


my $intFile=$libPath."/X_".$randStr.$chainAB."_int.pdb";
my $conFile=$libPath."/X_".$randStr.$chainAB."_con.lst";
if (! -s $intFile){
  $intFile="X_".$randStr.$chainAB."_int.pdb";
}
if (!-s $conFile){
  $conFile="X_".$randStr.$chainAB."_con.lst";
}
open (INTF, "<$intFile") || die "can't open intface pdb file";

my @resAArray;
my @resBArray;
my $resAIndex=0;
my $resBIndex=0;
my $querySeqA="";
my $querySeqB="";
while(my $line = <INTF>)
{
  if( $line =~ /ATOM/ )
  {
    my $ch = substr($line,21,1);
    my $resID = substr($line,22,5);
    chomp($resID);
    if( $ch eq $chainA)
    {
      $resAArray[$resAIndex]=$resID;
      my $type=changeAA(substr($line,17,3));
      $querySeqA=$querySeqA.$type;
      $resAIndex++;
    }
    if( $ch eq $chainB)
    {
      $resBArray[$resBIndex]=$resID;
      my $type=changeAA(substr($line,17,3));
                        $querySeqB=$querySeqB.$type;
      $resBIndex++;
    }
  }
}

my $ifResNum=$resAIndex + $resBIndex;

#print "$querySeqA\n$querySeqB\n";

open (OUT, ">$outputFile");

print OUT "Interface residues:\n";
print OUT "chain_ID_A: ".$chainA."\n";
print OUT "residue_A: ";

foreach my $resID (@resAArray)
{
  print OUT "$resID ";
}
print OUT "\n";

print OUT "chain_ID_B: ".$chainB."\n";
print OUT "residue_B: ";
foreach my $resID (@resBArray)
{
        print OUT "$resID ";
}
print OUT "\n";

print OUT "Alignments:\n";

print OUT "PDB   chA chB ifResNum alnResNum score ifSeqA ifSeqB\n";

printf OUT ("query  $chainA   $chainB  %3d  %3d  %5.3f  %s %s\n",$ifResNum,$ifResNum,1.000,$querySeqA,$querySeqB);


my $libIndex=0;
my $score=0.0;
my $libSeqA;
my $libSeqB;
my $libIfNum;
my $libAlnNum;
my $libChA;
my $libChB;
my $resAlignIndex=0;
my @libResA;
my @libResB;
my $libPDBID;

open(IALTMP, "<$tmpResult");
while(my $line = <IALTMP>)
{
  if( $line =~ /Structure 1:/)
  {
    $libChA=substr($line,41,1);
    $libChB=substr($line,43,1);
    $libPDBID=substr($line,13,6);
    $libIfNum=substr($line,46,3);
    next;
  }
  if( $line =~ /aligned residues/)
  {
    $libAlnNum=substr($line,30,3);
  }
  if( ($line =~ /IS-score/) || ($line =~ /TM-score/) )
  {
    $libIndex = $libIndex+1;
    $libSeqA="";
    $libSeqB="";
    $resAlignIndex=0;
    for(my $i=0;$i<$resAIndex;$i++)
    {
      $libResA[$i]='-';
    }
    for(my $i=0;$i<$resBIndex;$i++)
                {
                        $libResB[$i]='-';
                }

    chomp($score=(split /[\s+,]/, $line)[2]);
#    print "$score\n";
    next;
  }
        if($score < $is_cutoff)
        {
                next;
        }

  if( $line =~ / Index /)
  {
    $resAlignIndex=1;
    next;
  }
  if( $line =~ / Ansgtrom /)
  {
    $resAlignIndex=0;
  
    printf OUT ("%s $libChA   $libChB  %s  %s  %5.3f  ", $libPDBID, $libIfNum , $libAlnNum, $score);  
    for(my $i=0;$i<$resAIndex;$i++)
    {
      print OUT $libResA[$i];
    }
    print OUT " ";
    for(my $i=0;$i<$resBIndex;$i++)
                {
                        print OUT $libResB[$i];
                }
    print OUT "\n";
    next;
  }
  if($resAlignIndex > 0)
  {
    my $quResID=substr($line,29,5);
    my $libResType=changeAA(substr($line,18,3));
    if( $chainA eq substr($line,27,1))
    {
      my $idx=searchIndex($quResID, $resAIndex, @resAArray);
      $libChA=substr($line,9,1);
      $libResA[$idx]=$libResType;
    }
    if ( $chainB eq substr($line,27,1))
    {
      my $idx=searchIndex($quResID, $resBIndex, @resBArray);
      $libChB=substr($line,9,1);
      $libResB[$idx]=$libResType;
    }
  }  

}

sub changeAA {
  my $triName=$_[0];
  if($triName eq "ALA"){return 'A';}
  elsif($triName eq "CYS"){return 'C';}
  elsif($triName eq "ASP"){return 'D';}
  elsif($triName eq "GLU"){return 'E';}
  elsif($triName eq "PHE"){return 'F';}
  elsif($triName eq "GLY"){return 'G';}
  elsif($triName eq "HIS"){return 'H';}
  elsif($triName eq "ILE"){return 'I';}
  elsif($triName eq "LYS"){return 'K';}
  elsif($triName eq "LEU"){return 'L';}
  elsif($triName eq "MET"){return 'M';}
  elsif($triName eq "ASN"){return 'N';}
  elsif($triName eq "PRO"){return 'P';}
  elsif($triName eq "GLN"){return 'Q';}
  elsif($triName eq "ARG"){return 'R';}
  elsif($triName eq "SER"){return 'S';}
  elsif($triName eq "THR"){return 'T';}
  elsif($triName eq "VAL"){return 'V';}
  elsif($triName eq "TRP"){return 'W';}
  elsif($triName eq "TYR"){return 'Y';}
}

sub searchIndex {
  my $query=$_[0];
  my $num=$_[1];
  for(my $i=0;$i<$num;$i++)
  {
    if($query eq $_[$i+2])
    {  
      return $i;
    }
  }
  die "can't find residue ".$query;
}


`rm $randName`;
`rm $tmpResult`;
`rm $intFile`;
`rm $libPath/X_$randStr*`;
