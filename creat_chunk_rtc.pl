#! /usr/bin/perl
use File::Basename;
use FindBin '$Bin';
use Cwd qw(abs_path);

$chunk=$ARGV[0];
$outdir=$ARGV[1];
open IN,"$Bin/resolved-tool-contract.json";
$mark=0;
while (<IN>){
	$mark++;
	print "$_" if ($mark!=10 and $mark!=32);
	print "            \"$chunk\"\n" if ($mark==10);
	print "            \"$outdir\/ccs.consensusreadset.xml\"\n" if ($mark==32);
}
close IN;
