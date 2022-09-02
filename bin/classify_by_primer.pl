#! /usr/bin/perl
use Getopt::Long;
# use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Cwd qw(abs_path);
sub usage {
    print STDERR<< "USAGE";

Despriprion: BGI version's full-length transcript detection algorithm for PacBio official IsoSeq library construction protocol and BGI patented multi-isoforms in one ZMW library construction protocol.
Author: shizhuoxing\@bgi.com
Data: 20190419
Usage: perl $0 -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 8 -min_primerlen 15 -min_isolen 200 -outdir ./

Options:
	-blastm7*:		result of primer blast to ccs.fa in blast -outfmt 7 format
	-ccsfa*:		the ccs.fa you want to classify to get full-length transcript
	-umilen*:		the UMI length in your library, if set to 0 means nonUMI for library construction
	-min_primerlen*:	the minimum primer alignment length in ccs.fa
	-min_isolen*:		the minimum output's transcript length whithout polyA tail
	-outdir*:		output directory
	-help:			print this help
USAGE
    exit 0;
}

GetOptions(
	"blastm7:s" => \$blastm7,
	"ccsfa:s" => \$ccsfa,
	"umilen:s" => \$umilen,
	"min_primerlen:s" => \$min_primerlen,
	"min_isolen:s" => \$min_isolen,
	"outdir:s" => \$outdir,
	"help:s" => \$help
);
die &usage() if ((!defined $blastm7) or (!defined $ccsfa) or (!defined $umilen) or (!defined $min_primerlen) or (!defined $min_isolen) or (!defined $outdir) or (defined $help));

open IN,"$blastm7";
%m7=();
while(<IN>){
	chomp;
	@a=();@a=split(/\t/,$_);
	if($_=~/^m/ and $a[3]>=$min_primerlen){
		$m7{$a[0]}{$a[1]}{$a[6]}[0]=$a[3];
		$m7{$a[0]}{$a[1]}{$a[6]}[1]=$a[7];
		$m7{$a[0]}{$a[1]}{$a[6]}[2]=$a[8];
		$m7{$a[0]}{$a[1]}{$a[6]}[3]=$a[9];
		$m7{$a[0]}{$a[1]}{$a[6]}[4]="+" if($a[8]<$a[9]);
		$m7{$a[0]}{$a[1]}{$a[6]}[4]="-" if($a[8]>$a[9]);
	}
}
close IN;

open IN,"$ccsfa";
open FL,">$outdir/isoseq_flnc.fasta";
open NFL,">$outdir/isoseq_nfl.fasta";
open XLS,">$outdir/isoseq_flnc.polyAtail.xls";
# print XLS "SeqName\tUMI\tSeqLength\tPolyALength\tPolyAtail\n";
while(<IN>){
	chomp;
	if($_=~/^>/){
		$name=$_;
		$name=~s/^>//;
	}else{
		$seqlen{$name}=length($_);
		@a=();@a=split(/\//,$name);
		@key=();@key=keys %{$m7{$name}{primer_S}};
		@key2=();@key2=sort{$a<=>$b}@key;
		%primer=();
		$tname=$name;$tname=~s/ccs//;
		%sign=();$mm=0;$cc=@key2;
		foreach $k(@key2){
			$sign{$mm}=0;

			$primer{$mm}[0]=$m7{$name}{primer_S}{$k}[4];
			$tmpseq="N";$pas="N";
			if($m7{$name}{primer_S}{$k}[4] eq "+"){
				$tmpseq=substr($_,0,$k-9);
			}elsif($m7{$name}{primer_S}{$k}[4] eq "-"){
				$kt=$m7{$name}{primer_S}{$k}[1]+8;
				$tmpseq1=substr($_,$kt,) if($kt<$seqlen{$name});
				$tmpseq=reverse($tmpseq1);
				$tmpseq=~tr/ACGTacgt/TGCAtgca/;
			}
			$pa=denovo_poly($tmpseq) if($tmpseq ne "N");
			$pas=substr($tmpseq,$pa) if($tmpseq ne "N");
			$Alen=0;$Alen=$pas=~tr/A/A/;
			$primer{$mm}[1]="p3" if($Alen>=5);
			$primer{$mm}[1]="p5" if($Alen<5);
			$primer{$mm}[2]=length($pas);

			if($mm==0){
				$ff1=substr($_,0,$k);
				print NFL ">$tname" if($k>=$min_isolen);
				print NFL "0\_$k\_CCS\n$ff1\n" if($k>=$min_isolen);
			}elsif($mm>0){
				$p3p5="$primer{$mm-1}[1]$primer{$mm-1}[0]$primer{$mm}[1]$primer{$mm}[0]";
				
				$ff1=substr($_,$m7{$name}{primer_S}{$key2[$mm-1]}[1],$k-$m7{$name}{primer_S}{$key2[$mm-1]}[1]-1);
				if($p3p5 eq "p3-p5+" or $p3p5 eq "p5-p3+" and $sign{$mm-1}==0){
					if($p3p5=~/^p3/){
						$tmp=reverse($ff1);
						$tmp=~tr/ACGTacgt/TGCAtgca/;
						$umi=substr($tmp,-8);$polya=substr($tmp,-$primer{$mm-1}[2]-8,$primer{$mm-1}[2]);
						$fl=substr($tmp,0,-$primer{$mm-1}[2]-8);
					}elsif($p3p5=~/^p5/){
						$tmp=$ff1;
						$umi=substr($ff1,-8,8);$polya=substr($ff1,-$primer{$mm}[2]-8,$primer{$mm}[2]);
						$fl=substr($ff1,0,-$primer{$mm}[2]-8);
					}
					$fl=~s/^G+//;
					$fllen=length($fl);$palen=length($polya);
					print FL ">$tname$m7{$name}{primer_S}{$key2[$mm-1]}[1]\_$k\_CCS\n$fl\n" if($fllen>=$min_isolen);
					print XLS "$tname$m7{$name}{primer_S}{$key2[$mm-1]}[1]\_$k\_CCS\t$umi\t$fllen\t$palen\t$polya\n" if($fllen>=$min_isolen);
					$sign{$mm}=1;
				}else{
					print NFL ">$tname$m7{$name}{primer_S}{$key2[$mm-1]}[1]\_$k\_CCS\n$ff1\n" if(($k-$m7{$name}{primer_S}{$key2[$mm-1]}[1])>=$min_isolen);
				}

			}
			if($mm==$cc-1){
				$ff1=substr($_,$m7{$name}{primer_S}{$key2[$mm]}[1],);
				$ff1len=length($ff1);
				print NFL ">$tname$m7{$name}{primer_S}{$key2[$mm]}[1]\_$seqlen\_CCS\n$ff1\n" if($ff1len>=$min_isolen);
			}
			$mm++;
		}
	}
}
close IN;
close FL;
close NFL;
close XLS;

sub denovo_poly{
	$revseq=reverse($_[0]);$seqlen=length($revseq);
	$mark=0;%x=();
	while($revseq=~/(A+)/ig){
		$len=length $&;
		$pos=pos($revseq)-$len+1;
		$first=$pos if($mark==0);
		$end=$pos+$len;
		
		$x{$mark}[0]=$pos;$x{$mark}[1]=$end;
		$tail=substr($revseq,0,$end-1);
		$A=$tail=~tr/A/A/;$tAratio=0;$tAratio=sprintf("%.2f",$A/($end-$first)*100) if($end-1>0);
		
		$chunklen=$end-$x{$mark-1}[1] if($mark>0);
		$chunklen=$len if($mark==0);
		$chunk=substr($revseq,$x{$mark-1}[1]-1,$chunklen) if($mark>0);
		$chunk=substr($revseq,$pos-1,$chunklen) if($mark==0);
		$A=$chunk=~tr/A/A/;$cAratio=0;$cAratio=sprintf("%.2f",$A/$chunklen*100) if($chunklen>0);
		
		if($tAratio<80 or $cAratio<60){
			last;
		}
		$mark++;
	}
	$polystart=$seqlen-$x{$mark-1}[1]+1 if($mark>0);
	$polystart=$seqlen if($mark==0);
	return($polystart);
}
