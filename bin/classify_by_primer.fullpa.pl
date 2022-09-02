#! /usr/bin/perl
use Getopt::Long;
# use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Cwd qw(abs_path);
sub usage {
    print STDERR<< "USAGE";

Despriprion: BGI version's full-length transcript detection algorithm for BGI patented full-length polyA tail detection library construction protocol.
Author: shizhuoxing\@bgi.com
Data: 20190409
Usage: perl $0 -blastm7 mapped.m7 -ccsfa ccs.fa -umilen 6 -min_primerlen 19 -min_isolen 200 -outdir ./

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
%pp1=();
while(<IN>){
	chomp;
	@a=();@a=split(/\t/,$_);
	if($a[0]!~/^#/ and $a[3]>=$min_primerlen){
		$strand="+" if($a[8]<$a[9]);$strand="-" if($a[8]>$a[9]);
		$pp1{$a[0]}{$a[1]}{$a[6]}="$strand\t$a[3]\t$a[6]\t$a[7]\t$a[8]\t$a[9]" if($strand eq "+" and $a[8]<=5);
		$pp1{$a[0]}{$a[1]}{$a[6]}="$strand\t$a[3]\t$a[6]\t$a[7]\t$a[8]\t$a[9]" if($strand eq "-" and $a[9]<=5);
	}
}
close IN;

open IN,"$ccsfa";
open FL,">$outdir/isoseq_flnc.fasta";
open NFL,">$outdir/isoseq_nfl.fasta";
open POLYA,">$outdir/isoseq_flnc.polyAstat.xls";
while(<IN>){
	chomp;
	if($_=~/^>/){
		$name=$_;
		$name=~s/^>//;
	}else{
		$seqlen{$name}=length($_);
		@a=();@a=split(/\//,$name);$tname="$a[0]\/$a[1]";
		@key1=();@key1=sort keys %{$pp1{$name}{primer_U}};
		@key2=();@key2=sort{$a<=>$b}@key1;
		@key1=();@key1=sort keys %{$pp1{$name}{primer_S}};
		@key3=();@key3=sort{$a<=>$b}@key1;

		$m=0;$mx=@key2;%pp2=();
		foreach $k1(@key2){
			@tmp1=();@tmp1=split(/\t/,$pp1{$name}{primer_U}{$k1});
			$umi="NA";
			foreach $k2(@key3){
				@tmp2=();@tmp2=split(/\t/,$pp1{$name}{primer_S}{$k2});
				if($tmp1[0] eq "+" and $tmp2[0] eq "+" and abs($tmp1[2]-$tmp2[2])<=5){
					if($tmp2[1]>=32){
						$umi=substr($_,$tmp1[3],$umilen);
						$pp2{$m}[0]="p3+";$pp2{$m}[1]=$umi;$pp2{$m}[2]=$tmp2[2];$pp2{$m}[3]=$tmp2[3];
					}else{
						$pp2{$m}[0]="p5+";$pp2{$m}[1]=$umi;$pp2{$m}[2]=$tmp1[2];$pp2{$m}[3]=$tmp1[3];
					}
				}elsif($tmp1[0] eq "-" and $tmp2[0] eq "-" and abs($tmp1[3]-$tmp2[3])<=5){
					if($tmp2[1]>=32){
						$umi=reverse substr($_,$tmp1[2]-$umilen-1,$umilen);
						$umi=~tr/ATCGatcg/TAGCtacg/;
						$pp2{$m}[0]="p3-";$pp2{$m}[1]=$umi;$pp2{$m}[2]=$tmp2[2];$pp2{$m}[3]=$tmp2[3];
					}else{
						$pp2{$m}[0]="p5-";$pp2{$m}[1]=$umi;$pp2{$m}[2]=$tmp1[2];$pp2{$m}[3]=$tmp1[3];
					}
				}
			}
			
			if($m==0){
				$seq=substr($_,0,$pp2{$m}[2]-1) if(defined $pp2{$m}[2]);
				print NFL ">$tname\/0\_$pp2{$m}[2]\_CCS\n$seq\n" if(defined $pp2{$m}[2] and length($seq)>=$min_isolen);
			}elsif($m>0){
				if(($pp2{$m-1}[0] eq "p3+" and $pp2{$m}[0] eq "p5-") or ($pp2{$m-1}[0] eq "p5+" and $pp2{$m}[0] eq "p3-")){
					$seq=substr($_,$pp2{$m-1}[3],$pp2{$m}[2]-$pp2{$m-1}[3]-1) if($pp2{$m-1}[0] eq "p5+");
					$seq=reverse substr($_,$pp2{$m-1}[3],$pp2{$m}[2]-$pp2{$m-1}[3]-1) if($pp2{$m-1}[0] eq "p3+");
					$seq=~tr/ATCGatcg/TAGCtacg/ if($pp2{$m-1}[0] eq "p3+");
					$seq=~s/^G+//;
					$umi=$pp2{$m-1}[1] if($pp2{$m-1}[0] eq "p3+");$umi=$pp2{$m}[1] if($pp2{$m-1}[0] eq "p5+");
					
					$pa=denovo_poly($seq);
					$polyA="";$polyA=substr($seq,$pa,);$polyA=~s/G+$//;
					$Alen=0;$Alen=$polyA=~tr/A/A/;
					$Aratio=0;$Aratio=sprintf("%.2f",$Alen/length($polyA)) if(length($polyA)>0);
					$flseq="";$flseq=substr($seq,0,$pa) if($Aratio>=0.8 and length($polyA)>=3);
					$flseq=$seq if($Aratio<0.8 or length($polyA)<=3);
					$seqlen=length($flseq);
					$polyAlen=0;$polyAlen=length($polyA) if($Aratio>=0.8 and length($polyA)>=3);
					$finalpolya="NA";$finalpolya=$polyA if($Aratio>=0.8 and length($polyA)>=3);
					
					print FL ">$tname\/$pp2{$m-1}[3]\_$pp2{$m}[2]\_CCS\n$flseq\n" if($seqlen>=$min_isolen);
					print POLYA "$tname\/$pp2{$m-1}[3]\_$pp2{$m}[2]\_CCS\t$umi\t$seqlen\t$polyAlen\t$finalpolya\n" if($seqlen>=$min_isolen);
				}else{
					$seq=substr($_,$pp2{$m-1}[3],$pp2{$m}[2]-$pp2{$m-1}[3]-1);
					print NFL ">$tname\/$pp2{$m-1}[3]\_$pp2{$m}[2]\_CCS\n$seq\n" if(length($seq)>=$min_isolen);
				}
			}
			
			if($m==$mx-1){
				$seq=substr($_,$pp2{$m}[3],) if(defined $pp2{$m}[3]);
				print NFL ">$tname\/$pp2{$m}[3]\_$seqlen{$name}\_CCS\n$seq\n" if(defined $pp2{$m}[3] and length($seq)>=$min_isolen);
			}
			$m++;
		}
	}
}
close IN;
close FL;
close NFL;
close POLYA;

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
		
		if($tAratio<70 or $cAratio<50){
			last;
		}
		$mark++;
	}
	$polystart=$seqlen-$x{$mark-1}[1]+1 if($mark>0);
	$polystart=$seqlen if($mark==0);
	return($polystart);
}
