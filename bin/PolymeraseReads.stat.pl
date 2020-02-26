#! /usr/bin/perl
use Getopt::Long;
use FindBin '$Bin';
use List::Util qw/max min sum/;

open IN,"$ARGV[0]";
<IN>;
%polymer=();
while (<IN>){
	chomp;
	@a=();@a=split(/\t/,$_);
	@b=();@b=split(/\//,$a[0]);
	@c=();@c=split(/\_/,$b[2]);
        $polymer{$b[1]}[0]=$c[0] if(!defined $polymer{$b[1]});
        $polymer{$b[1]}[1]=$c[1] if(!defined $polymer{$b[1]});
	$polymer{$b[1]}[0]=$c[0] if($c[0]<$polymer{$b[1]}[0]);
	$polymer{$b[1]}[1]=$c[1] if($c[1]>$polymer{$b[1]}[1]);
}
close IN;

@key=keys %polymer;
@len=();
open OUT,">$ARGV[1]/polymerreads.len";
print OUT "ZMWS\tLength\n";
foreach $k(@key){
	$plen=$polymer{$k}[1]-$polymer{$k}[0];
	print OUT "$k\t$plen\n";
	push @len,$plen;
}
close OUT;

# $Rscript ||= "$Bin/Rscript";
$Rscript ||= "Rscript";
# $convert ||= "$Bin/convert";
 
open R,">$ARGV[1]/plot_polymerreads.R";
my $script=<<RSCRIPT;
library(ggplot2)
library(gridExtra)
library(grid)
rm(list=ls())
a<-read.csv("$ARGV[1]/polymerreads.len",header=TRUE,sep="\\t")
pdf("$ARGV[1]/PolymeraseReads.length.pdf",height=6,width=8)
ggplot(a,aes(x=Length))+
	geom_histogram(colour="white", fill="dodgerblue3",binwidth=1200)+
	labs(title="Polymerase Reads Length", x="Length", y="Number")+
	xlim(0,90000)+
	theme(plot.title=element_text(face='bold',size=15),axis.text=element_text(color='black'),panel.background = element_rect(fill='transparent'),panel.grid=element_line(color='grey'),panel.border=element_rect(fill='transparent',color='black'),axis.title=element_text(size=15))
dev.off()
RSCRIPT
print R $script;
close R;

# `$Rscript $ARGV[1]/plot_polymerreads.R && $convert -density 100 $ARGV[1]/plot_polymerreads.pdf $ARGV[1]/plot_polymerreads.png`;
`Rscript $ARGV[1]/plot_polymerreads.R && rm $ARGV[1]/plot_polymerreads.R`;

$tseq=@len;$tlen=sum(@len);$Nbase=sprintf("%.2f",sum(@len)/1000000000);$maxlen=max(@len);$meanlen=sprintf("%.2f",($tlen/$tseq));$minlen=min(@len);
open OUT,">$ARGV[1]/PolymeraseReads.stat.xls";
print OUT "ReadsType\tTotal Reads\tTotal Base(GB)\tMaxLength(bp)\tMeanLength(bp)\tN50 Length(bp)\n";
print OUT "PolymeraseReads\t$tseq\t$Nbase\t$maxlen\t$meanlen\t";
@len1=sort{$b<=>$a}@len;
$tmplen=0;
foreach $k(@len1){
	$tmplen+=$k;
	if($tmplen>=$tlen/2){
		print OUT "$k\n";
		exit;
	}
}
close OUT;
