#! /usr/bin/perl
use Getopt::Long;
use FindBin '$Bin';
use List::Util qw/max min sum/;

open IN,"$ARGV[0]";
@len=();
while (<IN>){
	chomp;
	@a=();@a=split(/\t/,$_);
	push @len,$a[1];
}
close IN;

# $Rscript ||= "$Bin/Rscript";
$Rscript ||= "Rscript";
$convert ||= "$Bin/convert";
 
open R,">$ARGV[1]/plot_subreads.R";
my $script=<<RSCRIPT;
library(ggplot2)
library(gridExtra)
library(grid)
rm(list=ls())
a<-read.csv("$ARGV[1]/subreads.len",header=TRUE,sep="\\t")
pdf("$ARGV[1]/SubReads.length.pdf",height=6,width=8)
ggplot(a,aes(x=Length))+
	geom_histogram(colour="white",fill="dodgerblue3",binwidth=150)+
	labs(title="SubReads Length",x="Length",y="Number")+
	xlim(0,9000)+
	theme(plot.title=element_text(face='bold',size=15),axis.text=element_text(color='black'),panel.background = element_rect(fill='transparent'),panel.grid=element_line(color='grey'),panel.border=element_rect(fill='transparent',color='black'),axis.title=element_text(size=15))
dev.off()
RSCRIPT
print R $script;
close R;

# `$Rscript $ARGV[1]/plot_subreads.R && $convert -density 100 $ARGV[1]/plot_subreads.pdf $ARGV[1]/plot_subreads.png`;
`Rscript $ARGV[1]/plot_subreads.R`;

$tseq=@len;$tlen=sum(@len);$Nbase=sprintf("%.2f",sum(@len)/1000000000);$maxlen=max(@len);$meanlen=sprintf("%.2f",($tlen/$tseq));$minlen=min(@len);
open OUT,">$ARGV[1]/SubReads.stat.xls";
print OUT "ReadsType\tTotal Reads\tTotal Base(GB)\tMaxLength(bp)\tMeanLength(bp)\tN50 Length(bp)\n";
print OUT "SubReads\t$tseq\t$Nbase\t$maxlen\t$meanlen\t";
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
