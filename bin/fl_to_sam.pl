#! /usr/bin/perl

open IN,"$ARGV[0]";
%zmw=();
while(<IN>){
	chomp;
	@a=();@a=split(/\t/,$_);
	@b=();@b=split(/\//,$a[0]);
	$zmw{$b[1]}=$_;
}
close IN;

open IN,"$ARGV[1]";
while(<IN>){
	chomp;
	if($_=~/^>/){
		$name=$_;
		$name=~s/^>//;
		@b=();@b=split(/\//,$name);
	}else{
		@a=();@a=split(/\t/,$zmw{$b[1]});
		$seqlen=length($_);
		print "$name\t4\t*\t0\t255\t*\t*\t0\t0\t$_\t";
		foreach ($i=0;$i<$seqlen;$i++){
			print "!";
		}
		print "\t$a[11]\t$a[12]\t$a[13]\t$a[14]\t$a[15]\t$a[16]\t$a[17]\n";
	}
}
close IN;
