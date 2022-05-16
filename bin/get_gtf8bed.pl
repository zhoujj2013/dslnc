#!/usr/bin/perl -w

use strict;

my ($bed_f,$gtf_f) = @ARGV;

my %index;
open IN,"$bed_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$index{$t[3]} = 1;
}
close IN;

open IN,"$gtf_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $trans_id = $1 if($t[8] =~ /transcript_id "(\S+)";/);
	if(exists $index{$trans_id}){
		print join "\t",@t;
		print "\n";
	}
}
close IN;
