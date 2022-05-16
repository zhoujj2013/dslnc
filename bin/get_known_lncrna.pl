#!/usr/bin/perl -w

use strict;

my $kl_f = shift;
my $f = shift;

my %kn;

open IN,"$kl_f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	$kn{$t[0]} = 1;
}
close IN;


open IN,"$f" || die $!;
while(<IN>){
	chomp;
	my @t = split /\t/;
	my $trans_id = $1 if($t[8] =~ /transcript_id "([^"]+)";/);
	if(exists $kn{$trans_id}){
		print join "\t",@t;
		print "\n";
	}
}
close IN;
