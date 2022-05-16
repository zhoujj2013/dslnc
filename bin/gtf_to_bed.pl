#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Data::Dumper;
use Cwd qw(abs_path);

&usage if @ARGV<1;

sub usage {
        my $usage = << "USAGE";

        This script create makefile for LncFunNet analysis.
        Author: zhoujj2013\@gmail.com 
        Usage: $0 config.cfg

USAGE
print "$usage";
exit(1);
};

my $f=shift;

open IN,"$f" || die $!;
while(<IN>){
	chomp;
	next if(/^#/);
	my @t = split /\t/;
	my $tran_id = $1 if($t[8] =~ /transcript_id "([^;]+)";/);
	if($t[2] eq "transcript"){
		print "$t[0]\t$t[3]\t$t[4]\t$tran_id\t.\t$t[6]\n";
	}
}
close IN;

