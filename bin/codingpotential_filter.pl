#!/usr/bin/perl -w
use strict;
use Math::BigFloat;
#my  $i = new Math::BigFloat '1.931533e-01';
my $f=shift;
open IN,"$f" || die $!;


my $index = 0;

while(<IN>){
    #chomp;
    #next if(/^#/);
    my @t = split /\t/;
    my  $i = new Math::BigFloat $t[2];
    #print join "\t",@t;
    #print "$i\n";
    if($index > 0){
        if($i > 0.8){
        print join "\t",@t;
        #print $t[2];
        }
    }
$index=$index+1;
}
close IN;
