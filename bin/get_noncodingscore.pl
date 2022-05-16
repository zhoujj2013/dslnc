#!/usr/bin/perl -w
use strict;

my $f=shift;
open IN,"$f" || die $!;
my $index = 0;


while(<IN>){
    my @t = split /\t/;
    #my @temp = split "/", $t;
    if($index == 0){
    print "Gene\tDSscore\n";
    }else{
        print "noncoding\t$t[-2]\n";
    }       
$index = $index + 1;
}

close IN;
