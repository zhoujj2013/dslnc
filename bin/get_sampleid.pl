#!/usr/bin/perl -w
use strict;

my $f=shift;
open IN,"$f" || die $!;
my $index = 0;


while(<IN>){
    my @t = split /\t/;
    #my @temp = split "/", $t;
    if($index == 0){
        for my $i(1..$#t){
            if($i >= 6){
                my @temp = split "/", $t[$i];
                my $final = $temp[-1];
                $final = (split(/\./,$final))[0];
                $t[$i] = $final;
            }
        }
    }
$index = $index + 1;
print join "\t",@t;
print "\n";
}

close IN;
