#!/usr/bin/perl -w

use strict;

my $fpkm_f = shift;


open IN,"$fpkm_f" || die $!;
my $index = 0;
#分组
my %mapOne = ();

while(<IN>){
   #记录每一列的值是否大于1,否就记为0,是则记为1
   my %mapTwo=();
   #每一组大于1的数量
   my %mapFour=();
   #每一组的总数量
   my %mapFive=();  
   my @t = split /\t/;
   if($index == 0){
     for my $i (1..$#t){
       my $item = $t[$i];
       #取疾病名称的前缀来分组
       $mapOne{$i}=(split(/_/,$item))[0];
     }
     print join "\t",@t;
   }else {
     for my $j (1..$#t){
       my $item = $t[$j];
       if( !exists($mapTwo{$j}) ){
         $mapTwo{$j} = 0;
       }
       #截取小数点前的数值
       my $tempItem = (split(/\./,$item))[0];
       #判断是否大于1
       if( $tempItem >= 1 ){
         if ( $item >1 ){
           $mapTwo{$j} = 1;
         }
       }
     } 
     my @mapTwoKeys = keys %mapTwo;
     for my $k (0..$#mapTwoKeys) {
        my $key = $mapTwoKeys[$k];
        my $group = $mapOne{$key};
        if ( !exists($mapFour{$group}) ) {
           $mapFour{$group} = 0;
        }      
        if( !exists($mapFive{$group}) ){
           $mapFive{$group} = 0;
        }
        $mapFour{$group} = $mapFour{$group} + $mapTwo{$key};
        $mapFive{$group} = $mapFive{$group} + 1;
     }
     my @mapFourKeys = keys %mapFour;
     for my $g (0..$#mapFourKeys){
        my $key = $mapFourKeys[$g];
        my $per = $mapFour{$key}/$mapFive{$key};
        if($per > 0.8){
           #输出结果
           print join "\t",@t;
           #结束循环
           last;
        }
     }
   }

   $index=$index+1;
}

close IN;

