#!/usr/bin/perl -w

use strict;

my $fpkm_f = shift;

open IN,"$fpkm_f" || die $!;
#第几行,第一行index=0
my $index = 0;
#分组
my %mapOne = ();

while(<IN>){
   my @currentRow = split /\t/;
   my %groups = ();
   if($index == 0){
     print $currentRow[0];
     for my $number (1..$#currentRow){
       my $item = $currentRow[$number];
       #取疾病名称的前缀来分组
       $mapOne{$number}=(split(/_/,$item))[0];
     }
   }else {
	 for my $number (1..$#currentRow){
		#把疾病前缀拿出来
		my $groupName = $mapOne{$number};
		#第$j列的数值
		my $lineData = $currentRow[$number];
        #是否满足条件,不满足为0,满足为1
		my $isMeetCondition = 0;
		if($lineData>1){
		    $isMeetCondition = 1;
		}
		if ( !exists($groups{$groupName}) ){
		   $groups{$groupName} = new GroupData($groupName, 1, $lineData, $isMeetCondition);
		}else{
		   my $groupData = $groups{$groupName};
		   my $currentTotalCount = $groupData->getTotalCount();
		   $groupData->setTotalCount($currentTotalCount+1);
		   
		   my $currentTotalSum = $groupData->getTotalSum();
		   $groupData->setTotalSum($currentTotalSum+$lineData);
		   
		   my $meetConditionCount = $groupData->getMeetConditionCount();
		   $groupData->setMeetConditionCount($meetConditionCount+$isMeetCondition);
		}
	 }  
	 
     my @groupNames = keys %groups;
	 if($index == 1){
		print "\t";
		print join "\t",(values @groupNames);
		print "\n";
	 }
	 print $currentRow[0];
	 for my $number (0..$#groupNames) {
	    my $groupName = $groupNames[$number];
		print "\t";
		#打印平均值
		print $groups{$groupName}->getTotalSum() / $groups{$groupName}->getTotalCount();
		#print "\t";
	 }
	 print "\n";
   }
   $index=$index+1;
}

close IN;


#定义一个对象
package GroupData;
sub new {
  my $class = shift;
  my $self = {
    #疾病分组名
    _groupName => shift,
    #该组疾病总样本量
    _totalCount => shift,
    #该组疾病样本总和
    _totalSum => shift,
    #该组疾病满足条件的样本数量,比如大于1的样本数量
    _meetConditionCount => shift,
  };
  bless $self, $class;
  return $self;
}
sub setGroupName {
	my ( $self, $groupName ) = @_;
	$self->{_groupName} = $groupName if defined($groupName);
	return $self->{_groupName};
}

sub getGroupName {
	my( $self ) = @_;
	return $self->{_groupName};
}

sub setTotalCount {
	my ( $self, $totalCount ) = @_;
	$self->{_totalCount} = $totalCount if defined($totalCount);
	return $self->{_totalCount};
}

sub getTotalCount {
	my( $self ) = @_;
	return $self->{_totalCount};
}

sub setTotalSum {
	my ( $self, $totalSum ) = @_;
	$self->{_totalSum} = $totalSum if defined($totalSum);
	return $self->{_totalSum};
}

sub getTotalSum {
	my( $self ) = @_;
	return $self->{_totalSum};
}

sub setMeetConditionCount {
	my ( $self, $meetConditionCount ) = @_;
	$self->{_meetConditionCount} = $meetConditionCount if defined($meetConditionCount);
	return $self->{_meetConditionCount};
}

sub getMeetConditionCount {
	my( $self ) = @_;
	return $self->{_meetConditionCount};
}
1;
