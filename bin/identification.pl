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

        This script is designed for mNGS analysis (species identification).
        Author: abehang\@sina.com zhoujj2013\@gmail.com 
        Usage: $0 <config.txt>

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
my %conf;
&load_conf($conf, \%conf);

my $all='all: ';
my $mk;$out/abinit.merged.gtf$out/abinit.merged.gtf
#my $fast5tofastq='';
#my $QC;
#my $kraken;

my $sample_f = abs_path($conf{SAMPLE});
my $out = abs_path($conf{OUTDIR});
my $thread = $conf{THREAD};

my $gtf = abs_path($conf{GTF});

#### write you things ###
#decide to add function in another file
#my %compare;
open IN,"$gtf" || die $!;
while(<IN>){
	chomp;
	next if(/^$/);
	next if(/^#/);
	my @t = split /\t/;
	my $id = shift @t;

	# determine nanopore
	my $type = "nanopore";
	$mk .= "perl \t$conf{FishInWinter.pl} $conf{Chr.lst} $out/lncRNA_identification/abinit.merged.gtf > abinit.merged.cl.gtf &&\\ \n";
    $fast5path = $t[0];
    #(这句话什么意思)
	## Fast5toFastq	start
	make_path "$out/lncRNA_identification";
	$mk .= "01rawdata.origin.finished: $fast5path\n";
    #(?)
	#$fast5tofastq .= "export PATH=\"/ifs/home/sunhang/miniconda3/envs/nanopore/bin\" &&\\ \n";
	$mk .= "\t$conf{STRINGTIE} --merge -G GTF -o abinit.merged.gtf -l SK gtf.lst &&\\ \n";
    $mk .= "perl \t$conf{FishInWinter.pl} $conf{Chr.lst} abinit.merged.gtf > abinit.merged.cl.gtf &&\\ \n";
    $mk .= "awk '$3 ~ /transcript/' abinit.merged.cl.gtf | perl -ne 'chomp; my @t = split /\t/; my $transid = $1 if($t[8] =~ /transcript_id "(\S+)";/); print "$t[0]\t$t[3]\t$t[4]\t$transid\t.\t$t[6]\n";' > abinit.merged.cl.bed &&\\ \n";
    

#lfast5tofastq .= "\t$conf{Rscript} $conf{MinINQC} -i $out/01_rawdata/$id/sequencing_summary.txt -o $out/01_rawdata/$id -p 2 &&\\n"
	$mk .= "\tcat $out/$id/01_rawdata/pass/*.fastq > $out/$id/01_rawdata/$id.raw.fastq &&\\ \n";
	my $rawfastq = "$out/$id/01_rawdata/$id.raw.fastq";
	$mk .= "\t$conf{nanoQC} $rawfastq -o $out/$id/01_rawdata &&\\";
	$mk .= "touch 01rawdata.origin.finished\n\n";
	$all .= "01rawdata.origin.finished ";
	## Fast5toFastq end
	
	## QC start
	make_path "$out/$id/02_QC";
	#my $Trimscore = 7;
	#my $Trimlenth = 1000;
	#my $Trimhead = 70;
	#my $Trimtail = 60;
	$mk .= "02QC.finished: 01rawdata.origin.finished\n"; 
	#$QC .= "export PATH=\"/ifs/home/sunhang/miniconda3/envs/nanopore/bin\" &&\\ \n";
	# $QC .= "conda activate nanopore &&\\ \n";
	$mk .= "\t$conf{NanoPlot} --summary $out/$id/01_rawdata/sequencing_summary.txt --loglength -o $out/$id/02_QC &&\\ \n";
	$mk .= "\t$conf{NanoFilt}  -q $conf{Trimscore} -l $conf{Trimlenth} --headcrop $conf{Trimhead} --tailcrop $conf{Trimtail} $rawfastq | gzip > $out/$id/02_QC/$id.clean.fastq.gz &&\\ \n";
	$mk .= "touch 02QC.finished\n\n";
	$all .= "02QC.finished ";
	my $cleanfastq = "$out/$id/02_QC/$id.clean.fastq.gz";
	
	## QC end ##
		
	##kraken2 start
	make_path "$out/$id/03_kraken";
	# database download
	# https://benlangmead.github.io/aws-indexes/k2
	$mk .= "03kraken.finished: 02QC.finished\n"; 
	$mk .= "\texport KRAKEN2_DB_PATH=\"$conf{KRAKEN2_DB_PATH}:\" && $conf{kraken2} --db $conf{kdb} --threads $thread --report $out/$id/03_kraken/$id.kraken.report --classified-out $out/$id/03_kraken/$id.cseqs.fq --unclassified-out $out/$id/03_kraken/$id.uncseqs.fq --use-names $cleanfastq > $out/$id/03_kraken/$id.spe.count 2> $out/$id/03_kraken/$id.kraken.log &&\\ \n";
	$mk .= "\tcut -f 3 $out/$id/03_kraken/$id.spe.count | sort | uniq -c | sort -k1nr > $out/$id/03_kraken/$id.spe.count.stat && perl $conf{BIN}/convert_whitespace2tab.pl $out/$id/03_kraken/$id.spe.count.stat > $out/$id/03_kraken/$id.spe.count.stat.rf &&\\ \n";
	$mk .= "\t$conf{BRACKEN} -d $conf{KRAKEN2_DB_PATH}/$conf{kdb} -i $out/$id/03_kraken/$id.kraken.report -o $out/$id/03_kraken/$id.bracken.report -w $out/$id/03_kraken/$id.bracken.report -r $conf{READLEN} -l S &&\\ \n";
	$mk .= "\t/ifs/home/sunhang/software/KrakenTools/kreport2mpa.py -r $out/$id/03_kraken/$id.kraken.report -o $out/$id/03_kraken/$id.S.report &&\\ \n";
	$mk .= "touch 03kraken.finished\n\n";
	$all .= "03kraken.finished ";
        ## kraken2 end
	
	##metaphlan start
	make_path "$out/$id/04_metaphlan";
	$mk .= "04metaphlan.finished: 02QC.finished\n";
	$mk .= "\t$conf{METAPHLAN} $cleanfastq --input_type fastq --bowtie2db $conf{METAPHLAN_DB_PATH} -x $conf{mdb} --nproc 6 -o $out/$id/04_metaphlan/$id.metaphlan3.out &&\\ \n";
	$mk .= "touch 04metaphlan.finisedh\n\n";
	$all .= "04metaphlan.finised";
	##metaphlan end
	
	make_path abs_path($conf{OUTDIR});
	open OUT, ">$out/$id/makefile";
	print OUT $all, "\n";
	print OUT $mk, "\n";
	close OUT;
	$all = "all: ";
	$mk = "";
}
close IN;

#make_path "$out/Shell"
#open OUT1, ">$out/01_rawdata/fast5trans.sh";
#print OUT1 "echo 'fast5tofastq start'", "\n";
#print OUT1 $fast5tofastq, "\n";
#print OUT1 "echo 'fast5tofastq finish'", "\n";
#close OUT1;

#open OUT2, ">$out/02_QC/QC.sh";
#print OUT2 "echo 'QC start'", "\n";
#print OUT2 $QC, "\n";
#print OUT2 "echo 'QC finish'", "\n";
#close OUT2;

#open OUT3, ">$out/03_kraken/kraken.sh";
#print OUT3 "echo 'kraken start'", "\n";
#print OUT3 $kraken, "\n";
#print OUT3 "echo 'kraken finish'", "\n";
#close OUT3;
#########################
sub load_conf
{
    my $conf_file=shift;
    my $conf_hash=shift; #hash ref
    open CONF, $conf_file || die "$!";
    while(<CONF>)
    {
        chomp;
        next unless $_ =~ /\S+/;
        next if $_ =~ /^#/;
        warn "$_\n";
        my @F = split"\t", $_;  #key->value
        $conf_hash->{$F[0]} = $F[1];
    }
    close CONF;
}
