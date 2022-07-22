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

        This script designed for disease-specificity analysis.
        Author: 2783313960@qq.com 07/05/2022
        Usage: $0 config.txt

USAGE
print "$usage";
exit(1);
};

my $conf=shift;
my %conf;
&load_conf($conf, \%conf);

my $all='all: ';
my $mk;


my $out = abs_path($conf{OUTDIR});
my $sample_f = abs_path($conf{SAMPLE});
my $thread = $conf{THREAD};
#### write you things ###
open IN,"$sample_f" || die $!;
while(<IN>){
	chomp;
	next if(/^#/);
	my @t = split /\t/;
	my $id = shift @t;
	
	my $type = "PE";
	my ($r1,$r2);
	if(scalar(@t) == 1){
		$type ="SE";
		$r1 = abs_path($t[0]);
		$r2 = $r1;
	}elsif(scalar(@t) == 2){
		$type = "PE";
		$r1 = abs_path($t[0]);
		$r2 = abs_path($t[1]);
	}
	
	#print "type: $type\n";
	my $encode_type = `perl $conf{BIN}/phredDetector.pl $r1`;
	my $filter_cutoff;
	    my $phred = "";
    if($encode_type == 33){
        $filter_cutoff = 38;
        $phred = " ";
    }elsif($encode_type == 64){
        $filter_cutoff = 69;
        $phred = " --phred64 ";
    }
	
	make_path "$out/1.standard/$id/00datafilter";
	$mk .= "00trim.finished: $r1 $r2\n";
	$mk .= "\t$conf{BIN}/trim_adaptor_PE_v3 $r1 $r2 $out/1.standard/$id/00datafilter/$id.Trim $conf{MINLEN} $filter_cutoff > $out/1.standard/$id/00datafilter/$id.trim.log 2>$out/1.standard/$id/00datafilter/$id.trim.err && touch 00trim.finished\n";
	$all .= "00trim.finished ";
	
	$mk .= "00fastqc.finished: $r1 $r2\n";
	$mk .= "\t$conf{FASTQC} --threads $thread -f fastq -o $out/1.standard/$id/00datafilter/ $r1 $r2 > /dev/null 2>/dev/null && touch 00fastqc.finished\n";
	$all .= "00fastqc.finished ";
	
	my $clear_r1 = "$out/1.standard/$id/00datafilter/$id.Trim.R1.fq";
	my $clear_r2 = "$out/1.standard/$id/00datafilter/$id.Trim.R2.fq";
	
		
	$mk .= "00rmdup.finished: 00trim.finished\n";
	$mk .= "\t$conf{BIN}/rmdup $clear_r1 $clear_r2 $out/1.standard/$id/00datafilter/$id.Trim.RD 1>$out/1.standard/$id/00datafilter/$id.rmdup.log 2>$out/1.standard/$id/00datafilter/$id.rmdup.err && touch 00rmdup.finished\n";
	$all .= "00rmdup.finished ";
	$clear_r1 = "$out/1.standard/$id/00datafilter/$id.Trim.RD.R1.fq";
	$clear_r2 = "$out/1.standard/$id/00datafilter/$id.Trim.RD.R2.fq";

	#my $geneset_index = "";
	#if(exists $conf{GENESETINDEX}){
	#	$geneset_index = "--transcriptome-index $conf{GENESETINDEX}";
	#}

	############################################################
	#outWigType          None
	#    string(s): type of signal output, e.g. "bedGraph" OR "bedGraph read1_5p". Requires sorted BAM: --outSAMtype BAM SortedByCoordinate .
	#                    1st word:
	#                    None       ... no signal output
	#                    bedGraph   ... bedGraph format
	#                    wiggle     ... wiggle format
	#                    2nd word:
	#                    read1_5p   ... signal from only 5' of the 1st read, useful for CAGE/RAMPAGE etc
	#                    read2      ... signal from only 2nd read
	#--outWigStrand Unstranded 
	#############################################################
	
	make_path "$out/1.standard/$id/01alignment";
	$mk .= "01align.finished: 00trim.finished 00rmdup.finished\n";
	if($type eq "SE"){
		$mk .= "\trm $clear_r2 && $conf{HISAT2} -x $conf{INDEX} -U $clear_r1 -p $thread -S $out/1.standard/$id/01alignment/$id.sam 1>$out/1.standard/$id/01alignment/$id.hisat2.log 2>$out/1.standard/$id/01alignment/$id.hisat2.err && touch 01align.finished\n";
	}elsif($type eq "PE"){
		#--library-type fr-firststrand
		# add --outSAMstrandField intronMotif for cufflinks, for correcting BAM record error: found spliced alignment without XS attribute
		$mk .= "\t$conf{HISAT2} -x $conf{INDEX} -1 $clear_r1 -2 $clear_r2 -p $thread -S $out/1.standard/$id/01alignment/$id.sam 1>$out/1.standard/$id/01alignment/$id.hisat2.log 2>$out/1.standard/$id/01alignment/$id.hisat2.err && touch 01align.finished\n";
	}
	$all .= "01align.finished ";
	
	#`rm $out/$id/01alignment/$id.bam` if(-f "$out/$id/01alignment/$id.bam");
	$mk .= "01linkbam.finished: 01align.finished\n";
	$mk .= "\t$conf{SAMTOOLS} view --threads $thread -b $out/1.standard/$id/01alignment/$id.sam > $out/1.standard/$id/01alignment/$id.raw.bam && $conf{SAMTOOLS} sort --threads $thread -o  $out/1.standard/$id/01alignment/$id.bam -T $id $out/1.standard/$id/01alignment/$id.raw.bam && $conf{SAMTOOLS} flagstat $out/1.standard/$id/01alignment/$id.bam > $out/1.standard/$id/01alignment/$id.bam.flagstat && touch 01linkbam.finished\n";
	$all .= "01linkbam.finished ";
	
	## quantification
	make_path "$out/1.standard/$id/02assembly";
	$mk .= "02assembly.finished: 01linkbam.finished\n";
	$mk .= "\t$conf{STRINGTIE} -o $out/1.standard/$id/02assembly/$id.stringtie.gtf -p $thread -G $conf{GTF} $out/1.standard/$id/01alignment/$id.bam  >$out/1.standard/$id/02assembly/$id.stringtie.log 2>$out/1.standard/$id/02assembly/$id.stringtie.err && touch 02assembly.finished\n";
	$all .= "02assembly.finished ";

	#my $spe;
	#my $ver;
	#if($conf{SPE} eq "human"){
	#	$spe = "hg";
	#	$ver = "hg19";
	#}elsif($conf{SPE} eq "mouse"){
	#	$spe = "mm";
	#	$ver = "mm9";
	#}
	#
	### visualization 
	#make_path "$out/$id/03visual";
	#$mk .= "03visual.finished: 02quantification.finished\n";
	#$mk .= "\t$conf{SAMTOOLS} view -h $out/$id/01alignment/$id.bam  > $out/$id/03visual/$id.sam && $conf{HOMER}/makeTagDirectory $out/$id/03visual/homer $out/$id/03visual/$id.sam -format sam -genome $ver -checkGC > $out/$id/03visual/makeTagDirectory.log 2>$out/$id/03visual/makeTagDirectory.err && rm $out/$id/03visual/$id.sam && $conf{HOMER}/makeUCSCfile $out/$id/03visual/homer -style rnaseq -o auto -strand both > $out/$id/03visual/makeUCSCfile.log 2>$out/$id/03visual/makeUCSCfile.err && gunzip -d -c $out/$id/03visual/homer/homer.ucsc.bedGraph.gz | grep -v \"^track\" > $out/$id/03visual/$id.bedGraph && $conf{BIN}/bedGraphToBigWig $out/$id/03visual/$id.bedGraph $conf{CHROMSIZE} $out/$id/03visual/$id.bw > $out/$id/03visual/$id.bw.log 2> $out/$id/03visual/$id.bw.err && touch 03visual.finished\n";
	#$all .= "03visual.finished ";
	
	make_path abs_path($conf{OUTDIR});
	open OUT, ">$out/1.standard/$id/makefile";
	print OUT $all, "\n";
	print OUT $mk, "\n";
	close OUT;
	$all= 'all: ';
	$mk = "";
}
close IN;

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


#lncRNA_inditification
make_path abs_path($conf{OUTDIR});
#open OUT, ">$out/1.standard/";
#my $gtf = "ls $out/1.standard/*/02assembly/*.stringtie.gtf > gtf.lst";
#my $run = "cut \-f 1 $conf{SAMPLE} \| while read line\; do echo \"cd \$line \&& make \&& cd \-";done > run.sh";
#open OUT, ">$out/run.sh";
#print OUT $run;
#close OUT;


make_path abs_path($conf{OUTDIR});
my $all1 = 'all: ';
#$all = 'all: ';
#$mk = "";
open OUT, ">$out/1.standard/makefile";
$all1 .= "get_gtf.finished";
my $mk1 .= "get_gtf.finished: $out/1.standard/*/02assembly.finished\n";
$mk1 .= "\tls $out/1.standard/*/02assembly/*.stringtie.gtf > gtf.lst && touch get_gtf.finished\n";
#$mk .= "\tget_gtf.finished";
#$all .= "get_gtf.finished";
print OUT $all1, "\n";
print OUT $mk1, "\n";
close;


make_path "$out/2.identification";
#make_path abs_path($conf{OUTDIR});
open OUT, ">$out/2.identification/makefile";
my $all2 = 'all: ';
#$all2 .= "identification.finished";

my $mk2 = "00stringtie.finished:$out/1.standard/get_gtf.finished\n";
$mk2 .= "\t$conf{STRINGTIE} --merge -G $conf{GTF} -o $out/2.identification/abinit.merged.gtf -l SK $out/1.standard/gtf.lst && touch 00stringtie.finished\n";
$all2 .= "00stringtie.finished ";

$mk2 .= "00gtftobed.finished:00stringtie.finished\n";
$mk2 .= "\tperl $conf{Gtf_to_bed_pl} $out/2.identification/abinit.merged.gtf > $out/2.identification/abinit.merged.bed && touch 00gtftobed.finished\n";
$all2 .= "00gtftobed.finished ";

$mk2 .= "00hg38tohg19.finished:00gtftobed.finished\n";
$mk2 .= "\t$conf{LiftOver} $out/2.identification/abinit.merged.bed $conf{Hg38ToHg19.over.chain} $out/2.identification/abinit.merged.hg19.bed unMapped && touch 00hg38tohg19.finished\n ";
$all2 .= "00hg38tohg19.finished ";



$mk2 .= "00rmknown.finished: 00hg38tohg19.finished\n";
$mk2 .= "\t$conf{Bedtools} intersect -a $out/2.identification/abinit.merged.hg19.bed -b $conf{GTF} -v -s -f 0.1 | egrep -v \"NM|NR\" > $out/2.identification/abinit.merged.cl.novo.bed && touch 00rmknown.finished\n";
$all2 .= "00rmknown.finished ";

$mk2 .= "00changegtf.finished: 00rmknown.finished\n";
$mk2 .= "\tperl $conf{Get_gtf_pl} $out/2.identification/abinit.merged.cl.novo.bed $out/2.identification/abinit.merged.cl.gtf > $out/2.identification/abinit.merged.cl.novo.gtf && touch 00changegtf.finished\n";
$all2 .= "00changegtf.finished ";

$mk2 .= "01iseerna.finished: 00changegtf.finished\n";
$mk2 .= "\t$conf{ISeeRNA} -c $conf{Iseerna_conf} -i $out/2.identification/abinit.merged.cl.novo.gtf -o $out/2.identification/abinit.iseerna.novo.out && touch 01iseerna.finished\n";
$all2 .= "01iseerna.finished ";

$mk2 .= "01make.finished: 01iseerna.finished\n";
$mk2 .= "\tcd $out/2.identification/abinit.iseerna.novo.out && make && touch 01make.finished\n";
$all2 .= "01make.finished ";

$mk2 .= "01getnovel.finished: 01make.finished\n";
#$mk2 .= "\tawk -F \'\\t\' \'{if((\$\$3>0.8)) print \$\$0}\' $out/2.identification/abinit.iseerna.novo.out/abinit.merged.cl.novo.result > $out/2.identification/abinit.iseerna.novo.out/abinit.merged.cl.novo.result.0.8.filtered && touch 01getnovel.finished\n";
$mk2 .= "\tperl $conf{Codingpotential_filter_pl} $out/2.identification/abinit.iseerna.novo.out/abinit.merged.cl.novo.result > $out/2.identification/abinit.iseerna.novo.out/abinit.merged.cl.novo.result.0.8.filtered && touch 01getnovel.finished\n";
$all2 .= "01getnovel.finished ";

$mk2 .= "01getnovelgtf.finished: 01getnovel.finished\n";
$mk2 .= "\tperl $conf{New_get_gtf_pl} $out/2.identification/abinit.iseerna.novo.out/abinit.merged.cl.novo.result.0.8.filtered $out/2.identification/abinit.merged.cl.gtf > $out/2.identification/abinit.new.gtf && touch 01getnovelgtf.finished\n";
$all2 .= "01getnovelgtf.finished ";

#接着这里
$mk2 .= "02changebed.finished: 01getnovelgtf.finished\n";
$mk2 .= "\tcd ..\n";
$mk2 .= "\tperl $conf{Generate_final_gtf_pl} $out/2.identification/abinit.new.gtf > $out/2.identification/novo.trans.gtf 2> $out/2.identification/novo.trans.bed  && touch 02changebed.finished\n";
$all2 .= "02changebed.finished ";

$mk2 .= "02rmIntergeniclnc.finished: 02changebed.finished\n";
$mk2 .= "\t$conf{Bedtools} intersect -a $out/2.identification/novo.trans.bed -b $conf{RefGene_bed} -v > $out/2.identification/novo.trans.non-ol.bed && touch 02rmIntergeniclnc.finished \n";
$all2 .= "02rmIntergeniclnc.finished ";

$mk2 .= "02getgtf.finished: 02rmIntergeniclnc.finished\n";
$mk2 .= "\tperl $conf{Get_gtf_pl} $out/2.identification/novo.trans.non-ol.bed $out/2.identification/abinit.new.gtf > $out/2.identification/novo.trans.non-ol.gtf && touch 02getgtf.finished\n";
$all2 .= "02getgtf.finished ";

$mk2 .= "02getfinalgtf.finished: 02getgtf.finished\n";
$mk2 .= "\tcat $conf{GTF} $out/2.identification/novo.trans.non-ol.gtf > $out/2.identification/all.final.gtf && touch 02getfinalgtf.finished\n";
$all2 .= "02getfinalgtf.finished ";

$mk2 .= "02getknownlnc.finished: 02getfinalgtf.finished\n";
$mk2 .= "\tless $conf{RefGene_anno} | grep lncRNA |cut -f 3 | sed \'s\/,\/\\n\/g\' | sort | uniq > $out/2.identification/lncRNA.trans.lst && touch 02getknownlnc.finished\n";
$all2 .= "02getknownlnc.finished ";

$mk2 .= "02getknownlncgtf.finished: 02getknownlnc.finished\n";
$mk2 .= "\tperl $conf{Get_known_lncrna_pl} $out/2.identification/lncRNA.trans.lst $out/2.identification/all.final.gtf > $out/2.identification/lncRNA.trans.gtf && touch 02getknownlncgtf.finished\n";
$all2 .= "02getknownlncgtf.finished ";

#$mk2 .= "02getnovellnc.finished: 02getknownlncgtf.finished\n";
#$mk2 .= "\tgrep \"SK\\.\" $out/2.identification/all.final.gtf | awk \'\$\$3 ~ \/transcript\/\' | awk \'{print \$\$10}\' | sed \'s\/\\\"\/\/g\' |  sed \'s/;\/\/g\' | sort | uniq > $out/2.identification/novo.lncrna.lst && touch 02getnovellnc.finished\n";
#$all2 .= "02getnovellnc.finished ";

$mk2 .= "02getnovellncgtf.finished:02getknownlncgtf.finished\n";
$mk2 .= "\tgrep \"SK\\.\" $out/2.identification/all.final.gtf > novo.rna.gtf && touch 02getnovellncgtf.finished\n";
$all2 .= "02getnovellncgtf.finished ";

$mk2 .= "02getmergegtf.finished:02getnovellncgtf.finished\n";
$mk2 .= "\tcat $out/2.identification/lncRNA.trans.gtf $out/2.identification/novo.rna.gtf > $out/2.identification/all.lncRNA.gtf && touch 02getmergegtf.finished\n";
$all2 .= "02getmergegtf.finished ";

$mk2 .= "02getalllnclst.finished:02getmergegtf.finished\n";
#$mk2 .= "\tperl -ne \'chomp; my \@t = split \/\\t/; my \$\$geneid = \$\$1 if(\$\$t[8] =~ /gene_id \"([^\"]+)\";/); print \"\$\$geneid\\n\";' $out/2.identification/all.lncRNA.gtf | sort |uniq > $out/2.identification/all.final.lncgeneid.lst && touch 02getalllnclst.finished\n";
$mk2 .= "\tperl $conf{Gtf_to_lnclst_pl} $out/2.identification/all.lncRNA.gtf | sort | uniq > $out/2.identification/all.final.lncgeneid.lst && touch 02getalllnclst.finished\n";
$all2 .= "02getalllnclst.finished ";

print OUT $all2, "\n";
print OUT $mk2, "\n";
close OUT;

#lncRNA quatification and filtering
make_path "$out/3.filter";
open OUT, ">$out/3.filter/makefile";
my $all3 = 'all: ';
#$all3 .= "filter.fnished";

#my $mk3 = "00getbam.finished: $out/2.identification/02getmergegtf.finished\n";
#$mk3 .= "\tbam \= \$\$(ls $out/1.standard/*/01alignment/*.bam | grep -v raw | xargs) && touch 00getbam.finished\n";
#my $bam = $out/standard/*/01alignment/*.bam | grep -v raw | xargs  
#$all3 .= "00getbam.finished ";

#bam=$(ls */01alignment/*.bam | grep -v raw | xargs)
my $mk3 = "00featurecount.finished: $out/2.identification/02getmergegtf.finished\n";
$mk3 .= "\t$conf{FeatureCounts} -T 6 -p -t exon -g gene_id -a $out/2.identification/all.final.gtf -o $out/3.filter/count.txt `ls $out/1.standard/*/01alignment/*.bam | grep -v raw | xargs` && touch 00featurecount.finished\n";
#$mk3 .= "\t$conf{FeatureCounts} -T 6 -p -t exon -g gene_id -a $out/2.identification/all.final.gtf -o $out/3.filter/count.txt \$bam && touch 00featurecount.fnished\n";
$all3 .= "00featurecount.finished ";

$mk3 .= "01rmheader.finished: 00featurecount.finished\n";
$mk3 .= "\tsed -i \'1d\' $out/3.filter/count.txt && touch 01rmheader.finished\n";
$all3 .= "01rmheader.finished ";

$mk3 .= "01getsampleid.finished: 01rmheader.finished\n";
$mk3 .= "\tperl $conf{Get_sampleid_pl} $out/3.filter/count.txt > $out/3.filter/count1.txt && touch 01getsampleid.finished\n";
$all3 .= "01getsampleid.finished ";

$mk3 .= "01counttofpkm.finished: 01getsampleid.finished\n";
$mk3 .= "\t$conf{Rscript} $conf{Count_to_fpkm_R} && touch 01counttofpkm.finished\n";
$all3 .= "01counttofpkm.finished ";

$mk3 .= "01getheader.finished: 01counttofpkm.finished\n";
$mk3 .= "\thead -1 $out/3.filter/fpkm_result.txt > $out/3.filter/header.txt && touch 01getheader.finished\n";
$all3 .= "01getheader.finished ";

$mk3 .= "02getalllncfpkm.finished: 01getheader.finished\n";
$mk3 .= "\tperl $conf{FishInWinter_pl} $out/2.identification/all.final.lncgeneid.lst $out/3.filter/fpkm_result.txt | cat $out/3.filter/header.txt - > $out/3.filter/all.lncRNA.fpkm.txt && touch 02getalllncfpkm.finished\n";
$all3 .= "02getalllncfpkm.finished ";

$mk3 .= "02lncfilter.finished: 02getalllncfpkm.finished\n";
$mk3 .= "\tperl $conf{Filter_lncrna_pl} $out/3.filter/all.lncRNA.fpkm.txt > $out/3.filter/all.lncRNA.filter.fpkm.txt && touch 02lncfilter.finished\n";
$all3 .= "02lncfilter.finished ";

$mk3 .= "03getpcglst.finished: 02lncfilter.finished\n";
$mk3 .= "\tless $conf{RefGene_anno} | grep \"coding gene\" | cut -f 1 | sed \'s\/,\/\\n\/g\' | sort | uniq > $out/3.filter/refgene.codingene.lst && touch 03getpcglst.finished\n";
$all3 .=  "03getpcglst.finished ";

$mk3 .= "03getpcgfpkm.finished: 03getpcglst.finished\n";
$mk3 .= "\tperl $conf{FishInWinter_pl} $out/3.filter/refgene.codingene.lst $out/3.filter/fpkm_result.txt | cat $out/3.filter/header.txt - > $out/3.filter/pcg.fpkm.txt && touch 03getpcgfpkm.finished\n";
$all3 .=  "03getpcgfpkm.finished ";


$mk3 .= "03pcgfilterfpkm.finished: 03getpcgfpkm.finished\n";
$mk3 .= "\tperl $conf{Filter_lncrna_pl} $out/3.filter/pcg.fpkm.txt > $out/3.filter/all.pcg.filter.fpkm.txt && touch 03pcgfilterfpkm.finished\n";
$all3 .=  "03pcgfilterfpkm.finished ";
print OUT $all3,"\n";
print OUT $mk3, "\n";
close OUT;

#lncRNA DS score
make_path "$out/4.ds_score";
open OUT, ">$out/4.ds_score/makefile";
my $all4 = 'all: ';
#my $all4 = "ds.fnished\n";
my $mk4 = "01getlncavg.finished: $out/3.filter/03pcgfilterfpkm.finished\n";
$mk4 .= "\tperl $conf{Get_average_pl} $out/3.filter/all.lncRNA.filter.fpkm.txt > $out/4.ds_score/all.lncRNA.filter.avg.fpkm.txt && touch 01getlncavg.finished\n";
$all4 .= "01getlncavg.finished ";

$mk4 .= "01lncds.finished: 01getlncavg.finished\n";
$mk4 .= "\tperl $conf{Ds_score_pl} $out/4.ds_score/all.lncRNA.filter.avg.fpkm.txt 2 n > $out/4.ds_score/all.lncRNA.filter.dsscore.txt && touch 01lncds.finished\n";
$all4 .= "01lncds.finished ";

$mk4 .= "02getpcgavg.finished: 01lncds.finished\n";
$mk4 .= "\tperl $conf{Get_average_pl} $out/3.filter/all.pcg.filter.fpkm.txt > $out/4.ds_score/all.pcg.filter.avg.fpkm.txt  && touch 02getpcgavg.finished\n";
$all4 .= "02getpcgavg.finished ";

$mk4 .= "02pcgds.finished: 02getpcgavg.finished\n";
$mk4 .= "\tperl $conf{Ds_score_pl} $out/4.ds_score/all.pcg.filter.avg.fpkm.txt 2 n > $out/4.ds_score/all.pcg.filter.dsscore.txt && touch 02pcgds.finished\n";
$all4 .= "02pcgds.finished ";

$mk4 .= "03getlncscore.finished: 02pcgds.finished\n";
$mk4 .= "\tperl $conf{Get_noncodingscore_pl} $out/4.ds_score/all.lncRNA.filter.dsscore.txt > $out/4.ds_score/noncodingscore.txt && touch 03getlncscore.finished\n";
$all4 .= "03getlncscore.finished ";

$mk4 .= "03getpcgscore.finished: 03getlncscore.finished\n";
$mk4 .= "\tperl $conf{Get_codingscore_pl} $out/4.ds_score/all.pcg.filter.dsscore.txt > $out/4.ds_score/codingscore.txt && touch 03getpcgscore.finished\n";
$all4 .= "03getpcgscore.finished ";

$mk4 .= "03rmheader.finished: 03getpcgscore.finished\n";
$mk4 .= "\t\tsed -i \'1d\' $out/4.ds_score/noncodingscore.txt && touch 03rmheader.finished\n";
$all4 .= "03rmheader.finished ";

$mk4 .= "03getallscore.finished: 03rmheader.finished\n";
$mk4 .= "\tcat $out/4.ds_score/codingscore.txt $out/4.ds_score/noncodingscore.txt > $out/4.ds_score/all.gene.filter.score.txt && touch 03getallscore.finished\n";
$all4 .= "03getallscore.finished ";

$mk4 .= "03destiny.finished: 03getallscore.finished\n";
$mk4 .= "\t$conf{Rscript} $conf{Destiny_R} && touch 03destiny.finished\n";
$all4 .= "03destiny.finished "; 

print OUT $all4,"\n";
print OUT $mk4,"\n";
close OUT;
