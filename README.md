Disease-specificity analysis using RNA-Seq (Workflow)
================

##### Liu Huiting (2202/05/07)

### Introduction

Not only does RNAseq have the ability to analyze differences in gene expression between samples, but can discover novel lncRNAs and genes expressed specifically. The dslnc.pl was aimed to identify novel lncRNAs and find disease-specific lncRNAs(including known and novel lncRNAs) in an organ system.

### Requirement

+ [FASTQ](Version 0.11.8)
+ [HISAT2] (v2.0.5)
+ [STRINGTIE] (v1.3.1 )
+ [SAMTOOLS] (v0.1.19)
+ [FeatureCounts] (v2.0.1)
+ [Bedtools] (https://github.com/arq5x/bedtools2)
+ [ISeeRNA] (http://www.myogenesisdb.org/iSeeRNA)
+ [Rscript] (r-3.6.0)
+ [GTF] (hg19) or (https://gitee.com/hui-tingzi/test/blob/master/data/refGene.rar)

### test data

```bash
Disease1: https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR8052751/SRR8052751.1
Disease2: https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-15/SRR8052708/SRR8052708.1
```

### Getting started

#### A.Input data: download and extract the relevant files  and write the samples.lst(the format should look something like this)

```bash
CTRL_002_healthy	../data/srr/SRR8052751_1.fastq	../data/srr/SRR8052751_2.fastq
AD_004_lesional	../data/srr/SRR8052708_1.fastq	../data/srr/SRR8052708_2.fastq
```

#### B.Edit configuration at bin/config.txt
Except for the software that needs to be installed(you should set the environment variable configuration), all the required scripts have been placed in https://gitee.com/hui-tingzi/test/blob/master/bin/, which can be downloaded and used directly.
Write the configure file, please refer to the format,please refer to the format here(https://gitee.com/hui-tingzi/test/blob/master/bin/config.txt).

```bash
OUTDIR  ./out
SAMPLE  ./samples.lst

#       parameter
READLEN 150
MINLEN  60
THREAD  24

#STANDTYPE      FR/FF/RF/RR
#       pro
BIN     bin
FASTQC  fastqc
HISAT2  hisat2
STRINGTIE       stringtie
SAMTOOLS        samtools
HOMER   /homer/bin

#       for     STAR
GTF     ../data/refGene.gtf
SPE     human
INDEX   ../data/hisat2_index/hg19
CHROMSIZE       ../data/hg19/hg19.chrom.size

#       for     dsanalysis
Chr     ../data/chr.lst
Bedtools        bedtools
ISeeRNA iSeeRNA
Iseerna_conf    ../data/iseerna.conf
FeatureCounts   featureCounts
Rscript Rscript
```

#### C.Make the appropriate directories as shown. 

```bash
├── bin
│   ├── auto_run.sh
│   ├── chr.lst
│   ├── codingpotential_filter.pl
│   ├── config.txt
│   ├── count_to_fpkm.R
│   ├── destiny.R
│   ├── dslnc.pl
│   ├── ds_score.pl
│   ├── fishInWinter.pl
│   ├── generate_final_gtf.pl
│   ├── get_average.pl
│   ├── get_codingscore.pl
│   ├── get_gtf8bed.pl
│   ├── get_known_lncrna.pl
│   ├── get_noncodingscore.pl
│   ├── get_sampleid.pl
│   ├── gtf_to_bed.pl
│   ├── gtf_to_lnclst.pl
│   ├── identification.pl
│   ├── iseerna.conf
│   ├── new.filter_lncrna.pl
│   ├── old_get_average.pl
│   ├── refGene.anno
│   └── refGene.bed
├── data
│   ├── chr.lst
│   ├── hg19.chrom.sizes
│   ├── hisat2_index
│   ├── iseerna.conf
│   ├── refGene.anno
│   ├── refGene.bed
│   └── refGene.gtf
└── test
    └── samples.lst
```


------------------------------------------------------------------------

### Workflow

##### Example data: If you would like to use example data for practicing the workflow, you could download human RNAseq data as below.

```bash
Disease1: https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos3/sra-pub-run-20/SRR8052751/SRR8052751.1
Disease2: https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-15/SRR8052708/SRR8052708.1
```

### Quick start

 We can issue all the commands step by step or use quick start as below. The result is the same.

#### Command

```bash
Quick Start
perl ../bin/dslnc.pl ../bin/config.txt 
sh ../bin/auto_run.sh
```

------------------------------------------------------------------------


### Step 1. RNAseq analysis

#### Command

```bash
perl ../bin/dslnc.pl ../bin/config.txt
sh out/1.standard/standard.sh
cd out/1.standard && make && cd -
```

#### Output

The first command will  generated a folder named out, which contain 4 subfolders , namely 1.standard, 2.identification, 3.filter, 4.ds_score desperately. There will be  a makefile in each subfolders that stores detailed commands. The second command will  automate the RNASeq workflow  for each sample one by one, and it might take a long time depending on your computer. The third command assumes RNAseq analysis  has been processesed correctly. A gtf.lst for gtf of each sample will be generated.

```bash
├── out
│   ├── 1.standard
│   │   ├── CTRL_001_health
│   │   │   ├── 00datafilter
│   │   │   ├── 01alignment
│   │   │   ├── 02assembly
│   │   │   └── makefile
│   │   ├── makefile
│   │   └── PCA_002_lesional
│   │       ├── 00datafilter
│   │       ├── 01alignment
│   │       ├── 02assembly
│   │       └── makefile
│   ├── 2.identification
│   │   └── makefile
│   ├── 3.filter
│   │   └── makefile
│   └── 4.ds_score
│       └── makefile
└── samples.lst
```
------------------------------------------------------------------------

### Step 2. Novel lncRNAs identification

#### Command

```bash
cd out/2.identification && make && cd -
```

#### Output

This run  will generate novel lncRNAs, which can be seen in test/out/2.identification/novo.lncrna.lst. This run will also generate a merged gtf for all the samples, which can be seen in test/out/2.identification/all.final.gtf.

```bash
├── 2.identification
│   ├── 00changebed.finished
│   ├── 00changegtf.finished
│   ├── 00getchr.finished
│   ├── 00rmknown.finished
│   ├── 00stringtie.finished
│   ├── 01getnovel.finished
│   ├── 01getnovelgtf.finished
│   ├── 01iseerna.finished
│   ├── 02changebed.finished
│   ├── 02getalllnclst.finished
│   ├── 02getfinalgtf.finished
│   ├── 02getgtf.finished
│   ├── 02getknownlnc.finished
│   ├── 02getknownlncgtf.finished
│   ├── 02getmergegtf.finished
│   ├── 02getnovellncgtf.finished
│   ├── 02rmIntergeniclnc.finished
│   ├── abinit.iseerna.novo.out
│   │   ├── 01make.finished
│   │   ├── abinit.merged.cl.novo.fa
│   │   ├── abinit.merged.cl.novo.feature
│   │   ├── abinit.merged.cl.novo.gtf
│   │   ├── abinit.merged.cl.novo.list
│   │   ├── abinit.merged.cl.novo.orf
│   │   ├── abinit.merged.cl.novo.predict
│   │   ├── abinit.merged.cl.novo.result
│   │   ├── abinit.merged.cl.novo.result.0.8.filtered
│   │   ├── abinit.merged.cl.novo.scaled
│   │   ├── abinit.merged.cl.novo.svm
│   │   ├── all_abinit.merged.cl.novo.consv
│   │   ├── consv
│   │   │   ├── chr10.consv
│   │   │   ├── chr10.dat
│   │   │   ├── chr10.list
│   │   │   ├── chr11.consv
│   │   │   ├── chr11.dat
│   │   │   ├── chr11.list
│   │   │   ├── chr12.consv
│   │   │   ├── chr12.dat
│   │   │   ├── chr12.list
│   │   │   ├── chr13.consv
│   │   │   ├── chr13.dat
│   │   │   ├── chr13.list
│   │   │   ├── chr14.consv
│   │   │   ├── chr14.dat
│   │   │   ├── chr14.list
│   │   │   ├── chr15.consv
│   │   │   ├── chr15.dat
│   │   │   ├── chr15.list
│   │   │   ├── chr16.consv
│   │   │   ├── chr16.dat
│   │   │   ├── chr16.list
│   │   │   ├── chr17.consv
│   │   │   ├── chr17.dat
│   │   │   ├── chr17.list
│   │   │   ├── chr18.consv
│   │   │   ├── chr18.dat
│   │   │   ├── chr18.list
│   │   │   ├── chr19.consv
│   │   │   ├── chr19.dat
│   │   │   ├── chr19.list
│   │   │   ├── chr1.consv
│   │   │   ├── chr1.dat
│   │   │   ├── chr1.list
│   │   │   ├── chr20.consv
│   │   │   ├── chr20.dat
│   │   │   ├── chr20.list
│   │   │   ├── chr21.consv
│   │   │   ├── chr21.dat
│   │   │   ├── chr21.list
│   │   │   ├── chr22.consv
│   │   │   ├── chr22.dat
│   │   │   ├── chr22.list
│   │   │   ├── chr2.consv
│   │   │   ├── chr2.dat
│   │   │   ├── chr2.list
│   │   │   ├── chr3.consv
│   │   │   ├── chr3.dat
│   │   │   ├── chr3.list
│   │   │   ├── chr4.consv
│   │   │   ├── chr4.dat
│   │   │   ├── chr4.list
│   │   │   ├── chr5.consv
│   │   │   ├── chr5.dat
│   │   │   ├── chr5.list
│   │   │   ├── chr6.consv
│   │   │   ├── chr6.dat
│   │   │   ├── chr6.list
│   │   │   ├── chr7.consv
│   │   │   ├── chr7.dat
│   │   │   ├── chr7.list
│   │   │   ├── chr8.consv
│   │   │   ├── chr8.dat
│   │   │   ├── chr8.list
│   │   │   ├── chr9.consv
│   │   │   ├── chr9.dat
│   │   │   ├── chr9.list
│   │   │   ├── chrM.consv
│   │   │   ├── chrM.dat
│   │   │   ├── chrM.list
│   │   │   ├── chrX.consv
│   │   │   ├── chrX.dat
│   │   │   ├── chrX.list
│   │   │   ├── chrY.consv
│   │   │   ├── chrY.dat
│   │   │   ├── chrY.list
│   │   │   └── Makefile
│   │   ├── iSeeRNA.conf
│   │   └── Makefile
│   ├── abinit.merged.cl.bed
│   ├── abinit.merged.cl.gtf
│   ├── abinit.merged.cl.novo.bed
│   ├── abinit.merged.cl.novo.gtf
│   ├── abinit.merged.gtf
│   ├── abinit.new.gtf
│   ├── all.final.gtf
│   ├── all.final.lncgeneid.lst
│   ├── all.lncRNA.gtf
│   ├── lncRNA.trans.gtf
│   ├── lncRNA.trans.lst
│   ├── makefile
│   ├── novo.rna.gtf
│   ├── novo.trans.bed
│   ├── novo.trans.gtf
│   ├── novo.trans.non-ol.bed
│   └── novo.trans.non-ol.gtf
```

------------------------------------------------------------------------

### Step 3. Expression filtering

#### Command

```bash
cd out/3.filter && make && cd -
```
#### Output

This run will get the expression levels of all genes (/test/out/3.filter/fpkm_result.txt) and the lncRNA and protein coding genes with low expression levels can be filtered out(/test/out/3.filter/all.lncRNA.filter.fpkm.txt and /test/out/3.filter/all.pcg.filter.fpkm.txt).

```bash
├── 3.filter
│   ├── 00featurecount.finished
│   ├── 01counttofpkm.finished
│   ├── 01getheader.finished
│   ├── 01getsampleid.finished
│   ├── 01rmheader.finished
│   ├── 02filter.finished
│   ├── 02getallfpkm.finished
│   ├── 02getalllncfpkm.finished
│   ├── 02lncfilter.finished
│   ├── 03getpcgfpkm.finished
│   ├── 03getpcglst.finished
│   ├── 03pcgfilterfpkm.finished
│   ├── all.lncRNA.filter.fpkm.txt
│   ├── all.lncRNA.fpkm.txt
│   ├── all.pcg.filter.fpkm.txt
│   ├── count.txt
│   ├── count.txt.summary
│   ├── fpkm_result.txt
│   ├── header.txt
│   ├── makefile
│   ├── pcg.fpkm.txt
│   └── refgene.codingene.lst
```

------------------------------------------------------------------------

### Step 4. Calculation of disease-specificity score

#### Command

```bash
cd out/4.ds_score && make && cd -
```

#### Output

The disease-specificity score of lncRNAs and protein coding genes will store scanning result in test/out/4.ds_score/all.lncRNA.filter.dsscore.txt and test/out/4.ds_score/all.pcg.filter.dsscore.txt desperately.The destiny.pdf will generated to show the comparion between lncRNA and protein coding genes.

```bash
└── 4.ds_score
    ├── 01getavg.finished
    ├── 01getlncavg.finished
    ├── 01lncds.finished
    ├── 02ds.finished
    ├── 02getpcgavg.finished
    ├── 02pcgds.finished
    ├── 03destiny.finished
    ├── 03getallscore.finished
    ├── 03getlncscore.finished
    ├── 03getpcgscore.finished
    ├── 03rmheader.finished
    ├── all.gene.filter.score.txt
    ├── all.lncRNA.filter.avg.fpkm.txt
    ├── all.lncRNA.filter.dsscore.txt
    ├── all.pcg.filter.avg.fpkm.txt
    ├── all.pcg.filter.dsscore.txt
    ├── codingscore.txt
    ├── destiny.pdf
    ├── makefile
    ├── noncodingscore.txt
    └── Rplots.pdf
```


The rusult of disease-specifity score.
```bash
#{id}tab{disease1}tab{disease2}tab{DS_scor _of_gene}tab{most_specific_disease}
LOC100133091            0.851892614542727               2.8187855877881 0.896895005858782	disease1	
LOC101928126            3.02669783420868                1.16114685818554        0.86991908642307        disease2
RALY-AS1                1.06407345733899                1.59204120009034        0.773022319339166	disease1
AS-PTPRE                2.0092520953125         1.42828260748092        0.75973106098206        disease2
```


This destiny plot is of DS score of lncRNAs (green) compared with protein coding genes (pink).
![destiny plot](README_files/destiny.png)


------------------------------------------------------------------------

### Citations:

When you use this pipwline, please site Liu et al. 2022 The long noncoding RNA atlas in healthy and diseased skin(manuscript in preparation).


