# fnumts
test
Finding NUMTs from long-read sequencing data.

Author: Chia-Hsuan Fan

Email: fanc@uab.edu

Draft date:

## Quick Start
```sh
git clone https://github.com/chiahsuanF/fnumts.git
cd fnumts/
```
```
./fnumts_main.py -b <index.bam> [option]
```

```sh
# Finding NUMTs with indexed bam. The reference using human hg38.
./fnumts_main.py -b <index.bam>

# Example usage:
./fnumts_main.py -b index_human_hg38.bam
```

```sh
# Finding NUMTs with indexed bam. The reference using other species' chrM indexed by BWA.
./fnumts_main.py -b <index.bam> -r <index_ref_chrM.fa>

# Example usage:
./fnumts_main.py -b index.bam -r path/directory/index_ref_chrM.fa
```

```sh
# Finding NUMTs with indexed bam. The reference using other species' chrM.fa indexed by BWA.
# Adjust distance, length, and support reads.
./fnumts_main.py -b <index.bam> -d <int> -l <int> -s <int> -r <index_ref_chrM.fa>

# Example usage:
./fnumts_main.py -b index.bam -d 20 -l 20 -s 10 -r path/directory/index_ref_chrM.fa
```

## Description
fnumts is a tool for finding insertion mitochondria DNA in nuclear genome. fnumts can detect NUMT length from 1 bp up to whole mitochondira DNA.      

The input file should be a long-read indexed BAM(PacBio Hifi and Oxford Nanopore) file. The input reference chrM needs to be same as aligned bam file reference, and needs to be duplicated and indexed by BWA. The output is a standard VCF file which passed filter. This program was test on a x86_64 Linux system with 40GB physical memory.

## Depencency
Dependencies for fnumts:

* python3
* pysam(version 0.22.0)
* argparse(version 1.4.0)
* pandas(version 2.1.1)
* BWA(version 0.7.17)
* bsalign(version 1.2.1)
* wtdbg2(version 2.5)

  
## Installation

### Conda
```
conda create -n fn python=3
conda activate fn

conda install bioconda::bwa
conda install bioconda::snpeff
conda install bioconda::bialign
pip install pysam
pip install argparse
pip install pandas

```
### Git
```
git clone https://github.com/chiahsuanF/fnumts.git
```
Next, add the directory to PATH:
```
export PATH=$PWD/fnumts/:$PATH
```
## General usage
| Option | Option | Description| Default|
|:---------|:-------------|:-------------|:-------------|
|-t |--type|Import datatype, hifi or ont|hifi|
|-b |--bam|Import aligned indexed BAM file|False|
|-c |--cram|Import aligned indexed CRAM file|False|
|-r |--reference|Import reference chromosome mitochondria. The reference needs to be indexed by BWA|userpath/fnumts/ref/hg38_chrM.fa|
|-o |--outpath|Output directory |./fnumts|
||--base_dir |Base directory for BAM file |../|
|-d |--distance|Maximun number of distance for merging intervals|30|
|-l |--length|Minmiun number of lenght for insertion sequence|10|
|-s |--supportread|Minimun number of support reads for insertion reads|3|

|-sr|Import reference for snpEff annotation|hg38|
## Use cases
fnumts required input as an aligned BAM file with header. Mapping raw data suggests using minimap2 and samtools. Input reference only needs reference chrM.fa, duplicated sequence and required indexed by using BWA.

### Mapping long-read raw data with minimap2 and samtools:
Download the reference FASTA file and sample's FASTQ file, then:
```
minimap2 -ax map-hifi ref.fa sample.fq > aligned.sam
samtools sort aligned.sam -o aligned.bam
samtools index -b aligned.bam
```
after this process, user will get <aligned.bam> file and <aligned.bam.bai> file.    

For more detail, please look at [Minimap2 GitHub pages](https://github.com/lh3/minimap2) and [SAMtools GitHub pages](https://github.com/samtools/samtools).

### Index reference chrM.fa required with BWA:
For fnumts, only needs to index reference chrM.fa.  
```
bwa index ref_chrM.fa
```
and user will see ref_chrM.fa file with five new files, ref_chrM.fa.amb, ref_chrM.fa.ann, ref_chrM.fa.bwt, ref_chrM.fa.pac, and ref_chrM.fa.sa. Next, copy the direct path of ref_chrM.fa as reference input index_ref_chrM.fa.   

For more detail, please look at [BWA GitHub pages](https://github.com/lh3/bwa).   

### Get chrM FASTA file
If user has trouble to get chrM fasta file, grep_chrM.py is a script for user to get only chrM sequence from reference FASTA file and save into FASTA standard format. The chrM.fa has duplicated sequences, so user can use this to do later analysis.
```
# Usage
./grep_chrM.py -i <input_fasta> [-o output.fa] [-c chrM_name_from_reference_fasta]

# Example
./grep_chrM.py -i GRCh38.p13.genome.fa
```
| Option | Option | Description| Default|
|:---------|:-------------|:-------------|:-------------|
|-i |--input_fasta|Input FASTA format file|False|
|-o |--output_fasta|Output FASTA format file|./chrM.fa|
|-c |--chromosome|Chromosome names from reference FASTA file to extract|chrM|


Download reference, please look at [UCSC](https://hgdownload2.soe.ucsc.edu/downloads.html) or [NCBI](https://www.ncbi.nlm.nih.gov/).

### fnumts can be using with full function:
```
./fnumts_main.py -b <index.bam> [-t <str>] [-d <int>] [-l <int>] [-s <int>] [-r <index_ref_chrM.fa>] [-o <str>]
```

### if user indexed BAM file is using human_hg38 as reference, and using defult options:
```
./fnumts_main.py -b <index.bam>
```
## Options of fnumts


### 1. -b, --bam, or -c, --cram,  input indexed bam file or cram file.
Option -b, --bam, or -c, --cram, is required. The index file, .bai or .crai need to be same location as .bam and .cram file. 
```
./fnumts_main.py -b <index.bam> [-r <path/directory/index_ref_chrM.fa>]
```
```
./fnumts_main.py -c <index.cram> [-r <path/directory/index_ref_chrM.fa>]
```
### 2. -r, --reference, input indexed reference chrM.fa.
Option -r, --refernece, is for user can input other human chrM reference or other species chrM reference. The refernce need to be same as input indexed BAM file refernce.
```
./fnumts_main.py -b <index.bam> [-r <path/directory/index_ref_chrM.fa>]
```
### 3. -t, --type, input datatype. hifi or ont.
Using this option to change the type you put in. PacBio HiFi or Oxford Nanopore(ONT). Default hifi.
```
./fnumts_main.py -b <index.bam> [-t <str>]
```
### 4. -d, --distance, maximun number of distance for merging intervals.
As this option, user can adjust the distance of each cluster range. Default 30.
```
./fnumts_main.py -b <index.bam> [-d <int>]
```
### 5. -l, --lenght, minmiun number of lenght for insertion sequence. 
For this option, user can adjust the lenght of the sequence. For example, if user set the number to 20, then fnumts will filter out the length of NUMTs less than 20bps. Default 10.
```
./fnumts_main.py -b <index.bam> [-l <int>]
```
### 6. -s, --supportread, minimun number of support reads for insertion reads.
-s, --supportread, is an option for user can adjust the support read number. For more accurate results, the number of support read needs to be higher. Default 3.
```
./fnumts_main.py -b <index.bam> [-s <int>]
```

## Output files
The output file directory includes:

| File name | Description| 
| :------------|:-------------|
|fnumts.vcf|VCF file of finding NUMTs. The chromosome, length, supportread, mapping quality, number of mismatch, start position, end position, type for each NUMT, mapping quality, and depth.|  

## Output files format
### fnumts.vcf
|Col|Type|Column Name|Description|
|-:|:--|:-------------|:-------------|
|1|str|CHROM|Target Chromosome.|
|2|int|POS|Target Position.|
|3|str|ID|Event ID|
|4|str|REF||
|5|str|ALT|The alternative allele is the allele found in studying.|
|6|int|QUAL|Event Score.|
|7|str|FILTER|PASS. If score lower than 40.00, shows LowQuality.|
|8|str|INFO|CHR,Chromosome where the NUMT sequence is located. LEN,Length of the NUMT sequence in base pairs. SUPREAD,Number of supporting reads for the NUMT sequence. NM,Number of mismatches in the NUMT sequence. MQ,Average mapping quality of the original support reads. chrMSPos,Start coordinate of the NUMT sequence in the mitochondrial genome. chrMEPos,End coordinate of the NUMT sequence in the mitochondrial genome. TYPE,Type of structural variant. METHOD,Detection method used for identifying the NUMT sequence. VAF,Variant allele frequency at the NUMT breakpoint.DEPTH,Total depth of coverage at the NUMT breakpoints|
|9|str|FORMAT|GT:DP:AD (Genomy Type, Total Reads, Support Reads)|
|10|str|SAMPLE_NAME|GenomeType:TotalReads:SupportReads|

#### QUAL system
Depend on the value of Number of missmatch, VAF, and Mapping Quality.


score = (0.3 * normalized_mismatches + 0.4 * normalized_vaf + 0.3 * normalized_mq) * 100


if score is lower than 40.00, FILTER will print LowQuality.

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##source=Fnumts
##fileDate=2024-09-30
##contig=<ID=chr21>
##INFO=<ID=CHR,Number=1,Type=String,Description="Chromosome where the NUMT sequence is located">
##INFO=<ID=LEN,Number=1,Type=Integer,Description="Length of the NUMT sequence in base pairs">
##INFO=<ID=SUPREAD,Number=1,Type=Integer,Description="Number of supporting reads for the NUMT sequence">
##INFO=<ID=NM,Number=1,Type=Integer,Description="Number of mismatches in the NUMT sequence">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality of the original support reads">
##INFO=<ID=chrMSPos,Number=1,Type=Integer,Description="Start coordinate of the NUMT sequence in the mitochondrial genome">
##INFO=<ID=chrMEPos,Number=1,Type=Integer,Description="End coordinate of the NUMT sequence in the mitochondrial genome">
##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=METHOD,Number=1,Type=String,Description="Detection method used for identifying the NUMT sequence">
##INFO=<ID=VAF,Number=1,Type=Float,Description="Variant allele frequency at the NUMT breakpoint">
##INFO=<ID=DEPTH,Number=1,Type=Integer,Description="Total depth of coverage at the NUMT breakpoints">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth of coverage at the NUMT breakpoints">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Number of supporting reads for the NUMT sequence">
##CommandLine=fnumts_main.py -b test.bam
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  test.bam
chr21   21708298        FN1     N       AACTGGCTTCAATCTACTTCTCCCGCCGCCGGGAAAAAAGGCGGGAGAAGCCCCGGCAGGTTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAATCACCTCATTTTCTCGTC  97.20   PASS    CHR=chr21;LEN=118;SUPREAD=25;NM=0;MQ=60;chrMSPos=5714;chrMEPos=5832;TYPE=NUMT;METHOD=Fnumt;VAF=0.93;DEPTH=27    GT:DP:AD        1/1:27:25
```

### fnumts.paf
|Col|Type|Description|
|--:|:---|:-------------|
|1|str|Query name|
|2|int|Query length|
|3|int|Query start|
|4|int|Query end|
|5|str|Relative strand: "+" or "-"|
|6|str|Target name|
|7|int|Target length|
|8|int|Target start|
|9|int|Target end|
|10|int|Number of mismatch|
|11|int|Alignment block length|


## Additional information
1. [pysam](https://pysam.readthedocs.io/en/latest/api.html)
2. [DeBreak](https://www.nature.com/articles/s41467-023-35996-1)
3. [argparse](https://docs.python.org/3/library/argparse.html)
4. [pandas](https://devdocs.io/pandas~0.25/)
5. [BWA](https://bio-bwa.sourceforge.net/bwa.shtml)
6. [bsalign](https://github.com/ruanjue/bsalign)
7. [wtdbg2](https://github.com/ruanjue/wtdbg2)
