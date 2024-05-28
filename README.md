# fnumts
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
./fnumts_main.py -bam <index.bam> [option]
```

```sh
# Finding NUMTs with indexed bam. The reference using human hg38.
./fnumts_main.py -bam <index.bam>

# Example usage:
./fnumts_main.py -bam index_human_hg38.bam
```

```sh
# Finding NUMTs with indexed bam. The reference using other species' chrM indexed by BWA.
./fnumts_main.py -bam <index.bam> -r <index_ref_chrM.fa>

# Example usage:
./fnumts_main.py -bam index.bam -r path/directory/index_ref_chrM.fa
```

```sh
# Finding NUMTs with indexed bam. The reference using other species' chrM.fa indexed by BWA.
# Adjust distance, length, and support reads.
./fnumts_main.py -bam <index.bam> -d <int> -l <int> -s <int> -r <index_ref_chrM.fa>

# Example usage:
./fnumts_main.py -bam index.bam -d 20 -l 20 -s 10 -r path/directory/index_ref_chrM.fa
```

## Description
fnumts is a tool for finding insertion mitochondria DNA in nuclear genome. fnumts can detect NUMT length from 1 bp up to whole mitochondira DNA.      

The input file should be a long-read indexed BAM(PacBio Hifi) file. The input reference chrM needs to be same as aligned bam file reference, and needs to be indexed by BWA. The output is a standard VCF file which passed filter. This program was test on a x86_64 Linux system with 40GB physical memory.

## Depencency
Dependencies for fnumts:

* python3
* pysam(version 0.22.0)
* argparse(version 1.4.0)
* pandas(version 2.1.1)
* BWA(version 0.7.17)
* bsalign(version 1.2.1)
* snpEff(version 5.2)(with java version 1.11 or later)
  
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
| Option | Description| Default|
| :------------|:-------------|:-------------|
|-bam|Import aligned BAM file|False|
|-o|Output directory |fnumts|
|-r|Import reference chromosome mitochondria. The reference needs to be indexed by BWA|userpath/fnumts/ref/chrM.fa|
|-sr|Import reference for snpEff annotation|hg38|
|-d|Maximun number of distance for merging intervals|10|
|-l|Minmiun number of lenght for insertion sequence|10|
|-s|Minimun number of support reads for insertion reads|5|

## Use cases
fnumts required input as an aligned BAM file with header. Mapping raw data suggests using minimap2 and samtools. Input reference only needs reference chrM.fa and required indexed by using BWA.

### Mipping long-read raw data with minimap2 and samtools:
Download the reference FASTA file and sample's FASTQ file, then:
```
minimap2 -ax map-hifi ref.fa sample.fq > aligned.sam
samtools sort aligned.sam -o aligned.bam
samtools index -b index.bam
```
after this process, user will get a index.bam file.    

For more detail, please look at [Minimap2 GitHub pages](https://github.com/lh3/minimap2) and [SAMtools GitHub pages](https://github.com/samtools/samtools).

### Index reference chrM.fa required with BWA:
For fnumts, only needs to index reference chrM.fa.  
```
bwa index ref_chrM.fa
```
and user will see ref_chrM.fa file with five new files, ref_chrM.fa.amb, ref_chrM.fa.ann, ref_chrM.fa.bwt, ref_chrM.fa.pac, and ref_chrM.fa.sa. Next, copy the direct path of ref_chrM.fa as reference input index_ref_chrM.fa.   

For more detail, please look at [BWA GitHub pages](https://github.com/lh3/bwa).   

### Get chrM FASTA file
If user has trouble to get chrM fasta file, grep_chrM.py is a script for user to get only chrM sequence from reference FASTA file and save into FASTA format.  
```
# Usage
./grep_chrM.py -i <input_fasta> [-o output.fa] [-c chrM_name_from_reference_fa]

# Example
./grep_chrM.py -i GRCh38.p13.genome.fa
```
Download reference, please look at [UCSC](https://hgdownload2.soe.ucsc.edu/downloads.html) or [NCBI](https://www.ncbi.nlm.nih.gov/).
### fnumts can be using with full function:
```
./fnumts_main.py -bam <index.bam> -d <int> -l <int> -s <int> -r <index_ref_chrM.fa>
```

### if userr indexed BAM file is using human_hg38 as reference, and don't want to adjust other options:
```
./fnumts_main.py -bam <index.bam>
```
## Options of fnumts

### 1. -r, --reference, input indexed reference chrM.fa.
Option -r, --refernece, is for user can input other human chrM reference or other species chrM reference. The refernce need to be same as input indexed BAM file refernce.
```
./fnumts_main.py -bam <index.bam> [-r <path/directory/index_ref_chrM.fa>]
```
### 2. -d, --distance, maximun number of distance for merging intervals.
As this option, user can adjust the distance of each cluster range. 
```
./fnumts_main.py -bam <index.bam> [-d <int>]
```
### 3. -l, --lenght, minmiun number of lenght for insertion sequence. 
For this option, user can adjust the lenght of the sequence. For example, if user set the number to 20, then fnumts will filter out the length of NUMTs less than 20bps.
```
./fnumts_main.py -bam <index.bam> [-l <int>]
```
### 4. -s, --supportread, minimun number of support reads for insertion reads.
-s, --supportread, is an option for user can adjust the support read number. For more accurate results, the number of support read needs to be higher.
```
./fnumts_main.py -bam <index.bam> [-s <int>]
```
### 5. -sr, --sref, Import reference for snpEff annotation.
the option -sr is for users to adjust the reference to do annotation on snpEff.
```
./ fnumts_main.py -bam <index.bam> [-sr reference]
```
For more detail, please look at [snpEff GitHub pages](https://pcingola.github.io/SnpEff/). 

## Output files
The output file directory includes:

| File name | Description| 
| :------------|:-------------|
|fnumts.ann.vcf|Standard VCF file of find NUMTs and annotation. The chromosome, length, supportread, mapping quality, number of mismatch, start position, end position, type, annotation for each NUMT.|
|fnumts.vcf|Temporery VCF file of finding NUMTs. The chromosome, length, supportread, mapping quality, number of mismatch, start position, end position, type for each NUMT.|  
|numt.fa|Standard FASTA format file that saves all the NUMT sequences.|
|fnumts.sam|Temporary SAM archive consistent with chrM.|
|clusters.ins.srot.sam|Temporary SAM archive of the cluster aligned with chrM.|
|poa_output|Temporary files during each cluster seeking consensus.|

## Output files format
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
7. [snpEff](https://pcingola.github.io/SnpEff/)
