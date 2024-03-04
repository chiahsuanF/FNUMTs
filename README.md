# fnumts

Author: Chia-Hsuan Fan

Email: fanc@uab.edu

Draft date:

## Quick Start

## Description

## Depencency
Dependencies for fnumts:

* python3
* pysam(version 0.22.0)
* argparse(version 1.4.0)
* pandas(version 2.1.1)
* bwa(version 0.7.17)
* bsalign(version 1.2.1)
* snpEff(version)
  
## Installation

### Conda
```
conda create -n fnumts
conda activate fnumts
conda install -c bioconda fnumts
```
### Git
```
git clone
```

## General usage

## Options of fnumts
| Option | Description| Default|
| :------------|:-------------|:-------------|
|-bam|Import aligned BAM file|False|
|-o|Output directory |fnumts|
|-r|Import reference chromosome mitochondria. The reference needs to be index by BWA|ref/chrM.fa|
|-d|Maximun distance for merging intervals|10|
|-l|Minmiun lenght for insertion sequence|10|
|-s|Minimun support reads for insertion reads|5|

## Output files
The output file directory includes:

| File name | Description| 
| :------------|:-------------|
|fnumts.ann.vcf|Standard VCF file of find NUMTs and annotation. The chromosome, length, supportread, mapping quality, number of missmatch, start position, end position, type, annotation for each NUMT.|
|fnumts.vcf|Temporery VCF file of finding NUMTs. The chromosome, length, supportread, mapping quality, number of missmatch, start position, end position, type for each NUMT.|  
|fnumts.sam|Temporary SAM archive consistent with chrM.|
|clusters.ins.srot.sam|Temporary SAM archive of the cluster aligned with chrM.|
|poa_output|Temporary files during each cluster seeking consensus.|
