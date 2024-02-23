#!/usr/bin/env python3

import os
import pandas as pd
import argparse
import sys
import pysam
from fnumts_cluster import *
from fnumts_bedtools import *
from fnumts_writevcf import *

def check_input_arguments(args):
    if not args.bam:
        print("Error: No input file given.")
        sys.exit(1)

def make_output_directory(output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

def parse_args():
    parser = argparse.ArgumentParser(description='finding NUMTs insertion')
    parser.add_argument('-bam', type=str, default=False)
    parser.add_argument('-o', '--outpath', type=str, default='./ffnumt_out/')
    parser.add_argument('-d', '--distance', type=int, default='10', help = 'Maximun distance for merging intervals')
    parser.add_argument('-l', '--lenght', type=int, default='10', help = 'Minmiun lenght for insertion sequence')
    parser.add_argument('-s', '--supportread', type=int, default='5', help = 'Minimun support reads for insertion reads')
    return parser.parse_args()
args = parse_args()
writepath = args.outpath

# cluster collecting 
check_input_arguments(args)
make_output_directory(writepath)
chr_list_file = '/data/scratch/fanc/fnumts/chr_list'
with open(chr_list_file, 'r') as file:  
    chr_list = file.read().splitlines()
chr_list2 = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",  "chr10", "chr11", "chr12", "chr13", 
             "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]


logfile = open(writepath + 'login_note', 'w')
logfile.write('Start collecting clusters...\n')
logfile.close()


process_bam_file(args.bam, writepath)
sort_and_filter_cigar_data(args, writepath)
bed_input = pd.read_csv(writepath + '53_cut_3.srt.txt', sep='\t', header=None)
bed = bed_input.values.tolist()
merged_output = merge_intervals(bed, int(args.distance))
output1 = pd.DataFrame(merged_output)
output2 = output1[output1[3] >= args.supportread]

output2.columns = ["chr", "startp", "endp", "depth"]
output3 = output2[output2["chr"].isin(chr_list2)].sort_values(['chr', 'startp'])
output3.to_csv(writepath + '54_sdtin.txt', sep='\t', header=None, index=False)

a_intervals = read_intervals(writepath + '54_sdtin.txt')
b_intervals = read_intervals(writepath + '53_cut_3.srt.txt')
intersecting_intervals = bed_intersect(a_intervals, b_intervals)

result_df = pd.DataFrame(intersecting_intervals)
result_df.to_csv(writepath + '55_bedtools_f.txt', sep='\t', header=None, index=None)

input_file = writepath +'55_bedtools_f.txt'
output_file = writepath +'58_clusters.ins.sorted.txt'
processed_data = read_and_process_data(input_file, chr_list)
save_sorted_data(processed_data, output_file)

input_file_58 = writepath +'58_clusters.ins.sorted.txt'
output_file_59 = writepath +'59_pre_filename.fna'
pre_data = read_data(input_file_58)
processed_data = add_header_column(pre_data)
save_to_fna(processed_data, output_file_59)
os.chdir(writepath)
os.system("/share/apps/rc/software/BWA/0.7.17-GCC-8.3.0/bin/bwa mem -t4 /data/scratch/fanc/fnumts/ref/chrM.fa 59_pre_filename.fna > clusters.ins.srot.sam")

samfile1 = pysam.AlignmentFile("clusters.ins.srot.sam", "r")
iter = samfile1.fetch()
with open ('60_.txt', 'w') as pbs1:
    for read1 in iter:
        if read1.reference_name == 'chrM':
            pbs1.write(read1.query_name + '\n')
pbs1.close()


logfile = open('login_note', 'a')
logfile.write('End collecting clusters...\n')
logfile.write('Start calling bsalign poa ...\n')
logfile.close()


# bsalign
pbsw_file = '60_.txt'
read2_file = '58_clusters.ins.sorted.txt'
preprocess_and_write(pbsw_file, read2_file)

os.system("ls *_*.fa > fastalist")
fastalist=open("fastalist",'r').read().split('\n')[:-1]
counter = 0
cmd = '/home/fanc/path/bsalign'
for k in fastalist:
    counter = counter +1
    os.system(cmd + ' poa '+ k +' -o ' +  k[:-3] + ".align.fa ")
    os.system("sed 's/cns_seq/"  +k[:-3] + "/' < " + k[:-2] + "align.fa > "+ k[:-3] + ".rename.fa")
    
os.system("cat *.rename.fa > temp_final.align.fa")

pre_data_2 = pd.read_csv('temp_final.align.fa', header = None)
pre_data_2.columns = ["col1"]

with open("final.align.fa", 'w') as svfile:
    for l in range(0, len(pre_data_2),2):
        m = l +1
        len1 = len(pre_data_2["col1"][m])
        svfile.write('{}_{}\n{}\n'.format(pre_data_2["col1"][l], str(len1), pre_data_2["col1"][m]))
svfile.close


logfile = open('login_note', 'a')
logfile.write('End calling bsalign poa ...\n')
logfile.close()


os.system("/share/apps/rc/software/BWA/0.7.17-GCC-8.3.0/bin/bwa mem -k13 -T0 /data/scratch/fanc/fnumts/ref/chrM.fa final.align.fa > fnumts.sam")

logfile = open('login_note', 'a')
logfile.write('Start writing VCF file and variant annotation ...\n')
logfile.close()

# snpeff
samfile2 = pysam.AlignmentFile("fnumts.sam", "r")
iter2 = samfile2.fetch()

with open ('70_.txt', 'w') as snp1:
    for read2 in iter2:
        Ms = int(read2.reference_start)
        nm_tag = None
        for tag, value in read2.get_tags():
            if tag == "NM":
                nm_tag = value
                break
        if nm_tag is not None:  
            snp1.write(read2.query_name + '\t' + read2.query_sequence + '\t' + str(read2.mapping_quality) + '\t' + str(Ms +1) + '\t' + str(nm_tag) +'\n')

snp2 = pd.read_csv('70_.txt', sep = '\t', header = None)
snp2 = snp2.replace('_', ' ', regex=True)
snp2.to_csv( "71_.txt", sep = "\t", header=None, index=0)

snp3 = pd.read_csv('71_.txt', sep="\t| ", engine='python', header=None)
snp3.columns = ["chr", "start_pos", "depth", "lenght", "seq", "MQ", "M_start", "NM"]

with open('fnumts.vcf', 'w') as snp4:
    vcf_rows = snp3.apply(generate_vcf_row, axis=1)
    snp4.write("\n".join(vcf_rows))

os.mkdir('poa_output')
os.system('mv -t poa_output *.rename.fa')
os.system('rm *.fa *.txt 59_pre_filename.fna fastalist')
os.system("java -Xmx8g -jar /data/scratch/fanc/app/snpEff/snpEff.jar  -v hg38 fnumts.vcf > fnumts.ann.vcf")


logfile = open('login_note', 'a')
logfile.write('End writing VCF file and variant annotation ...\n')
logfile.close()
