#!/usr/bin/env python3

import os, argparse, sys


def check_input_arguments(args):
    if not args.bam and not args.cram:
        print("Error: No input file given.")
        sys.exit(1)

def make_output_directory(output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

def parse_args():
    parser = argparse.ArgumentParser(description='Finding NUMTs from long-read sequencing data', usage="fnumt.py [-h] -b <index.bam> [options]")
    
    parser.add_argument('-t', '--type', type=str, choices=['hifi', 'ont'], default='hifi', help="[STR] Input datatype, HIFI or ONT")
    parser.add_argument('-ts', "--threads", type=int, default=4, help="[INT] Number of CPU threads to use")
    parser.add_argument('-b', '--bam', type=str, help="[STR] Input sorted bam, index required")#
    parser.add_argument('-c', '--cram', type=str, help="[STR] Input sorted cram file, index required")
    parser.add_argument('-mr', '--reference', type=str, default='/scratch/fanc/fnumts/ref/d_chrM.fa', help="[STR] Import reference chromosome mitochondria. The reference needs to be indexed by BWA")
    parser.add_argument('-dr', '--deletion_reference', type=str, default='/scratch/fanc/fnumts/ref/GRCh38.p13.genome.fa', help="[STR] Import reference chromosome genome")
    parser.add_argument('-o', '--outpath', type=str, default='./fnumt_out/', help='[STR] Name of output directory')
    parser.add_argument('-ov', '--output_vcfname', type=str, default='fnumts', help="[STR] The name of the output VCF file")
    parser.add_argument('-bd',"--base_dir", type=str, default='../', help="[STR] Base directory for BAM file")
    parser.add_argument('-cr',"--cram_reference", type=str, default='/scratch/fanc/fnumts/ref/1KG_ONT_VIENNA_hg38.fa', help="[STR] Import cram reference fasta file. Cram file required.")
    parser.add_argument('-d', '--distance', type=int, default=30, help='[INT] Maximum number of distance for merging intervals') #30
    parser.add_argument('-ml', '--max_length', type=int, default=16569, help='[INT] Maximum number of deletion length. it should be less than chrM total length. default 16568(hg38 chrM)') #16568
    parser.add_argument('-l', '--length', type=int, default=10, help='[INT] Minimum length for insertion sequence')
    parser.add_argument('-s', '--support_read', type=int, default=3, help='[INT] Minimum number of support reads for insertion reads')#3
    args = parser.parse_args()

    return args
