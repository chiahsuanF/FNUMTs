#!/usr/bin/env python3

import os
import argparse
import sys

def check_input_arguments(args):
    if not args.bam:
        print("Error: No input file given.")
        sys.exit(1)

def make_output_directory(output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

def parse_args():
    parser = argparse.ArgumentParser(description='Finding NUMTs from long-read sequencing data', usage="fnumt.py [-h] -bam <index.bam> [options]")
    parser.add_argument('-bam', type=str, default=False, help = "input sorted bam. index required")
    parser.add_argument('-o', '--outpath', type=str, default='./4test/')
    parser.add_argument('-r', '--reference', type=str, default='/data/scratch/fanc/fnumts/ref/chrM.fa', help = "Import reference chromosome mitochondria. The reference needs to be indexed by BWA")
    parser.add_argument('-sr', '--sref', type=str, default='hg38', help = "Import reference for snpEff annotation")    
    parser.add_argument('-d', '--distance', type=int, default='10', help = 'Maximun number of distance for merging intervals')
    parser.add_argument('-l', '--lenght', type=int, default='10', help = 'Minmiun number of lenght for insertion sequence')
    parser.add_argument('-s', '--supportread', type=int, default='5', help = 'Minimun number of support reads for insertion reads')
    return parser.parse_args()
