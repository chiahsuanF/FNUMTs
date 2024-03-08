#!/usr/bin/env python3

import argparse

def extract_chromosomes(input_fasta, output_fasta, chromosome_names):
    selected_chromosomes = {}
    current_chromosome = None
    with open(input_fasta, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                current_chromosome = line.strip()
                for name in chromosome_names:
                    if name in current_chromosome:
                        selected_chromosomes[current_chromosome] = []
                        break
                else:
                    current_chromosome = None
            else:
                if current_chromosome:
                    selected_chromosomes[current_chromosome].append(line.strip())

    with open(output_fasta, 'w') as outfile:
        for chromosome, sequence_lines in selected_chromosomes.items():
            outfile.write(chromosome + '\n')
            outfile.write('\n'.join(sequence_lines) + '\n')
        

parser = argparse.ArgumentParser(description='Extract chromosome sequences from a FASTA file', usage="grep_chrM.py [-h] -i <input_fasta> [options]")
parser.add_argument('-i', '--input_fasta', type=str, default=False, help='Input FASTA file')
parser.add_argument('-o', '--output_fasta', type=str, default='./chrM.fa', help='Output FASTA file')
parser.add_argument('-c', '--chromosomes', nargs='+', default=['chrM'], help='chrM names from reference FASTA file to extract')
args = parser.parse_args()

extract_chromosomes(args.input_fasta, args.output_fasta, args.chromosomes)
