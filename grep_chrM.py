#!/usr/bin/env python3

import argparse

# Function to extract chromosomes and combine every two lines into one
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
    
    # Process sequences to merge every two lines into one
    def merge_sequences(sequence_lines):
        merged_sequences = []
        # Iterate in pairs, combining two lines into one
        for i in range(0, len(sequence_lines), 2):
            if i + 1 < len(sequence_lines):
                merged_sequences.append(sequence_lines[i] + sequence_lines[i + 1])
            else:
                # If there's an odd line at the end, just add it as is
                merged_sequences.append(sequence_lines[i])
        return merged_sequences

    with open(output_fasta, 'w') as outfile:
        for chromosome, sequence_lines in selected_chromosomes.items():
            outfile.write(chromosome + '\n')
            # Merge every two lines and write them to the file
            merged_lines = merge_sequences(sequence_lines)
            outfile.write('\n'.join(merged_lines) + '\n')

# Argument parser
parser = argparse.ArgumentParser(description='Extract chromosome sequences from a FASTA file', usage="grep_chrM.py [-h] -i <input_fasta> [options]")
parser.add_argument('-i', '--input_fasta', type=str, required=True, help='Input FASTA file')
parser.add_argument('-o', '--output_fasta', type=str, default='./chrM.fa', help='Output FASTA file')
parser.add_argument('-c', '--chromosomes', nargs='+', default=['chrM'], help='Chromosome names from reference FASTA file to extract')

# Parse arguments
args = parser.parse_args()

# Call the function to extract and process chromosomes
extract_chromosomes(args.input_fasta, args.output_fasta, args.chromosomes)
