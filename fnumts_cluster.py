#!/usr/bin/env python3

import pysam

def write_ins_details(ct, line, s, r, length):
    ins_start = s + 1
    ins_end = s + 1 + length
    ins_quality = line.mapping_quality
    ins_seq = line.query_sequence[r:r+length]
    ct.write(f'{line.reference_name}\t{ins_start}\t{ins_end}\t{ins_quality}\tINS\t{length}\t{ins_seq}\n')

def process_soft_clips(ct, line):
    try:
        sa_tag = line.get_tag('SA')
        if 'chrM' in sa_tag:
            cigar = line.cigarstring
            if 'S' in cigar:
                start = line.reference_start
                quality = line.mapping_quality
                seq_length = line.query_length
                seq = line.query_sequence
                ct.write(f'{line.reference_name}\t{start+1}\t{start+seq_length}\t{quality}\tINS\t{seq_length}\t{seq}\n')
    except KeyError:
        pass

def process_cigar_operations(ct, line):
    s = line.reference_start
    r = 0
    for operation, length in line.cigartuples:
        if operation == 4:  # S
            r += length
        elif operation == 0:  # M
            s += length
            r += length
        elif operation == 1:  # I
            write_ins_details(ct, line, s, r, length)
            r += length
        elif operation == 2:  # D
            s += length

def process_bam_file(bam_filename, writepath):
    bamfile = pysam.AlignmentFile(bam_filename, 'rb')
    cigar_file = f'{writepath}51_cigar.txt'
    
    with open(cigar_file, 'w') as ct:
        for line in bamfile.fetch():
            if line.flag & 256:
                continue
            process_soft_clips(ct, line)
            process_cigar_operations(ct, line)
            
    bamfile.close()
