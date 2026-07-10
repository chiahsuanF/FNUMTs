#!/usr/bin/env python3

import os, math
import pandas as pd
from multiprocessing import Pool, cpu_count


def sort_and_filter_cigar_data(args, writepath):
    chr_list2 = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
        
    sort_test = pd.read_csv(writepath+'51_cigar.txt', sep="\t", header=None, 
                            names=["chr", "start_pos", "end", "MQ", "INS", "length", "seq"])
    k = sort_test[sort_test.length >= int(args.length)]
    k.chr = k.chr.astype("category")
    k.chr = k.chr.cat.set_categories(chr_list2)
    newsort = k.sort_values(['chr','start_pos'])
    newsort["end_pos"] = newsort["start_pos"] + 1
    newsort = newsort.dropna(how='any',axis=0)
    newsort[['chr', 'start_pos', 'end_pos', 'INS', 'MQ', 'length', 'seq']].to_csv(writepath+'53_cut_3.srt.txt', sep='\t', header=False, index=False)

def merge_intervals(bed_list, max_distance):
    if not bed_list:
        return []

    bed_list = sorted(bed_list, key=lambda x: x[1])
    merged = []
    current_interval = bed_list[0]
    count = 1  

    for i in range(1, len(bed_list)):
        if bed_list[i][0] == current_interval[0] and bed_list[i][1] <= current_interval[2] + max_distance:
            current_interval = (current_interval[0], current_interval[1], max(current_interval[2], bed_list[i][2]))
            count += 1
        else:
            merged.append((current_interval[0], current_interval[1], current_interval[2], count))
            current_interval = bed_list[i]
            count = 1
    merged.append((current_interval[0], current_interval[1], current_interval[2], count))
    return merged

def read_intervals(filename):
    intervals = pd.read_csv(filename, sep='\t', header=None).values
    return sorted(intervals, key=lambda x: (x[0], x[1], x[2]))

def bed_intersect(a_intervals, b_intervals):
    intersecting_intervals = []
    i, j = 0, 0

    while i < len(a_intervals) and j < len(b_intervals):
        a_interval = a_intervals[i]
        b_interval = b_intervals[j]
        
        if a_interval[0] < b_interval[0] :
            i += 1
        elif a_interval[0] > b_interval[0] :
            j += 1       
        elif (a_interval[0] == b_interval[0])  and  (a_interval[1] < b_interval[2]) and ((a_interval[2]) >= (b_interval[1])):
            intersecting_intervals.append(list(a_interval) + list(b_interval))
            j += 1             
        elif (a_interval[0] == b_interval[0]) and ((a_interval[2] ) < (b_interval[1])):
            i += 1
        else:
            j += 1
            
    return intersecting_intervals
