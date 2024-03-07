#!/usr/bin/env python3

import pandas as pd

def sort_and_filter_cigar_data(args, writepath):
    chr_list2 = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
        
    sort_test = pd.read_csv(writepath+'51_cigar.txt', sep="\t", header=None, 
                            names=["chr", "start_pos", "end", "MQ", "INS", "lenght", "seq"])
    k = sort_test[sort_test.lenght >= int(args.lenght)]
    k.chr = k.chr.astype("category")
    k.chr = k.chr.cat.set_categories(chr_list2)
    newsort = k.sort_values(['chr','start_pos'])
    newsort["end_pos"] = newsort["start_pos"] + 1
    newsort = newsort.dropna(how='any',axis=0)
    newsort[['chr', 'start_pos', 'end_pos', 'INS', 'MQ', 'lenght', 'seq']].to_csv(writepath+'53_cut_3.srt.txt', sep='\t', header=False, index=False)
