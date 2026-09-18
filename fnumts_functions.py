#!/usr/bin/env python3

import os, time
import pandas as pd
from multiprocessing import Pool, cpu_count
 
def sort_and_filter_cigar_data(args, writepath):
    chr_list2 = [
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
        "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
        "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"
    ]

    input_file = os.path.join(writepath, '51_cigar.txt')
    output_file = os.path.join(writepath, '53_cut_3.srt.txt')

    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Missing input file: {input_file}")

    sort_test = pd.read_csv(
        input_file, sep="\t", header=None,
        names=["chr", "start_pos", "end", "MQ", "INS", "length", "seq"]
    )

    k = sort_test[sort_test["length"] >= int(args.length)].copy()
    k["chr"] = k["chr"].astype("category")
    k["chr"] = k["chr"].cat.set_categories(chr_list2)

    newsort = k.sort_values(["chr", "start_pos"]).dropna(how='any', axis=0)
    newsort["end_pos"] = newsort["start_pos"] + 1

    newsort[["chr", "start_pos", "end_pos", "INS", "MQ", "length", "seq"]].to_csv(
        output_file, sep='\t', header=False, index=False
    )
