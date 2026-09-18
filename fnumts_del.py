#!/usr/bin/env python3

import os 
import pandas as pd
from Bio import SeqIO
from fnumts_sortcigar import *

def count_del(args, writepath):
    import tempfile
    chunksize = 1_000_000
    chr_list2 = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    inpath = os.path.join(writepath, '51_cigar_del.txt')
    outpath = os.path.join(writepath, '53_cut_3_del.txt')
    min_len = int(args.length)

    # temp files per chromosome
    tmpfiles = {c: tempfile.NamedTemporaryFile(delete=False, mode="w")
                for c in chr_list2}

    reader = pd.read_csv(
        inpath, sep='\t', header=None, chunksize=chunksize,
        names=["chr", "start_pos", "end", "MQ", "DEL", "length"]
    )

    for chunk in reader:
        chunk = chunk[chunk["length"] >= min_len]

        for c in chr_list2:
            sub = chunk[chunk["chr"] == c]
            if len(sub) > 0:
                sub.to_csv(tmpfiles[c].name, sep="\t",
                           header=False, index=False, mode="a")

    for f in tmpfiles.values():
        f.close()

    # final sort + concat
    with open(outpath, "w") as outfile:
        for c in chr_list2:
            tf = tmpfiles[c].name
            if os.path.getsize(tf) == 0:
                continue

            df = pd.read_csv(tf, sep="\t", header=None,
                             names=["chr","start_pos","end","MQ","DEL","length"])
            df = df.sort_values(["start_pos"])
            df.to_csv(outfile, sep="\t", header=False, index=False)
            os.remove(tf)

def extract_and_sort_columns_del(df):
    selected_columns = ["chr2", "sp2", "ep2", "depth" , "DEL", "MQ", "len", "avg_MQ"]
    sorted_df = df[selected_columns].sort_values(by=["chr2", "sp2", "len"])
    return sorted_df


def read_and_process_data_del(input_file, chr_list):
    df = pd.read_csv(input_file, sep="\t", header=None) #MQ
    df.columns = ["chr", "sp1", "ep1", "depth", "chr2", "sp2", "ep2", "MQ", "DEL", "len"]
    
    # Group by chr and sp1, and calculate the average MQ
    df['avg_MQ'] = df.groupby(['chr', 'sp1'])['MQ'].transform('mean').astype(int)
    
    idx = (
        df.groupby(['chr', 'sp1'])['sp2']
        .apply(lambda x: x.value_counts().idxmax())   # most frequent sp2
        .reset_index()
    )
    
    df_merged = pd.merge(df, idx, on=['chr', 'sp1', 'sp2'], how='inner')
    df_merged = df_merged.drop_duplicates(subset=['chr', 'sp1', 'sp2'])
    
    sorted_df = extract_and_sort_columns_del(df_merged)
    sorted_df["chr2"] = sorted_df["chr2"].astype("category")
    sorted_df["chr2"] = sorted_df["chr2"].cat.set_categories(chr_list)

    return sorted_df

def sort_and_uniq(input_file, output_file, uniq_cols=None):
    df = pd.read_csv(input_file, sep="\t", header=None)
    
    df[1] = pd.to_numeric(df[1], errors='coerce')
    df[2] = pd.to_numeric(df[2], errors='coerce')
    df_sorted = df.sort_values(by=[0, 1, 2])

    df_unique = df_sorted.drop_duplicates(subset=uniq_cols)
    df_unique.to_csv(output_file, sep="\t", index=False, header=False)
 
def extract_subsequences_from_bed(fasta_file, bed_file, output_file):
    ref_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    with open(bed_file) as bed, open(output_file, "w") as out_fa:
        for idx, line in enumerate(bed, start=1):  # line index starts at 1
            if line.strip() == "" or line.startswith("#"):
                continue

            fields = line.strip().split()

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            depth = int(fields[3])
            ntype = fields[4]
            mq = int(fields[5])            
            length = int(fields[6])  # column 9

            if chrom not in ref_dict:
                print(f"Warning: {chrom} not found in reference.")
                continue

            seq_record = ref_dict[chrom]
            subseq = seq_record.seq[start:end]

            # Optional sanity check
            expected_len = end - start
            if len(subseq) != expected_len:
                print(f"Warning: {chrom}:{start}-{end} expected {expected_len} bp but got {len(subseq)} bp")
            del_len = int(length -1)
            # Header with length and line index
            header = f"{chrom}_{start}_{depth}_{ntype}_{mq}"
            out_fa.write(f">{header}\n{str(subseq)}\n")
