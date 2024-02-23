#!/usr/bin/env python3

import pandas as pd
import pysam

def process_bam_file(bam_filename, writepath):
    bamfile = pysam.AlignmentFile(bam_filename, 'rb')
    cigar_file = f'{writepath}51_cigar.txt'
    with open(cigar_file, 'w') as ct:
        for line in bamfile.fetch():
            if line.flag & 256:
                continue
            s = line.reference_start
            r = 0
            # M=0, I=1, D=2, S=4, H=5
            for operation, lenght in line.cigartuples:
                if operation == 4:  # S
                    r += lenght
                elif operation == 0:  # M
                    s += lenght
                    r += lenght
                elif operation == 1:  # I
                    ct.write(f'{line.query_name}\t{line.flag}\t{line.reference_name}\t{s+1}\t{s+1+lenght}\t{line.mapping_quality}\tINS\t{lenght}\t{line.query_sequence[r:r+lenght]}\n')
                    r += lenght
                elif operation == 2:  # D
                    s += lenght
    ct.close()
    
def sort_and_filter_cigar_data(args, writepath):
    chr_list2 = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
        
    sort_test = pd.read_csv(writepath+'51_cigar.txt', sep="\t", header=None, 
                            names=["name", "flag", "chr", "start_pos", "end", "MQ", "INS", "lenght", "seq"])
    k = sort_test[sort_test.lenght >= int(args.lenght)]
    k.chr = k.chr.astype("category")
    k.chr = k.chr.cat.set_categories(chr_list2)
    newsort = k.sort_values(['chr','start_pos'])
    newsort["end_pos"] = newsort["start_pos"] + 1
    newsort = newsort.dropna(how='any',axis=0)
    newsort[['chr', 'start_pos', 'end_pos', 'INS', 'MQ', 'lenght', 'seq']].to_csv(writepath+'53_cut_3.srt.txt', sep='\t', header=False, index=False)

def preprocess_and_write(pbsw_file, read2_file):
    pbsw1 = pd.read_csv(pbsw_file, sep='\t', header=None)
    pbsw1 = pbsw1.replace('_', ' ', regex=True)
    pbsw1.to_csv("61_.txt", sep="\t", header=None, index=0)

    pbsw1 = pd.read_csv('61_.txt', sep="\t| ", engine='python', header=None)
    pbsw1.columns = ['col1', 'col2', 'col3', 'col4', 'col5']

    read2 = pd.read_csv(read2_file, sep='\t', header=None)
    read2.columns = ['col1', 'col2', 'depth', 'MQ', 'len', 'seq']

    open_files = {}

    for i in range(len(pbsw1)):
        chr_val = pbsw1.loc[i, "col1"]
        start_pos = pbsw1.loc[i, "col2"]
        depth = pbsw1.loc[i, "col3"]
        
        file_name = f"{chr_val}_{start_pos}_{depth}.fa"

        if file_name in open_files:
            svfile = open_files[file_name]
        else:
            svfile = open(file_name, 'a')
            open_files[file_name] = svfile

        seq = read2[(read2['col1'] == chr_val) & (read2['col2'] == start_pos)]['seq'].values
        seq = seq[0] if len(seq) > 0 else ""  # If there's a match, use the sequence, otherwise, use an empty string

        svfile.write(f'>{chr_val}_{start_pos}_{i}\n{seq}\n')

    for svfile in open_files.values():
        svfile.close()
        
def read_and_process_data(input_file, chr_list):
    df = pd.read_csv(input_file, sep="\t", header=None)
    df.columns = ["chr", "sp1", "ep1", "depth", "chr2", "sp2", "ep2", "INS", "MQ", "len", "seq"]

    sorted_df = extract_and_sort_columns(df)
    sorted_df["chr"] = sorted_df["chr"].astype("category")
    sorted_df["chr"] = sorted_df["chr"].cat.set_categories(chr_list)

    return sorted_df

def extract_and_sort_columns(df):
    selected_columns = ["chr", "sp1", "depth", "INS", "len", "seq"]
    sorted_df = df[selected_columns].sort_values(by=["chr", "sp1", "len"])
    return sorted_df

def save_sorted_data(sorted_df, output_file):
    sorted_df.to_csv(output_file, sep='\t', header=None, index=False)

def read_data(input_file):
    pre_data = pd.read_csv(input_file, sep="\t", header=None)
    pre_data.columns = ["chr", "start_pos", "depths", "INS", "lenght", "seq"]
    return pre_data

def add_header_column(data):
    data['header'] = data.apply(generate_header, axis=1)
    return data[['header', 'seq']]

def generate_header(row):
    return f'>{row["chr"]}_{row["start_pos"]}_{row["depths"]}_INS_{row["lenght"]}'

def save_to_fna(data, output_file):
    data.to_csv(output_file, sep='\n', header=False, index=False)
