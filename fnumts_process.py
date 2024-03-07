#!/usr/bin/env python3

import pandas as pd

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
  
def read_data(input_file):
    pre_data = pd.read_csv(input_file, sep="\t", header=None)
    pre_data.columns = ["chr", "start_pos", "depths", "INS", "lenght", "seq"]
    return pre_data

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
