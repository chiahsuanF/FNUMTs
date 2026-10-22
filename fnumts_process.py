#!/usr/bin/env python3

import os, subprocess
import pandas as pd
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm


def extract_and_sort_columns(df):
    selected_columns = ["chr", "sp1", "depth", "INS",  "avg_MQ", "len", "seq"]
    sorted_df = df[selected_columns].sort_values(by=["chr", "sp1", "len"])
    return sorted_df

def read_and_process_data(input_file, chr_list):
    df = pd.read_csv(input_file, sep="\t", header=None)
    df.columns = ["chr", "sp1", "ep1", "depth", "chr2", "sp2", "ep2", "INS", "MQ", "len", "seq"]

    # 🔥 recompute support reads (this is what you want)
    df['support_reads'] = df.groupby(['chr', 'sp1', 'ep1']).transform('size')

    # replace old depth
    df['depth'] = df['support_reads']

    # average MQ (your existing logic)
    df['avg_MQ'] = df.groupby(['chr', 'sp1'])['MQ'].transform('mean').astype(int)

    sorted_df = extract_and_sort_columns(df)

    sorted_df["chr"] = sorted_df["chr"].astype("category")
    sorted_df["chr"] = sorted_df["chr"].cat.set_categories(chr_list)

    return sorted_df

def read_data(input_file):
    pre_data = pd.read_csv(input_file, sep="\t", header=None)
    pre_data.columns = ["chr", "start_pos", "depths", "INS", "MQ", "lenght", "seq"]
    return pre_data

def preprocess_and_write_wt(pbsw_file, read2_file):


    if os.stat(pbsw_file).st_size == 0:
        print(f"{pbsw_file} is empty. Skipping preprocessing.")
        return []

    pbsw1 = pd.read_csv(pbsw_file, sep='\t', header=None).replace('_', ' ', regex=True)
    pbsw1.to_csv("61_SC.txt", sep="\t", header=None, index=False)
    pbsw1 = pd.read_csv('61_SC.txt', sep="\t| ", engine='python', header=None)
    pbsw1.columns = ['col1', 'col2', 'col3', 'col4', 'col5', 'col6']

    read2 = pd.read_csv(read2_file, sep='\t', header=None)
    read2.columns = ['col1', 'col2', 'depth', 'INS', 'len', 'MQ', 'seq']

    pbsw_index = 0
    read2_index = 0

    seq_dict = {}   # file_name -> list of fasta records

    while pbsw_index < len(pbsw1) and read2_index < len(read2):
        pbsw_line = pbsw1.iloc[pbsw_index]
        read2_line = read2.iloc[read2_index]

        chr_val_pbsw = pbsw_line["col1"]
        start_pos_pbsw = pbsw_line["col2"]
        depth_pbsw = pbsw_line["col3"]
        type_pbsw = pbsw_line["col4"]
        mq_pbsw = pbsw_line["col5"]

        chr_val_read2 = read2_line["col1"]
        start_pos_read2 = read2_line["col2"]
        type2 = read2_line["INS"]
        seq_read2 = read2_line["seq"]

        if chr_val_pbsw == chr_val_read2 and start_pos_pbsw == start_pos_read2 and type_pbsw == type2:
            file_name = f"{chr_val_pbsw}_{start_pos_pbsw}_{depth_pbsw}_{type_pbsw}_{mq_pbsw}.wt.fa"

            header = f">{chr_val_pbsw}_{start_pos_pbsw}_{pbsw_index}_{type_pbsw}_{mq_pbsw}"
            record = f"{header}\n{seq_read2}\n"

            seq_dict.setdefault(file_name, []).append(record)

            pbsw_index += 1
            read2_index += 1

        elif chr_val_pbsw < chr_val_read2 or (
            chr_val_pbsw == chr_val_read2 and start_pos_pbsw < start_pos_read2
        ):
            pbsw_index += 1
        else:
            read2_index += 1

    output_files = []
    for file_name, records in seq_dict.items():
        if len(records) >= 2:
            with open(file_name, "w") as f:
                f.writelines(records)
            output_files.append(file_name)

    return output_files



def preprocess_and_write_bs(pbsw_file, read2_file):
    if os.stat(pbsw_file).st_size == 0:
        print(f"{pbsw_file} is empty. Skipping preprocessing.")
        return []

    pbsw1 = pd.read_csv(pbsw_file, sep='\t', header=None).replace('_', ' ', regex=True)
    pbsw1.to_csv("61_INS.txt", sep="\t", header=None, index=False)
    pbsw1 = pd.read_csv('61_INS.txt', sep="\t| ", engine='python', header=None)
    pbsw1.columns = ['col1', 'col2', 'col3', 'col4', 'col5', 'col6']

    read2 = pd.read_csv(read2_file, sep='\t', header=None)
    read2.columns = ['col1', 'col2', 'depth', 'INS', 'len', 'MQ', 'seq']

    open_files = {}
    pbsw_index = 0
    read2_index = 0
    output_files = []

    while pbsw_index < len(pbsw1) and read2_index < len(read2):
        pbsw_line = pbsw1.iloc[pbsw_index]
        read2_line = read2.iloc[read2_index]

        chr_val_pbsw = pbsw_line["col1"]
        start_pos_pbsw = pbsw_line["col2"]
        depth_pbsw = pbsw_line["col3"]
        btyep = pbsw_line["col4"]
        mq_pbsw = pbsw_line["col5"]

        chr_val_read2 = read2_line["col1"]
        start_pos_read2 = read2_line["col2"]
        type2 = read2_line["INS"]
        seq_read2 = read2_line["seq"]

        if chr_val_pbsw == chr_val_read2 and start_pos_pbsw == start_pos_read2 and btyep == type2:
            file_name = f"{chr_val_pbsw}_{start_pos_pbsw}_{depth_pbsw}_{btyep}_{mq_pbsw}.bs.fa"
            if file_name not in open_files:
                open_files[file_name] = open(file_name, 'a')
                output_files.append(file_name)

            open_files[file_name].write(f'>{chr_val_pbsw}_{start_pos_pbsw}_{pbsw_index}_{mq_pbsw}\n{seq_read2}\n')

            pbsw_index += 1
            read2_index += 1
        elif chr_val_pbsw < chr_val_read2 or (chr_val_pbsw == chr_val_read2 and start_pos_pbsw < start_pos_read2):
            pbsw_index += 1
        else:
            read2_index += 1

    for f in open_files.values():
        f.close()

    return output_files

def run_cmd(cmd):
    """Run a shell command safely."""
    subprocess.run(cmd, shell=True, check=True)

def run_wtdbg2(file, threads):
    prefix = file[:-3]  # remove .fa
    lay = f"{prefix}.wtdbg.ctg.lay.gz"

    try:
        subprocess.run(
            f"wtdbg2 -t {threads} -q -i {file} -fo {prefix}.wtdbg",
            shell=True,
            check=True
        )

        if not os.path.exists(lay) or os.path.getsize(lay) < 1000:
            print(f"[SKIP] Invalid layout: {lay}")
            return prefix, "NO_LAYOUT"

        subprocess.run(
            f"wtpoa-cns -t 1 -q -i {lay} -fo {prefix}.align.fa",
            shell=True,
            check=True
        )

        subprocess.run(
            f"sed -E 's/ctg[0-9]+ len=[0-9]+/{prefix[:-3]}/' "
            f"{prefix}.align.fa > {prefix}.rename.fa",
            shell=True,
            check=True
        )

        return prefix, "OK"

    except subprocess.CalledProcessError as e:
        print(f"[FAIL] wtpoa-cns failed for {prefix}")
        return prefix, "FAILED"

    except Exception as e:
        print(f"[ERROR] Unexpected error for {prefix}: {e}")
        return prefix, "ERROR"

def run_bsalign(file, threads):
    if not os.path.exists(file) or os.stat(file).st_size == 0:
        print(f"Skipping empty/missing file: {file}")
        return
    prefix = file[:-3]
    run_cmd(f"bsalign poa {file} -o {prefix}.align.fa")
    run_cmd(f"sed -E 's/cns_seq/{prefix[:-3]}/' {prefix}.align.fa > {prefix[:-9]}.rename.fa")

def run_parallel_task(file_list, task_func, threads):
    with Pool(processes=threads) as pool:
        list(tqdm(pool.imap_unordered(partial(task_func, threads=threads), file_list),
                  total=len(file_list),
                  desc="Processing consensus"))
