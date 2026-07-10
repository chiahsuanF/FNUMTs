#!/usr/bin/env python3

import os, psutil, glob, time, resource, argparse, re, sys, subprocess, pysam, traceback
import pandas as pd
from multiprocessing import Pool, cpu_count

import fnumts_check, fnumts_cluster, fnumts_bedtools, fnumts_del, fnumts_writevcf, fnumts_sortcigar, fnumts_process, fnumts_calculate, fnumts_generate


t = time.time()

args = parse_args()
writepath = args.outpath

# cluster collecting 
check_input_arguments(args)
make_output_directory(writepath)

chr_list2 = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

logfile_path = os.path.join(writepath, 'login_note')
with open(logfile_path, 'w') as logfile:
    logfile.write('Hello, welcome to the world of fnumts \n----------------------------------\n')
    # Record memory usage at start
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / (1024 ** 2)
    logfile.write(f'Initial memory used: {mem_mb:.2f} MB\n')
    logfile.write('----------------------------------\n\nStart collecting clusters ...\n')

step_start = time.time()

if args.bam:
    run_parallel(args.bam, args.outpath, args, is_cram=False)
elif args.cram:
    run_parallel(args.cram, args.outpath, args, is_cram=True, ref=args.cram_reference)
else:
    print("Please provide either a BAM or CRAM file.")

sort_and_filter_cigar_data(args, writepath)
count_del(args, writepath)
run_parallel_merge(writepath, args, chr_list2, threads=args.threads)

run_parallel_intersections(writepath, threads=args.threads)

input_55_ins = os.path.join(writepath, '55_bedtools_f.txt')
output_58_ins = os.path.join(writepath, '58_clusters.ins.sorted.txt')
processed_data = read_and_process_data(input_55_ins, chr_list2)
save_sorted_data(processed_data, output_58_ins)

input_55_del = os.path.join(writepath, '55_bedtools_f_del.txt')
output_59_del = os.path.join(writepath, '59_pre_filename_del.bed')
processed_data = read_and_process_data_del(input_55_del, chr_list2)
save_sorted_data(processed_data, output_59_del)


input_file_58 = output_58_ins
output_file_59 = os.path.join(writepath, '59_pre_filename.fna')
pre_data = read_data(input_file_58)
processed_data = add_header_column(pre_data)
save_to_fna(processed_data, output_file_59)


os.chdir(writepath)
os.system(f"bwa mem -T0 -t {args.threads} {args.reference} 59_pre_filename.fna > clusters.ins.srot.sam")

process = psutil.Process()
mem_used = process.memory_info().rss / (1024 ** 2)
elapsed = time.time() - step_start

with open('login_note', 'a') as logfile:
    logfile.write('End collecting clusters ...\n')
    logfile.write(f'## Time used: {elapsed:.2f} sec\n')
    logfile.write(f'## Memory used: {mem_used:.2f} MB\n')
    logfile.write(f'## Total runtime so far: {time.time() - t:.2f} sec\n\n')
    logfile.write('----------------------------------\n\nStart consensus ...\n')

step_2 = time.time()

ins_file = '60_INS.txt'
sc_file = '60_SC.txt'

with pysam.AlignmentFile("clusters.ins.srot.sam", "r") as samfile:
    with open(ins_file, 'w') as ins_outfile, open(sc_file, 'w') as sc_outfile:
        for read in samfile.fetch():
            if read.tid != -1 and samfile.getrname(read.tid) == 'chrM':
                line = f'{read.query_name}\n'
                if 'INS' in read.query_name:
                    ins_outfile.write(line)
                elif 'SC' in read.query_name:
                    sc_outfile.write(line)

pbsw1_file_wt = '60_SC.txt'
pbsw1_file_bs = '60_INS.txt'
read2_file = '58_clusters.ins.sorted.txt'

fastalist2 = preprocess_and_write_wt(pbsw1_file_wt, read2_file)
fastalist1 = preprocess_and_write_bs(pbsw1_file_bs, read2_file)

print("\n=== Running wtdbg2 / wtpoa-cns (WT) ===")
run_parallel_task(fastalist2, run_wtdbg2, args.threads)

print("\n=== Running bsalign (BS) ===")
run_parallel_task(fastalist1, run_bsalign, args.threads)

print("\nAll consensus jobs finished.\n")

process = psutil.Process()
mem_used = process.memory_info().rss / (1024 ** 2)
elapsed = time.time() - step_2

with open('login_note', 'a') as logfile:
    logfile.write('End consensus ...\n')
    logfile.write(f'## Time used: {elapsed:.2f} sec\n')
    logfile.write(f'## Memory used: {mem_used:.2f} MB\n')
    logfile.write(f'## Total runtime so far: {time.time() - t:.2f} sec\n\n')
    logfile.write('----------------------------------\n\nStart deletion mapping ...\n')

step_3 = time.time()

fasta = args.deletion_reference
bed = "59_pre_filename_del.bed"
output = "60_del.fna"

extract_subsequences_from_bed(fasta, bed, output)

os.system(f"bwa mem -T0 -t {args.threads} {args.reference} 60_del.fna > clusters.del.srot.sam")

with pysam.AlignmentFile("clusters.del.srot.sam", "r") as samfile:
    with open('61_DEL.fa', 'w') as del_outfile:
        for read in samfile.fetch():
            if read.tid != -1 and samfile.getrname(read.tid) == 'chrM':
                line = f'>{read.query_name}\n{read.query_sequence}\n'
                del_outfile.write(line)

process = psutil.Process()
mem_used = process.memory_info().rss / (1024 ** 2)
elapsed = time.time() - step_3

with open('login_note', 'a') as logfile:
    logfile.write('End deletion mapping ...\n')
    logfile.write(f'## Time used: {elapsed:.2f} sec\n')
    logfile.write(f'## Memory used: {mem_used:.2f} MB\n')
    logfile.write(f'## Total runtime so far: {time.time() - t:.2f} sec\n\n')
    logfile.write('----------------------------------\n\nStart double check candidates ...\n')
step_4 = time.time()
os.system("cat *.rename.fa 61_DEL.fa > temp_final.align.fa")

os.system("bwa mem -T0 -t4 -k20 " + args.reference + " temp_final.align.fa > fnumts.sam")

filter_best_primary_alignments("fnumts.sam", "fnumts2.sam")

sam_file = "fnumts2.sam"
paf_file = "fnumts.paf"
vcf_tmp = "fnumts.vcf.tmp"
snp_file = "70_.txt"
processed_snp_file = "71_.txt"

passing_reads = sam_to_paf(sam_file, paf_file, args.type)
with open(snp_file, "w") as snp1:
    for read in passing_reads:
        process_read(read, snp1)
snp_df = pd.read_csv(snp_file, sep="\t", header=None)
snp_df = snp_df.replace("_", " ", regex=True)
snp_df.to_csv(processed_snp_file, sep="\t", header=None, index=False)

snp_df2 = pd.read_csv(
    processed_snp_file,
    sep=r"\s+",
    header=None,
    engine="python",
    skip_blank_lines=True,
    on_bad_lines="skip"
)
snp_df2.columns = ["chr", "start_pos", "depth", "type", "aveMQ",
                   "seq", "MQ", "M_start", "NM"]
with open(vcf_tmp, "w") as vcf_out:
    for idx, row in snp_df2.iterrows():
        if row["chr"] != "chrM":
            vcf_out.write(generate_vcf_row(row, idx + 1) + "\n")
if args.bam:
    bam_filename = os.path.join(args.base_dir, args.bam)
elif args.cram:
    bam_filename = os.path.join(args.base_dir, args.cram)
else:
    raise ValueError("Provide either -b <BAM> or -c <CRAM>.")

if not os.path.exists(bam_filename):
    raise FileNotFoundError(f"Missing file: {bam_filename}")

print(f"Using BAM/CRAM file: {bam_filename}")
output_vcf_cf = "fnumts.cf.vcf.tmp"
process_vcf_and_filter(vcf_tmp, bam_filename, output_vcf_cf, args.type, reference=args.cram_reference)
fa_df = pd.read_csv(output_vcf_cf, sep="\t", header=None)

with open(args.output_vcfname + ".fasta", "w") as fasta_out:
    for _, row in fa_df.iterrows():
        info = row[7]

        if "TYPE=INS" in info:
            seq = row[4]
        elif "TYPE=DEL" in info:
            seq = row[3]
        else:
            seq = row[4]

        fasta_out.write(f">{row[0]}_{row[1]}\n{seq}\n")
process_vcf(output_vcf_cf, "fnumt.filter", args.type)
process = psutil.Process()
mem_used = process.memory_info().rss / (1024 ** 2)
elapsed = time.time() - step_3

with open("login_note", "a") as log:
    log.write("End double check candidates ...\n")
    log.write(f"## Time used: {elapsed:.2f} sec\n")
    log.write(f"## Memory used: {mem_used:.2f} MB\n")
    log.write(f"## Total runtime so far: {time.time() - t:.2f} sec\n\n")
    log.write("----------------------------------\n\nStart generating VCF file ...\n")
df = pd.read_csv("fnumt.filter", sep="\t", header=None)
df[0] = pd.Categorical(df[0], categories=chr_list2, ordered=True)
df = df.sort_values(by=[0, 1])

# drop old ID
df = df.drop(columns=[2])

# create new ID
df.insert(2, "ID", [f"FN{i}" for i in range(1, len(df) + 1)])
df.to_csv("fnumt.filter.sorted", sep="\t", index=False, header=False)
with open("fnumt.filter.sorted") as fh:
    sorted_text = fh.read()

updated_text = add_format_and_sample_data(sorted_text)
chromosomes = extract_chromosomes(sorted_text)
header_obj = vcf_header(chromosomes, bam_filename)

with open(args.output_vcfname + ".vcf", "w") as out_vcf:
    out_vcf.write(str(header_obj))
    out_vcf.write(updated_text)

os.system('rm -r poa_output *.wtdbg.*')

process = psutil.Process()
mem_used = process.memory_info().rss / (1024 ** 2)
elapsed = time.time() - step_4

with open('login_note', 'a') as logfile:
    logfile.write('End generating VCF file ...\n')
    logfile.write(f'## Time used: {elapsed:.2f} sec\n')
    logfile.write(f'## Memory used: {mem_used:.2f} MB\n')
    logfile.write(f'## Total runtime so far: {time.time() - t:.2f} sec\n\n')
    logfile.write('----------------------------------\n    Thank you for using fnumts\n Fun fact, workout for 30 minutes\n')
    logfile.write(' everyday is good for your health!\n----------------------------------\n\n')
    logfile.write('==================================\n############ Summary #############\n==================================\n')

    total_time = time.time() - t

    current_mem_mb = process.memory_info().rss / (1024 ** 2)
    peak_mem_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mem_mb = peak_mem_kb / 1024
    total_mem_mb = peak_mem_mb + current_mem_mb

    logfile.write(f'Total time: {total_time:.2f} sec\n')
    logfile.write(f'Peak memory used: {peak_mem_mb:.2f} MB\n')
    logfile.write(f'Total memory used (approx.): {total_mem_mb:.2f} MB\n')
    logfile.write('==================================\n##################################\n==================================\n')

