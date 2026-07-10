#!/usr/bin/env python3

import os, psutil, glob, time, resource, argparse, re, sys, subprocess, pysam, traceback
import pandas as pd
from multiprocessing import Pool, cpu_count

import fnumts_check, fnumts_cluster, fnumts_bedtools, fnumts_del, fnumts_writevcf, fnumts_sortcigar, fnumts_process, fnumts_calculate, fnumts_generate

# ---------------------------------------------------------
# Step 1: Setup
# ---------------------------------------------------------
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

# ---------------------------------------------------------
# Step 2: Process BAM or CRAM
# ---------------------------------------------------------

step_start = time.time()

if args.bam:
    run_parallel(args.bam, args.outpath, args, is_cram=False)
elif args.cram:
    run_parallel(args.cram, args.outpath, args, is_cram=True, ref=args.cram_reference)
else:
    print("Please provide either a BAM or CRAM file.")

""" old, but works
if args.bam:
    process_bam_file(args.bam, args.outpath, args)
elif args.cram:
    process_cram_file(args.cram, args.outpath, args)
else:
    print("Please provide either a BAM or CRAM file.")
"""
# ---------------------------------------------------------
# Step 3: Sort and Filter
# ---------------------------------------------------------
sort_and_filter_cigar_data(args, writepath)
count_del(args, writepath)
run_parallel_merge(writepath, args, chr_list2, threads=args.threads)
# ---------------------------------------------------------
# Step 4: Insertions
# ---------------------------------------------------------
# main_merge_and_intersect(writepath, chr_list2, args)
# ins
"""#
input_53_srt = os.path.join(writepath, '53_cut_3.srt.txt')
output_54_sdtin = os.path.join(writepath, '54_sdtin.txt')
process_intervals(input_53_srt, output_54_sdtin, chr_list2, max_distance=int(args.distance), min_support=args.support_read, n_workers=args.threads)
"""
"""
#####
input_53_srt = os.path.join(writepath, '53_cut_3.srt.txt')
output_54_sdtin = os.path.join(writepath, '54_sdtin.txt')

bed_input = pd.read_csv(input_53_srt, sep='\t', header=None)
bed = bed_input.values.tolist()
merged_output = merge_intervals(bed, int(args.distance))
output1 = pd.DataFrame(merged_output)
output2 = output1[output1[3] >= args.support_read]
output2.columns = ["chr", "startp", "endp", "depth"]
output3 = output2[output2["chr"].isin(chr_list2)].sort_values(['chr', 'startp'])
output3.to_csv(output_54_sdtin, sep='\t', header=None, index=False)
"""
# ---------------------------------------------------------
# Step 5: Deletions
# ---------------------------------------------------------
# del
""" #
input_53_del = os.path.join(writepath, '53_cut_3_del.txt')
output_54_sdtin_del = os.path.join(writepath, '54_sdtin_del.txt')
process_intervals(input_53_del, output_54_sdtin_del, chr_list2, max_distance=int(args.distance), min_support=args.support_read, n_workers=args.threads)
"""
"""
input_53_del = os.path.join(writepath, '53_cut_3_del.txt')
output_54_sdtin_del = os.path.join(writepath, '54_sdtin_del.txt')

bed_input = pd.read_csv(input_53_del, sep='\t', header=None)
bed = bed_input.values.tolist()
merged_output = merge_intervals(bed, int(args.distance))
output1 = pd.DataFrame(merged_output)
output2 = output1[output1[3] >= args.support_read]
output2.columns = ["chr", "startp", "endp", "depth"]
output3 = output2[output2["chr"].isin(chr_list2)].sort_values(['chr', 'startp'])
output3.to_csv(output_54_sdtin_del, sep='\t', header=None, index=False)
"""
# ---------------------------------------------------------
# Step 6: Intersections (bedtools-like)
# ---------------------------------------------------------
run_parallel_intersections(writepath, threads=args.threads)
"""
# insertions
a_intervals = read_intervals(output_54_sdtin)
b_intervals = read_intervals(input_53_srt)
intersecting_intervals = bed_intersect(a_intervals, b_intervals)
pd.DataFrame(intersecting_intervals).to_csv(
    os.path.join(writepath, '55_bedtools_f.txt'), sep='\t', header=None, index=None
)

# deletions
a_intervals = read_intervals(output_54_sdtin_del)
b_intervals = read_intervals(input_53_del)
intersecting_intervals = bed_intersect(a_intervals, b_intervals)
pd.DataFrame(intersecting_intervals).to_csv(
    os.path.join(writepath, '55_bedtools_f_del.txt'), sep='\t', header=None, index=None
)
"""
# ---------------------------------------------------------
# Step 7: Cluster Sorting
# ---------------------------------------------------------

# insertions
input_55_ins = os.path.join(writepath, '55_bedtools_f.txt')
output_58_ins = os.path.join(writepath, '58_clusters.ins.sorted.txt')
processed_data = read_and_process_data(input_55_ins, chr_list2)
save_sorted_data(processed_data, output_58_ins)

# deletions
input_55_del = os.path.join(writepath, '55_bedtools_f_del.txt')
output_59_del = os.path.join(writepath, '59_pre_filename_del.bed')
processed_data = read_and_process_data_del(input_55_del, chr_list2)
save_sorted_data(processed_data, output_59_del)


# ---------------------------------------------------------
# Step 8: Sequence Extraction
# ---------------------------------------------------------

# insertions
input_file_58 = output_58_ins
output_file_59 = os.path.join(writepath, '59_pre_filename.fna')
pre_data = read_data(input_file_58)
processed_data = add_header_column(pre_data)
save_to_fna(processed_data, output_file_59)
"""
# deletions (unique positions)
sort_and_uniq(
    output_58_del,
    os.path.join(writepath, '59_pre_filename_del.bed'),
    uniq_cols=[0, 1, 2]
)
"""

os.chdir(writepath)
os.system(f"bwa mem -T0 -t {args.threads} {args.reference} 59_pre_filename.fna > clusters.ins.srot.sam")
# os.system("bwa mem -T0 -t4 " + args.reference + " 59_pre_filename.fna > clusters.ins.srot.sam")
# os.system("/share/apps/rc/software/BWA/0.7.17-GCC-8.3.0/bin/bwa mem -T20 -t4 " + args.reference + " 59_pre_filename.fna > clusters.ins.srot.sam")

"""
samfile1 = pysam.AlignmentFile("clusters.ins.srot.sam", "r")
with open('60_.txt', 'w') as pbs1:
    for read1 in samfile1.fetch():
        if read1.tid != -1 and samfile1.getrname(read1.tid) == 'chrM':
            pbs1.write(read1.query_name + '\n')
"""

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
# ins
#60_ins 60_sc

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

# wtdbg2
"""
pbsw1_file_wt = '60_SC.txt'
pbsw1_file_bs = '60_INS.txt'
read2_file = '58_clusters.ins.sorted.txt'
preprocess_and_write_wt(pbsw1_file_wt, read2_file)
preprocess_and_write_bs(pbsw1_file_bs, read2_file)
"""
pbsw1_file_wt = '60_SC.txt'
pbsw1_file_bs = '60_INS.txt'
read2_file = '58_clusters.ins.sorted.txt'

# Generate WT/BS files

fastalist2 = preprocess_and_write_wt(pbsw1_file_wt, read2_file)
fastalist1 = preprocess_and_write_bs(pbsw1_file_bs, read2_file)
"""
fastalist2 = preprocess_and_write_parallel(pbsw1_file_wt, read2_file, threads=args.threads)
fastalist1 = preprocess_and_write_parallel(pbsw1_file_bs, read2_file, threads=args.threads)
"""
# Run consensus in parallel
print("\n=== Running wtdbg2 / wtpoa-cns (WT) ===")
run_parallel_task(fastalist2, run_wtdbg2, args.threads)

print("\n=== Running bsalign (BS) ===")
run_parallel_task(fastalist1, run_bsalign, args.threads)

print("\nAll consensus jobs finished.\n")


""" # old consuses
#old, but correct
os.system("ls *_SC*.wt.fa > fastalist2") 
fastalist2=open("fastalist2",'r').read().split('\n')[:-1]
counter = 0
for k in fastalist2:
    counter = counter +1
    #os.system(cmd + ' poa '+ k +' -o ' +  k[:-3] + ".align.fa ")
    os.system("wtdbg2 -q -i " + k + " -fo " + k[:-3] + ".wtdbg")
    os.system("wtpoa-cns -q -i " + k[:-3] + ".wtdbg.ctg.lay.gz  -fo "+ k[:-3] + ".align.fa")
    os.system("sed 's/ctg[0-9]\+ len=[0-9]\+/"  +k[:-6] + "/' " + k[:-2] + "align.fa > "+ k[:-3] + ".rename.fa")

# bsalign

#old, but correct
os.system("ls *_INS*.bs.fa > fastalist1")
fastalist1=open("fastalist1",'r').read().split('\n')[:-1]
counter = 0
# cmd = '/home/fanc/path/bsalign'
for k in fastalist1:
    counter = counter +1
    os.system('bsalign poa '+ k +' -o ' +  k[:-3] + ".align.fa ")
    os.system("sed 's/cns_seq/"  +k[:-6] + "/' < " + k[:-2] + "align.fa > "+ k[:-9] + ".rename.fa")
"""
####
process = psutil.Process()
mem_used = process.memory_info().rss / (1024 ** 2)
elapsed = time.time() - step_2

with open('login_note', 'a') as logfile:
    logfile.write('End consensus ...\n')
    logfile.write(f'## Time used: {elapsed:.2f} sec\n')
    logfile.write(f'## Memory used: {mem_used:.2f} MB\n')
    logfile.write(f'## Total runtime so far: {time.time() - t:.2f} sec\n\n')
    logfile.write('----------------------------------\n\nStart deletion mapping ...\n')
####

step_3 = time.time()

# del
fasta = args.deletion_reference
bed = "59_pre_filename_del.bed"
output = "60_del.fna"

extract_subsequences_from_bed(fasta, bed, output)

os.system(f"bwa mem -T0 -t {args.threads} {args.reference} 60_del.fna > clusters.del.srot.sam")
# os.system("bwa mem -T0 -t4 " + args.reference + " 60_del.fna > clusters.del.srot.sam")

# del_file = '61_DEL.fa'

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
# cat del.fa and temp_final.align.fa
os.system("cat *.rename.fa 61_DEL.fa > temp_final.align.fa")
# os.system("cat *.rename.fa > temp_final.align.fa")

#
"""
def convert_fasta_to_single_line(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        seq = ''
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    outfile.write(seq + '\n')
                outfile.write(line + '\n')
                seq = ''
            else:
                seq += line
        if seq:
            outfile.write(seq + '\n')

# Input and output file paths
input_fasta = 'temp_final.align.fa'
output_fasta = 'temp_final_signal.align.fa'

convert_fasta_to_single_line(input_fasta, output_fasta)


pre_data_2 = pd.read_csv('temp_final_signal.align.fa', header = None)
pre_data_2.columns = ["col1"]

with open("final.align.fa", 'w') as svfile:
    for l in range(0, len(pre_data_2),2):
        m = l +1
        len1 = len(pre_data_2["col1"][m])
        svfile.write('{}_{}\n{}\n'.format(pre_data_2["col1"][l], str(len1), pre_data_2["col1"][m]))
svfile.close
"""
"""
if args.type == 'hifi': 
    #os.system("bwa mem -T0 -t4 -k20 " + args.reference2 + " final.align.fa > fnumts.sam")
    os.system("bwa mem -T0 -k20 -t4 " + args.reference1 + " final.align.fa > fnumts.sam")
    #os.system("/share/apps/rc/software/BWA/0.7.17-GCC-8.3.0/bin/bwa mem -T20 -t4 " + args.reference + " final.align.fa > fnumts.sam")
elif args.type == 'ont':
    os.system("bwa mem -T0 -k20 -t4 " + args.reference1 + " final.align.fa > fnumts.sam") 
    #os.system("/share/apps/rc/software/BWA/0.7.17-GCC-8.3.0/bin/bwa mem -T20 -B0 " + args.reference + " final.align.fa > fnumts.sam")    
elif args.type == 'cram':
    os.system("bwa mem -T10 -B3 -t4 " + args.reference2 + " final.align.fa > fnumts.sam")
   """ 

#os.system("bwa mem -T0 -t4 -k20 " + args.reference + " temp_final.align.fa > fnumts.sam")
os.system("bwa mem -T0 -t4 -k20 " + args.reference + " temp_final.align.fa > fnumts.sam")

filter_best_primary_alignments("fnumts.sam", "fnumts2.sam")

"""# 

samfile = pysam.AlignmentFile("fnumts.sam", "r")
iter2 = samfile.fetch()
####

with open('70_.txt', 'w') as snp1:
    for read2 in iter2:
        process_read(read2, snp1)

####

snp2 = pd.read_csv('70_.txt', sep = '\t', header = None)
snp2 = snp2.replace('_', ' ', regex=True)
snp2.to_csv( "71_.txt", sep = "\t", header=None, index=0)

snp3 = pd.read_csv('71_.txt', sep="\t| ", engine='python', header=None)
snp3.columns = ["chr", "start_pos", "depth", "lenght", "seq", "MQ", "M_start", "NM"]



with open('fnumts.vcf.tmp', 'w') as snp4:
    for id_value, (index, row) in enumerate(snp3.iterrows(), start=1):
        if row["chr"] != "chrM":  # Check if chr is not chrM
            vcf_row = generate_vcf_row(row, id_value)
            snp4.write(vcf_row)
            snp4.write("\n")
"""

# ---------------------------------------------
# File paths
# ---------------------------------------------
sam_file = "fnumts2.sam"
paf_file = "fnumts.paf"
vcf_tmp = "fnumts.vcf.tmp"
snp_file = "70_.txt"
processed_snp_file = "71_.txt"

# ---------------------------------------------
# Step 1–2: Convert SAM → PAF and get passing reads
# ---------------------------------------------
passing_reads = sam_to_paf(sam_file, paf_file, args.type)

# ---------------------------------------------
# Step 3: Write passing reads to temporary SNP file
# ---------------------------------------------
with open(snp_file, "w") as snp1:
    for read in passing_reads:
        process_read(read, snp1)

# ---------------------------------------------
# Step 4: Pandas preprocessing
# ---------------------------------------------
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


# ---------------------------------------------
# Step 5: Generate VCF entries
# ---------------------------------------------
with open(vcf_tmp, "w") as vcf_out:
    for idx, row in snp_df2.iterrows():
        if row["chr"] != "chrM":
            vcf_out.write(generate_vcf_row(row, idx + 1) + "\n")

# ---------------------------------------------
# Step 6: Select BAM/CRAM file
# ---------------------------------------------
if args.bam:
    bam_filename = os.path.join(args.base_dir, args.bam)
elif args.cram:
    bam_filename = os.path.join(args.base_dir, args.cram)
else:
    raise ValueError("Provide either -b <BAM> or -c <CRAM>.")

if not os.path.exists(bam_filename):
    raise FileNotFoundError(f"Missing file: {bam_filename}")

print(f"Using BAM/CRAM file: {bam_filename}")

# ---------------------------------------------
# Step 7: Filter VCF using alignments
# ---------------------------------------------
output_vcf_cf = "fnumts.cf.vcf.tmp"
process_vcf_and_filter(vcf_tmp, bam_filename, output_vcf_cf, args.type, reference=args.cram_reference)

# ---------------------------------------------
# Step 8: Create FASTA from VCF records
# ---------------------------------------------
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

# ---------------------------------------------
# Step 9: Final filtering
# ---------------------------------------------
process_vcf(output_vcf_cf, "fnumt.filter", args.type)

# ---------------------------------------------
# Step 10: Logging
# ---------------------------------------------
process = psutil.Process()
mem_used = process.memory_info().rss / (1024 ** 2)
elapsed = time.time() - step_3

with open("login_note", "a") as log:
    log.write("End double check candidates ...\n")
    log.write(f"## Time used: {elapsed:.2f} sec\n")
    log.write(f"## Memory used: {mem_used:.2f} MB\n")
    log.write(f"## Total runtime so far: {time.time() - t:.2f} sec\n\n")
    log.write("----------------------------------\n\nStart generating VCF file ...\n")

# ---------------------------------------------
# Step 11: VCF sorting & reassign FN IDs
# ---------------------------------------------
df = pd.read_csv("fnumt.filter", sep="\t", header=None)
df[0] = pd.Categorical(df[0], categories=chr_list2, ordered=True)
df = df.sort_values(by=[0, 1])

# drop old ID
df = df.drop(columns=[2])

# create new ID
df.insert(2, "ID", [f"FN{i}" for i in range(1, len(df) + 1)])
df.to_csv("fnumt.filter.sorted", sep="\t", index=False, header=False)

# ---------------------------------------------
# Step 12: Add FORMAT/sample columns
# ---------------------------------------------
with open("fnumt.filter.sorted") as fh:
    sorted_text = fh.read()

updated_text = add_format_and_sample_data(sorted_text)

# ---------------------------------------------
# Step 13: Build final VCF header & write file
# ---------------------------------------------
chromosomes = extract_chromosomes(sorted_text)
header_obj = vcf_header(chromosomes, bam_filename)

with open(args.output_vcfname + ".vcf", "w") as out_vcf:
    out_vcf.write(str(header_obj))
    out_vcf.write(updated_text)


##############################
"""
# Example VCF file data processing
with open('fnumt.filter', 'r') as vcf_file:
    file_data = vcf_file.read()

chromosomes = extract_chromosomes(file_data)  # Extract chromosomes dynamically from the VCF file
if args.bam:
    bam_filename = os.path.join(args.base_dir, args.bam) 
    vcf_header_obj = vcf_header(chromosomes, args.bam) 
elif args.cram:
    bam_filename = os.path.join(args.base_dir, args.cram) 
    vcf_header_obj = vcf_header(chromosomes, args.cram) 
#vcf_header_obj = vcf_header(chromosomes, args.bam)      # Create VCF header with contig info

# Write the header and data to a new VCF file
with open('fnumts.vcf', 'w') as new_file:
    new_file.write(str(vcf_header_obj))  
    new_file.write(file_data) 
"""

########################################
# Add header and file data to a new file

#os.system('rm fnumts.vcf.tmp fnumts.cf.vcf.tmp fnumt.filter* fnumts.paf *.txt *.fa 59_pre_filename* *.wtdbg.* clusters.*.srot.sam 60_del.fna fnumts2.sam fnumts.sam')
# os.system('rm fnumts.vcf.tmp fnumts.cf.vcf.tmp fnumt.filter* fnumts.paf *.txt *.fa 59_pre_filename*  fastalist* *.wtdbg.* clusters.*.srot.sam 60_del.fna fnumts2.sam')
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

    # Current memory
    current_mem_mb = process.memory_info().rss / (1024 ** 2)
    # Peak memory
    peak_mem_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_mem_mb = peak_mem_kb / 1024

    # Approx. total memory used (sum-like measure)
    total_mem_mb = peak_mem_mb + current_mem_mb

    logfile.write(f'Total time: {total_time:.2f} sec\n')
    #logfile.write(f'Current memory used: {current_mem_mb:.2f} MB\n')
    logfile.write(f'Peak memory used: {peak_mem_mb:.2f} MB\n')
    logfile.write(f'Total memory used (approx.): {total_mem_mb:.2f} MB\n')
    logfile.write('==================================\n##################################\n==================================\n')

