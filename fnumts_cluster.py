#!/usr/bin/env python3

import pysam, os, multiprocessing

TARGET_CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# -----------------------------
# Helper functions
# -----------------------------

def write_ins(ct, line, s, r, length):
    seq = line.query_sequence[r:r+length]
    ct.write(f"{line.reference_name}\t{s+1}\t{s+1+length}\t{line.mapping_quality}\tINS\t{length}\t{seq}\n")

def write_del(cd, line, s, length):
    cd.write(f"{line.reference_name}\t{s+1}\t{s+length}\t{line.mapping_quality}\tDEL\t{length}\n")

def calculate_end(line):
    pos = line.reference_start

    if line.cigartuples:
        for op, length in line.cigartuples:
            if op in (0, 2, 7, 8):  # M, D, =, X
                pos += length

    return pos

def extract_soft_clips(line):
    r = 0
    scs = []

    if line.cigartuples:
        for op, length in line.cigartuples:
            if op == 4:  # S
                scs.append((r, length))
                r += length
            elif op in (0, 1, 7, 8):  # M, I, =, X
                r += length
            # D does not consume query
    return scs

def write_soft(ct, line, scs):
    end_pos = calculate_end(line)

    for r, length in scs:
        if r == 0:
            start = line.reference_start + 1
        else:
            start = end_pos + 1

        seq = line.query_sequence[r:r + length]
        ct.write(
            f"{line.reference_name}\t"
            f"{start}\t{start + length - 1}\t"
            f"{line.mapping_quality}\tSC\t{length}\t{seq}\n"
        )
        
def process_line(ct, cd, line, args):
    s = line.reference_start  # reference position (0-based)
    r = 0                     # query position

    if not line.cigartuples:
        return

    for op, length in line.cigartuples:

        if op == 4:  # S
            r += length

        elif op in (0, 7, 8):  # M, =, X
            s += length
            r += length

        elif op == 1:  # I
            if length <= args.max_length:
                write_ins(ct, line, s, r, length)
            r += length

        elif op == 2:  # D
            if length <= args.max_length:
                write_del(cd, line, s, length)
            s += length

# -----------------------------
# Worker function for multiprocessing
# -----------------------------

def worker(task):
    file_path, contig, out_prefix, args, is_cram, ref = task
    cigar_file = f"{out_prefix}_{contig}.cigar.txt"
    del_file = f"{out_prefix}_{contig}.del.txt"

    mode = "rc" if is_cram else "rb"
    aln = pysam.AlignmentFile(file_path, mode, reference_filename=ref)

    with open(cigar_file, "w") as ct, open(del_file, "w") as cd:
        for line in aln.fetch(contig):
            if line.flag & 256:  # skip secondary
                continue
            scs = extract_soft_clips(line)
            write_soft(ct, line, scs)
            process_line(ct, cd, line, args)

    aln.close()
    return cigar_file, del_file

# -----------------------------
# Main function to run parallel processing
# -----------------------------

def run_parallel(file_path, outpath, args, is_cram=False, ref=None):
    out_prefix = os.path.join(outpath, "51")
    tasks = []

    mode = "rc" if is_cram else "rb"
    aln = pysam.AlignmentFile(file_path, mode, reference_filename=ref)

    for contig in aln.references:
        if contig not in TARGET_CHROMOSOMES:
            continue
        tasks.append((file_path, contig, out_prefix, args, is_cram, ref))
    aln.close()

    # Run one CPU per chromosome (or limit to threads)
    with multiprocessing.Pool(processes=min(args.threads, len(tasks))) as pool:
        results = pool.map(worker, tasks)

    # Merge temporary files
    cigar_out = os.path.join(outpath, "51_cigar.txt")
    del_out = os.path.join(outpath, "51_cigar_del.txt")
    with open(cigar_out, "w") as ct, open(del_out, "w") as cd:
        for cig, de in results:
            with open(cig) as f:
                ct.write(f.read())
            with open(de) as f:
                cd.write(f.read())
            os.remove(cig)
            os.remove(de)

    print(f"FINISHED multiprocessing {'CRAM' if is_cram else 'BAM'}.")
