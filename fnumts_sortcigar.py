#!/usr/bin/env python3

import pandas as pd
from multiprocessing import Pool, cpu_count
import os, time, traceback, csv, tempfile, fnumts_bedtools

def log(msg):
    t = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{t}] {msg}")

def process_chr(tf_name_c):
    """
    Load a chromosome temp file, sort by start_pos, compute end_pos,
    and return rows as list of lists. Deletes temp file afterward.
    """
    tf_name, chrom = tf_name_c
    if os.path.getsize(tf_name) == 0:
        return []

    df = pd.read_csv(
        tf_name, sep="\t", header=None,
        names=["chr", "start_pos", "end", "MQ", "INS", "length", "seq"]
    )
    df = df.sort_values("start_pos")
    df["end_pos"] = df["start_pos"] + 1
    result = df[["chr","start_pos","end_pos","INS","MQ","length","seq"]].values.tolist()
    os.remove(tf_name)
    return result

# -------------------------
# Main function
# -------------------------
def sort_and_filter_cigar_data(args, writepath):
    """
    Hybrid memory-efficient version.
    - Small files: read fully in memory, filter, sort, write.
    - Large files: chunked read, per-chromosome temp files, parallel sort.
    """
    csv.field_size_limit(100000000)

    inpath = os.path.join(writepath, '51_cigar.txt')
    outpath = os.path.join(writepath, '53_cut_3.srt.txt')

    # Detect chromosomes present in the file
    first_column = pd.read_csv(inpath, sep="\t", header=None, usecols=[0], dtype=str)
    chr_list2 = sorted(first_column[0].unique(), key=lambda x:
                    (not x[3:].isdigit(), int(x[3:]) if x[3:].isdigit() else 1e9))
    """
    chr_list2 = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

    """

    if not os.path.exists(inpath):
        raise FileNotFoundError(inpath)

    file_size_gb = os.path.getsize(inpath) / 1e9
    min_len = int(args.length)

    # -------------------------
    # LIGHTWEIGHT version for small files
    # -------------------------
    if file_size_gb <= 2:  # threshold in GB, adjust as needed
        df = pd.read_csv(
            inpath, sep='\t', header=None,
            names=["chr", "start_pos", "end", "MQ", "INS", "length", "seq"],
            dtype={"chr": str},
            quoting=csv.QUOTE_NONE,
            engine='c'
        )
        df = df[df["length"] >= min_len].copy()
        df["chr"] = pd.Categorical(df["chr"], categories=chr_list2, ordered=True)
        df = df.sort_values(["chr", "start_pos"])
        df["end_pos"] = df["start_pos"] + 1
        df[["chr","start_pos","end_pos","INS","MQ","length","seq"]].to_csv(
            outpath, sep="\t", header=False, index=False
        )
        return outpath

    # -------------------------
    # CHUNKED + per-chromosome + multi-CPU for large files
    # -------------------------
    tmpfiles = {c: tempfile.NamedTemporaryFile(delete=False, mode="w") for c in chr_list2}
    chunksize = 1_000_000

    reader = pd.read_csv(
        inpath, sep='\t', header=None,
        names=["chr", "start_pos", "end", "MQ", "INS", "length", "seq"],
        dtype={"chr": str},
        quoting=csv.QUOTE_NONE,
        engine='c',
        chunksize=chunksize
    )

    for chunk in reader:
        chunk = chunk[chunk["length"] >= min_len]
        for c in chr_list2:
            sub = chunk[chunk["chr"] == c]
            if len(sub) > 0:
                sub.to_csv(tmpfiles[c].name, sep="\t", header=False, index=False, mode="a")

    for f in tmpfiles.values():
        f.close()

    # -------------------------
    # Parallel sorting per chromosome
    # -------------------------
    cpu_num = min(args.threads, len(chr_list2))
    with Pool(cpu_num) as pool:
        results = pool.map(process_chr, [(tmpfiles[c].name, c) for c in chr_list2])

    # -------------------------
    # Write final output
    # -------------------------
    with open(outpath, "w") as outfile:
        for res in results:
            if len(res) > 0:
                for row in res:
                    outfile.write("\t".join(map(str, row)) + "\n")

    return outpath

# wtdbg2
def run_parallel_task(task_list, func, max_threads=None):
    if not task_list:
        print("[INFO] No tasks, skipping")
        return [], 0, 0

    if max_threads is None:
        max_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", cpu_count()))

    n_workers = min(max_threads, len(task_list))
    print(f"[INFO] Running {n_workers} workers")

    start = time.time()
    with Pool(processes=n_workers) as p:
        results = p.map(func, task_list)
    runtime = time.time() - start

    print(f"[INFO] Completed in {runtime:.2f} sec")
    return results, n_workers, runtime

def run_wtdbg2(k):
    prefix = k.replace(".wt.fa", "")
    os.system(f"wtdbg2 -t 4 -q -i {k} -fo {prefix}.wtdbg")
    os.system(f"wtpoa-cns -t 4 -q -i {prefix}.wtdbg.ctg.lay.gz -fo {prefix}.align.fa")
    os.system(f"sed 's/ctg[0-9]\\+ len=[0-9]\\+/{prefix}/' {prefix}.align.fa > {prefix}.rename.fa")

# bsalign
def run_bsalign(k):
    prefix = k.replace(".bs.fa", "")
    out_align = f"{prefix}.align.fa"
    out_rename = f"{prefix}.rename.fa"

    # Step 1: bsalign poa
    os.system(f"bsalign poa {k} -o {out_align}")

    # Step 2: rename header
    os.system(f"sed 's/cns_seq/{prefix}/' {out_align} > {out_rename}")

    return out_rename

def process_bed_file(file_path, output_path, max_distance, support_read, chr_filter):
    """
    Read bed-like file, merge intervals, filter by support and chromosome list, write result.
    Returns the output_path on success.
    """
    # defensive: check input exists
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"process_bed_file: input does not exist: {file_path}")

    bed_input = pd.read_csv(file_path, sep='\t', header=None, engine='python', comment='#')
    bed = bed_input.values.tolist()

    # merge_intervals should be provided by your fnumts_bedtools module
    merged_output = merge_intervals(bed, int(max_distance))
    output1 = pd.DataFrame(merged_output)
    if output1.shape[1] < 4:
        # unexpected output from merge_intervals
        raise ValueError("merge_intervals returned unexpected shape: " + repr(output1.shape))

    output2 = output1[output1[3] >= int(support_read)].copy()
    output2.columns = ["chr", "startp", "endp", "depth"]
    output3 = output2[output2["chr"].isin(chr_filter)].sort_values(['chr', 'startp'])
    output3.to_csv(output_path, sep='\t', header=False, index=False)
    return output_path

def run_parallel_merge(writepath, args, chr_list2, threads=None):
    """
    Build tasks for merging intervals and run them in multiprocessing.
    Skip tasks whose input files are missing and log warnings instead of crashing.
    """
    tasks = []

    input_53_srt = os.path.join(writepath, '53_cut_3.srt.txt')
    output_54_sdtin = os.path.join(writepath, '54_sdtin.txt')
    if os.path.exists(input_53_srt):
        tasks.append((input_53_srt, output_54_sdtin, int(args.distance), args.support_read, chr_list2))
    else:
        log(f"[WARN] missing standard intervals input, skipping: {input_53_srt}")

    input_53_del = os.path.join(writepath, '53_cut_3_del.txt')
    output_54_sdtin_del = os.path.join(writepath, '54_sdtin_del.txt')
    if os.path.exists(input_53_del):
        tasks.append((input_53_del, output_54_sdtin_del, int(args.distance), args.support_read, chr_list2))
    else:
        log(f"[WARN] missing deletion intervals input, skipping: {input_53_del}")

    if not tasks:
        log("[INFO] No merge tasks to run (all inputs missing).")
        return []

    threads = threads or cpu_count()
    nprocs = min(max(1, threads), len(tasks))
    log(f"[INFO] Running merge with {nprocs} worker(s) on {len(tasks)} task(s)")

    outputs = []
    # Use apply_async and capture exceptions inside worker to avoid silent crash
    with Pool(processes=nprocs) as pool:
        async_results = []
        for t in tasks:
            async_results.append(pool.apply_async(_worker_safe_process_bed_file, t))

        for ar in async_results:
            try:
                res = ar.get()
                outputs.append(res)
                log(f"[OK] merge task finished: {res}")
            except Exception as e:
                # log traceback and continue
                log(f"[ERROR] merge task failed: {e}")
                tb = traceback.format_exc()
                log(tb)

    return outputs

def _worker_safe_process_bed_file(file_path, output_path, max_distance, support_read, chr_filter):
    """
    Wrapper executed in worker to call process_bed_file and surface errors.
    """
    try:
        return process_bed_file(file_path, output_path, max_distance, support_read, chr_filter)
    except Exception as e:
        # Convert to explicit exception so parent process sees it
        raise RuntimeError(f"process_bed_file failed for {file_path} -> {output_path}: {str(e)}")

def process_intersection(a_file, b_file, out_file):
    if not os.path.exists(a_file):
        raise FileNotFoundError(f"Intersection a_file missing: {a_file}")
    if not os.path.exists(b_file):
        raise FileNotFoundError(f"Intersection b_file missing: {b_file}")

    a_intervals = read_intervals(a_file)
    b_intervals = read_intervals(b_file)
    intersecting_intervals = bed_intersect(a_intervals, b_intervals)
    pd.DataFrame(intersecting_intervals).to_csv(out_file, sep='\t', header=None, index=False)
    return out_file

def run_parallel_intersections(writepath, threads=None):
    tasks = []

    input_53_srt = os.path.join(writepath, '53_cut_3.srt.txt')
    output_54_sdtin = os.path.join(writepath, '54_sdtin.txt')
    output_55_bedtools = os.path.join(writepath, '55_bedtools_f.txt')
    if os.path.exists(output_54_sdtin) and os.path.exists(input_53_srt):
        tasks.append((output_54_sdtin, input_53_srt, output_55_bedtools))
    else:
        log(f"[WARN] skipping standard intersection because required files missing: {output_54_sdtin} or {input_53_srt}")

    input_53_del = os.path.join(writepath, '53_cut_3_del.txt')
    output_54_sdtin_del = os.path.join(writepath, '54_sdtin_del.txt')
    output_55_bedtools_del = os.path.join(writepath, '55_bedtools_f_del.txt')
    if os.path.exists(output_54_sdtin_del) and os.path.exists(input_53_del):
        tasks.append((output_54_sdtin_del, input_53_del, output_55_bedtools_del))
    else:
        log(f"[WARN] skipping deletion intersection because required files missing: {output_54_sdtin_del} or {input_53_del}")

    if not tasks:
        log("[INFO] No intersection tasks to run.")
        return []

    threads = threads or cpu_count()
    nprocs = min(max(1, threads), len(tasks))
    log(f"[INFO] Running intersections with {nprocs} worker(s) on {len(tasks)} task(s)")

    outputs = []
    with Pool(processes=nprocs) as pool:
        async_results = []
        for (a, b, out) in tasks:
            async_results.append(pool.apply_async(_worker_safe_intersection, (a, b, out)))

        for ar in async_results:
            try:
                res = ar.get()
                outputs.append(res)
                log(f"[OK] intersection finished: {res}")
            except Exception:
                log(f"[ERROR] intersection failed:\n{traceback.format_exc()}")

    return outputs

def _worker_safe_intersection(a_file, b_file, out_file):
    try:
        return process_intersection(a_file, b_file, out_file)
    except Exception as e:
        raise RuntimeError(f"process_intersection failed for {a_file} & {b_file}: {str(e)}")
