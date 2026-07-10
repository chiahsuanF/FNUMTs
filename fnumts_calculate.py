#!/usr/bin/env python3

import pysam
from fnumts_check import parse_args
from collections import defaultdict

#from io import StringIO


def sam_to_paf(sam_file, paf_file, seq_type):
    passing_reads = []
    with pysam.AlignmentFile(sam_file, "r") as samfile, open(paf_file, 'w') as paffile:
        for read in samfile:
            if read.is_unmapped:
                continue

            query_name = read.query_name
            query_length = read.query_length
            query_start = int(read.query_alignment_start) + 1
            query_end = int(read.query_alignment_end) + 1
            strand = '-' if read.is_reverse else '+'
            target_name = read.reference_name
            target_length = int(samfile.get_reference_length(target_name))//2
            target_start = read.reference_start
            target_end = read.reference_end
            num_mismatches = read.get_tag("NM")  # Number of mismatches
            alignment_block_length = read.query_alignment_length

            paf_line = [
                query_name,
                str(query_length),
                str(query_start),
                str(query_end),
                strand,
                target_name,
                str(target_length),
                str(target_start),
                str(target_end),
                str(alignment_block_length),
                str(num_mismatches)  # Mismatches
            ]

            paffile.write('\t'.join(paf_line) + '\n')

            # Check if the read passes the criteria for further processing
            if target_name == 'chrM':
                query_coverage = (query_end - query_start) / query_length * 100
                if (seq_type == 'hifi' and query_coverage > 45.00) or (seq_type == 'ont' and query_coverage > 39.00):
                    matches = (query_end - query_start - num_mismatches)
                    match_ratio = matches / (query_end - query_start) * 100
                    if (seq_type == 'hifi' and match_ratio > 79.99) or (seq_type == 'ont' and match_ratio > 69.99):
                        passing_reads.append(read)
    return passing_reads

def generate_vcf_row(row, id_value):
    seq_length = len(row["seq"])
    chrMEPos = row["M_start"] + seq_length

    # Adjust mitochondrial coordinates
    adj_M_start = row["M_start"]
    if adj_M_start > 16569:
        adj_M_start -= 16568

    adj_chrMEPos = chrMEPos
    if chrMEPos > 16569:
        adj_chrMEPos -= 16568

    # Normalize and classify variant type
    type_str = str(row["type"]).strip().upper()
    if type_str in ["INS", "SC"]:
        vcf_type = "INS"
        ref_allele = "N"
        alt_allele = row["seq"]
    elif type_str == "DEL":
        vcf_type = "DEL"
        ref_allele = row["seq"]
        alt_allele = "N"
    else:
        vcf_type = "XXXX"
        ref_allele = "N"
        alt_allele = row["seq"]

    # Build VCF row
    vcf_data = [
        row["chr"],                 # CHROM
        row["start_pos"],           # POS
        f"FN{id_value}",            # ID
        ref_allele,                 # REF
        alt_allele,                 # ALT
        ".",                        # QUAL
        "PASS",                     # FILTER
        f"CHR={row['chr']};LEN={seq_length};SUPREAD={row['depth']};"
        f"NM={row['NM']};MQ={row['aveMQ']};"
        f"chrMSPos={adj_M_start};chrMEPos={adj_chrMEPos};"
        f"TYPE={vcf_type};METHOD=Fnumt"  # INFO
    ]

    return "\t".join(map(str, vcf_data))


def process_read(read, snp1):
    Ms = read.reference_start if read.reference_start else 0
    nm_tag = next((value for tag, value in read.get_tags() if tag == "NM"), "NA")
    write_chrM_read(snp1, read, Ms, nm_tag)

def write_chrM_read(snp1, read, pos, nm_tag):
    snp1.write(f"{read.query_name}\t{read.query_sequence}\t{read.mapping_quality}\t{pos + 1}\t{nm_tag}\n")


##########

def get_depth(bam_filename, chrom, pos,
              min_base_quality=0,
              min_mapping_quality=0,
              reference=None):
    try:
        # Open BAM or CRAM correctly
        if bam_filename.endswith(".cram"):
            if reference is None:
                raise ValueError("CRAM requires reference FASTA")
            bamfile = pysam.AlignmentFile(
                bam_filename, "rc", reference_filename=reference
            )
        else:
            bamfile = pysam.AlignmentFile(bam_filename, "rb")

        # Fix chr / non-chr mismatch
        if chrom not in bamfile.references:
            chrom_alt = chrom.replace("chr", "") if chrom.startswith("chr") else "chr" + chrom
            if chrom_alt in bamfile.references:
                chrom = chrom_alt
            else:
                bamfile.close()
                return 0

        depth = 0
        for col in bamfile.pileup(
            chrom,
            pos - 1,
            pos,
            truncate=True,
            stepper="all",
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
            multiple_iterators=False
        ):
            if col.pos == pos - 1:
                depth = col.nsegments

        bamfile.close()
        return depth

    except Exception:
        # pysam pileup iteration error → treat as zero depth
        return 0

def process_vcf_and_filter(input_vcf, bam_filename, output_vcf, seq_type, reference=None):
    if seq_type not in ("hifi", "ont"):
        raise ValueError("Unsupported sequencing type. Use 'hifi' or 'ont'.")

    threshold = 0.15

    with open(input_vcf) as vcf_in, open(output_vcf, "w") as vcf_out:
        for line in vcf_in:
            if line.startswith("#"):
                vcf_out.write(line)
                continue

            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = int(fields[1])

            info_dict = parse_info(fields[7])
            support_read = float(info_dict.get("SUPREAD", 0))
            vcf_type = info_dict.get("TYPE", "INS")

            depth = get_depth(
                bam_filename,
                chrom,
                pos,
                reference=reference
            )

            total_depth = depth + support_read if vcf_type == "DEL" else depth
            if total_depth <= 0:
                continue

            vaf = support_read / total_depth
            if vaf <= threshold:
                continue

            info_dict["VAF"] = f"{min(vaf, 1.0):.2f}"
            info_dict["DEPTH"] = str(int(total_depth))
            fields[7] = format_info(info_dict)

            vcf_out.write("\t".join(fields) + "\n")
            
def parse_info(info):
    d = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = True
    return d


def format_info(d):
    return ";".join(
        k if v is True else f"{k}={v}"
        for k, v in d.items()
    )



def filter_best_primary_alignments(input_sam: str, output_sam: str) -> None:
    samfile = pysam.AlignmentFile(input_sam, "r")
    out = pysam.AlignmentFile(output_sam, "w", template=samfile)

    best_alignment = {}

    for read in samfile.fetch(until_eof=True):
        # Skip secondary or supplementary alignments
        if read.is_secondary or read.is_supplementary:
            continue

        name = read.query_name
        score = read.get_tag("AS") if read.has_tag("AS") else 0

        # Keep read if no previous or current has better score
        if name not in best_alignment or score > best_alignment[name].get_tag("AS"):
            best_alignment[name] = read

    # Write best alignments to output
    for read in best_alignment.values():
        out.write(read)

    out.close()
    samfile.close()
