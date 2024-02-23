#!/usr/bin/env python3

import pysam

def process_alignment_file(input_path, output_path):
    samfile = pysam.AlignmentFile(input_path, "r")
    iter_records = samfile.fetch()

    with open(output_path, 'w') as output_file:
        process_records(iter_records, output_file)

def process_records(records, output_file):
    for record in records:
        process_record(record, output_file)

def process_record(record, output_file):
    reference_start = int(record.reference_start)
    nm_tag = get_nm_tag(record)

    if nm_tag is not None:
        output_file.write(f"{record.query_name}\t{record.query_sequence}\t{record.mapping_quality}\t{reference_start + 1}\t{nm_tag}\n")

def get_nm_tag(record):
    for tag, value in record.get_tags():
        if tag == "NM":
            return value
    return None

def generate_vcf_row(row):
    chrMEPos = row["M_start"] + row["lenght"]
    vcf_data = [
        row["chr"],
        row["start_pos"],
        ".",
        "N",
        row["seq"],
        ".",
        "PASS",
        f"CHR={row['chr']};LEN={row['lenght']};SUPREAD={row['depth']};MAPQ={row['MQ']};NM={row['NM']};chrMSPos={row['M_start']};chrMEPos={chrMEPos};TYPE=INS;METHOD=Fnumt"
    ]
    return "\t".join(map(str, vcf_data))
