#!/usr/bin/env python3

import pysam, datetime, sys, os, argparse

def extract_chromosomes(file_data):
    # Assuming chromosomes are in the #CHROM column
    chromosomes = set()
    for line in file_data.splitlines():
        if not line.startswith("#"):
            fields = line.split('\t')
            chromosomes.add(fields[0])  # Assuming the chromosome is in the first column
    return list(chromosomes)

def vcf_header(chromosomes, bampath):
    filedate = str(datetime.datetime.now()).split(' ')[0]
    vcfheader = pysam.VariantHeader()
    
    sample_name = os.path.basename(bampath)

    bamfile = pysam.AlignmentFile(bampath, 'r')
    vcfheader.add_sample(sample_name)
    
    vcfheader.add_meta('source', 'Fnumts')
    vcfheader.add_meta('fileDate', filedate)

    for chrom in chromosomes:
        vcfheader.add_line(f'##contig=<ID={chrom}>')

    vcfheader.add_meta('INFO', items=[('ID', 'CHR'), ('Number', 1), ('Type', 'String'), ('Description', 'Chromosome where the NUMT sequence is located')])
    vcfheader.add_meta('INFO', items=[('ID', 'LEN'), ('Number', 1), ('Type', 'Integer'), ('Description', 'Length of the NUMT sequence in base pairs')])
    vcfheader.add_meta('INFO', items=[('ID', 'SUPREAD'), ('Number', 1), ('Type', 'Integer'), ('Description', 'Number of supporting reads for the NUMT sequence')])
    vcfheader.add_meta('INFO', items=[('ID', 'NM'), ('Number', 1), ('Type', 'Integer'), ('Description', 'Number of mismatches in the NUMT sequence')])
    vcfheader.add_meta('INFO', items=[('ID', 'MQ'), ('Number', 1), ('Type', 'Integer'), ('Description', 'Average mapping quality of the original support reads')])
    vcfheader.add_meta('INFO', items=[('ID', 'chrMSPos'), ('Number', 1), ('Type', 'Integer'), ('Description', 'Start coordinate of the NUMT sequence in the mitochondrial genome')])
    vcfheader.add_meta('INFO', items=[('ID', 'chrMEPos'), ('Number', 1), ('Type', 'Integer'), ('Description', 'End coordinate of the NUMT sequence in the mitochondrial genome')])
    vcfheader.add_meta('INFO', items=[('ID', 'TYPE'), ('Number', 1), ('Type', 'String'), ('Description', 'Structural variant type: INS indicates an insertion relative to the reference; DEL indicates a deletion relative to the reference')])
    vcfheader.add_meta('INFO', items=[('ID', 'METHOD'), ('Number', 1), ('Type', 'String'), ('Description', 'Detection method used for identifying the NUMT sequence')])
    vcfheader.add_meta('INFO', items=[('ID', 'VAF'), ('Number', 1), ('Type', 'Float'), ('Description', 'Variant allele frequency at the NUMT breakpoint')])
    vcfheader.add_meta('INFO', items=[('ID', 'DEPTH'), ('Number', 1), ('Type', 'Integer'), ('Description', 'Total depth of coverage at the NUMT breakpoints')])
    vcfheader.add_meta('FORMAT', items=[('ID',"GT"), ('Number',1), ('Type','String'),('Description','Genotype')])
    vcfheader.add_meta('FORMAT', items=[('ID',"DP"), ('Number',1), ('Type','Integer'),('Description','Total depth of coverage at the NUMT breakpoints')])
    vcfheader.add_meta('FORMAT', items=[('ID',"AD"), ('Number',1), ('Type','Integer'),('Description','Number of supporting reads for the NUMT sequence')])
    vcfheader.add_meta('FILTER', items=[('ID',"LowQuality"), ('Number',1), ('Type','String'),('Description','Low quality variant')])    
    # Add command line info
    vcfheader.add_meta('CommandLine', ' '.join(sys.argv))

    return vcfheader

def add_format_and_sample_data(file_data):
    new_lines = []
    
    for line in file_data.splitlines():
        # Skip the header line
        if line.startswith("#CHROM"):
            new_lines.append(line + '\tFORMAT\tSAMPLE')
            continue
        
        fields = line.split('\t')
        info_field = fields[7]  # INFO field is the 8th column
        info_dict = {kv.split('=')[0]: kv.split('=')[1] for kv in info_field.split(';')}
        vaf = float(info_dict.get('VAF', 0))  # Extract VAF
        depth = int(float(info_dict.get('DEPTH', '0')))
        supread = info_dict.get('SUPREAD', '0')  # Extract SUPREAD

        if vaf < 0.3:
            gt = '0/0'
        elif 0.3 <= vaf <= 0.69:
            gt = '0/1'
        else:
            gt = '1/1'
        
        format_value = 'GT:DP:AD'
        sample_value = f'{gt}:{depth}:{supread}'
        
        new_line = line + f'\t{format_value}\t{sample_value}'
        new_lines.append(new_line)
    
    return '\n'.join(new_lines)
