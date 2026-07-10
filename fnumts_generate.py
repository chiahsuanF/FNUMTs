#!/usr/bin/env python3



def save_sorted_data(sorted_df, output_file):
    sorted_df.to_csv(output_file, sep='\t', header=None, index=False)

def add_header_column(data):
    data['header'] = data.apply(generate_header, axis=1)
    return data[['header', 'seq']]

def generate_header(row):
    return f'>{row["chr"]}_{row["start_pos"]}_{row["depths"]}_{row["INS"]}_{row["MQ"]}_{row["lenght"]}'

def save_to_fna(data, output_file):
    data.to_csv(output_file, sep='\n', header=False, index=False)

def calculate_score(mismatches, vaf, mq, seq_type, max_mismatches=50, max_mq=60, weights=(0.3, 0.4, 0.3)):
    def normalize_mismatches(mismatches, max_mismatches):
        return 1 - min(mismatches / max_mismatches, 1)

    def normalize_vaf(vaf):
        return min(vaf, 1)
    
    def normalize_mq(mq, max_mq):
        return min(mq / max_mq, 1)

    normalized_mismatches = normalize_mismatches(mismatches, max_mismatches)
    normalized_vaf = normalize_vaf(vaf)  # Ensure VAF is capped at 1
    normalized_mq = normalize_mq(mq, max_mq)
    
    score = (
        (weights[0] * normalized_mismatches +
        weights[1] * normalized_vaf +
         weights[2] * normalized_mq) * 100
    )
    return score
    
def process_vcf(input_file, output_file, seq_type):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
                continue
            
            fields = line.strip().split('\t')
            chr_value = fields[0]
            pos = int(fields[1])
            id_value = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            filter_status = fields[6]
            info = fields[7]
            
            info_dict = dict(item.split('=') for item in info.split(';'))
            support_read = float(info_dict.get('SUPREAD', 0))
            mismatches = int(info_dict.get('NM', 0))
            vaf = float(info_dict.get('VAF', 0))
            mq = int(info_dict.get('MQ', 0))

            score = calculate_score(mismatches, vaf, mq, seq_type)

            if score < 40.00:
                filter_status = 'LowQuality'
            else:
                filter_status = 'PASS'
            
            qual = f"{score:.2f}"

            info_dict['VAF'] = f'{vaf:.2f}'
            updated_info = ';'.join([f'{key}={value}' for key, value in info_dict.items()])

            updated_line = '\t'.join([
                chr_value,
                str(pos),
                id_value,
                ref,
                alt,
                qual,
                filter_status,
                updated_info
            ])
            outfile.write(updated_line + '\n')
