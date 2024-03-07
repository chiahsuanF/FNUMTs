#!/usr/bin/env python3


#72
def save_sorted_data(sorted_df, output_file):
    sorted_df.to_csv(output_file, sep='\t', header=None, index=False)

#77
def add_header_column(data):
    data['header'] = data.apply(generate_header, axis=1)
    return data[['header', 'seq']]

def generate_header(row):
    return f'>{row["chr"]}_{row["start_pos"]}_{row["depths"]}_INS_{row["lenght"]}'

#78
def save_to_fna(data, output_file):
    data.to_csv(output_file, sep='\n', header=False, index=False)
