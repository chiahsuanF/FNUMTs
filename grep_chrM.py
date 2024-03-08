def extract_chromosomes(input_fasta, output_fasta, chromosome_names):
    selected_chromosomes = []
    with open(input_fasta, 'r') as infile:
        sequence = ''
        for line in infile:
            if line.startswith('>'):
                for name in chromosome_names:
                    if name in line:
                        selected_chromosomes.append(line.strip())
                        break
                else:
                    if sequence:
                        selected_chromosomes.append(sequence)
                        sequence = ''
                    continue
            else:
                sequence += line.strip()
        if sequence:  
            selected_chromosomes.append(sequence)

    # Write selected chromosome sequences to output file
    with open(output_fasta, 'w') as outfile:
        outfile.write('\n'.join(selected_chromosomes))

# Example Usage
input_fasta = 'input.fasta'
output_fasta = 'chromosomesM.fasta'
chromosome_names = ['chrM']  # The chromosome M name needs to be the same as the fasta file
extract_chromosomes(input_fasta, output_fasta, chromosome_names)
