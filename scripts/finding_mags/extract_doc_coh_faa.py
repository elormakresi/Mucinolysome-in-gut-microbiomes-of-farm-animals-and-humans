#######################################
# This script extracts the coh/doc sequences that after mmseq2 and hmmer
# input is protein sequences of putative genomes and tsv results from hmmer
#######################################
import os

tsv_dir = './genomes-human-gutMAGSnewfaa/putative_minus_5496/coh_doc_results/selected_cols/coh_doc_combined'
faa_dir = './genomes-human-gutMAGSnewfaa/putative_minus_5496'

def read_protein_ids(tsv_file):
    protein_ids = []
    with open(tsv_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            protein_ids.append(parts[0])
    return protein_ids

def index_faa_file(faa_file):
    index = {}
    pos = 0
    with open(faa_file, 'r') as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                protein_id = line.split()[0][1:].split(' ')[0]
                index[protein_id] = pos
            pos = file.tell()
            line = file.readline()
    return index

def extract_sequences(faa_file, protein_ids):
    index = index_faa_file(faa_file)
    sequences = {}
    with open(faa_file, 'r') as file:
        for protein_id in protein_ids:
            if protein_id in index:
                file.seek(index[protein_id])
                header = file.readline().strip()
                seq = []
                while True:
                    pos = file.tell()
                    line = file.readline()
                    if line.startswith('>'):
                        file.seek(pos)
                        break
                    if line.strip() and not line.startswith('>'):
                        seq.append(line.strip())
                seq = ''.join(seq)
                if seq.endswith('*'):
                    seq = seq[:-1]  # Remove the trailing asterisk. After i used prodigal, some of the faa seqs had asterisks
                sequences[protein_id] = seq
    return sequences

# Process each TSV file
for tsv_filename in os.listdir(tsv_dir):
    if tsv_filename.endswith('.tsv'):
        genome_id = tsv_filename.split('_')[0]
        tsv_file_path = os.path.join(tsv_dir, tsv_filename)
        faa_file_path = os.path.join(faa_dir, f'{genome_id}.faa')
        
        if os.path.exists(faa_file_path):
            protein_ids = set(read_protein_ids(tsv_file_path))  # Use a set for faster lookups
            sequences = extract_sequences(faa_file_path, protein_ids)
            
            output_file = os.path.join(tsv_dir, f'{genome_id}_extracted_coh_doc_proteins.faa')
            with open(output_file, 'w') as out_file:
                for pid, seq in sequences.items():
                    out_file.write(f'>{pid}\n{seq}\n')
            print(f'Extracted proteins for {genome_id} written to {output_file}')
        else:
            print(f'No FAA file found for {genome_id}')

