#######################################
# to parse hmmer results using coverage cutoffs
# input is hmmer tsv results
#######################################

import os

def parse_files_for_coverage(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".tsv"):  # Modify this if needed to match your file types
            file_path = os.path.join(directory, filename)
            output_filename = filename.replace(".tsv", "_Coverage_parsed.tsv")
            output_path = os.path.join(directory, output_filename)

            with open(file_path, 'r') as infile, open(output_path, 'w') as outfile:
                for line in infile:
                    if line.startswith('#') or line.strip() == '':
                        continue
                    parts = line.split()
                    try:
                        qlen = int(parts[5])  # Query length
                        start = int(parts[15])  # Start position
                        end = int(parts[16])  # End position
                        coverage = (end - start) / qlen
                        if coverage >= 0.6:
                            outfile.write(line)
                    except IndexError:
                        print(f"Skipping malformed line in {filename}: {line}")
                    except ValueError:
                        print(f"Data conversion error in {filename}: {line}")

directory = "./faa_files_putative_mucins/mucosal_ref_genomes/coh_doc_results"
parse_files_for_coverage(directory)
