#!/bin/bash

# run prokka to obtain gff3 since panaroo prefers its gff version (when i tried the prodigal gff it failed if i remember correctly.
prokka --outdir "$temp_dir" --prefix "$(basename "${fna_file%.*}")" --force --gffver 3 --genus Unknown --species Unknown --strain "$(basename "${fna_file%.*}")" --usegenus --quiet "$fna_file"

# run panaroo for gene alignment
panaroo -i ./gff_files_prokka/*.gff -o panaroo_prokka_gff_putative_mucins --clean-mode strict -a core -t 32 --aligner mafft --core_threshold 1

# run for tree inference and as input into itol
FastTree -nt ./panaroo_prokka_gff_putative_mucins/core_gene_alignment_filtered.aln > putative_mucin_core_gene_alignment_fasttree.nwk
