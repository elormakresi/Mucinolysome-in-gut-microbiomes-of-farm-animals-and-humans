#!/bin/bash

# run gtdbtk to assign taxonomic classification
gtdbtk classify_wf --genome_dir ./putative_mucins_fna/all_fna_seqs --out_dir ./gtdb_classification --skip_ani_screen --cpus 32 --pplacer_cpus 32
