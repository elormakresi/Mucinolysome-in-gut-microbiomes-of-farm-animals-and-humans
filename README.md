## Tools and databases

- [Pfam v 37.0](https://www.ebi.ac.uk/interpro/entry/pfam/#table) — Protein domain family database. 
- [HMMER v3.4](https://github.com/EddyRivasLab/hmmer) — Biological sequence analysis using profile HMMs.
- [MMseqs2 v16.747c6](https://github.com/soedinglab/MMseqs2) — Ultra fast and sensitive sequence search and clustering suite.
- [DIAMOND v2.1.10.164](https://github.com/bbuchfink/diamond) — Sequence aligner for protein and translated DNA searches.
- [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) — Assigning objective taxonomic classifications to bacterial and archaeal genomes based on the Genome Database Taxonomy.
- [CheckM2 v1.0.2](https://github.com/chklovski/CheckM2) — Rapid assessment of genome bin quality using machine learning.
- [Panaroo](https://github.com/gtonkinhill/panaroo) - Pipeline for pangenome investigation. 
- [FastTree v2.2.0](https://github.com/morgannprice/fasttree) — Inference of approximately-maximum-likelihood phylogenetic trees from alignments of nucleotide or protein sequences.
- [iTOL v7](https://itol.embl.de/) — Display, annotation and management of phylogenetic and other trees.
- [FastANI](https://github.com/ParBLiSS/FastANI) — Fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI).
- [Bowtie2 v2.3](https://github.com/BenLangmead/bowtie2) — Fast and sensitive gapped read aligner.
- [Samtools v1.19](https://github.com/samtools/samtools) — Parse and manipulate alignments in the SAM/BAM format.
- [inStrain v1.3](https://github.com/MrOlm/inStrain) — Analysis of co-occurring genome populations from metagenomes.
- [CoverM v0.6](https://github.com/wwood/CoverM) — Fast DNA read coverage and relative abundance calculator focused on metagenomics applications.
- [BBTools](https://github.com/bbushnell/BBTools) — Fast, multithreaded bioinformatics tools for DNA/RNA analysis.
- [Kraken2 v2.1](https://github.com/DerrickWood/kraken2) — Metagenomic taxonomic sequence classification system.
- [Seqtk v1.2](https://github.com/lh3/seqtk) — Fast and lightweight tool for processing sequences in the FASTA or FASTQ format.
- [SeqKit](https://github.com/shenwei356/seqkit) — A cross-platform and ultrafast toolkit for FASTA/Q file manipulation.
- [FastQC v0.12](https://github.com/s-andrews/FastQC) — A Quality Control application for FastQ files.
- [MultiQC v1.8](https://github.com/MultiQC/MultiQC) — Aggregate bioinformatics results across many samples into a single report.
- [Proteinortho](http://legacy.bioinf.uni-leipzig.de/Software/proteinortho/) — Detection orthologous genes within different species.
- [Barrnap v0.9](https://github.com/tseemann/barrnap) — BAsic Rapid Ribosomal RNA Predictor.
- [Alphafold3](https://github.com/google-deepmind/alphafold3) — Accurate structure prediction of biomolecular interactions.
- [Folddock](https://gitlab.com/ElofssonLab/FoldDock) — Prediction of protein-protein interactions.
- [Trimmomatic v0.39](https://github.com/usadellab/Trimmomatic) — A Java-based processing and trimming tool for Illumina NGS sequencing data.
- [Hisat2 v2.2](https://github.com/DaehwanKimLab/hisat2) — Fast and sensitive alignment program for mapping next-generation sequencing reads.
- [Gffread](https://github.com/gpertea/gffread) — GFF/GTF utility providing format conversions, filtering, FASTA sequence extraction and more.
- [Subread v2.0](https://github.com/ShiLab-Bioinformatics/subread) — Tool kit for processing next-gen sequencing data.
- [DESeq2](https://github.com/thelovelab/DESeq2) — Differential gene expression analysis of RNA sequencing data.
- [clusterProfiler v4.12.6](https://github.com/YuLab-SMU/clusterProfiler) — Statistical analysis and visualization of functional profiles for genes and gene clusters.
- []()
## MATERIALS AND METHODS

### Search for mucinolysome-encoding MAGs in UHGG

- [run_hmm_cohesin.sh](/scripts/finding_mags/run_hmm_cohesin.sh)
- [run_hmm_dockerin.sh](/scripts/finding_mags/run_hmm_dockerin.sh)
- [parse_hmmer_coverage.py](/scripts/finding_mags/parse_hmmer_coverage.py)
- [extract_doc_coh_faa.py](/scripts/finding_mags/extract_doc_coh_faa.py)
- [run_dbcan.sh](/scripts/finding_mags/run_dbcan.sh)
### Expand the search for mucinolysome-encoding MAGs in other sources

- [mmseqs.sh](/scripts/finding_mags/mmseqs.sh)
- [get_peptidase.sh](/scripts/finding_mags/get_peptidase.sh)
- [parse_sulfatlas_hmmscan_results.sh](/scripts/finding_mags/parse_sulfatlas_hmmscan_results.sh)
### MAG taxonomic classification, quality check and phylogenetic analysis

- [gtdb.sh](/scripts/finding_mags/gtdb.sh)
- [mag_quality.sh](/scripts/finding_mags/mag_quality.sh)
- [core_gene_aligment.sh](/scripts/finding_mags/core_gene_aligment.sh)
### MAG relative abundance and prevalence by read mapping of microbiome samples

- [get_reads.sh](/scripts/read_mapping_original_and_sample_environments/get_reads.sh)
- [bbduk.sh](/scripts/read_mapping_original_and_sample_environments/bbduk.sh)
- [kraken_seqtk.sh](/scripts/read_mapping_original_and_sample_environments/kraken_seqtk.sh)
- [subsampling.sh](/scripts/read_mapping_original_and_sample_environments/subsampling.sh)
- [read_mapping_original.sh](/scripts/read_mapping_original_and_sample_environments/read_mapping_original.sh) — read mapping against original samples
- [read_mapping_whole_genome.sh](/scripts/read_mapping_original_and_sample_environments/read_mapping_whole_genome.sh) — whole genome read mapping
- [read_mapping_cazyme_scaf_cds.sh](/scripts/read_mapping_original_and_sample_environments/read_mapping_cazyme_scaf_cds.sh) — read mapping against CDS
- [read_mapping_16s.sh](/scripts/read_mapping_original_and_sample_environments/read_mapping_16s.sh) — read mapping against 16s
- [barrnap.sh](/scripts/read_mapping_original_and_sample_environments/barrnap.sh)
### Protein structure-based protein-protein interactions

- [slurm_generate_AF3.py](/scripts/structure_prediction/slurm_generate_AF3.py)
- [ppi_metrics_from_pdb.py](/scripts/structure_prediction/ppi_metrics_from_pdb.py)
- [get_other_metrics_from_json.py](/scripts/structure_prediction/get_other_metrics_from_json.py)
- [make_final_data_file_4_heatmaps.py](/scripts/structure_prediction/make_final_data_file_4_heatmaps.py)
- [cper_combined_images_all_6.py](/scripts/structure_prediction/cper_combined_images_all_6.py)
### RNA-seq differential gene expression analysis

- [fastqc_01.sh](/scripts/RNA_seq/fastqc_01.sh)
- [trimmomatic_02.sh](/scripts/RNA_seq/trimmomatic_02.sh)
- [fastqc_01-post_trim.sh](/scripts/RNA_seq/fastqc_01-post_trim.sh)
- [hisat_samtools.sh](/scripts/RNA_seq/hisat_samtools.sh)
- [read_counting.sh](/scripts/RNA_seq/read_counting.sh)
- [deseq2_w_heatmap_w_tpm_rpkm.R](/scripts/RNA_seq/deseq2_w_heatmap_w_tpm_rpkm.R)
