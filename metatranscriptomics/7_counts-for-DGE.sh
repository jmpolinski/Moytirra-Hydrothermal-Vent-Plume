#!/bin/bash 

## Use align & estimate script that comes with Trinity to generate count data for DGE
export PATH=$PATH:/data/app/RSEM-1.3.3/:/data/app/bowtie2/
/data/app/trinityrnaseq-v2.13.2/util/align_and_estimate_abundance.pl --transcripts prokka+transD_transcripts.fa --est_method RSEM --aln_method bowtie2 --samples_file quant-samples.txt --seqType fq --thread_count 12 --prep_reference --coordsort_bam

# get a text file listing all results files
ls ./C*/*isoforms.results > results.files

# use trinity script to convert to count matrix
/data/app/trinityrnaseq-v2.13.2/util/abundance_estimates_to_matrix.pl --est_method RSEM --quant_files results.files --gene_trans_map none --name_sample_by_basedir


