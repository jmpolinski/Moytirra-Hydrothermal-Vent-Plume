#!/bin/bash

# set the path for wehre Virsorter2 resides
export PATH=/data/app/conda/bin/:$PATH

# run VirSorter2 classification with 5 threads (-j) and virsort/prokkatrans as output directory (-w) with a minimum contig length of 1000 bp
virsorter run -w virsort/prokkatrans -i ../assemblies/prokka+transD_transcripts.fa -j 5 --min-length 1000

# run CheckV to look for non-viral sequences (requires prodigal)
export PATH=/data/app/conda/bin/:$PATH
export PATH=data/app/Prodigal/:$PATH

checkv end_to_end final-vifral-combined.fa checkv -t 5 -d /data/app/Virsorter2/checkv/checkv-db-v1.0

# Combine CheckV output
cat checkv/proviruses.fna checkv/viruses.fan > checkv/combined.fna

