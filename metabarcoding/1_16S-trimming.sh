#!/bin/bash

#working in directory with raw fastq files

gunzip *gz
ls *R1* > files
sed -i -e 's/_L001_R1_001.fastq//g' files

# exces length trim (too much overlap results in poor merging)
for i in `cat files`; do /data/app/fastx-toolkit/fastx_trimmer -Q33 -f 15 -l 200 -i ${i}_L001_R1_001.fastq -o ${i}_R1-trim.fq; done
for i in `cat files`; do /data/app/fastx-toolkit/fastx_trimmer -Q33 -f 15 -l 200 -i ${i}_L001_R2_001.fastq -o ${i}_R2-trim.fq; done

# quality trim
for i in `cat files`; do java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 5 -phred33 ${i}_R1-trim.fq ${i}_R2-trim.fq ${i}_R1_qc.fq ${i}_R1.orphan ${i}_R2_qc.fq ${i}_R2.orphan TRAILING:15 MINLEN:150; done

# tidy up directory
mkdir unpaired raw
mv *fastq raw/
mv *orphan unpaired/
