#!/bin/bash

# working in /data/prj/ecosystem-diverity/OceanX/YEP1_2021/metatranscriptome/data/

for i in `cat samples`; do \
java -jar /data/app/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 8 -phred33 ./merged/${i}_all_R1.fastq ./merged${i}_all_R2.fastq \
${i}_R1-trim.fq ${i}_R1-unpaired.fq ${i}_R2-trim.fq ${i}_R2-unpaired.fq \
ILLUMINACLIP:/data/app/Trimmomatic-0.38/adapters/NexteraPE-PE.fa:2:10:30:8 SLIDINGWINDOW:5:20 HEADCROP:12 MINLEN:50; \
done

# move unpaired fastq files
mkdir unpaired-trimmed && mv *unpaired.fq unpaired-trimmed
mkdir trimmed && mv *trim.fq trimmed
