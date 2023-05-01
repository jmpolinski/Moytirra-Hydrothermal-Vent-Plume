#!/bin/bash

# working in /data/prj/ecosystem-diverity/OceanX/YEP1_2021/metatranscriptome/data/

# create sample list
ls ././metaT_test-run/C*R1*fastq > samples
sed -i -e 's/_S*.*//g' samples

for i in `cat samples`; do \
cat ./metaT_test-run/${i}_*R1*fastq ./run1/${i}_*R1*fastq ./run2/${i}_*R1*fastq > ./merged/${i}_all__R1.fastq && \
cat ./metaT_test-run/${i}_*R2*fastq ./run1/${i}_*R2*fastq ./run2/${i}_*R2*fastq > ./merged/${i}_all_R2.fastq; 
done
