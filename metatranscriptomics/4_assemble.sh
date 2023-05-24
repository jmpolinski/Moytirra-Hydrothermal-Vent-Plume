#!/bin/bash

# assembling metatranscriptomes with trinity and IDBA

## TRINITY
### working in /data/prj/ecosystem-diverity/OceanX/YEP1_2021/metatranscriptome/assemblies
### files "samples" is list of samples IDs that are in file names

export PATH=$PATH:/data/app/trinityrnaseq-v2.13.2/
export TRINITY_HOME=/data/app/trinityrnaseq-v2.13.2/
for i in `cat samples`; do \
Trinity --seqType fq --samples_file ${i}.files --SS_lib_type RF --CPU 24 --max_memory 75G --output trinity_${i}; done

### map reads back to check assembly quality
for i in `cat assembled`; \
do /data/app/bowtie2-2.4.5/bowtie2-build ../trinity_${i}.Trinity.fasta ${i}-trinity &&  \
bowtie2 -p 10 -q --no-unal -k 20 -x ${i}-trinity -1 ../${i}F1_QC_R1.fq -2 ../${i}F1_QC_R2.fq 2>${i}_align_stats.txt| samtools view -@10 -Sb -o ${i}_bowtie2.bam; done



## IDBA 
### file "idba.assembler.sh" included as separate script 
for i in `cat samples`; do ./idba-assembler.sh ${i}; done

