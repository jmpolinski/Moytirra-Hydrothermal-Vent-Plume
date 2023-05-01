#!/bin/bash

# working in /data/prj/ecosystem-diverity/OceanX/YEP1_2021/metatranscriptome/data

#prepare directory & download reference genome
mkdir bmtagger && cd bmtagger && mkdir tmp
wget https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz
gunzip GRCh38_latest_genomic.fna.gz
export PATH=$PATH:/data/app/bmtools/bmtagger/

# create index/databbase files for BMTagger
# make index for bmfilter:
bmtool -d GRCh38_latest_genomic.fna -o human_GRCh38.bitmask -A 0 -w 18
# make index for srprism:
srprism mkindex -i GRCh38_latest_genomic.fna -o human_GRCh38.srprism -M 7168
# make blast database:
makeblastdb -in GRCh38_latest_genomic.fna -dbtype nucl -out human_GRCh38

# run BMTagger in a loop
mv ../samples .
for i in `cat samples`; 
do bmtagger.sh -b human_GRCh38.bitmask -x human_GRCh38.srprism -T tmp -q 1 -1 ../${i}_all-trim_R1.fq -2 ../${i}_all-trim_R2.fq -o ${i}-bmtagger; 
done

# remove sequences tagged as human
# bmtagger output only gives header up to the space character -- remove what's after in the fastq files so can use the list of contamination sequences from bmtagger to subset -- then subset
## read 1
for i in `cat samples`; 
do sed -i -e 's/ 1:N:0*.*//g' ${i}_all-trim_R1.fq && \
/data/app/seqkit grep -v -n -f bmtagger/${i}-bmtagger ${i}_all-trim_R1.fq -o ${i}_QC_R1.fq; 
done
## read 2
for i in `cat samples`; 
do sed -i -e 's/ 2:N:0*.*//g' ${i}_all-trim_R2.fq && \
/data/app/seqkit grep -v -n -f bmtagger/${i}-bmtagger ${i}_all-trim_R2.fq -o ${i}_QC_R2.fq; 
done

