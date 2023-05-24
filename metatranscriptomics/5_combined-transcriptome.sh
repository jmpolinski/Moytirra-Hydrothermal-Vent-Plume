#!/bin/bash

#Combine each sample's Trinity and IDBA-UD assemblies and use CD-HIT to get a consensus transcript set:
##pwd = /data/prj/ecosystem-diverity/OceanX/YEP1_2021/metatranscriptome/assemblies/

# first get Trinity longest isoform:
for i in `cat all-samples`; do /data/app/trinityrnaseq-v2.13.2/util/misc/get_longest_isoform_seq_per_trinity_gene.pl trinity_${i}.Trinity.fasta > trinity_${i}.longest.fasta; done

# combine Trinity and IDBA assemblies:
mkdir cdHit_transcript-sets && cd cdHit_transcript-sets
for i in `cat all-samples`; do cat ../trinity_${i}.longest.fasta ../${i}_idba/${i}_idba-ud/contig.fa > ./${i}_trinity+idba.fa; done

# run cd-hit on each sample's combined set:
for i in `cat all-samples`; do /data/app/cd-hit-otu-illumina-0.0.1/cd-hit/cd-hit -i ${i}_trinity+idba.fa -o ${i}_cdhit-90.fa -c 0.9 -T 8 -b 50; done

# combine and run ch-hit for final consensus set
cat *_cdhit-90.fa > combined_cdhit.fa
/data/app/cd-hit-otu-illumina-0.0.1/cd-hit/cd-hit -i ${i}_trinity+idba.fa -o /data/app/cd-hit-otu-illumina-0.0.1/cd-hit/cd-hit -i combined_cdhit.fa -o combined_cdhit90-FINAL.fa -c 0.9 -T 8 -b 50





# PROKKA to identify bacterial genes
/data/app/prokka/bin/prokka --outdir cdhit90_prokka ./combined_cdhit90-FINAL.fa


# Subset contigs that do not have a PROKKA annotation/hit:
## get headers for contigs with a prokka hit(s)
cd cdhit90_prokka
grep -B 1 -v ">Feature" PROKKA_04062022.tbl > feature-w-hits
grep ">" feature-w-hits > tmp
rm feature-w-hits && mv tmp feature-w-hits
sed -i -e 's/>Feature //g' feature-w-hits
## use prokka hit ID list to subset sequences not in the list
cd ..
/data/app/seqkit grep -v -f cdhit90_prokka/feature-w-hits combined_cdhit90-FINAL.fa > contigs_no-prokka.fa


# TransDecoder to identify eukaryotic genes/proteins:
## extract open reading frames
/data/app/TransDecoder-v5.5.0/TransDecoder.LongOrfs -t contigs_no-prokka.fa --output_dir cdhit90_transDecoder
## get HMM and blastP results for the longest ORF predictions:
cd cdhit90_transDecoder
/data/app/Diamond/diamond blastp -d /data/app/databases/refseq_2022-01/refseq-diamond.dmnd -q ./longest_orfs.pep -o refseq.outfmt6 -p 8 -f 6 -e 1e-5
/data/app/hmmer-3.3/bin/hmmscan --domtblout hmm-results.txt /data/app/databases/HMM-databases/Pfam-A.hmm longest_orfs.pep
## get final protein set:
cd ..
/data/app/TransDecoder-v5.5.0/TransDecoder.Predict -t contigs_no-prokka.fa --retain_pfam_hits ./cdhit90_transDecoder/hmm-results.txt --retain_blastp_hits ./cdhit90_transDecoder/refseq.outfmt6 --output_dir cdhit90_transDecoder

# Combine predicted transcripts from PROKKA and TransDecoder for downstream analysis:
##pwd = /data/prj/ecosystem-diverity/OceanX/YEP1_2021/metatranscriptome/assemblies
cat ./cdhit90_prokka/PROKKA_04062022.ffn ./cdhit90_transDecoder/contigs_no-prokka.fa.transdecoder.cds > prokka+transD_transcripts.fa
