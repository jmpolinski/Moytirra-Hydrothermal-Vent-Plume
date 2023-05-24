# Functional Characterization via Metatranscriptomics

## Library Preparation

RNA quantity was measured with the Qubit HS RNA assay and quality was checked on an Agilent Fragment Analyzer with the HS RNA kit. Illumina TruSeq Stranded Total RNA Prep kit with Ribo-Zero Plus was used to prepare libraries following manufacturer's SOP.  Samples were pooled and sequenced on an Illumina NovaSeq (multiple runs).  

## Bioinformatic Analysis

1. Separate run files for each sample were combined

2. Trimmomatic (v0.38) was used for quality-based read trimming and to remove adapter sequence 

3. BMTagger and human reference genome GRCh38 were used to remove human contamination  

4. Assemble transcriptomes for each sample with Trinity and IDBA  

5. Create a combined transcriptome with CD-Hit and identify coding sequences with Prokka & TransDecoder  

6. Identify viral sequences in transcript set from step 5

7. Generate normalized count table with Trinity scripts (align_and_estimate_abundance.pl & abundance_estimates_to_matrix.pl)

8. Differential gene expression analysis
