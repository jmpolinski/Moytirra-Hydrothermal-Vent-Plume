# Prokaryotic (16S) & Eukaryotic (18S) Community Characterization

Metabarcoding using prokaryote- and eukaryote-specific primers was conducted in order to characterize and identify differences in community composition across the transect.  

## Primers

### Prokaryotic Communities

A ~291bp section of the V4 region of the 16S SSU rRNA gene was targeted using the following primers:  
- 16S-V4_515F: FL + GTGYCAGCMGCCGCGGTAA [(Parada et al. 2016)](https://ami-journals.onlinelibrary.wiley.com/doi/epdf/10.1111/1462-2920.13023)
- 16S-V4_806F: RL + GGACTACNVGGGTWTCTAAT [(Apprill et al. 2015)](https://www.int-res.com/abstracts/ame/v75/n2/p129-137/)

### Eukaryotic Communities

Primers from [Blaxter et al. 1998](https://www.nature.com/articles/32160) were used to target overall eukaryotic diversity.
- 18S-SSU_FO4: FL + GCTTGTCTCAAAGATTAAGCC 
- 18S-SSU_R22: RL + GCCTGCTGCTGCCTTCCTTGGA  
  
Both 16S and 18S primers included linking seqeunces (FL/RL) for addition of Illumina-compatible dual-index sequencing adapters. Sequencing was conducted on a MiSeq using v3 reagents (2x300).  

## Bioinformatic analysis

1. Trimming & QC of raw fastq with fastx-toolkit and Trimmomatic (v0.38)  

2. ASV identificiation with USEARCH (v11): 
- 16S: ```/data/app/usearch_pipeline.sh OXR21-16S```  
- 18S: ```/data/app/usearch_pipeline.sh OXR21-18S```

3. Classify ASV sequences with Silva v138.1 (mothur classify.seqs)  

4. Visualization & analysis in R
