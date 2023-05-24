###-------- load packages you need --------###
library(DESeq2)
library(DEVis)

###-------- load data process data --------###
dta <- read.csv("input/RSEM.counts.matrix.csv", 
                header=TRUE, row.names = 1)

# remove rows with all zeros
dta <- dta[rowSums(dta[])>0,]

###-------- define data structure --------###
samples <- data.frame(row.names=c("C1B1F1","C1B1F2","C1B1F3","C1B2F1","C1B2F2",
                                  "C1B2F3","C1B3F1","C1B3F2","C1B3F3","C1B4F1",
                                  "C1B4F2","C1B4F3","C1B5F1","C1B5F2","C1B5F3",
                                  "C1B6F1","C1B6F2","C1B6F3","C1B7F1","C1B7F2",
                                  "C1B7F3","C1B8F1","C1B8F2","C1B8F3","C1B9F1",
                                  "C1B9F2","C1B9F3","C1B10F1","C1B10F2","C1B10F3",
                                  "C1B11F1","C1B11F2","C1B11F3","C2B2F1","C2B2F2",
                                  "C2B2F3","C2B3F1","C2B3F2","C2B3F3","C2B4F1",
                                  "C2B4F2","C2B4F3","C2B5F1","C2B5F2","C2B5F3",
                                  "C1B12F1","C1B12F2","C1B12F3","C2B1F1","C2B1F2",
                                  "C2B1F3"), 
                      condition=as.factor(c(rep("C1B1",3),rep("C1B2",3),rep("C1B3",3),
                                            rep("C1B4",3),rep("C1B5",3),rep("C1B6",3),
                                            rep("C1B7",3),rep("C1B8",3),rep("C1B9",3),
                                            rep("C1B10",3),rep("C1B11",3),rep("C2B2",3),
                                            rep("C2B3",3),rep("C2B4",3),rep("C2B5",3),
                                            rep("Out",6))))
###-------- create deseq object and process --------###
cds <- DESeqDataSetFromMatrix(countData = dta, colData=samples, design=~condition)
cds1 <- DESeq(cds)

# conduct comparisons
c1b1 <- results(cds1, contrast=c("condition","C1B1","Out"))
c1b2 <- results(cds1, contrast=c("condition","C1B2","Out"))
c1b3 <- results(cds1, contrast=c("condition","C1B3","Out"))
c1b4 <- results(cds1, contrast=c("condition","C1B4","Out"))
c1b5 <- results(cds1, contrast=c("condition","C1B5","Out"))
c1b6 <- results(cds1, contrast=c("condition","C1B6","Out"))
c1b7 <- results(cds1, contrast=c("condition","C1B7","Out"))
c1b8 <- results(cds1, contrast=c("condition","C1B8","Out"))
c1b9 <- results(cds1, contrast=c("condition","C1B9","Out"))
c1b10 <- results(cds1, contrast=c("condition","C1B10","Out"))
c1b11 <- results(cds1, contrast=c("condition","C1B11","Out"))
c2b2 <- results(cds1, contrast=c("condition","C2B2","Out"))
c2b3 <- results(cds1, contrast=c("condition","C2B3","Out"))
c2b4 <- results(cds1, contrast=c("condition","C2B4","Out"))
c2b5 <- results(cds1, contrast=c("condition","C2B5","Out"))


###-------- make merged table of results --------###
# using DEVis to make a combined table of results showing only significant DEGs, padj <0.05, log2fdc>2
# define the base directory for this analysis.  
base_dir <- "C:/Users/Matt/OneDrive - Gloucester Marine Genomics Institute/Desktop/rdata/oxr/results/"

# generate the directory structure for results.
create_dir_struct(base_dir)

# specify input directories.
cnt_dir        <- "results/counts/"
tgt_dir        <- "results/targets/"
init_data_paths(cnt_dir, tgt_dir)

# set significance cutoffs for p-values and log fold-change.
init_cutoffs(p_signif=0.05, lfc_cut=2)

# make a list of all contrasts
results_list <- list(c1b1, c1b2, c1b3, c1b4, c1b5, c1b6, c1b7, c1b8, c1b9, c1b10,
                     c1b11, c2b2, c2b3, c2b4, c2b5)

# aggregate differentially expressed genes across all contrasts
master_df <- create_master_res(results_list, filename="master_DE_list.txt", method="union", 
                               lfc_filter = TRUE)

###-------- filter individual results by padj <0.05 and save --------###
# filter by only significant results
c1b1Sig <- subset(c1b1, padj < 0.05)
c1b2Sig <- subset(c1b2, padj < 0.05)
c1b3Sig <- subset(c1b3, padj < 0.05)
c1b4Sig <- subset(c1b4, padj < 0.05)
c1b5Sig <- subset(c1b5, padj < 0.05)
c1b6Sig <- subset(c1b6, padj < 0.05)
c1b7Sig <- subset(c1b7, padj < 0.05)
c1b8Sig <- subset(c1b8, padj < 0.05)
c1b9Sig <- subset(c1b9, padj < 0.05)
c1b10Sig <- subset(c1b10, padj < 0.05)
c1b11Sig <- subset(c1b11, padj < 0.05)
c2b2Sig <- subset(c2b2, padj < 0.05)
c2b3Sig <- subset(c2b3, padj < 0.05)
c2b4Sig <- subset(c2b4, padj < 0.05)
c2b5Sig <- subset(c2b5, padj < 0.05)

# save output
write.csv(as.data.frame(c1b1Sig), file="results/c1b1results.csv")
write.csv(as.data.frame(c1b2Sig), file="results/c1b2results.csv")
write.csv(as.data.frame(c1b3Sig), file="results/c1b3results.csv")
write.csv(as.data.frame(c1b4Sig), file="results/c1b4results.csv")
write.csv(as.data.frame(c1b5Sig), file="results/c1b5results.csv")
write.csv(as.data.frame(c1b6Sig), file="results/c1b6results.csv")
write.csv(as.data.frame(c1b7Sig), file="results/c1b7results.csv")
write.csv(as.data.frame(c1b8Sig), file="results/c1b8results.csv")
write.csv(as.data.frame(c1b9Sig), file="results/c1b9results.csv")
write.csv(as.data.frame(c1b10Sig), file="results/c1b10results.csv")
write.csv(as.data.frame(c1b11Sig), file="results/c1b11results.csv")
write.csv(as.data.frame(c2b2Sig), file="results/c2b2results.csv")
write.csv(as.data.frame(c2b3Sig), file="results/c2b3results.csv")
write.csv(as.data.frame(c2b4Sig), file="results/c2b4results.csv")
write.csv(as.data.frame(c2b5Sig), file="results/c2b5results.csv")
