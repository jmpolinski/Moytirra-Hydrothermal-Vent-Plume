###-------- load packages you need --------###
library("ggplot2")
library("vegan")
library("dplyr")
library("scales")
library("grid")
library("reshape2")
library("phyloseq")
library("dada2")
library("microbiome")
library("knitr")
library("RColorBrewer")
library("stringr")

###-------- load metabarcoding count data --------###
otutab <- read.csv("input/OXR-16S_zotutab.csv", check.names = F)
rownames(otutab)<- otutab[,1]
otutab <- otutab[,2:52]
otutab <- otu_table(otutab, taxa_are_rows = T)

###-------- load taxonomy classifications --------###
tax <- read.csv("input/OXR21-16S_taxonomy_tree-updated_combined.csv", header=F)
rownames(tax)<- taxa_names(otutab)
tax <- tax[,2:8]
tax <- as.matrix(tax)
tax <- tax_table(tax)
taxa_names(tax) <- taxa_names(otutab)

###-------- load sample metadata --------###
meta <- read.csv("input/oxr-meta.csv")
meta$replicate <- as.factor(meta$replicate)
meta$bottle <- factor(meta$bottle, 
                         levels = c("C1B1", "C1B2", "C1B3", "C1B4", "C1B5", "C1B6",
                                  "C1B7", "C1B8", "C1B9", "C1B10", "C1B11", "C1B12",
                                  "C2B1", "C2B2", "C2B3", "C2B4", "C2B5"))
meta$SampleID <- factor(meta$SampleID, levels = meta$SampleID)
meta <- sample_data(meta)
rownames(meta) <- meta$SampleID

###-------- create phyloseq object and curate --------###
physeq <- merge_phyloseq(otutab, tax, meta)

# fix taxonomy column names
colnames(tax_table(physeq))
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", 
                                 "Family", "Genus", "Taxa")

# remove unknowns and chloroplasts/cyanobacteria
physeq <- subset_taxa(physeq, Kingdom!="unknown")
physeq <- subset_taxa(physeq, Phylum!="Cyanobacteria")

# observe read counts per sample
sample_sums(physeq)


###-------- generate alpha diversity measures --------###
alpha <- microbiome::alpha(physeq)
write.csv(alpha, "results/OXR-16S_alpha-diversity.csv")


###-------- Create plots for manuscript --------###
# prepare colors for plots
div_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
div_col = unlist(mapply(brewer.pal, div_col_pals$maxcolors, rownames(div_col_pals)))
qual_col = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_pal = c(div_col, qual_col)

# summary violin plot
p <- subset_taxa(physeq) %>%
  tax_glom(taxrank = "Kingdom") %>%                     #agglomerate to phylum level
  transform_sample_counts(function(x){x/sum(x)}) %>%   #transform to relative abundance
  psmelt()
write.csv(p, "results/OXR-16S_king.csv")

mean_ci <- function(x){
  m = mean(x)
  se = sd(x)/sqrt(length(x))
  ql = qnorm(1-0.025)
  c('y'=m, 'ymin'=m-ql*se, 'ymax'=m+ql*se)
}

ggplot(p, aes(x=Kingdom, fill=Kingdom, y=Abundance))+
  geom_violin()+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size = 12))+
  stat_summary(fun.data = mean_ci,
               geom="pointrange",
               fatten = 2,
               size=1)+
  scale_y_continuous('Ave Relative Read Abundnace', breaks = seq(0, 1, by=0.05))+
  scale_fill_manual(values = c("#0072B2", "#E69F00"))+
  theme(legend.position = 'none')

# Archaea bar plot, merged reps (Figure 4)
df <- as.data.frame(lapply(sample_data(physeq),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq)
sample_data(physeq) <- sample_data(df)

phybottle <- merge_samples(physeq, "bottle") #merge replicates
a <- subset_taxa(phybottle, Kingdom=="Archaea") %>%
  tax_glom(taxrank = "Taxa") %>%                     #agglomerate to phylum level
  transform_sample_counts(function(x){x/sum(x)}) %>%   #transform to relative abundance
  psmelt()

colourCount = length(unique(a$Taxa))
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

facetlab <- c("South Plume", "Above Vent", "North Plume", "Out")
facet_labeller <- function(variable,value){
  return(facetlab[value])
}

ggplot(arrange(a, Taxa), aes(x = factor(Sample, level = c('C1B1', 'C1B2', 'C1B3',
                                                          'C1B4', 'C1B5', 'C1B6',
                                                          'C1B7', 'C1B8', 'C1B9',
                                                          'C1B10', 'C1B11', 'C2B2',
                                                          'C2B3', 'C2B4', 'C2B5',
                                                          'C1B12', 'C2B1')),
                             y=Abundance, fill=Taxa)) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12, hjust = 0.5), 
        legend.text=element_text(size=12),
        legend.position = "top") +
  ylab("Relative Read Abundance") +
  theme(axis.text = element_text(size = 12))+
  scale_fill_manual(values = getPalette(colourCount)) +
  scale_y_continuous(limits = c(0,1.001), expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  guides(fill=guide_legend(ncol = 2))+
  facet_grid(.~type, scales = "free_x",space="free_x", 
             labeller = facet_labeller)+
  theme(strip.text.x = element_text(size = 12))

# Bacteria bar plot, merged reps (Figure 5)
phybottle <- merge_samples(physeq, "bottle") #merge replicates
b <- subset_taxa(phybottle, Kingdom=="Bacteria") %>%
  tax_glom(taxrank = "Taxa") %>%                     #agglomerate to phylum level
  transform_sample_counts(function(x){x/sum(x)}) %>%   #transform to relative abundance
  psmelt()

colourCount = length(unique(b$Taxa))
getPalette = colorRampPalette(brewer.pal(11, "Spectral"))

ggplot(arrange(b, Taxa), aes(x = factor(Sample, level = c('C1B1', 'C1B2', 'C1B3',
                                                          'C1B4', 'C1B5', 'C1B6',
                                                          'C1B7', 'C1B8', 'C1B9',
                                                          'C1B10', 'C1B11', 'C2B2',
                                                          'C2B3', 'C2B4', 'C2B5',
                                                          'C1B12', 'C2B1'))
                             , y=Abundance, fill=Taxa)) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=12, hjust = 0.5), 
        legend.text=element_text(size=12),
        legend.position = "top") +
  ylab("Relative Read Abundance") +
  theme(axis.text = element_text(size = 12))+
  scale_fill_manual(values = col_pal) +
  scale_y_continuous(limits = c(0,1.001), expand = c(0,0))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  guides(fill=guide_legend(ncol = 3))+
  facet_grid(.~type, scales = "free_x",space="free_x", 
             labeller = facet_labeller)+
  theme(strip.text.x = element_text(size = 12))

# PCoA of prokaryotic communities across the sample set (Supplementary Figure 2)
# first rarefy to even depth
rare <- rarefy_even_depth(physeq)
rare
sample_sums(rare)

# create the ordination
pcoa <- ordinate(physeq=rare, method = "PCoA", distance = "bray")

# plot the ordination
plot_ordination(physeq=rare, ordination = pcoa, color = "bottle", shape= "location", 
                title = "PCoA of 16S Data") +
  geom_point(size = 5)+
  scale_shape_manual(values=c(2, 1, 9, 0))+
  theme_bw()+
  theme(text = element_text(size = 14))+
  scale_color_manual(values = qual_col[10:35])+
  guides(shape = guide_legend(title = "Location"))+
  guides(color = guide_legend(title = "Niskin", ncol = 2))


#------- Use DESeq2 to test for significant abundance patterns -------#
library(DESeq2)
diagdds = phyloseq_to_deseq2(physeq, ~ type)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
head(sigtab)
write.csv(sigtab, "results/OXR-16S_deseq2.csv")

