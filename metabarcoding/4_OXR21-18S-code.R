###-------- load packages you need --------###
library("ggplot2")
library("vegan")
library("dplyr")
library("scales")
library("grid")
library("reshape2")
library("phyloseq")
library("RColorBrewer")
library("microbiome")
library("knitr")


###-------- load metabarcoding count data --------###
otutab <- read.csv("input/OXR-18S_zotutab.csv", check.names = F)
rownames(otutab)<- otutab[,1]
otutab <- otutab[,2:51]
otutab <- otu_table(otutab, taxa_are_rows = T)


###-------- load taxonomy classifications --------###
tax <- read.csv("input/OXR-18S_taxa_updated_merge.csv", header=F)
colnames(tax) <- c("otuID", "Kingdom", "Supergroup", "Division", "Class", "Order", "Family", "Genus", "Species", "Taxa")
head(tax)

# remove header
rownames(tax)<- taxa_names(otutab)
tax <- tax[,2:10]
tax <- as.matrix(tax)
tax <- tax_table(tax)
taxa_names(tax) <- taxa_names(otutab)


###-------- load sample metadata --------###
meta <- read.csv("input/oxr-meta.csv")
meta$replicate <- as.factor(meta$replicate)
meta$bottle <- factor(meta$bottle, levels = c("C1B1", "C1B2", "C1B3", "C1B4", "C1B5", "C1B6",
                                              "C1B7", "C1B8", "C1B9", "C1B10", "C1B11", "C1B12",
                                              "C2B1", "C2B2", "C2B3", "C2B4", "C2B5"))
meta$SampleID <- factor(meta$SampleID, levels = meta$SampleID)
meta <- sample_data(meta)
rownames(meta) <- meta$SampleID


###-------- create phyloseq object and curate --------###
physeq <- merge_phyloseq(otutab, tax, meta)

# remove unknowns
physeq <- subset_taxa(physeq, Kingdom!="unknown")

# observe read counts per sample
sample_sums(physeq)


###-------- generate alpha diversity measures --------###
alpha <- microbiome::alpha(physeq)
write.csv(alpha, "results/OXR-18S_alpha-diversity.csv")


###-------- Create plots for manuscript --------###
# prepare colors for plots
div_col_pals = brewer.pal.info[brewer.pal.info$category == 'div',]
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
div_col = unlist(mapply(brewer.pal, div_col_pals$maxcolors, rownames(div_col_pals)))
qual_col = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_pal = c(div_col, qual_col)

# Summary violin plot (supplementary figure 3)
mean_ci <- function(x){
  m = mean(x)
  se = sd(x)/sqrt(length(x))
  ql = qnorm(1-0.025)
  c('y'=m, 'ymin'=m-ql*se, 'ymax'=m+ql*se)
}

# generate object by Kingdom
p <- subset_taxa(physeq) %>%
  tax_glom(taxrank = "Kingdom") %>%                     #agglomerate to phylum level
  transform_sample_counts(function(x){x/sum(x)}) %>%   #transform to relative abundance
  psmelt()
write.csv(p, "results/OXR-18S_king.csv")

# merge in with archaea and bacteria, then reimport
p <- read.csv("input/OXR-16S_18S_king.csv")

ggplot(p, aes(x=Kingdom, fill=Kingdom, y=Abundance))+
  geom_violin()+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size = 12))+
  stat_summary(fun.data = mean_ci,
               geom="pointrange",
               fatten = 2,
               size=1)+
  scale_y_continuous('Ave Relative Read Abundnace', breaks = seq(0, 1, by=0.1))+
  scale_fill_manual(values = c("#0072B2", "#E69F00", "#999999"))+
  theme(legend.position = 'none')

# eukaryotic bar plot, merged rep (Figure 6)
phybottle <- merge_samples(physeq, "bottle") #merge replicates
a <- subset_taxa(phybottle) %>%
  tax_glom(taxrank = "Taxa") %>%                     #agglomerate to phylum level
  transform_sample_counts(function(x){x/sum(x)}) %>%   #transform to relative abundance
  psmelt()
write.csv(a, "results/OXR-18S_taxa.csv")

colourCount = length(unique(a$Taxa))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))


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
  guides(fill=guide_legend(ncol = 3))+
  facet_grid(.~type, scales = "free_x",space="free_x", 
             labeller = facet_labeller)+
  theme(strip.text.x = element_text(size = 12))

# metazoa bar plot, merged rep (supplementary figure 4)
phybottle <- merge_samples(physeq, "bottle") #merge replicates
metaz <- subset_taxa(phybottle, Taxa=="Opisthokonta-Metazoa") %>%
  tax_glom(taxrank = "Class") %>%                     #agglomerate to phylum level
  transform_sample_counts(function(x){x/sum(x)}) %>%   #transform to relative abundance
  psmelt()
write.csv(metaz, "results/OXR-18S_metaz.csv")

colourCount = length(unique(metaz$Class))
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

facetlab <- c("South Plume", "Above Vent", "North Plume", "Out")
facet_labeller <- function(variable,value){
  return(facetlab[value])
}

ggplot(arrange(metaz, Class), aes(x = factor(Sample, level = c('C1B1', 'C1B2', 'C1B3',
                                                          'C1B4', 'C1B5', 'C1B6',
                                                          'C1B7', 'C1B8', 'C1B9',
                                                          'C1B10', 'C1B11', 'C2B2',
                                                          'C2B3', 'C2B4', 'C2B5',
                                                          'C1B12', 'C2B1')),
                             y=Abundance, fill=Class)) +
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
  guides(fill=guide_legend(ncol = 3))+
  facet_grid(.~type, scales = "free_x",space="free_x", 
             labeller = facet_labeller)+
  theme(strip.text.x = element_text(size = 12))

# radiolarian bar plot, merged rep (supplementary figure5)
phybottle <- merge_samples(physeq, "bottle") #merge replicates
radio <- subset_taxa(phybottle, Taxa=="Rhizaria-Radiolaria") %>%
  tax_glom(taxrank = "Class") %>%                     #agglomerate to phylum level
  transform_sample_counts(function(x){x/sum(x)}) %>%   #transform to relative abundance
  psmelt()
write.csv(radio, "results/OXR-18S_radio.csv")

colourCount = length(unique(metaz$Class))
getPalette = colorRampPalette(brewer.pal(8, "Set3"))

facetlab <- c("South Plume", "Above Vent", "North Plume", "Out")
facet_labeller <- function(variable,value){
  return(facetlab[value])
}

ggplot(arrange(radio, Class), aes(x = factor(Sample, level = c('C1B1', 'C1B2', 'C1B3',
                                                               'C1B4', 'C1B5', 'C1B6',
                                                               'C1B7', 'C1B8', 'C1B9',
                                                               'C1B10', 'C1B11', 'C2B2',
                                                               'C2B3', 'C2B4', 'C2B5',
                                                               'C1B12', 'C2B1')),
                                  y=Abundance, fill=Class)) +
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
  guides(fill=guide_legend(ncol = 3))+
  facet_grid(.~type, scales = "free_x",space="free_x", 
             labeller = facet_labeller)+
  theme(strip.text.x = element_text(size = 12))

# cilliate bar plot, merged rep (supplementary figure 6)
phybottle <- merge_samples(physeq, "bottle") #merge replicates
cilia <- subset_taxa(phybottle, Taxa=="Alveolata-Ciliates") %>%
  tax_glom(taxrank = "Genus") %>%                     #agglomerate to phylum level
  transform_sample_counts(function(x){x/sum(x)}) %>%   #transform to relative abundance
  psmelt()
write.csv(cilia, "results/OXR-18S_cilia.csv")

colourCount = length(unique(metaz$Class))
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

facetlab <- c("South Plume", "Above Vent", "North Plume", "Out")
facet_labeller <- function(variable,value){
  return(facetlab[value])
}

ggplot(arrange(cilia, Genus), aes(x = factor(Sample, level = c('C1B1', 'C1B2', 'C1B3',
                                                               'C1B4', 'C1B5', 'C1B6',
                                                               'C1B7', 'C1B8', 'C1B9',
                                                               'C1B10', 'C1B11', 'C2B2',
                                                               'C2B3', 'C2B4', 'C2B5',
                                                               'C1B12', 'C2B1')),
                                  y=Abundance, fill=Genus)) +
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
  guides(fill=guide_legend(ncol = 3))+
  facet_grid(.~type, scales = "free_x",space="free_x", 
             labeller = facet_labeller)+
  theme(strip.text.x = element_text(size = 12))

# PCoA of eukaryotic communities across the sample set
# first rarefy to even depth
rare <- rarefy_even_depth(physeq)
rare
sample_sums(rare)

# create the ordination
pcoa <- ordinate(physeq = rare, method = "PCoA", distance = "bray")

# plot the ordination
plot_ordination(physeq = rare, ordination = pcoa, color = "bottle", shape = "location", 
                title = "PCoA of 18S Data") +
  geom_point(aes(color = bottle), size = 5)+
  scale_shape_manual(values=c(2, 1, 9, 0))+
  theme_bw()+
  theme(text = element_text(size = 14))+
  scale_color_manual(values = qual_col[10:35])+
  guides(shape = guide_legend(title = "Location"))+
  guides(color = guide_legend(title = "Niskin", ncol = 2))

