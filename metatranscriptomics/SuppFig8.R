###-------- load packages you need --------###
library(pheatmap)

###-------- load data --------###
dta <- read.csv("input/pathtoko_all.csv", header=TRUE, row.names = 1)
geneName <- read.csv("input/pathtokoNames.csv", header = TRUE, row.names = 1)
turb <- read.csv("input/turbid_withOut.csv", header = TRUE, row.names = 1)

###-------- set colors and other formatting --------###
ann_colors = list(turbidity = c("firebrick", "yellow", "darkgreen"),
                  class = c("Amino acid metabolism"="#EB6662",
                           "Carbohydrate metabolism"=  "#F7B172",
                           "Energy metabolism"= "#F7D37E",
                           "Nucleotide metabolism" = "#82C881",
                           "Translation"  = "#1D8F94",
                           "Xenobiotics biodegradation and metabolism" = "#203D85"))

# set breaks to center around zero and colors
col <- colorRampPalette(c("white","#3AA2F0","#005797"))(256)
myBreaks <- c(seq(min(dta), 0, length.out=ceiling(256/2) + 1), 
                  seq(max(dta)/256, max(dta), length.out=floor(256/2)))

###-------- generate pheatmap plot --------###
pheatmap(dta, 
         cluster_cols=FALSE,
         cluster_rows = FALSE,
         color = col,
         gaps_col = c(6,8,15),
         fontsize_row = 10, 
         annotation_col=turb, 
         annotation_colors=ann_colors,
         annotation_row = geneName,
         show_rownames = TRUE)

