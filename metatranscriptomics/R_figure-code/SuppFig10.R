###-------- load packages you need --------###
library(pheatmap)

###-------- load data --------###
dta <- read.csv("input/nitrotde.csv", header=TRUE, row.names = 1)
sd <- read.csv("input/nitrotsd.csv", header=TRUE, row.names = 1)
turb <- read.csv("input/turbid.csv", header = TRUE, row.names = 1)

###-------- set colors and other formatting --------###
ann_colors = list(turbidity = c("firebrick", "yellow", "darkgreen"))

# set breaks to center around zero and colors
col <- colorRampPalette(c("dodgerblue","white","red"))(256)
myBreaks <- c(seq(min(dta), 0, length.out=ceiling(256/2) + 1), 
                  seq(max(dta)/256, max(dta), length.out=floor(256/2)))

###-------- generate pheatmap plot --------###
pheatmap(dta, 
         breaks = myBreaks,
         color = col,
         #gaps_row = c(5,12),
         gaps_col = c(6,8),
         cluster_cols=FALSE, 
         cluster_rows = FALSE,
         fontsize_row = 10, 
         annotation_col=turb,
         annotation_colors=ann_colors,
         #annotation_row = geneName,
         show_rownames = TRUE,
         #cellwidth = 15,
         #cellheight = 5,
         display_numbers = matrix(ifelse(sd < 0.05, "*", ""),nrow(sd))
         )

