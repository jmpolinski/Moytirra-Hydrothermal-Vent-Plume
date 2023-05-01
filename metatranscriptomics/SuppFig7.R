###-------- load packages you need --------###
library(FactoMineR)
library(factoextra)

###-------- import data for pca, data must be transformed --------###
dta <- read.csv("input/pathtoko_pca.csv", header=TRUE, row.names = 1)

###-------- calculate PCA --------###
res.pca <- PCA(dta[,-155], graph = FALSE)

###-------- generate plot --------###
fviz_pca_ind(res.pca,
             geom.ind = "point", # show points only (but not "text")
             pointshape = 21,
             pointsize = 2.5,
             fill.ind = dta$location,
             col.ind = dta$location, # color by groups
             palette = c("firebrick1", "green4", "dodgerblue", "gold"),
             addEllipses = TRUE, ellipse.type = "confidence", ellipse.level=0.95,
             mean.point = FALSE,
             legend.title = "Location")+theme_bw()
