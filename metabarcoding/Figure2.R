###-------- load packages you need --------###
library(corrplot)
library(qgraph)

###-------- load data --------###
cordta <- read.table(file="input/cormat.txt",
                     sep="\t", row.names = 1, header=TRUE)

###-------- generate Spearman Rank correlation and test for significance --------###

# load function to get p-values
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# calculate spearman rank correlation
res2 <-cor(cordta, method = c("spearman"))
pmat <- cor.mtest(cordta)

# adjust p-values for multiple testing with Benjamin-Hochberge method
pAdj <- p.adjust(pmat, method = "BH")

# convert padj values back into matrix for plotting
resAdj <- matrix(pAdj, ncol = dim(pmat))

# manually curate labels
colnames(resAdj) <- c("Temp", "Salinity", "Depth", "Turbidity",
                      "18S Richness", "16S Richness", "PO4", "NO23", "TDP", "TDN", "Si",
                      "SO4", "DOC", "DIC", "DON", "DOP", "break18", "break16")
rownames(resAdj) <- c("Temp", "Salinity", "Depth", "Turbidity",
                      "18S Richness", "16S Richness", "PO4", "NO23", "TDP", "TDN", "Si",
                      "SO4", "DOC", "DIC", "DON", "DOP", "break18", "break16")
colnames(res2) <- c("Temp", "Salinity", "Depth", "Turbidity",
                    "18S Richness", "16S Richness", "PO4", "NO23", "TDP", "TDN", "Si",
                    "SO4", "DOC", "DIC", "DON", "DOP", "break18", "break16")
rownames(res2) <- c("Temp", "Salinity", "Depth", "Turbidity",
                    "18S Richness", "16S Richness", "PO4", "NO23", "TDP", "TDN", "Si",
                    "SO4", "DOC", "DIC", "DON", "DOP", "break18", "break16")

# generate correlation plot
corrplot(res2, type="upper", tl.col="black", diag=FALSE,
         p.mat = resAdj, sig.level = c(.001, .01, .05),
         pch.cex = 1.3, insig = "label_sig", pch.col = "white", tl.cex=1.5)
