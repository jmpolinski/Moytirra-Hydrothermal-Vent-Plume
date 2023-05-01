###-------- load packages you need --------###
library(pheatmap)
library(ggplot2)
library(tidyr)
library(dplyr)

###-------- load data --------###
dta <- read.csv("input/sulfurde.csv", header=TRUE, row.names = 1)
sd <- read.csv("input/sulfursd.csv", header=TRUE, row.names = 1)
turb <- read.csv("input/turbid.csv", header = TRUE, row.names = 1)
geneName <- read.csv("input/sulfurNames.csv", header = TRUE, row.names = 1)
dt <- read.csv(file="input/sulfurPerc.csv", header = TRUE)

###-------- set colors and other formatting --------###
ann_colors = list(turbidity = c("firebrick", "yellow", "darkgreen"),
                  substrate = c("Sulfide->Sulfur"="brown",
                                "Sulfite->Trisulfide"="gold",
                                "Sulfate->APS"=  "darkorchid",
                                "APS->Sulfite"= "darkgreen",
                                "Thiosulfate->Sulfate" = "darkturquoise",
                                "APS->PAPS"  = "gray",
                                "PAPS->Sulfite" = "deepskyblue4",
                                "Sulfite->Sulfide" = "lightsalmon"))

# set breaks to center around zero and colors
col <- colorRampPalette(c("dodgerblue","white","red"))(256)
myBreaks <- c(seq(min(dta), 0, length.out=ceiling(256/2) + 1), 
                  seq(max(dta)/256, max(dta), length.out=floor(256/2)))

###-------- generate pheatmap plot --------###
pheatmap(dta, 
         breaks = myBreaks,
         color = col,
         gaps_row = c(1,6,12),
         gaps_col = c(6,8),
         cluster_cols=FALSE, 
         cluster_rows = FALSE,
         fontsize_row = 10, 
         annotation_col=turb,
         annotation_colors=ann_colors,
         annotation_row = geneName,
         show_rownames = TRUE,
         cellwidth = 15,
         cellheight = 15,
         display_numbers = matrix(ifelse(sd < 0.05, "*", ""),nrow(sd)))


###-------- generate supplementary figure 9 --------###
s.perc <- data.frame(dt[,c(1,2,3:19)])
colnames(s.perc) <- c("Gene","Function","C1B1","C1B2","C1B3","C1B4","C1B5","C1B6","C1B7","C1B8",
                      "C1B9","C1B10","C1B11","C2B2","C2B3","C2B4","C2B5","C1B12","C2B1")
s.perc <- s.perc %>% 
  gather(key="Bottle", value="value", -Gene, -Function) %>%
  mutate(Gene=factor(Gene)) %>%
  mutate(value=as.numeric(value))

s.perc$Bottle <-factor(s.perc$Bottle,  
                       c(levels = "C1B1","C1B2","C1B3","C1B4","C1B5","C1B6","C1B7","C1B8",
                         "C1B9","C1B10","C1B11","C2B2","C2B3","C2B4","C2B5","C1B12","C2B1")) 

s.perc$Gene <-as.factor(s.perc$Gene) 

ggplot(data=s.perc, aes(x=Bottle, y=Gene)) + 
  geom_point(aes(size=value, color=Function)) +  
  scale_size_continuous(range=c(0,10)) + 
  scale_color_manual(values=c("#ff8aff", "#bdbf00"))+
  theme_classic() + 
  theme(text = element_text(size = 12)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
