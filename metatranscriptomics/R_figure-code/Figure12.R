###-------- load packages you need --------###
library(ggplot2)

###-------- load data --------###
dta <- read.csv("input/growthmat.csv", header=TRUE)

###-------- set colors and other formatting --------###
cbPalette <- c("#F0E442", "#0072B2", "#D55E00")

###-------- generate plot --------###
ggplot(dta, aes(x=factor(bottle, level=c('C1B1', 'C1B2', 'C1B3', 'C1B4','C1B5',
                                         'C1B6', 'C1B7', 'C1B8', 'C1B9', 'C1B10',
                                         'C1B11', 'C2B2', 'C2B3', 'C2B4', 'C2B5')),
                y=LFC)) +
  geom_point(size=4, aes(color = kegg, shape = sig))+
  theme_bw() +
  theme(text = element_text(size=12))+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(.~location,scales="free",space="free")#+
  scale_y_continuous(expand = c(0,0))

  