###-------- load packages you need --------###
library(ggplot2)
library(dplyr)
library(gridExtra)

###-------- load data --------###
df <- read.csv(file="input/shannon_all_16s.csv", header = TRUE) 
df2 <- read.csv(file="input/shannon_all_18s.csv", header = TRUE) 

###-------- Genearte plots --------###
# plot, comparison by distance 16S
pro <- ggplot(df, aes(factor(p2, level= c('Above Plume', 'Below Plume', 'South of Vent',
                                   'Vent', 'North of Vent')), value))+
  geom_boxplot(aes(fill = factor(location)))+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values = c("firebrick1", "dodgerblue", "gold"))+
  labs(y = "Shannon Index")+ 
  ggtitle("16S Diversity")+
  theme(text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# plot, comparison by distance 18S
euk <- ggplot(df2, aes(factor(p2, level= c('Above Plume', 'Below Plume', 'South of Vent',
                                    'Vent', 'North of Vent')), value))+
  geom_boxplot(aes(fill = factor(location)))+
  theme_classic()+
  theme(legend.position="none")+
  scale_fill_manual(values = c("firebrick1", "dodgerblue", "gold"))+
  labs(y = "Shannon Index")+ 
  ggtitle("18S Diversity")+
  theme(text = element_text(size=14))+
  theme(plot.title = element_text(size=14))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# plot both stacked on top of each other
grid.arrange(pro, euk, ncol = 1, heights = c(1, 1))


###-------- test for significance --------###
# Kruskal-Wallis test by biotype
# generate summary stats for 16S
group_by(df, p2) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
    
  )

# perform kruskal test, 16S data
kruskal.test(value ~ p2, data = df)
pairwise.wilcox.test(df$value, df$p2,
                     p.adjust.method = "BH")

# generate summary stats for 18S
group_by(df2, p2) %>%
  summarise(
    count = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
  )

# perform kruskal test, 18S data
kruskal.test(value ~ p2, data = df2)
pairwise.wilcox.test(df2$value, df2$p2,
                     p.adjust.method = "BH")
