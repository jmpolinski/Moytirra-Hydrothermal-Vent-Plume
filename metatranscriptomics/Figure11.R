###-------- load packages you need --------###
library(vegan)
library(ggplot2)
library(gridExtra)

###-------- load metatranscriptomic count data --------###
# load normalized data
df <- read.csv("input/normcount.csv", 
               header=TRUE, row.names =1)

# remove rows with all zeros
df1 <- df[rowSums(df[])>0,]

# transform rows to columns
df1 <- as.data.frame(t(df1))


###-------- conduct bray-curtis analysis --------###
bray <- vegdist(df1, method = "bray")
bray <- as.matrix(bray)

# convert to data frame and save
y <- as.data.frame.table(as.matrix(bray))|>
  transform(Var1 = as.character(Var1), Var2 = as.character(Var2)) |>
  subset(Var1<Var2)
write.table(y,'results/brayfun.csv',sep = ",")


###-------- load metabarcoding count data --------###
at <- read.csv("input/OXRall_zotutab.csv",  header=TRUE, row.names = 1)

# remove rows with all zeros
at <- at[rowSums(at[])>0,]

# transform rows to columns
af <- as.data.frame(t(at))

# square-root transformation
sqaf <- sqrt(af)

###-------- conduct bray-curtis analysis --------###
bray <- vegdist(sqaf, method = "bray")
bray <- as.matrix(bray)
write.table(bray,'results/braydivmat2.csv',sep = ",")

# convert to data frame and save
x <- as.data.frame.table(as.matrix(bray))|>
  transform(Var1 = as.character(Var1), Var2 = as.character(Var2)) |>
  subset(Var1<Var2)
write.table(x,'results/braydiv.csv',sep = ",")

# NOTE: a new .csv file was then created from these results by manually extracting comparisons
# between any sample and samples C1B7 and C1B8
# Results were then annotated for distance from Moytirra (C1B7 or C1B8) and designated as 
# "upstream" or "downstream" depending on whether north or south of the vent, respectively

###-------- plot results --------###
# load data
dta <- read.csv("results/bray_df3.csv", header=TRUE)

# stacked boxplots with lm trend line away from vent site

f <- ggplot(dta, aes(factor(Distance), fun)) +
  geom_boxplot() +
  geom_jitter()+
  geom_smooth(aes(group=group), method="lm")+
  labs(x="Distance from Moytirra (m)",
       y="Functional dissimilarity")+
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

d <- ggplot(dta, aes(factor(Distance), div)) +
  geom_boxplot() +
  geom_jitter()+
  geom_smooth(aes(group=group), method="lm")+
  labs(x="Distance from Moytirra (m)",
       y="Diversity dissimilarity")+
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))


grid.arrange(f, d, ncol = 1, heights = c(1, 1))

###-------- generate linear regression stats for each trend line --------###
###-------- stats were placed manually into figure using Inkscape --------###
# stats on linear regression - diversity
dfu <-dta[ which(dta$group=='Upstream'), ]
dfud.lm <- lm(formula = div ~ Distance, data = dfu) 
summary(dfud.lm)

dfd <-dta[ which(dta$group=='Downstream'), ]
dfdd.lm <- lm(formula = div ~ Distance, data = dfd) 
summary(dfdd.lm)

# stats on linear regression - function
dfu <-dta[ which(dta$group=='Upstream'), ]
dfuf.lm <- lm(formula = fun ~ Distance, data = dfu) 
summary(dfuf.lm)

dfd <-dta[ which(dta$group=='Downstream'), ]
dfdf.lm <- lm(formula = fun ~ Distance, data = dfd) 
summary(dfdf.lm)

