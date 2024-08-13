setwd("~/Documents/IndividualProjects/Nicole/HETH/het/")

library(tidyverse)

##Color palette
##Color palette
cols <- viridis(10)
light <- lightness(cols, scalefac(1.3))
dark <- cols
m.cols <- light[3]
f.cols <- light[5]
a.cols <- alpha(c(f.cols,m.cols),alpha=0.7)

##Read in individuals and get file names
inds <- read.delim("../Data/all.219.inds.txt",header=F)
head(inds)
inds$IDs <- substr(inds$V1,2,7)
files <- paste(inds[,1],".ml",sep="")

#Get metadata
meta <- read.csv("../HETH_meta.csv")
ordermeta <- meta[match(inds$IDs,meta$Sample),]
ordermeta$het <- NA

#Extract heterozygosity from files
for (i in 1:length(files)) {
  a <- scan(files[i])
  ordermeta$het[i] <- a[2]/sum(a)
}
het <- ggplot(ordermeta,aes(x=Year,y=het,color=Sex)) + geom_jitter(size=2) + theme_classic() +
  scale_color_manual(labels=c("Female","Male"),values=a.cols) +
  ylab("Heterozygosity")
saveRDS(het,file="hetplot.rds")

m1 <- lm(het~Year*Sex,data=ordermeta)
m2 <- lm(het~Year+Sex,data=ordermeta)
m3 <- lm(het~Year,data=ordermeta)
anova(m2,m3)
summary(m3)
