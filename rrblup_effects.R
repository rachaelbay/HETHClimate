setwd("~/Documents/IndividualProjects/Nicole/HETH/cands")

library(tidyverse)
library(rrBLUP)

##Color palette
cols <- viridis(10)
light <- lightness(cols, scalefac(1.3))
dark <- cols
m.cols <- light[3]
f.cols <- light[5]
a.cols <- alpha(c(f.cols,m.cols),alpha=0.7)

##Read in dosage information (from beagle)
dos <- read.delim("cutoff5/impute.HETH.contemp.female.tarsPCcovar.lrt1.beagle.gz.dose.gz", sep=" ")
dos.t <- t(dos[,4:ncol(dos)]) - 1

#Get metadata
meta <- read.csv("../GWAS/HETH.219.sample",sep=" ")[-1,]
meta[,2:ncol(meta)] <- apply(meta[,2:ncol(meta)],2,as.numeric)
meta$Sex <- as.factor(meta$Sex)

##Sex sexes (2 is male)
in.sex=1
out.sex=2
morph="Bill"
rel.morph=paste("rel",morph,sep="")
morph.ind <- which(names(meta)==morph)
rel.morph.ind <- which(names(meta)==rel.morph)

#Get effect size primary analysis: Sex "2" is male
sub.dos <- dos.t[(meta$Sex==in.sex & meta$Year>=2005),]
sub.meta <- meta[(meta$Sex==in.sex & meta$Year>=2005),]
morph.beta <- mixed.solve(sub.meta[,morph.ind],sub.dos,K=NULL,
                          SE=F,return.Hinv = F)

#Get effect size other sex
sex.dos <- dos.t[(meta$Sex==out.sex & meta$Year>=2005),]
sex.meta <- meta[(meta$Sex==out.sex & meta$Year>=2005),]
sex.beta <- mixed.solve(sex.meta[,morph.ind],sex.dos,K=NULL, 
                          SE=F,return.Hinv = F)

#Get effect size relative morphology
rel.sub.dos <- dos.t[(meta$Sex==in.sex & meta$Year>=2005),]
rel.sub.meta <- meta[(meta$Sex==in.sex & meta$Year>=2005),]
rel.morph.beta <- mixed.solve(rel.sub.meta[,rel.morph.ind],rel.sub.dos,K=NULL,
                          SE=F,return.Hinv = F,method="ML")

#Get effect size relative morphology other sex
rel.sex.dos <- dos.t[(meta$Sex==out.sex & meta$Year>=2005),]
rel.sex.meta <- meta[(meta$Sex==out.sex & meta$Year>=2005),]
rel.sex.beta <- mixed.solve(rel.sex.meta[,rel.morph.ind],sex.dos,K=NULL, 
                        SE=F,return.Hinv = F)

#Get effect size temporal analysis - same sex
temp.dos <- dos.t[(meta$Sex==in.sex),]
dos.meta <- meta[(meta$Sex==in.sex),]
temp.beta <- mixed.solve(dos.meta$Year,temp.dos,K=NULL,
                         SE=F,return.Hinv = F)

#Get effect size temporal - other sex
st.dos <- dos.t[(meta$Sex==out.sex),]
st.meta <- meta[(meta$Sex==out.sex),]
st.beta <- mixed.solve(st.meta$Year,st.dos,K=NULL,
                         SE=F,return.Hinv = F)


frame <- data.frame(morph=morph.beta$u,rel.morph=rel.morph.beta$u,
                    temp=temp.beta$u,sex.temp=st.beta$u,
                    sex.morph=sex.beta$u,rel.sex.morp=rel.sex.beta$u)
saveRDS(frame,"FemaleBill.rds")

ggplot(frame,aes(x=morph,y=temp)) + geom_point() + theme_classic() + stat_smooth(method="lm")

dim(frame)
cor.test(frame$morph,frame$temp,method="spearman") #morphology v. time
cor.test(frame$m,frame$sm, method="spearman") #sex comparison
