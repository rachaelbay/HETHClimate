setwd("~/Documents/IndividualProjects/Nicole/HETH/GWAS/")

library(vroom)
library(tidyverse)
library(fastman)
library(cowplot)
library(shades)
library(viridis)

m="filtLRT/sub.HETH.contemp.male.PCcovar.lrt3.gz"
f="filtLRT/sub.HETH.contemp.female.PCcovar.lrt3.gz"

##Read in data for males
m.lrt <- vroom(m,delim="\t")
f.lrt <- vroom(f,delim="\t")

##Filter out P=O
m.sub <- m.lrt %>% filter(P!=0)
f.sub <- f.lrt %>% filter(P!=0)

##plot data
m.plot <- m.sub %>% filter(P<0.05)
f.plot <- f.sub %>% filter(P<0.05)

##Replace chromosome names
chr <- read.delim("sequence_report.tsv")
chr$RefSeq.seq.accession <- sub("_","-",chr$RefSeq.seq.accession)
m.merge <- merge(m.plot,chr,by.x="Chromosome",by.y="RefSeq.seq.accession")
m.merge$Chromosome.name <- factor(m.merge$Chromosome.name,levels=c(order(as.numeric(as.character(unique(m.merge$Chromosome.name)))), "scaffolds"))
m.merge$Chromosome.name[is.na(m.merge$Chromosome.name)] <- "scaffolds"
f.merge <- merge(f.plot,chr,by.x="Chromosome",by.y="RefSeq.seq.accession")
f.merge$Chromosome.name <- factor(f.merge$Chromosome.name,levels=c(order(as.numeric(as.character(unique(f.merge$Chromosome.name)))), "scaffolds"))
f.merge$Chromosome.name[is.na(f.merge$Chromosome.name)] <- "scaffolds"


cols <- viridis(10)
dark <- cols
light <- lightness(cols, scalefac(1.6))
m.cols <- c(dark[3],light[3])
f.cols <- c(dark[5],light[5])

m.man <- fastman_gg(m.merge,chr="Chromosome.name",bp="Position",p="P",speedup=T,col=m.cols,
           suggestiveline = NA,genomewideline=-log10(1e-5),ylim=c(0,15),xlab="",cex=2) +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + ggtitle("Male")
f.man <- fastman_gg(f.merge,chr="Chromosome.name",bp="Position",p="P",speedup=T,col=f.cols,
           suggestiveline = NA,genomewideline=-log10(1e-5),ylim=c(0,15),xlab="",cex=2) +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + ggtitle("Female")
m.qq <- fastqq_gg(m.sub$P, speedup=T,fix_zero=FALSE,cex=2)
f.qq <- fastqq_gg(f.sub$P, speedup=T,fix_zero=FALSE,cex=2)

jpeg("../Figures/New_07.24/TarsusLengthGWAS.jpg",width=2700,height=1500,units="px",res=250)
plot_grid(m.man,m.qq[[1]],f.man,f.qq[[1]],nrow=2,labels=c("A","B","C","D"),rel_widths = rep(c(3,1),2))
dev.off()

