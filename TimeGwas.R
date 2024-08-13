setwd("~/Documents/IndividualProjects/Nicole/HETH/GWAS/")

library(vroom)
library(tidyverse)
library(fastman)
library(cowplot)
library(shades)

y="filtLRT/sub.HETH.all.female.Year.PCCovar.lrt0.gz"
bt="filtLRT/sub.HETH.all.female.Year.PCCovar.lrt1.gz"
wt="filtLRT/sub.HETH.all.female.Year.PCCovar.lrt3.gz"
wp="filtLRT/sub.HETH.all.female.Year.PCCovar.lrt4.gz"

##Read in data for males
y.lrt <- vroom(y,delim="\t")
bt.lrt <- vroom(bt,delim="\t")
wt.lrt <- vroom(wt,delim="\t")
wp.lrt <- vroom(wp,delim="\t")

##Filter out P=O
y.sub <- y.lrt %>% filter(P!=0)
bt.sub <- bt.lrt %>% filter(P!=0)
wt.sub <- wt.lrt %>% filter(P!=0)
wp.sub <- wp.lrt %>% filter(P!=0)

##plot data
y.plot <- y.sub %>% filter(P<0.05)
bt.plot <- bt.sub %>% filter(P<0.05)
wt.plot <- wt.sub %>% filter(P<0.05)
wp.plot <- wp.sub %>% filter(P<0.05)

##Replace chromosome names
chr <- read.delim("sequence_report.tsv")
chr$RefSeq.seq.accession <- sub("_","-",chr$RefSeq.seq.accession)
y.merge <- merge(y.plot,chr,by.x="Chromosome",by.y="RefSeq.seq.accession")
y.merge$Chromosome.name <- factor(y.merge$Chromosome.name,levels=c(order(as.numeric(as.character(unique(y.merge$Chromosome.name)))), "scaffolds"))
y.merge$Chromosome.name[is.na(y.merge$Chromosome.name)] <- "scaffolds"
bt.merge <- merge(bt.plot,chr,by.x="Chromosome",by.y="RefSeq.seq.accession")
bt.merge$Chromosome.name <- factor(bt.merge$Chromosome.name,levels=c(order(as.numeric(as.character(unique(bt.merge$Chromosome.name)))), "scaffolds"))
bt.merge$Chromosome.name[is.na(bt.merge$Chromosome.name)] <- "scaffolds"
wt.merge <- merge(wt.plot,chr,by.x="Chromosome",by.y="RefSeq.seq.accession")
wt.merge$Chromosome.name <- factor(wt.merge$Chromosome.name,levels=c(order(as.numeric(as.character(unique(wt.merge$Chromosome.name)))), "scaffolds"))
wt.merge$Chromosome.name[is.na(wt.merge$Chromosome.name)] <- "scaffolds"
wp.merge <- merge(wp.plot,chr,by.x="Chromosome",by.y="RefSeq.seq.accession")
wp.merge$Chromosome.name <- factor(wp.merge$Chromosome.name,levels=c(order(as.numeric(as.character(unique(wp.merge$Chromosome.name)))), "scaffolds"))
wp.merge$Chromosome.name[is.na(wp.merge$Chromosome.name)] <- "scaffolds"


cols <- viridis(10)
dark <- cols
light <- lightness(cols, scalefac(1.6))
m.cols <- c(dark[3],light[3])
f.cols <- c(dark[5],light[5])

y.man <- fastman_gg(y.merge,chr="Chromosome.name",bp="Position",p="P",speedup=T,col=f.cols,
                    suggestiveline = NA,genomewideline=-log10(1e-5),ylim=c(0,15),xlab="",cex=2) +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + ggtitle("Year")
bt.man <- fastman_gg(bt.merge,chr="Chromosome.name",bp="Position",p="P",speedup=T,col=f.cols,
                    suggestiveline = NA,genomewideline=-log10(1e-5),ylim=c(0,15),xlab="",cex=2) +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + ggtitle("Breeding Temperature")
wt.man <- fastman_gg(wt.merge,chr="Chromosome.name",bp="Position",p="P",speedup=T,col=f.cols,
                    suggestiveline = NA,genomewideline=-log10(1e-5),ylim=c(0,15),cex=2) +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + ggtitle("Wintering Temperature")
wp.man <- fastman_gg(wp.merge,chr="Chromosome.name",bp="Position",p="P",speedup=T,col=f.cols,
                     suggestiveline = NA,genomewideline=-log10(1e-5),ylim=c(0,15),cex=2) +
  theme(plot.margin = unit(c(0.2,0.5,0,0), "cm")) + ggtitle("Wintering Precipitation")

y.qq <- fastqq_gg(y.sub$P, speedup=T,fix_zero=FALSE,cex=2)
bt.qq <- fastqq_gg(bt.sub$P, speedup=T,fix_zero=FALSE,cex=2)
wt.qq <- fastqq_gg(wt.sub$P, speedup=T,fix_zero=FALSE,cex=2)
wp.qq <- fastqq_gg(wp.sub$P, speedup=T,fix_zero=FALSE,cex=2)



jpeg("Figures/FemaleTimeGWAS.jpg",width=2700,height=2600,units="px",res=250)
plot_grid(y.man,y.qq[[1]],bt.man,bt.qq[[1]],wt.man,wt.qq[[1]],wp.man,wp.qq[[1]],nrow=4,labels=c("A","B","C","D","E","F","G","H"),rel_widths = rep(c(3,1),4))
dev.off()

