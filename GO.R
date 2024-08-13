setwd("~/Documents/IndividualProjects/Nicole/HETH/annot/")

library(biomaRt)
library(topGO)
library(tidyverse)
library(stringr)
library(vroom)

options(timeout = 30000)

########

ld <- vroom("HETH_0.9_coding.txt")
ld <- separate(ld,col="annotation",into=c("ID","Dbxref","Name","gbkey","gene","gene_biotype"),sep=";")
ld <- separate(ld,col="gene",into=c("g","gene"),sep="=")
ld$chr <- paste(gsub("chr","NC-",ld$chromosome),".1",sep="")
ld$chr <- gsub("scaff","NW-",ld$chr)
ld$snp <- paste(ld$chr,ld$region_start,sep="_")

###Grab zebra finch annotations
###Get annotations from biomaRt
# ensembl = useEnsembl(biomart="ensembl")
# zf = useEnsembl(biomart="ensembl", dataset="tguttata_gene_ensembl")
# listAttributes(zf)
# zf_genes <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "description",
#                                     "chromosome_name", "start_position", "end_position",
#                                     "external_gene_name", "gene_biotype", "go_id", "name_1006", "definition_1006",
#                                     "namespace_1003","ensembl_transcript_id"),
#                      values = T,
#                      mart = zf)
# head(zf_genes)
# saveRDS(zf_genes,"ZFgenes.rds")

###Grab chicken annotations
# ch = useEnsembl(biomart="ensembl",dataset="ggallus_gene_ensembl")
# ch_genes <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id", "description",
#                                     "chromosome_name", "start_position", "end_position",
#                                     "external_gene_name", "gene_biotype", "go_id", "name_1006", "definition_1006",
#                                     "namespace_1003","ensembl_transcript_id"),
#                      values = T,
#                      mart = ch)
# saveRDS(ch_genes,"CHgenes.rds")

ch_genes <- readRDS("CHgenes.rds")
zf_genes <- readRDS("ZFgenes.rds")

###Add go terms
genes <- na.omit(unique(ld$gene))
gene_data <-data.frame(gene=genes,
                       zdesc=zf_genes$description[match(genes,zf_genes$external_gene_name)],
                       cdesc=ch_genes$description[match(genes,ch_genes$external_gene_name)])
gene_data$cdesc[gene_data$cdesc==""] <- NA
length(grep("LOC",gene_data$gene)) # Not really annotated
length(which(!is.na(gene_data$cdesc) | !is.na(gene_data$zdesc))) # Annotated with some gene name

gos <- data.frame(gene_id=genes,desc=NA,GO=NA)
for (g in genes) {
  if(!is.na(gene_data$zdesc[gene_data$gene==g])) {
    sub <- subset(zf_genes,external_gene_name==g)
    gostring <- paste(sub$go_id,collapse=",")
    gos$desc[gos$gene_id==g] <- sub$description[1]
  }
  else if (!is.na(gene_data$cdesc[gene_data$gene==g])) {
    sub <- subset(ch_genes,external_gene_name==g)
    gostring <- paste(sub$go_id,collapse=",")   
    gos$desc[gos$gene_id==g] <- sub$description[1]
  }
  else {gostring <- NA}
  gos$GO[gos$gene_id==g] <- gostring
}
write.table(gos[,c("gene_id","GO")],"GOmap.txt",quote=F,row.names=F,col.names=F,sep="\t")


###################################
parname = "FemaleBill"
file="../cands/cutoff5/HETH.contemp.female.tarsPCcovar.lrt1.sites"
sites <- read.delim(file,header=F)
dim(sites)

####Pull out genes linked to significant snps
candSNP <- sites
candGenes <- unique(ld$gene[ld$snp%in%candSNP[,1]])
candGeneDesc <- sapply(candGenes,function(x) gos$desc[gos$gene_id==x])
candGeneDescShort <- gsub(" \\[.*\\]","",candGeneDesc)
geneFrame <- data.frame(GeneID=candGenes,Description=candGeneDescShort)
write.csv(geneFrame,paste("output/genelists/",parname,".genes.csv",sep=""),row.names=F)

###For a given column, do GO enrichment for significant SNPs
geneID2GO <- readMappings(file="GOmap.txt",sep="\t",IDsep=",")
allgenes <- unique(ld$gene) #This is the 'background'
myIG <- factor(as.numeric(allgenes%in%candGenes))
names(myIG) <- allgenes


##BP
GOdata <- new("topGOdata",ontology="BP",allGenes=myIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>=5,]
write.csv(filt,paste("output/",parname,".GO_BP.csv",sep=""))

##MF
GOdata <- new("topGOdata",ontology="MF",allGenes=myIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>=5,]
write.csv(filt,paste("output/",parname,".GO_MF.csv",sep=""))

##CC
GOdata <- new("topGOdata",ontology="CC",allGenes=myIG, nodeSize=10,
              annotationFun=annFUN.gene2GO,gene2GO=geneID2GO)
test.stat <- new("classicCount",testStatistic=GOFisherTest,name="Fisher test",cutoff=0.01)
resultFisher <- getSigGroups(GOdata,test.stat)
res <- GenTable(GOdata,classic=resultFisher,topNodes=length(resultFisher@score),numChar=100)
filt <- res[res$classic<0.05 & res$Significant>=5,]
write.csv(filt,paste("output/",parname,".GO_CC.csv",sep=""))
