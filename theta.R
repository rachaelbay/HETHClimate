setwd("~/Documents/IndividualProjects/Nicole/HETH/theta/")

library(tidyverse)
library(patchwork)
library(shades)

###################
## Make sample files for each year > 5 individuals
#####################

##Read in individuals and get file names
# inds <- read.delim("./Data/all.219.inds.txt",header=F)
# head(inds)
# inds$IDs <- substr(inds$V1,2,7)
# files <- paste(inds[,1],".ml",sep="")
# 
# #Get metadata
# meta <- read.csv("./HETH_meta.csv")
# ordermeta <- meta[match(inds$IDs,meta$Sample),]
# ordermeta$Zname <- inds[,1]
# table(ordermeta$Year)
# write.table(ordermeta$Zname[ordermeta$Year==2014],
#             file="theta/samples.2014.txt",
#             row.names=F,col.names=F,quote=F)


##################

file.list <- list.files(pattern="*.pestPG")

# Create a function to read a file and add a column with the file name
read_and_add_filename <- function(file) {
  data <- read.delim(file, header = TRUE)  # Adjust read function based on your file format
  data <- data %>% mutate(filename = file)
  return(data)
}

# Use map_dfr to apply the function to each file and concatenate the results
combined_data <- map_dfr(file.list, read_and_add_filename)
combined_data <- separate(combined_data,col="filename",into=c("Year","thetas","pestPG"),sep="\\.")
combined_data$Year <- as.numeric(combined_data$Year)
thet <- ggplot(na.omit(combined_data),aes(x=Year,y=tW/nSites)) + geom_jitter(alpha=0.2,col="gray") + 
  theme_classic() + stat_summary(fun.data="mean_cl_normal",geom="crossbar") + ylab("Theta")
div <- ggplot(na.omit(combined_data),aes(x=Year,y=tP/nSites)) + geom_jitter(alpha=0.2,col="gray") + 
  theme_classic() + stat_summary(fun.data="mean_cl_normal",geom="crossbar") + ylab("Nucleotide Diversity")


saveRDS(thet,"thetaplot.rds")
saveRDS(div,"divplot.rds")

lm <- lm(tP/nSites~Year,data=na.omit(combined_data))
hist(residuals(lm))
summary(lm)


###
# Genetic diversity plot
het <- readRDS("../het/hetplot.rds")

jpeg("../Figures/Diversity.jpg", width=3000,height=800, units="px",res=250)
#pdf("../Figures/Diversity.pdf",width=7,height=3)
het + thet + div+plot_annotation(tag_levels='A')
dev.off()
