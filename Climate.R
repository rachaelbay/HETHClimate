setwd("~/Documents/IndividualProjects/Nicole/HETH/Climate/")

library(ebirdst)
library(tidyverse)
library(sf)
library(sfheaders)
library(terra)
library(daymetr)
library(rnaturalearth)
library(patchwork)
library(viridisLite)
library(ggstar)

sf_use_s2(FALSE)
#####Just get full species ranges - OLD - now using genoscape files
##First get the range boundaries from eBird
#set_ebirdst_access_key("c2lgknbh7u5l") #Works until July 2024

# ebirdst_runs %>% filter(species_code=="herthr") %>% glimpse()
# #ebirdst_download_status(species="herthr",download_all=T)
# 
# ranges <- load_ranges("herthr",resolution="27km",path = ebirdst_data_dir())
# br_full <- filter(ranges,season=="breeding")
# wi_full <- filter(ranges,season=="nonbreeding")
# plot(br_full)
# plot(wi_full)


###Using genoscape shapefile
breeding <- st_read("../shapefile/HETH.genoscape_brick.WGS84_3.5no_ovlp2.shp") #This is the eastern taiga population only
br_full <- st_read("../shapefile/HETHbreedingrange.shp")
#plot(br_full)
#br_union <- st_union(br_full)
#br_buffered <- st_buffer(br_union, dist = 0.1)
#br_simplified <- st_simplify(br_buffered, dTolerance = 0.2)
#br_border <- st_boundary(br_simplified)
#saveRDS(br_border,"../shapefile/BreedingBorder.rds")
br_border <- readRDS("../shapefile/BreedingBorder.rds")

# Plot the external border
ggplot() +
  geom_sf(data = br_border) +
  theme_minimal()



################
####Breeding range
##################

##Make the grid (1 degree)
# x.inc=1
# y.inc=0.5
# g.br=expand.grid(
#   x=seq(st_bbox(breeding)$xmin,st_bbox(breeding)$xmax,by=x.inc),
#   y=seq(st_bbox(breeding)$ymin,st_bbox(breeding)$ymax,by=y.inc)
# )
# g.br_sfc <- sfc_point(as.matrix(g.br)) %>%
#   st_set_crs(st_crs(breeding))
# 
# ##Trim the grid
# br <- g.br_sfc[unlist(st_contains(breeding,g.br_sfc))]
# plot(br,pch=19,cex=0.5,col="red")
# br_coords <- do.call(rbind, st_geometry(br)) %>%
#   as_tibble() %>% setNames(c("lon","lat"))
# br.frame <- data.frame(site=paste("br",1:nrow(br_coords),sep=""),
#                        latitude=br_coords$lat,
#                        longitude=br_coords$lon)
# write.csv(br.frame,"BreedingPoints.csv",row.names=F)
# 
# ##Pull climate data
# breed <- download_daymet_batch(file_location = 'BreedingPoints.csv',
#                                start = 1981,
#                                end = 2016,
#                                internal = T,
#                                simplify=F)
# saveRDS(breed,"BreedClimate.rds")
breed <- readRDS("BreedClimate.rds")
br.frame <- read.csv("BreedingPoints.csv")
b.cor.frame <- data.frame(br.frame,tmin=NA,tmax=NA,precip=NA)

##correlations
for (i in 1:length(breed)) {
  if (length(breed[[i]])==1) {
    b.cor.frame$tmin[i] <- NA
    b.cor.frame$tmax[i] <- NA
    b.cor.frame$precip[i] <- NA
  }
  else {
    breed.tr <- breed[[i]]$data
    ## take only June and July (days 152-212)
    breed.jj <- breed.tr %>% filter(yday>=152 & yday<=212)
    breed.agg <- breed.jj %>% group_by(year) %>%
      summarise(precip=sum(prcp..mm.day.),
                tmax=quantile(tmax..deg.c.,0.95),
                tmin=quantile(tmin..deg.c.,0.05))
    b.cor.frame$tmin[i] <- cor(breed.agg$year,breed.agg$tmin,method="spearman")
    b.cor.frame$tmax[i] <- cor(breed.agg$year,breed.agg$tmax,method="spearman")
    b.cor.frame$precip[i] <- cor(breed.agg$year,breed.agg$precip,method="spearman")
  }
}

##Plot
ggplot(na.omit(b.cor.frame), aes(x=longitude,y=latitude,col=tmax)) + geom_point(size=1) +
  theme_classic() + scale_color_gradient2(midpoint=0)



################
####Wintering range
##################
wi_shape <- st_read("../shapefile/HETH.winterEcoregions.simple.sf.WGS84.eBird.shp")
wi_boundary <- st_boundary(st_union(wi_shape))
wintering <- wi_shape

##Make the grid (1 degree)
# increment=1
# g.wi=expand.grid(
#   x=seq(st_bbox(wintering)$xmin,st_bbox(wintering)$xmax,by=increment),
#   y=seq(st_bbox(wintering)$ymin,st_bbox(wintering)$ymax,by=increment/2)
# )
# g.wi_sfc <- sfc_point(as.matrix(g.wi)) %>%
#   st_set_crs(st_crs(wintering))
# 
# ##Trim the grid
# wi <- g.wi_sfc[unlist(st_contains(wintering,g.wi_sfc))]
# plot(wi,pch=19,cex=0.5,col="red")
# wi_coords <- do.call(rbind, st_geometry(wi)) %>%
#   as_tibble() %>% setNames(c("lon","lat"))
# wi.frame <- data.frame(site=paste("wi",1:nrow(wi_coords),sep=""),
#                        latitude=wi_coords$lat,
#                        longitude=wi_coords$lon)
# write.csv(wi.frame,"WinteringPoints.csv",row.names=F)
# 
# ##Pull climate data
# winter <- download_daymet_batch(file_location = 'WinteringPoints.csv',
#                                 start = 1981,
#                                 end = 2016,
#                                 internal = T,
#                                 simplify=F)
# saveRDS(winter,"WintClimate.rds")
winter <- readRDS("WintClimate.rds")
wi.frame <- read.csv("WinteringPoints.csv")
w.cor.frame <- data.frame(wi.frame,tmin=NA,tmax=NA,precip=NA)

##correlations
for (i in 1:length(winter)) {
  if (length(winter[[i]])==1) {
    w.cor.frame$tmin[i] <- NA
    w.cor.frame$tmax[i] <- NA
    w.cor.frame$precip[i] <- NA
  }
  else {
    wint.tr <- winter[[i]]$data
    ## take only January February (days 1-60)
    wint.jj <- wint.tr %>% filter(yday>=1 & yday<=60)
    wint.agg <- wint.jj %>% group_by(year) %>%
      summarise(precip=sum(prcp..mm.day.),
                tmax=quantile(tmax..deg.c.,0.95),
                tmin=quantile(tmin..deg.c.,0.05))
    w.cor.frame$tmin[i] <- cor(wint.agg$year,wint.agg$tmin,method="spearman")
    w.cor.frame$tmax[i] <- cor(wint.agg$year,wint.agg$tmax,method="spearman")
    w.cor.frame$precip[i] <- cor(wint.agg$year,wint.agg$precip,method="spearman")
  }
}

##Plot
ggplot(na.omit(w.cor.frame), aes(x=longitude,y=latitude,col=tmax)) + geom_point(size=2) +
  theme_classic() + scale_color_gradient2(midpoint=0)


#######################
### Plot
#######################
br_box <- st_bbox(br_full)
wi_box <- st_bbox(wi_shape)
xmin <- min(c(br_box[1],wi_box[1]))-2
xmax <- max(c(br_box[3],wi_box[3])) +2
ymin <- min(c(br_box[2],wi_box[2]))-2
ymax <- max(c(br_box[4],wi_box[4]))+2
sf_use_s2(FALSE)

world <- ne_countries(scale = "medium", returnclass = "sf")
world_crop <- st_crop(world, xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax)
plot(world_crop)


cols <- viridis(10,option="C")[c(3,7)]

map1 <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1')+
  geom_sf(data=br_border,lwd=0.1) +
  geom_sf(data=br_simplified,fill="grey95") +
  geom_point(data = na.omit(b.cor.frame), aes(x=longitude, y=latitude,color=tmin), size=0.05, shape=15)+
  scale_color_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  geom_sf(data = world_crop, fill = NA) +
  theme(legend.position = 'right',axis.text=element_text(size=6)) +
  scale_x_continuous(labels=~.x) + 
  scale_y_continuous(labels=~.x)

hist1 <- ggplot(data=na.omit(b.cor.frame),aes(x=tmin)) + geom_histogram(aes(fill=..x..),col="grey30") +
  xlim(-1,1) + theme_classic() +   scale_fill_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  theme(legend.position="none",axis.text=element_text(size=6))

map2 <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1')+
  geom_sf(data=br_border,lwd=0.1) +
  geom_sf(data=br_simplified,fill="grey95") +
  geom_point(data = na.omit(b.cor.frame), aes(x=longitude, y=latitude,color=tmax), size=0.05, shape=15)+
  scale_color_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  geom_sf(data = world_crop, fill = NA) +
  theme(legend.position = 'right',axis.text=element_text(size=6)) +
  scale_x_continuous(labels=~.x) + 
  scale_y_continuous(labels=~.x)

hist2 <- ggplot(data=na.omit(b.cor.frame),aes(x=tmax)) + geom_histogram(aes(fill=..x..),col="grey30") +
    xlim(-1,1) + theme_classic() +   scale_fill_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  theme(legend.position="none",axis.text=element_text(size=6))

map3 <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1')+
  geom_sf(data=br_border,lwd=0.1) +
  geom_sf(data=br_simplified,fill="grey95") +
  geom_point(data = na.omit(b.cor.frame), aes(x=longitude, y=latitude,color=precip), size=0.05, shape=15)+
  scale_color_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  geom_sf(data = world_crop, fill = NA) +
  theme(legend.position = 'right',axis.text=element_text(size=6)) +
  scale_x_continuous(labels=~.x) + 
  scale_y_continuous(labels=~.x)

hist3 <- ggplot(data=na.omit(b.cor.frame),aes(x=precip)) + geom_histogram(aes(fill=..x..),col="grey30") +
    xlim(-1,1) + theme_classic() +   scale_fill_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  theme(legend.position="none",axis.text=element_text(size=6))
  
map4 <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1')+
  geom_sf(data=wi_shape,fill='grey95',lwd=0) +
  geom_sf(data=wi_boundary,lwd=0.1) +
  geom_point(data = na.omit(w.cor.frame), aes(x=longitude, y=latitude,color=tmin), size=0.05, shape=15)+
  scale_color_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  geom_sf(data = world_crop, fill = NA) +
  theme(legend.position = 'right',axis.text=element_text(size=6)) +
  scale_x_continuous(labels=~.x) + 
  scale_y_continuous(labels=~.x)

hist4 <- ggplot(data=na.omit(w.cor.frame),aes(x=tmin)) + geom_histogram(aes(fill=..x..),col="grey30") +
    xlim(-1,1) + theme_classic() +   scale_fill_gradient2(midpoint=0,low=cols[1],high=cols[2])+
  theme(legend.position="none",axis.text=element_text(size=6))
  
map5 <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1')+
  geom_sf(data=wi_shape,fill='grey95',lwd=0) +
  geom_sf(data=wi_boundary,lwd=0.1) +
  geom_point(data = na.omit(w.cor.frame), aes(x=longitude, y=latitude,color=tmax), size=0.05, shape=15)+
  scale_color_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  geom_sf(data = world_crop, fill = NA) +
  theme(legend.position = 'right',axis.text=element_text(size=6)) +
  scale_x_continuous(labels=~.x) + 
  scale_y_continuous(labels=~.x)

hist5 <- ggplot(data=na.omit(w.cor.frame),aes(x=tmax)) + geom_histogram(aes(fill=..x..),col="grey30") +
    xlim(-1,1) + theme_classic() +   scale_fill_gradient2(midpoint=0,low=cols[1],high=cols[2])  +
  theme(legend.position="none",axis.text=element_text(size=6))
  
map6 <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1')+
  geom_sf(data=wi_shape,fill='grey95',lwd=0) +
  geom_sf(data=wi_boundary,lwd=0.1) +
  geom_point(data = na.omit(w.cor.frame), aes(x=longitude, y=latitude,color=precip), size=0.05, shape=15)+
  scale_color_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  geom_sf(data = world_crop, fill = NA) +
  theme(legend.position = 'right',axis.text=element_text(size=6)) +
  scale_x_continuous(labels=~.x) + 
  scale_y_continuous(labels=~.x)

hist6 <- ggplot(data=na.omit(w.cor.frame),aes(x=precip)) + geom_histogram(aes(fill=..x..),col="grey30") +
    xlim(-1,1) + theme_classic() +   scale_fill_gradient2(midpoint=0,low=cols[1],high=cols[2])+
  theme(legend.position="none",axis.text=element_text(size=6))

pdf(file="../Figures/New_07.24/AllClimate.pdf",width=8.5,height=7)
#jpeg(filename="../Figures/New_07.24/AllClimate.jpeg",width=700,height=800,units="px")
map1 + map2 + map3 + hist1 + hist2 + hist3 + map4 + map5 + map6 + hist4 + hist5 + hist6 +
  plot_annotation(tag_levels='A') +
  theme(plot.tag=element_text(hjust=0,vjust=0)) +
  plot_layout(guides="collect",nrow=4,heights=c(3,1,3,1)) &
  scale_color_gradient2(limits=c(-1,1),name="Trend (rho)",low=cols[1],high=cols[2])
dev.off()



#######################
## Just tmin in one plot for main figure
#####################

b <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1')+
  geom_sf(data=br_simplified,fill="grey95") +
  geom_point(data = na.omit(b.cor.frame), aes(x=longitude, y=latitude,color=tmin), size=0.2, shape=15)+
  scale_color_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  geom_sf(data = world_crop, fill = NA) +
  geom_sf(data=br_border,lwd=0.01) +
  geom_star(aes(x=-87.62,y=41.87),size=3, fill="darkgoldenrod1",starstroke=0.3) +
  theme(legend.position = 'right',axis.text=element_text(size=6)) +
  scale_x_continuous(labels=~.x) + 
  scale_y_continuous(labels=~.x) + ggtitle("Breeding")+
  theme(plot.title = element_text(hjust = 0.5))
# 
w <- ggplot(data = world) +
  theme_bw()+
  geom_sf(data = world_crop, fill = 'antiquewhite1')+
  geom_sf(data=wi_shape,fill='grey95',lwd=0) +
  geom_point(data = na.omit(w.cor.frame), aes(x=longitude, y=latitude,color=tmin), size=0.2, shape=15)+
  scale_color_gradient2(midpoint=0,low=cols[1],high=cols[2]) +
  geom_sf(data = world_crop, fill = NA) +
  geom_sf(data=wi_boundary,lwd=0.01) +
  geom_star(aes(x=-87.62,y=41.87),size=3, fill="darkgoldenrod1",starstroke=0.3) +
  theme(legend.position = 'right',axis.text=element_text(size=6)) +
  scale_x_continuous(labels=~.x) + 
  scale_y_continuous(labels=~.x) + ggtitle("Wintering")+
  theme(plot.title = element_text(hjust = 0.5))

pdf(file="../Figures/New_07.24/Tmin_maps.pdf",width=7,height=4)
b + w + plot_layout(guides="collect",nrow=1) &
  scale_color_gradient2(limits=c(-1,1),name="Tmin Trend",low=cols[1],high=cols[2])
dev.off()
