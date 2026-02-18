# R script for plotting fine-scale PCAs, first one described in depth. Also explores the east results in depth.

# load envs
library("adegenet")
library("ade4")
library("vcfR")
library("scales")
library("parallel")
library("viridis")
library("wesanderson")
library("tidyverse")
library("fuzzyjoin")
setwd('C:/Documents/Projects/Greenbul/finescale')

#formula for approximate p-values from "cell" in PCA
# p = 1 - exp( -[ c[cell]^2 /2 )
# 1.5 = 67%
# 2.5 = 95%

# popmaps, assumes popmap is sorted
map_east <- read.csv('east_pops.csv', sep = ',')
map_west <- read.csv('west_pops.csv', sep = ',')
map_central <- read.csv('central_pops.csv', sep = ',')

# east ####
# colors
colors_east <- c("coral4", "coral3", "chocolate2", "chocolate1", "gold4", "gold3", "gold2", "gold1")
east_vcf <- read.vcfR("east22_noZ_90_10k.recode.vcf") # load VCF
east_gl <- vcfR2genlight(east_vcf) # convert to genlight
pop(east_gl) <- map_east$Site # popmap
pca_east <- glPca(east_gl, n.cores=4, nf=4) # run the PCA
scatter(pca_east, cex=.25) # scatter plot to troubleshoot
# labeled PCA
s.class(pca_east$scores[,1:2], pop(east_gl), col=colors_east, clab=1, cell=2, cstar=0)
# unlabeled
s.class(pca_east$scores[,1:2], pop(east_gl), col=colors_east, clab=0, cell=2)
# PC 1 and 3
s.class(pca_east$scores[,c(1,3)], pop(east_gl), col=colors_east, clab=1, cell=2, cstar=0)
# PC 1 and 4
s.class(pca_east$scores[,c(1,4)], pop(east_gl), col=colors_east, clab=1, cell=2)
# barplot of loadings
barplot(pca_east$eig/sum(pca_east$eig), main="eigenvalues", col=heat.colors(length(pca_east$eig)))
# write out east scores, not needed
east_scores <- cbind(pca_east$scores[,1:4], east_gl$pop) %>% as.data.frame()
east_scores$names <- rownames(east_scores)
write.csv(east_scores, file="east_scores.csv", quote = F)


# central ####
central_vcf <- read.vcfR("central22_noZ_90_10k.recode.vcf")
central_gl <- vcfR2genlight(central_vcf)
pop(central_gl) <- map_central$Site
pca_central <- glPca(central_gl, n.cores=4, nf=4)
scatter(pca_central, cex=.25)
s.class(pca_central$scores[,1:2], pop(central_gl), col=magma(8,begin=.1,end=.9), clab=1, cell=2)
s.class(pca_central$scores[,c(1,3)], pop(central_gl), col=magma(8,begin=.1,end=.9), clab=1, cell=2, cstar=0)
s.class(pca_central$scores[,c(1,4)], pop(central_gl), col=magma(8,begin=.1,end=.9), clab=1, cell=2)
barplot(pca_central$eig/sum(pca_central$eig), main="eigenvalues", col=heat.colors(length(pca_central$eig)))

# west ####
west_vcf <- read.vcfR("west5_noZ_90_10k.recode.vcf")
west_gl <- vcfR2genlight(west_vcf)
pop(west_gl) <- map_west$Site
pca_west <- glPca(west_gl, n.cores=4, nf=4)
scatter(pca_west, cex=.25)
s.class(pca_west$scores[,1:2], pop(west_gl), col=magma(8,begin=.1,end=.9), clab=1, cell=2)
s.class(pca_west$scores[,c(1,3)], pop(west_gl), col=magma(8,begin=.1,end=.9), clab=1, cell=2, cstar=0)
s.class(pca_west$scores[,c(1,4)], pop(west_gl), col=magma(8,begin=.1,end=.9), clab=1, cell=2)
barplot(pca_west$eig/sum(pca_west$eig), main="eigenvalues", col=heat.colors(length(pca_west$eig)))

# east in depth stuff ####
# read in scores
east_scores <- read.csv("east_scores.csv", header=T)
# join with popmap
combined <- inner_join(east_scores, map_east, by = join_by(names==Sample))

# map libraries
library(elevatr)
library(rgeoboundaries)
library(sf)
library(raster)
library(ggplot2)
library(viridis)
library(ggnewscale)

# read in coordinates
coords <- read.csv("pop_coords.csv") %>% filter(!is.na(lon))
# convert to shape file
coord_sf <- st_as_sf(coords, coords=c("lon","lat"))
st_crs(coord_sf) <- 4326
# set bounds
x<-c(24.8,34); y<-c(-5,4.8)
bounds <- data.frame(x,y)
bounds <- st_as_sf(bounds, coords=c("x", "y"))
st_crs(bounds) <- 4326
# base map
map <- map_data("world")
# download elevation raster, convert to df
elevation <- get_elev_raster(locations = bounds, z = 6, clip = "locations")
elev_df <- as.data.frame(elevation, xy=TRUE)
colnames(elev_df)[3] <- "elevation"
elev_df <- elev_df[complete.cases(elev_df), ]

# map colors
colors_map <- c("coral4", "coral3", "chocolate2", "chocolate1", "gold4", "gold3", "gold2", "gold1", "red1")

# plot it all together
ggplot() +
  geom_raster(data = elev_df, aes(x = x, y = y, fill = elevation)) +
  scale_fill_gradient2(low = "white", mid = "gray20", high="black", midpoint=2500) +
  geom_path(data = map, aes(x = long, y = lat, group = group),
            color = "black", linewidth = 0.7) +
  new_scale_fill() +
  geom_point(data = filter(coords), aes(x=lon, y=lat, fill=as.factor(Site)),color="black", size=4, shape=21) +
  scale_fill_manual(values=colors_map) +
  geom_sf(data = bounds, color = "white", fill = NA) +
  coord_sf() +
  xlim(c(x)) + ylim(c(y)) +
  theme_void()

# unzoomed
coords <- read.csv("pop_coords.csv") %>% filter(!is.na(lon))
coord_sf <- st_as_sf(coords, coords=c("lon","lat"))
st_crs(coord_sf) <- 4326
x<-c(-16.7,52); y<-c(-12,12)
bounds <- data.frame(x,y)
bounds <- st_as_sf(bounds, coords=c("x", "y"))
map <- map_data("world")
st_crs(bounds) <- 4326
elevation <- get_elev_raster(locations = bounds, z = 3, clip = "locations")
elev_df <- as.data.frame(elevation, xy=TRUE)
colnames(elev_df)[3] <- "elevation"
elev_df <- elev_df %>%
  mutate(elevation_floor = case_when(elevation>0 ~ elevation,
                                     .default = 0))
elev_df <- elev_df[complete.cases(elev_df), ]


ggplot() +
  geom_raster(data = elev_df, aes(x = x, y = y, fill = elevation_floor)) +
  scale_fill_gradient2(low = "white", mid = "gray20", high="black", midpoint=2500) +
  geom_path(data = map, aes(x = long, y = lat, group = group),
            color = "black", linewidth = 0.7) +
  new_scale_fill() +
  geom_sf(data = bounds, color = "white", fill = NA) +
  coord_sf() +
  xlim(c(x)) + ylim(c(y)) +
  theme_void()

