##########
# By: Ethan Gyllenhaal
# Updated 17Feb2026
#
# R script used for plotting EEMS results for greenbul dataset.
# Uses reemsplots2 to generate plots, then adds maps and sampling points under them.
# Take 4 chains of EEMS output and a csv mapping coordinates to populations
########

# install for reemsplots2
#library("devtools")
#install_github("dipetkov/reemsplots2")

# load libraries
library(reemsplots2)
# plotting libraries
library(rworldmap)
library(rworldxtra)
library(ggplot2)
library(broom)
library(gridExtra)
library(ggnewscale)
library(sf)
library(tidyverse)

setwd("C:/Documents/Projects/Greenbul/eems")

# make vector of directory names for eems outputs
mcmcpath <- c("greenbul49_noZ_90_10k_chain1", "greenbul49_noZ_90_10k_chain2", "greenbul49_noZ_90_10k_chain3", "greenbul49_noZ_90_10k_chain4", "greenbul49_noZ_90_10k_chain5", "greenbul49_noZ_90_10k_chain6", "greenbul49_noZ_90_10k_chain7", "greenbul49_noZ_90_10k_chain8", "greenbul49_noZ_90_10k_chain9", "greenbul49_noZ_90_10k_chain10", "greenbul49_noZ_90_10k_chain11", "greenbul49_noZ_90_10k_chain12", "greenbul49_noZ_90_10k_chain13", "greenbul49_noZ_90_10k_chain14", "greenbul49_noZ_90_10k_chain15", "greenbul49_noZ_90_10k_chain16", "greenbul49_noZ_90_10k_chain17", "greenbul49_noZ_90_10k_chain18", "greenbul49_noZ_90_10k_chain19", "greenbul49_noZ_90_10k_chain20")
#mcmcpath <- c("greenbul49_noZ_90_10k_chain4") # check indiv chains
# run the main plotting function, look at a couple
plots <- make_eems_plots(mcmcpath, longlat = TRUE)
# checking out plots
plots$qrates01 # diversity surface
plots$mrates01 # migration surface
plots$pilogl01 # MCMC trace, seems like it was level after burnin

# load in coordinates
coords <- read.csv("pop_coords.csv") %>%
  mutate(mismatch = if_else(pop==mtDNA, 0, 1)) %>%
  rowwise() %>%
  mutate(lon_jitter=lon+runif(1,-0.8,0.8)) %>%
  mutate(lat_jitter=lat+runif(1,-0.8,0.8))

# download map, use broom's tidy function to make it neat, per documentation of eems plotting package
#map <- getMap(resolution="high")
map <- map_data("world")

# make migration rate plot
mrate <- plots$mrates01 +
  geom_path(data = map, aes(x = long, y = lat, group = group),
            color = "#888888", linewidth = 0.5) +
  coord_quickmap() +
  new_scale_fill() +
  geom_point(data = filter(coords, mismatch==0), aes(x=lon_jitter, y=lat_jitter, fill=pop), size=5, shape=21) +
  scale_fill_manual(values=c("#E6480B", "#E6D839","#FE9E00", "#4688FB")) +
  geom_point(data = filter(coords, mismatch==1), aes(x=lon_jitter, y=lat_jitter, fill=pop), size=5, shape=21) +
  geom_point(data = filter(coords, mismatch==1), aes(x=lon_jitter, y=lat_jitter, fill=mtDNA), size=3, shape=21)
mrate



# make diversity surface plot
qrate <- plots$qrates01 +
  geom_path(data = map, aes(x = long, y = lat, group = group),
            color = "#888888", size = 0.5) +
  coord_quickmap() +
  new_scale_fill() +
  geom_point(data = filter(coords, mismatch==0), aes(x=lon_jitter, y=lat_jitter, fill=pop), size=5, shape=21) +
  scale_fill_manual(values=c("#E6480B", "#E6D839","#FE9E00", "#4688FB")) +
  geom_point(data = filter(coords, mismatch==1), aes(x=lon_jitter, y=lat_jitter, fill=pop), size=5, shape=21) +
  geom_point(data = filter(coords, mismatch==1), aes(x=lon_jitter, y=lat_jitter, fill=mtDNA), size=3, shape=21)
qrate
# plot together
grid.arrange(mrate, qrate, nrow=2)

