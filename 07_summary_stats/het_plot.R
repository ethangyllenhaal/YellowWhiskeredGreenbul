# load packages
library(tidyverse)
library(viridis)
library(boot)
library(ggbeeswarm)
library(stringr)

# set working directory
setwd("C:/Documents/Projects/Greenbul/winstats")

# read in focal Tajima's D data
taj <- read.csv("window_noZ100k_fine_Tajima.txt", sep="\t") %>%
  filter(pop1 != "Outgroup") %>%
  mutate(weight = calculated_stat*number_sites)

# function for making bootstrapped estimates of weighted mean genome-wide Tajima's D
taj_bootstrap <- function(data, indicies, pop){
  boot_data <- data[indicies,]
  return(sum(filter(boot_data, pop1==pop)$weight, na.rm=T) /
    sum(filter(boot_data, pop1==pop)$number_sites, na.rm=T))
}

# make dataframe for bootstraps
bootstrap_df <- data.frame(WeightedTajD=double(), Population=character())
# loop through populations and get bootstraps
for(pop in c("West", "Central", "C_discordant", "East", "Rwenzori", "Northeast")){
  boot_temp <- boot(data=taj, statistic = taj_bootstrap, pop=pop, R=50)
  df <- data.frame(boot_temp$t) %>%
    mutate(p=pop)
  colnames(df) <- c("WeightedTajD", "Population")
  bootstrap_df <- rbind(bootstrap_df, df)
}

# point plot of bootstrapped estimates (almost no variance)
ggplot(data=bootstrap_df, aes(x=Population, fill=Population, y=WeightedTajD))+
  stat_summary(fun = "mean", position = position_dodge(0.9), shape=23, size=2) +
  ylim(c(-2,0)) +
  scale_x_discrete(limits=c("West", "Central", "C_discordant", "East", "Rwenzori", "Northeast")) +
  scale_fill_manual(values=c("#E6480B", "#E6480B", "#E6D839", "#E6D839", "#E6D839", "#4688FB")) +
  theme_bw()

# violin plot for windowed stat
ggplot(data=taj, aes(x=pop1, fill=pop1, y=calculated_stat))+
  geom_violin() +
  scale_x_discrete(limits=c("West", "Central", "C_discordant", "East", "Rwenzori", "Northeast")) +
  stat_summary(fun = "mean", position = position_dodge(0.9), shape=16) + 
  scale_fill_manual(values=c("#E6480B", "#E6480B", "#E6D839", "#E6D839", "#E6D839", "#4688FB")) +
  theme_bw()

# heterozygostiy section ####
# populations
west <- c("E_latirostris__SLE_KambuiHills__KU19685", "E_latirostris__SLE_KambuiHills__KU19697", "E_latirostris__SLE_KambuiHills__KU19776", "E_latirostris__GHA_AssinFoso__FMNH396575", "E_latirostris__GHA_AssinFoso__FMNH396576")
central <- c("E_latirostris__GNQ_Moka__KU32866", "E_latirostris__GNQ_Moka__KU32961", "E_latirostris__GNQ_Moka__KU32908", "E_latirostris__GNQ_MonteAlen__KU8325", "E_latirostris__GNQ_MonteAlen__KU8326", "E_latirostris__GNQ_MonteAlen__KU8332", "E_latirostris__CMR_Nzimdipnenkum__FMNH511209", "E_latirostris__CMR_Nzimdipnenkum__FMNH511210", "E_latirostris__GAB_Doussala__FMNH396331", "E_latirostris__GAB_Doussala__FMNH396330", "E_latirostris__GAB_Doussala__FMNH389319", "E_latirostris__COD_Luki__FMNH490024", "E_latirostris__COD_Luki__FMNH490025", "E_latirostris__COD_Luki__FMNH490026", "E_latirostris__CAF_Dzanga__FMNH429448", "E_latirostris__CAF_Dzanga__FMNH429445", "E_latirostris__CAF_Dzanga__FMNH429442", "E_latirostris__COD_Boende__FMNH490182", "E_latirostris__COD_Boende__FMNH490193", "E_latirostris__COD_Boende__FMNH490213", "E_latirostris__COD_Tschambi__FMNH473401", "E_latirostris__COD_Tschambi__FMNH473402")
east <- c("E_latirostris__COD_KahuziBiega__FMNH481070", "E_latirostris__COD_KahuziBiega__FMNH481071", "E_latirostris__COD_KahuziBiega__FMNH481072", "E_latirostris__COD_Kangala__FMNH481083", "E_latirostris__COD_Kangala__FMNH481084", "E_latirostris__COD_Kangala__FMNH481085", "E_latirostris__COD_Idjwi__FMNH429729", "E_latirostris__BDI_Kibira__FMNH357966", "E_latirostris__BDI_Kibira__FMNH346297", "E_latirostris__UGA_Nteko__FMNH384909", "E_latirostris__UGA_Nteko__FMNH384910", "E_latirostris__UGA_Nteko__FMNH384911", "E_latirostris__UGA_Echuya__FMNH384894", "E_latirostris__UGA_Echuya__FMNH384895", "E_latirostris__UGA_Echuya__FMNH384896", "E_latirostris__UGA_Rwenzori__FMNH355307", "E_latirostris__UGA_Rwenzori__FMNH355308", "E_latirostris__UGA_Rwenzori__FMNH355309", "E_latirostris__UGA_Rwenzori__FMNH355310", "E_latirostris__UGA_Rwenzori__FMNH355311")
northeast <- c("E_latirostris__UGA_AgoroAgu__FMNH489311", "E_latirostris__UGA_AgoroAgu__FMNH489312")

# make big summary dataframe with weighted means
full_het_summarize <- read.csv("window_noZ100k_heterozygosity.txt", sep='\t') %>%
  mutate(weight = calculated_stat*number_sites)%>%
  group_by(pop1) %>%
  summarise(sum_weight = sum(weight, na.rm=T), sum_sites = sum(number_sites, na.rm=T)) %>%
  mutate(mean_het = sum_weight/sum_sites)  %>%
  mutate(cluster = case_when(pop1 %in% west ~ "1west",
                             pop1 %in% central ~ "2central",
                             pop1 %in% northeast ~ "4northeast",
                             .default = "3east"))

# violin plot of heterozygosity values
ggplot(data=full_het_summarize, aes(x=cluster, y=mean_het, fill=cluster)) +
  geom_violin() + geom_beeswarm() +
  scale_fill_manual(values=c("#4688FB", "#E6480B", "#E6D839","#FE9E00")) +
  theme_bw()

# ANOVA and post-hoc comparisons
summary(aov(data=full_het_summarize, mean_het ~ cluster))
TukeyHSD(aov(data=full_het_summarize, mean_het ~ cluster))

# means for each cluster
mean(filter(full_het_summarize, cluster=="2central")$mean_het)
mean(filter(full_het_summarize, cluster=="3east")$mean_het)


