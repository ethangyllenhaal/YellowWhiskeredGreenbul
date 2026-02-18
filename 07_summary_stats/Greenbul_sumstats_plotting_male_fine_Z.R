options(scipen=999)
library(Hmisc)
setwd("C:/Documents/Projects/Greenbul/winstats")
dxy <- read.table("window_male100k_fine_dxy.txt", header=T)
fst <- read.table("window_male100k_fine_fst.txt", header=T)
pi <- read.table("window_male100k_fine_pi.txt", header=T)
pi$calculated_stat[is.na(pi$calculated_stat)] <- 0
tajima <- read.table("window_male100k_fine_tajima.txt", header=T)
# combine all
x <- rbind(dxy, fst, pi, tajima) %>%
  filter(chr=="NC_089174.1_RagTag")

# determine order of windows
scaffolds <- c("NC_089174.1_RagTag")
scaffold_names <- c("Z")
window_order <- c()
for(a in 1:length(scaffolds)) {
  a_rep <- x[x$chr == scaffolds[a], ]
  a_rep <- a_rep[match(unique(a_rep$start), a_rep$start),]
  a_rep <- a_rep[order(a_rep$start),]
  window_order <- rbind(window_order, a_rep)
}
window_order <- window_order[,4:6]

stat_order <- rbind(c("East", "Central", "Dxy"),
                    c("East", "C_discordant", "Dxy"),
                    c("Central", "C_discordant", "Dxy"),
                    c("East", "Central", "Fst"),
                    c("East", "C_discordant", "Fst"),
                    c("Central", "C_discordant", "Fst"),
                    c("East", "none", "pi"),
                    c("Central", "none", "pi"),
                    c("C_discordant", "none", "pi"),
                    c("East", "none", "Tajima_D"),
                    c("Central", "none", "Tajima_D"),
                    c("C_discordant", "none", "Tajima_D"))


# calculate mean values for each stat of interest
for(a in 1:12) {
  if(a >= 1 & a <= 3) { #dxy
    a_rep <- x[x$pop1 == stat_order[a,1] & x$pop2 == stat_order[a,2] & x$stat == stat_order[a,3],]
    a_rep <- sum(a_rep$calculated_stat * a_rep$number_sites) / sum(a_rep$number_sites)
    writeLines(paste(stat_order[a,1], stat_order[a,2], stat_order[a,3], a_rep, sep="\t\t"))
  } else if(a >= 4 & a <= 6) { #fst
    a_rep <- x[x$pop1 == stat_order[a,1] & x$pop2 == stat_order[a,2] & x$stat == stat_order[a,3],]
    a_rep <- sum(a_rep$calculated_stat * a_rep$number_variable_sites) / sum(a_rep$number_variable_sites)
    writeLines(paste(stat_order[a,1], stat_order[a,2], stat_order[a,3], a_rep, sep="\t\t"))
  } else {
    a_rep <- x[x$pop1 == stat_order[a,1] & x$pop2 == stat_order[a,2] & x$stat == stat_order[a,3],]
    a_rep <- sum(a_rep$calculated_stat * a_rep$number_sites, na.rm = T) / sum(a_rep$number_sites, na.rm = T)
    writeLines(paste(stat_order[a,1], stat_order[a,2], stat_order[a,3], a_rep, sep="\t\t"))
  }
}

# loop for each window to summarize
stat_adds <- list()
for(a in 1:nrow(window_order)) {
  if(a %% 1000 == 0) { print(a) }
  a_rep <- x[x$chr == window_order$chr[a] & x$start == window_order$start[a],]
  stat_adds_rep <- c()
  # for each stat
  for(b in 1:nrow(stat_order)) {
    b_rep <- a_rep[a_rep[,1] == stat_order[b,1] & a_rep[,2] == stat_order[b,2] & a_rep[,3] == stat_order[b,3],]
    if(nrow(b_rep) == 1) {
      stat_adds_rep <- c(stat_adds_rep, b_rep[,9])
    } else {
      stat_adds_rep <- c(stat_adds_rep, NA)
    }
  }
  stat_adds[[a]] <- stat_adds_rep
}

#unlist stat_adds
stat_adds2 <- do.call(rbind.data.frame, stat_adds)

colnames(stat_adds2) <- c("dxy_East_Central", "dxy_East_CoreC", "dxy_Central_CoreC", 
                          "fst_East_Central", "fst_East_CoreC", "fst_Central_CoreC", 
                          "pi_East", "pi_Central", "pi_CoreC", 
                          "TajD_East", "TajD_Central", "TajD_CoreC")
# combine data frames
x <- cbind(window_order, stat_adds2) %>%
  select(-end)

x <- na.omit(x)
rownames(x) <- seq(from=1, to=nrow(x), by=1)


##################################
#Subsetting data to plot

kept_stats <- c("chr", "start", "dxy_East_Central", "dxy_East_CoreC", "dxy_Central_CoreC", 
                "fst_East_Central", "fst_East_CoreC", "fst_Central_CoreC", 
                "pi_East", "pi_Central", "pi_CoreC", 
                "TajD_East", "TajD_Central", "TajD_CoreC")
x_subset <- x[kept_stats]

##################################

# windows
windows <- as.numeric(rownames(x_subset))
# sliding window size
window_size <- 10
# set up row numbers for line plots
line_rows <- seq(from=1,to=nrow(x_subset), by=1)[seq(from=1,to=nrow(x_subset), by=1) %% window_size == 0]

# what are the unique chromosomes?
chr <- unique(x_subset$chr)
scaffold <- chr

# define population order and naming conventions
population_order <- c("East", "Central", "CoreC")
population_names <- c("East", "Central", "CoreC")

# define population comparisons order and naming conventions
pop_split_order <- c("East_Central", "East_CoreC", "Central_CoreC")
pop_split_names <- c("East x Central", "East x CoreC", "Central x CoreC")

# ylim boundaries for each stat.
pi_ylim <- c(0, 0.04)
fst_ylim <- c(0, 0.8)
dxy_ylim <- c(0,0.04)
tajd_ylim <- c(-2.5, 0.5)

#########################################################################
#########################################################################
#########################################################################
# Plot Dxy
#########################################################################
#########################################################################
#########################################################################

# make the plotting polygons
# How this works:
# 1. Creates a list of chromosomes (32 in this case)
# 2. Performs a for loop that uses "a" as whichever # chromosome we're referring to (1-32). 
        # a1 pulls however many windows are in each chromosome (for chr30, it's 3056, but for chr2, it's 2335, and for chr30, it's only 24)
        # a2 pulls which window # each chromosome ends on.
        # chr_polygons[a] creates matrices that ultimately spawn the "polygons" or rectangles that say where each chromosome starts and ends on the plot. 
        # For instance, the rectangle for Chr 1 starts at window #1, ends at window #3056, and spans 0 dxy to 0.012 dxy. The first column is the x values, the second column is the y values. 

chr_polygons <- list()
for(a in 1:length(chr)) {
  a1 <- windows[x_subset$chr == chr[a]]
  a2 <- a1[length(a1)]
  a1 <- a1[1]
  chr_polygons[[a]] <- rbind(c(a1, dxy_ylim[1]), c(a2, dxy_ylim[1]), c(a2, dxy_ylim[2]), c(a1, dxy_ylim[2]), c(a1, dxy_ylim[1]))
} 

# loop for each plot
# How this works:
# 1. First line creates a pdf of certain size dimensions
# 2. Second line creates a matrix of different plots based on however many things are are in "pop_split_order". Here, they're all just stacked on top of each other. Only 1 column. 
# 3. Third line determines the margins or something
# 4. Performs a for loop that uses "b" as whichever population comparison we're referring to. 
    # The first part just plots the chromosome windows. 
    # The second part creates b_rep for each comparison (overwriting it after each plot is made). b_rep is a simplified dataframe containing the chromosome ID, start position of the window, stat for that window, and window #. 
    # Then it plots the points and adds a title. 
    # The third part and fourth parts deal with sliding means or something. I'm not entirely sure. 

pdf("_window_DXY_male_fine.pdf", width=10, height=10)
par(mfrow=c(length(pop_split_order),1))
par(mar=c(2, 2, 0.1, 0.1))
for(b in 1:length(pop_split_order)) {
  # plot the chromosomes
  plot(c(-1,-1), ylim=c(dxy_ylim[1], dxy_ylim[2]), xlim=c(1, length(windows)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="", xlab="")
  odd <- 0
  for(a in 1:length(chr_polygons)) {
    if(odd == 1) {
      polygon(chr_polygons[[a]], col="snow2", border="white")
      odd <- 0	
    } else {
      odd <- 1
    }
  }
  
  # plot the stat
  # subset this individual or comparison
  b_rep <- data.frame(chromosome=as.character(x_subset$chr), start=as.numeric(x_subset$start), stat=as.numeric(x_subset[,colnames(x_subset) == paste("dxy_", pop_split_order[b], sep="")]), window=as.numeric(rownames(x)))
  # plot
  points(rownames(b_rep), b_rep$stat, pch=19, cex=0.05, col="gray65")
  title(paste(pop_split_names[b], sep=""), adj=0.01, line=-1)
  
  # plot sliding mean line plots
  line_x_axis <- line_rows
  line_y_axis <- c()
  line_scaffold <- c()
  for(e in line_rows) {
    line_y_axis <- c(line_y_axis, mean(b_rep$stat[b_rep$window %in% (b_rep$window[e] - (floor(window_size / 2))):(b_rep$window[e] + (floor(window_size / 2))) & b_rep$chromosome == b_rep$chromosome[e]][1:window_size]))
    line_scaffold <- c(line_scaffold, b_rep$chromosome[e])
  }
  line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
  # plot each scaffold at a time (so lines don't connect between scaffolds)
  for(e in 1:length(unique(scaffold))) {
    a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[e],]
    lines(a_rep[,1:2], lwd=0.8, col="gray29")
  }
  
}
dev.off()


#########################################################################
#########################################################################
#########################################################################
# Plot Fst
#########################################################################
#########################################################################
#########################################################################

# make the plotting polygons
chr_polygons <- list()
for(a in 1:length(chr)) {
  a1 <- windows[x_subset$chr == chr[a]]
  a2 <- a1[length(a1)]
  a1 <- a1[1]
  chr_polygons[[a]] <- rbind(c(a1, fst_ylim[1]), c(a2, fst_ylim[1]), c(a2, fst_ylim[2]), c(a1, fst_ylim[2]), c(a1, fst_ylim[1]))
}

# loop for each plot
pdf("_window_FST_male_fine.pdf", width=10, height=10)
par(mfrow=c(length(pop_split_order),1))
par(mar=c(2, 2, 0.1, 0.1))
for(b in 1:length(pop_split_order)) {
  # plot the chromosomes
  plot(c(-1,-1), ylim=c(fst_ylim[1], fst_ylim[2]), xlim=c(1, length(windows)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="", xlab="")
  odd <- 0
  for(a in 1:length(chr_polygons)) {
    if(odd == 1) {
      polygon(chr_polygons[[a]], col="snow2", border="white")
      odd <- 0	
    } else {
      odd <- 1
    }
  }
  
  # plot the stat
  # subset this individual or comparison
  b_rep <- data.frame(chromosome=as.character(x_subset$chr), start=as.numeric(x_subset$start), stat=as.numeric(x_subset[,colnames(x_subset) == paste("fst_", pop_split_order[b], sep="")]), window=as.numeric(rownames(x_subset)))
  
  # plot
  points(rownames(b_rep), b_rep$stat, pch=19, cex=0.05, col="gray65")
  title(paste(pop_split_names[b], sep=""), adj=0.01, line=-1)
  
  # plot sliding mean line plots
  line_x_axis <- line_rows
  line_y_axis <- c()
  line_scaffold <- c()
  for(e in line_rows) {
    line_y_axis <- c(line_y_axis, mean(b_rep$stat[b_rep$window %in% (b_rep$window[e] - (floor(window_size / 2))):(b_rep$window[e] + (floor(window_size / 2))) & b_rep$chromosome == b_rep$chromosome[e]][1:window_size]))
    line_scaffold <- c(line_scaffold, b_rep$chromosome[e])
  }
  line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
  # plot each scaffold at a time (so lines don't connect between scaffolds)
  for(e in 1:length(unique(scaffold))) {
    a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[e],]
    lines(a_rep[,1:2], lwd=0.8, col="gray29")
  }
}
dev.off()

# This paper should help explain the relationship between Fst and Dxy: https://onlinelibrary.wiley.com/doi/pdf/10.1111/evo.14234 
# And see this paper too: https://www.pnas.org/doi/pdf/10.1073/pnas.1713288114 
# Fst is influenced by within-population variation. Dxy is JUST variation between populations. Within-population variation does not elevate Dxy.
# So with low Dxy and high Fst, there hasn't been a lot of variation between populations but there are a lot of ancestral polymorphisms. 
# Dxy only accounts for present variation. Fst accounts for past variation in ancestral populations. 

# "DXY is impacted by levels of ancestral diversity but not current levels of diversity, and may therefore be reduced in FST peaks that 
# have experienced selection within a common ancestor and then recurrently in the descendent lineages"


# "analyses of Dxy will only discriminate between some of the alternative factors that cause genomic islands
# of divergence identified by FST. These problems can be alleviated by (i) reconstructing the past demography of the lineages concerned and (ii) comparing lineages
# with contrasted gene flow and times of divergence along the speciation continuum. This is especially critical for plants in which interspecific gene flow occurs more frequently
# than in animals."

#########################################################################
#########################################################################
#########################################################################
# Plot pi
#########################################################################
#########################################################################
#########################################################################

# make the plotting polygons
chr_polygons <- list()
for(a in 1:length(chr)) {
  a1 <- windows[x_subset$chr == chr[a]]
  a2 <- a1[length(a1)]
  a1 <- a1[1]
  chr_polygons[[a]] <- rbind(c(a1, pi_ylim[1]), c(a2, pi_ylim[1]), c(a2, pi_ylim[2]), c(a1, pi_ylim[2]), c(a1, pi_ylim[1]))
}

# loop for each plot
pdf("_window_pi_male_fine.pdf", width=10, height=10)
par(mfrow=c(length(population_order),1))
par(mar=c(2, 2, 0.1, 0.1))
for(b in 1:length(population_order)) {
  # plot the chromosomes
  plot(c(-1,-1), ylim=c(pi_ylim[1],pi_ylim[2]), xlim=c(1, length(windows)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="", xlab="")
  odd <- 0
  for(a in 1:length(chr_polygons)) {
    if(odd == 1) {
      polygon(chr_polygons[[a]], col="snow2", border="white")
      odd <- 0	
    } else {
      odd <- 1
    }
  }
  
  # plot the stat
  # subset this individual or comparison
  b_rep <- data.frame(chromosome=as.character(x_subset$chr), start=as.numeric(x_subset$start), stat=as.numeric(x_subset[,colnames(x_subset) == paste("pi_", population_order[b], sep="")]), window=as.numeric(rownames(x_subset)))
  
  # plot
  points(rownames(b_rep), b_rep$stat, pch=19, cex=0.05, col="gray65")
  title(paste(population_names[b], sep=""), adj=0.01, line=-1)
  
  # plot sliding mean line plots
  line_x_axis <- line_rows
  line_y_axis <- c()
  line_scaffold <- c()
  for(e in line_rows) {
    line_y_axis <- c(line_y_axis, mean(b_rep$stat[b_rep$window %in% (b_rep$window[e] - (floor(window_size / 2))):(b_rep$window[e] + (floor(window_size / 2))) & b_rep$chromosome == b_rep$chromosome[e]][1:window_size]))
    line_scaffold <- c(line_scaffold, b_rep$chromosome[e])
  }
  line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
  # plot each scaffold at a time (so lines don't connect between scaffolds)
  for(e in 1:length(unique(scaffold))) {
    a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[e],]
    lines(a_rep[,1:2], lwd=0.8, col="gray29")
  }
}
dev.off()


#########################################################################
#########################################################################
#########################################################################
# Plot Tajima D
#########################################################################
#########################################################################
#########################################################################

# make the plotting polygons
chr_polygons <- list()
for(a in 1:length(chr)) {
  a1 <- windows[x_subset$chr == chr[a]]
  a2 <- a1[length(a1)]
  a1 <- a1[1]
  chr_polygons[[a]] <- rbind(c(a1, tajd_ylim[1]), c(a2, tajd_ylim[1]), c(a2, tajd_ylim[2]), c(a1, tajd_ylim[2]), c(a1, tajd_ylim[1]))
}

# loop for each plot
pdf("_window_TajimaD_male_fine.pdf", width=10, height=10)
par(mfrow=c(length(population_order),1))
par(mar=c(2, 2, 0.1, 0.1))
for(b in 1:length(population_order)) {
  # plot the chromosomes
  plot(c(-1,-1), ylim=c(tajd_ylim[1], tajd_ylim[2]), xlim=c(1, length(windows)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="", xlab="")
  odd <- 0
  for(a in 1:length(chr_polygons)) {
    if(odd == 1) {
      polygon(chr_polygons[[a]], col="snow2", border="white")
      odd <- 0	
    } else {
      odd <- 1
    }
  }
  
  # plot the stat
  # subset this individual or comparison
  b_rep <- data.frame(chromosome=as.character(x_subset$chr), start=as.numeric(x_subset$start), stat=as.numeric(x_subset[,colnames(x_subset) == paste("TajD_", population_order[b], sep="")]), window=as.numeric(rownames(x_subset)))
  # plot
  points(rownames(b_rep), b_rep$stat, pch=19, cex=0.05, col="gray65")
  title(paste(population_names[b], sep=""), adj=0.01, line=-1)
  
  # plot sliding mean line plots
  line_x_axis <- line_rows
  line_y_axis <- c()
  line_scaffold <- c()
  for(e in line_rows) {
    line_y_axis <- c(line_y_axis, mean(b_rep$stat[b_rep$window %in% (b_rep$window[e] - (floor(window_size / 2))):(b_rep$window[e] + (floor(window_size / 2))) & b_rep$chromosome == b_rep$chromosome[e]][1:window_size]))
    line_scaffold <- c(line_scaffold, b_rep$chromosome[e])
  }
  line_plotting <- data.frame(line_x_axis=as.numeric(line_x_axis), line_y_axis=as.numeric(unlist(line_y_axis)), line_scaffold=as.character(unlist(line_scaffold)))
  # plot each scaffold at a time (so lines don't connect between scaffolds)
  for(e in 1:length(unique(scaffold))) {
    a_rep <- line_plotting[line_plotting[,3] == unique(scaffold)[e],]
    lines(a_rep[,1:2], lwd=0.8, col="gray29")
  }
}
dev.off()

