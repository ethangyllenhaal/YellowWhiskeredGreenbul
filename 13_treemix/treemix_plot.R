# the source is the path to the plotting_funcs.R file
# (see TreemixManual for more info: 
## http://gensoft.pasteur.fr/docs/treemix/1.12/treemix_manual_10_1_2012.pdf)
# the input for "plot_tree" is the path to the folder your treemix output
## with the stem of the output file names (e.g. sympos for me)
# includes the 
source("C:/Documents/Projects/Greenbul/treemix/plotting_funcs.R")
setwd("C:/Documents/Projects/Greenbul/treemix")

par(mfrow=c(1,3), mar=c(0,2,0,0), oma=c(1,2,1,0))

plot_tree("out/greenbul51_90_jk0") # supp
plot_tree("out/greenbul51_90_jk1") # supp
plot_tree("out/greenbul51_90_jk2") # main
