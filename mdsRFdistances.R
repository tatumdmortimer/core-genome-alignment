#!/usr/bin/Rscript --vanilla

# This script performs multidimensional analysis on RF distances output by   
# RAxML and plots them
# Usage: mdsRFdistance.R <filename>

# Load required libraries
library(ggplot2)
library(reshape2)

# Read data from file provided by user
args <- commandArgs(trailingOnly = TRUE)
data <- read.table(args[1])
d <- data[,c(1,2,4)]
m <- acast(d, V2~V1, value.var="V4")
dmatrix <- as.dist(m)
fit <- cmdscale(dmatrix, eig=TRUE, k=2)
x <- fit$points[,1]
y <- fit$points[,2]


ExportPlot <- function(gplot, filename, width=2, height=1.5) {
    # Export plot in PDF and EPS.
    # Notice that A4: width=11.69, height=8.27
    ggsave(paste(filename, '.pdf', sep=""), gplot, width=width, height=height)
    postscript(file=paste(filename,'.eps', sep=""), width=width, height=height)
    print(gplot)
    dev.off()
    png(file=paste(filename,'.png',sep=""), width=width*100, height=height*100)
    print(gplot)
    dev.off()
}


p <- qplot(x,y) + geom_point() + theme_bw()
ExportPlot(p, "mds_RF", width=3, height=3)
