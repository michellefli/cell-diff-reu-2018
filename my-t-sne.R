# Michelle Li
# t-SNE code from van der Maaten/jdonaldson

# last updated: July 1, 2018
# successfully ran Nestorowa and Paul data, but not Grover data
# added option to select cluster(s) to highlight in tsne plots

# load required packages: tsne (can download from CRAN) and rgl for 3D plot
library(tsne)
library(rgl)


#### PREPARING DATA ----------------------------------------------

## NOTES ON OUR DATA:
## Nestorowa data: 1645 cells x 4290 genes 
## Paul data: 8716 genes x 2730 cells (need transpose for t-SNE)
## Grover data: 46175 genes x 135 cells (need transpose for t-SNE)

## import data and convert into matrix
rawdata <- read.delim("~/GitHub/cell-diff-reu-2018/data/coordinates_gene_counts_flow_cytometry.txt", row.names=1)
# rawdata <- read.csv("~/GitHub/cell-diff-reu-2018/data/grover_expression.txt", row.names=1, sep="")
cleandata <- na.omit(rawdata) # remove NAs
data <- cleandata[,5:13] # choose flow cytometry columns
# data <- cleandata[,14:ncol(cleandata)] # choose gene exp columns
data <- as.matrix(data) # convert to matrix

## load .RData (Paul)
# load("~/GitHub/cell-diff-reu-2018/data/Paul_Cell_MARSseq_GSE72857.RData") # paul data

dataname <- "Nestorowa Data (Gene Exp)" # set name

## create factor for grouping
rownames <- row.names(data)
groups <- as.character(rownames)
groups[grepl("HSPC",groups)] <- "HSPC" # replace all labels with hspc
groups[grepl("LT",groups)] <- "LT.HSC" # replace with lt.hsc
groups[grepl("Prog",groups)] <- "PROG" # same for prog
groupsf <- factor(groups)
# groupsf <- factor(cluster.id) # paul data

nclusters <- nlevels(groupsf) # number of groups to cluster
colors <- rainbow(nclusters)


#### T-SNE -------------------------------------------------

# for t-SNE: rows=cells, columns=genes
# Nestorowa: use original
# Paul and Grover: use transpose
data <- t(data)

# Note: can take very, very long with large datasets
# 200-300 max iterations tends to work well (default is 1000), 100 is too little
ydata3d <- tsne(data, k=3, max_iter = 200) # ydata = matrix(rnorm(k * n),n)


#### PLOTTING -------------------------------------------------

par(mar = c(5,5,5,5), xpd = "TRUE") 

## plot 2d
plot(x=ydata3d[,1], y=ydata3d[,2],
     col = colors[groupsf], 
     pch=20,
     main = paste("2D t-SNE Plot of", dataname)) 
legend("topright", inset = c(0,0), 
       legend = levels(groupsf), 
       pch = 20, 
       col = colors, 
       ncol=3, cex=0.5) 

## plot 3d 
plot3d(x=ydata3d[,1], y=ydata3d[,2], z=ydata3d[,3], 
       col=colors[groupsf], 
       pch=20,
       main = paste("3D t-SNE Plot of", dataname))
par3d(windowRect=c(0,0,1000,1000))
legend3d("topright", inset=c(0.1,0.1),
         legend = levels(groupsf),  
         pch = 20, 
         col = colors, 
         ncol=3, cex=1)


## select cluster(s) to highlight 
clusterindex <- 10 # select index of the cluster you want to visualize
color1 <- colors
color1[-clusterindex] <- "black"

## plot 2d select clusters
plot(x=ydata3d[,1], y=ydata3d[,2], 
     col = color1[groupsf], 
     pch=20,
     main = paste("2D t-SNE Plot of Cluster", clusterindex, dataname, "(tsne)")) 

## plot 3d select clusters
plot3d(x=ydata3d[,1], y=ydata3d[,2], z=ydata3d[,3], 
       col = color1[groupsf],
       pch=20,
       main = paste("3D t-SNE Plot of Cluster", clusterindex, dataname, "(tsne)"))



