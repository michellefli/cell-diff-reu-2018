# Michelle Li
# t-SNE code from van der Maaten/jdonaldson

# last updated: June 27, 2018 
# ran Nestorowa flow cytometry and gene exp
# ran Paul data, takes forever, very poor plot
# 2D and 3D plot

# load required packages: tsne (can download from CRAN) and rgl for 3D plot
library(tsne)
library(rgl)


#### PREPARING DATA ----------------------------------------------

## import data and convert into matrix
rawdata <- read.delim("~/GitHub/cell-diff-reu-2018/data/coordinates_gene_counts_flow_cytometry.txt", row.names=1)
# rawdata <- read.csv("~/GitHub/cell-diff-reu-2018/data/grover_expression.txt", row.names=1, sep="")
cleandata <- na.omit(rawdata) # remove NAs
data <- cleandata[,14:ncol(cleandata)] # choose gene exp rows
X <- as.matrix(data) # convert to matrix

## load .RData (Paul)
# load("~/GitHub/cell-diff-reu-2018/data/Paul_Cell_MARSseq_GSE72857.RData") # paul data
# X <- data

dataname <- "Nestorowa Data (Gene Exp)" # set name

## create factor for grouping
rownames <- row.names(X)
groups <- as.character(rownames)
groups[grepl("HSPC",groups)] <- "HSPC" # replace all labels with hspc
groups[grepl("LT",groups)] <- "LT.HSC" # replace with lt.hsc
groups[grepl("Prog",groups)] <- "PROG" # same for prog
groupsf <- factor(groups)
# groupsf <- factor(cluster.id) # paul data

nclusters <- nlevels(groupsf) # number of groups to cluster
colors <- rainbow(nclusters)


#### T-SNE -------------------------------------------------

# Note: can take very, very long with large datasets
# 300 max iterations tends to work well (default is 1000)
ydata <- tsne(X, max_iter = 300) # ydata = matrix(rnorm(k * n),n) 
ydata3d <- tsne(X, k=3, max_iter = 300)


#### PLOTTING -------------------------------------------------

par(mar = c(5,5,5,5), xpd = "TRUE") 

## plot 2d
plot(ydata, 
     col = colors[groupsf], 
     pch=20,
     main = paste("2D t-SNE Plot of ", dataname)) 
legend("topright", inset = c(0,0), 
       legend = levels(groupsf), 
       pch = 20, 
       col = colors, 
       ncol=3, cex=0.5) 

## plot 3d 
plot3d(ydata3d, 
       col=colors[groupsf], 
       pch=20,
       main = paste("3D t-SNE Plot of ", dataname))
legend3d("topright", inset=c(0.1,0.1),
         legend = levels(groupsf),  
         pch = 20, 
         col = colors, 
         ncol=3, cex=1)

