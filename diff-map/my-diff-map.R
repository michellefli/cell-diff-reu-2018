# Michelle Li
# diffusion mapping using destiny 

# updated June 27, 2018
# Got destiny package working with Nestorowa data and paul data
# Added interactive 3D plot

# requires destiny package (Github download), rgl for 3D plot
library(destiny)
library(Biobase)
library(rgl)


#### PREPARING DATA ----------------------------------------------

## read in data and convert to exp set (Nestorowa, non-.RData file)
coordinates_gene_counts_flow_cytometry <- read.delim("~/Downloads/coordinates_gene_counts_flow_cytometry.txt", row.names=1)
rawdata <- coordinates_gene_counts_flow_cytometry
cleandata <- na.omit(rawdata) # remove rows with NA
data <- cleandata[,14:ncol(cleandata)] # choose relevant flow cytometry data
dataES <- as.ExpressionSet(data) # convert to expression set

## load .RData file and convert to exp set
load("/Users/Michelle/Downloads/data_Paul_Cell2015/Paul_Cell_MARSseq_GSE72857.RData")
dataES <- ExpressionSet(data)

## note on Expression Set: expression set is formatted so you have cells as the columns and genes as the rows
## it seems that somehow it automatically knows which is which -- when you put in data as itself or its transpose, it comes out the same

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
colors <- rainbow(nclusters) # choose colors

## normalization -- unsure if we need this?
normalizations <- colMeans(exprs(dataES))
normdata <- dataES
exprs(normdata) <- exprs(normdata) - normalizations


#### DIFFUSION MAP -------------------------------------------------

## diffusion map (un-normalized data)
dm <- DiffusionMap(dataES)
# diffusion map and plot using normalized data
dm1 <- DiffusionMap(normdata)


#### PLOTTING ------------------------------------------------------

par(mar = c(5,5,5,5), xpd = "TRUE")

## plot 2d
plot(dm, col = colors[groupsf], pch=20, main = paste("Diffusion Map, Un-Normalized, ", dataname))
plot(dm1, col = colors[groupsf], pch=20, main = paste("Diffusion Map, Normalized, ", dataname))
legend("bottomright", inset=c(0,-0.2), legend=levels(groupsf), col=colors, pch=20, ncol=1, cex=0.75) 

## plot 3d
plot3d(eigenvectors(dm)[,1:3], col = colors[groupsf], pch=20, main = paste("Diffusion Map, Un-Normalized, ", dataname))
plot3d(eigenvectors(dm1)[,1:3], col = colors[groupsf], pch=20, main = paste("Diffusion Map, Normalized, ", dataname))
legend3d("topright", legend = levels(groupsf), inset=c(0.1,0.1), pch = 20, col = colors, ncol=1, cex=1)

## sigmas
sigmas <- find_sigmas(dataES) # normalized
plot(sigmas)
sigmas1 <- find_sigmas(normdata) # unnormalized 
plot(sigmas1)

