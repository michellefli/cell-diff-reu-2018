# Michelle Li
# diffusion mapping using destiny package

# updated June 28, 2018
# Got destiny package working with Nestorowa data and paul data
# Added interactive 3D plot

# required packages
library(destiny)
library(Biobase)
library(rgl)


#### PREPARING DATA ----------------------------------------------

## NOTES ON DATA:
## nestorowa data: 1645 cells x 4290 genes (need transpose for diff mapping)
## for Paul data: 8716 genes x 2730 cells
## Grover data: 46175 genes x 135 cells

## read in data and convert to matrix (Nestorowa, Grover, non-.RData file)
# rawdata <- read.delim("~/GitHub/cell-diff-reu-2018/data/coordinates_gene_counts_flow_cytometry.txt", row.names=1)
rawdata <- read.csv("~/GitHub/cell-diff-reu-2018/data/grover_expression.txt", row.names=1, sep="")
cleandata <- na.omit(rawdata) # remove rows with NA
# data <- cleandata[,14:ncol(cleandata)] # choose relevant flow cytometry data
data <- cleandata
data <- as.matrix(data)

## load .RData file 
# load("~/GitHub/cell-diff-reu-2018/data/Paul_Cell_MARSseq_GSE72857.RData")

dataname <- "Grover Data" # set name

## convert to expression set
## NOTE ON EXPRESSION SET DIMENSIONS: 
## cells = columns and genes = rows
# data <- t(data) # may need to transpose to have correct rows/columns
dataES <- ExpressionSet(assayData=data) # you must include the 'assayData=' for this to come out correctly

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
plot(x=eigenvectors(dm)[,1], y=eigenvectors(dm)[,2],
     xlab="DC1", ylab="DC2",
     # col = colors[groupsf], 
     pch=20, 
     main = paste("Diffusion Map (destiny), Un-Normalized, ", dataname))
legend("bottomright", inset=c(0,0), 
       legend=levels(groupsf), 
       col=colors, 
       pch=20, 
       ncol=3, cex=0.3) 

## plot 3d
plot3d(eigenvectors(dm)[,1:3], 
       # col = colors[groupsf], 
       pch=20, 
       main = paste("Diffusion Map (destiny), Un-Normalized, ", dataname))
par3d(windowRect=c(0,0,1000,1000))
legend3d("bottomright", inset=c(0.1,0.1),
         legend = levels(groupsf),
         pch = 20, 
         col = colors, 
         ncol=3, cex=1)

## sigmas
sigmas <- find_sigmas(dataES) # unnormalized
plot(sigmas)
sigmas1 <- find_sigmas(normdata) # normalized 
plot(sigmas1)

