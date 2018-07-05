# Michelle Li
# pca

# last updated: July 3, 2018
# tried 2 PCA packages, working with Nestorowa and Paul data
# added handwritten version

library(rgl)


#### PREPARING DATA ----------------------------------------------

## NOTES ON OUR DATA:
## Nestorowa data: 1645 cells x 4290 genes (need transpose for PCA)
## Paul data: 8716 genes x 2730 cells
## Grover data: 46175 genes x 135 cells
## Russ data: 12674 genes x 101 cells

## IMPORTANT NOTE ON DIMENSIONS FOR YOUR DATA MATRIX:
## whatever you want plotted in PCA, whether the cells or the genes, put that in the COLUMNS of your data matrix
## e.g. to get PCA that plots for each CELL, make sure you have cells as columns, genes as rows

## read in data and convert to matrix (Nestorowa, Grover, non-.RData file)
rawdata <- read.delim("~/GitHub/cell-diff-reu-2018/data/coordinates_gene_counts_flow_cytometry.txt", row.names=1)
# rawdata <- read.csv("~/GitHub/cell-diff-reu-2018/data/grover_expression.txt", row.names=1, sep="")
cleandata <- na.omit(rawdata) # remove rows with NA
data <- cleandata[,14:ncol(cleandata)] # choose relevant gene exp data
data <- cleandata
data <- as.matrix(data)

## load data (Paul, .RData file)
load("~/GitHub/cell-diff-reu-2018/data/Paul_Cell_MARSseq_GSE72857.RData")

dataname <- "Paul Data" # enter title of data

## create factor for grouping
# rownames <- row.names(data)
# groups <- as.character(rownames)
# groups[grepl("HSPC",groups)] <- "HSPC" # replace all labels with hspc
# groups[grepl("LT",groups)] <- "LT.HSC" # replace with lt.hsc
# groups[grepl("Prog",groups)] <- "PROG" # same for prog
# groupsf <- factor(groups)
groupsf <- factor(cluster.id) # paul data

nclusters <- nlevels(groupsf) # number of groups to cluster
colors <- rainbow(nclusters) # choose colors


## TOY DATA
data <- matrix(c(1,2,3,2,4,6,4,8,12,3,6,9,5,10,15,6,12,18),
               nrow=6, ncol=3) # rows>columns
dataname <- "Toy Data"



################### Handwritten PCA ##########################

## refer to above to choose whether you need to transpose your data or not
# data <- t(data) # may need to take transpose

## center data
M1centered <- scale(data, center=TRUE, scale=FALSE) # subtracts column means from each column


#### EIG/SVD ------------------------------------------------

## testing covariance
C <- cov(M1centered)
eigC <- eigen(C)



## eig
# C <- cov(M1centered) # will return a cov matrix, ncol x ncol
# eig <- eigen(C) 
# eigvalues <- eig$values[1:3]
# V <- eig$vectors[,1:3] # V=columns are the eig vectors of Gram matrix
# V <- -V # by default, eig vectors in R point in the negative direction. Multiply by -1 to use the positive-pointing vector
# U1eig <- M1centered %*% V # U: columns are the eig vectors of original covariance matrix
# Ueig <- scale(U1eig, center=FALSE, scale=TRUE) # normalize U

# eig: U (eig vectors of original covariance matrix) = M1centered * V (eig vectors of fake covariance matrix (Gram Matrix))
# U from eig is a scale of the U from svd
# unclear on the exact relationship between U from eig and U from svd

## svd
svd <- svd(M1centered)
Usvd <- svd$u # dim: nrowx3
Ssvd <- svd$d
Vsvd <- svd$v # dim: ncolx3

## Calculate PC scores
weights <- t(M1centered) %*% Usvd # columns of scores are the principal components, dim: ncolx3


#### PLOTTING -------------------------------------------------

par(mar = c(5,5,5,5), xpd = "TRUE") # add extra space to margins for legend 

## plot weights 2d 
plot(x=weights[,1], y=weights[,2], 
     xlab="PC1 weights", ylab="PC2 weights", 
     # col = colors[groupsf],
     pch=20, 
     main = paste("PCA Weights (Handwritten),", dataname))
# legend("bottomright", inset = c(0,0), 
#        legend = levels(groupsf), 
#        pch = 20, 
#        col = colors, 
#        ncol=3, cex=0.4)

## plot weights 3d 
plot3d(x=weights[,1], y=weights[,2], z=weights[,3], 
       xlab="PC1 weights", ylab="PC2 weights", zlab="PC3 weights", 
       # col = colors[groupsf], 
       pch=20, 
       main = paste("PCA Weights (Handwritten),", dataname))
par3d(windowRect=c(0,0,1000,1000))
# legend3d("bottomright", inset=c(0.2,0.2), 
#          legend = levels(groupsf), 
#          pch = 20, 
#          col = colors, 
#          ncol=3, cex=1)
rglwidget()


## plot PCs 2d (eigenvectors of original covariance matrix)
plot(x=Usvd[,1], y=Usvd[,2], 
     xlab="PC1", ylab="PC2", 
     # col = colors[groupsf], 
     pch=20, 
     main = paste("PCA (Handwritten),", dataname))
legend("bottomleft", inset = c(0,0), 
       legend = levels(groupsf), 
       pch = 20, 
       col = colors, 
       ncol=3, cex=0.4)

## plot PCs 3d
plot3d(x=Usvd[,1], y=Usvd[,2], z=Usvd[,3], 
       xlab="PC1", ylab="PC2", zlab="PC3", 
       # col = colors[groupsf], 
       pch=20, 
       main = paste("PCA (Handwritten),", dataname))
par3d(windowRect=c(0,0,1000,1000))
legend3d("bottomright", inset=c(0.2,0.2),
         legend = levels(groupsf),
         pch = 20, 
         col = colors, 
         ncol=3, cex=1)
rglwidget()


################### PCA Packages #######################
## Note: I haven't looked at this in a long time, might not work/make sense anymore

#### PCA with function prcomp, uses svd -----------------------

## do pca
pca1 <- prcomp(data, scale = TRUE)
# pca1t <- prcomp(t(datamat), scale=TRUE)
loadings1 <- pca1$rotation # columns are the eigenvectors ncol x ncol (# cells)
x1 <- pca1$x # x = rotated data (centered and scaled if requested) %*% rotation matrix

A <- loadings1 %*% t(x1)
B <- x1 %*% loadings1

## plot
par(mar = c(5,5,5,5), xpd = "TRUE") # add extra space to margins for legend 
plot(x=loadings1[,1], y=loadings1[,2], 
     # col = colors[groupsf], 
     pch=20, 
     main = paste("PCA (prcomp), ", dataname))

plot(x=x1[1,], y=x1[2,], 
     # col = colors[groupsf], 
     pch=20, 
     main = paste("PCA (prcomp),", dataname))

plot(x=A[,1], y=A[,2], 
     # col = colors[groupsf], 
     pch=20, 
     main = paste("PCA (prcomp), ", dataname))

plot(x=B[,1], y=B[,2], 
     # col = colors[groupsf], 
     pch=20, 
     main = paste("PCA (prcomp), ", dataname))


plot(x=loadings1[,1], y=loadings1[,2], col = colors[groupsf], pch=20, main = paste("PCA (prcomp), ", dataname))

legend("topright", inset=c(0,0), legend=levels(groupsf), col=colors,pch=20, ncol=1, cex=0.5)


#### PCA with function PCA() from package "FactoMineR", uses eig -----------------------

library(FactoMineR)

# takes in a dataframe X with n rows (individuals) and p columns (numeric variables)
data <- as.data.frame(t(data))
pca2 <- PCA(data)
coord <- pca2$ind$coord

## plot 2d
plot(pca2, choix="ind")
plot(x=coord[,1], y=coord[,2], habillage="ind", 
     # col.hab=colors[groupsf], 
     label="none",
     title=paste("2D PCA (FactoMineR), ", dataname))
legend("topright", inset = c(0,0), legend = levels(groupsf), pch = 20, col = colors, ncol=3, cex=0.6) # add legend

## plot 3d 
plot3d(x=coord[,1], y=coord[,2], z=coord[,3], xlab="PC1", ylab="PC2", zlab="PC3", col=colors[groupsf], main = paste("3D PCA (FactoMineR), ", dataname))
par3d(windowRect=c(0,0,1000,1000))
legend3d("topright", legend = levels(groupsf), inset=c(0.1,0.1), pch = 20, col = colors, ncol=1, cex=1) # add legend

