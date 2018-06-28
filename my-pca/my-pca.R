# Michelle Li
# pca

# last updated: June 28, 2018
# tried 2 PCA packages, working with Nestorowa and Paul data
# added handwritten version

library(rgl)


#### PREPARING DATA ----------------------------------------------

## IMPORTANT NOTE ON DIMENSIONS FOR YOUR DATA MATRIX:
## whatever you want plotted in PCA, whether the cells or the genes, put that in the COLUMNS of your data matrix
## e.g. to get PCA that plots for each CELL, make sure you have cells as columns, genes as rows
## nestorowa data: 1645 cells, 4290 genes. 1645 x 4290. want to use transpose
## for Paul data: use regular data matrix (genes x cells) 

## read in data (Nestorowa, .txt file)
coordinates_gene_counts_flow_cytometry <- read.delim("~/Downloads/coordinates_gene_counts_flow_cytometry.txt", row.names=1)
rawdata <- coordinates_gene_counts_flow_cytometry
cleandata <- na.omit(rawdata) # remove rows with NA
data <- cleandata[,14:ncol(cleandata)] # choose relevant data
data <- as.matrix(data)

## load data (Paul, .RData file)
load("/Users/Michelle/Downloads/data_Paul_Cell2015/Paul_Cell_MARSseq_GSE72857.RData")

dataname <- "Nestorowa Data" # enter title of data

## refer to above to choose whether you need to transpose your data or not
X <- data
# X <- t(data) # may need to take transpose

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


################### Handwritten PCA ##########################

## center data
M1centered <- scale(X, center=TRUE, scale=FALSE) # subtracts column means from each column


#### EIG/SVD ------------------------------------------------

## eig
C <- cov(M1centered) # will return a cov matrix, ncol x ncol
eig <- eigen(C) 
V <- eig$vectors # V=columns are the eig vectors of Gram matrix
V <- -V # by default, eig vectors in R point in the negative direction. Multiply by -1 to use the positive-pointing vector
U1 <- M1centered %*% V # U: columns are the eig vectors of original covariance matrix
U <- scale(U1, center=FALSE, scale=TRUE) # normalize U

# eig: U (eig vectors of original covariance matrix) = M1centered * V (eig vectors of fake covariance matrix (Gram Matrix))
# U from eig is a scale of the U from svd
# unclear on the exact relationship between U from eig and U from svd

## svd
svd <- svd(M1centered, nu=3, nv=0)
U <- svd$u # dim: nrowx3

## Calculate PC scores
scores <- t(M1centered) %*% U # columns of scores are the principal components, dim: ncolx3


#### PLOTTING -------------------------------------------------

par(mar = c(5,5,5,5), xpd = "TRUE") # add extra space to margins for legend 

## plot scores 2d
plot(x=scores[,1], y=scores[,2], 
     xlab="PC1 scores", ylab="PC2 scores", 
     col = colors[groupsf], 
     pch=20, 
     main = paste("PCA Scores (Handwritten),", dataname))
legend("bottomright", inset = c(0,0), 
       legend = levels(groupsf), 
       pch = 20, 
       col = colors, 
       ncol=3, cex=0.4)

## plot scores 3d
plot3d(x=scores[,1], y=scores[,2], z=scores[,3], 
       xlab="PC1 scores", ylab="PC2 scores", zlab="PC3 scores", 
       col = colors[groupsf], 
       pch=20, 
       main = paste("PCA Scores (Handwritten),", dataname))
par3d(windowRect=c(0,0,1000,1000))
legend3d("bottomright", inset=c(0.2,0.2), 
         legend = levels(groupsf), 
         pch = 20, 
         col = colors, 
         ncol=3, cex=1)

## plot PCs 2d (eigenvectors of original covariance matrix)
plot(x=U[,1], y=U[,2], 
     xlab="PC1", ylab="PC2", 
     col = colors[groupsf], 
     pch=20, 
     main = paste("PCA (Handwritten),", dataname))
legend("bottomleft", inset = c(0,0), 
       legend = levels(groupsf), 
       pch = 20, 
       col = colors, 
       ncol=3, cex=0.4)

## plot PCs 3d
plot3d(x=U[,1], y=U[,2], z=U[,3], 
       xlab="PC1", ylab="PC2", zlab="PC3", 
       col = colors[groupsf], 
       pch=20, 
       main = paste("PCA (Handwritten),", dataname))
par3d(windowRect=c(0,0,1000,1000))
legend3d("bottomright", inset=c(0.2,0.2),
         legend = levels(groupsf),
         pch = 20, 
         col = colors, 
         ncol=3, cex=1)


################### PCA Packages #######################
## Note: I haven't looked at this in a long time, might not work/make sense anymore

#### PCA with function prcomp, uses svd -----------------------

## do pca
pca1 <- prcomp(X, scale = TRUE)
# pca1t <- prcomp(t(datamat), scale=TRUE)
loadings1 <- pca1$rotation # columns are the eigenvectors
scores1 <- pca1$x # PCs (scores)

A <- loadings1%*%scores1

## plot
par(mar = c(5,5,5,5), xpd = "TRUE") # add extra space to margins for legend 
plot(x=scores1[1,], y=scores1[2,], col = colors[groupsf], pch=20, main = paste("PCA (prcomp), ", dataname))

plot(x=A[,1], y=A[,2], col = colors[groupsf], pch=20, main = paste("PCA (prcomp), ", dataname))


plot(x=loadings1[,1], y=loadings1[,2], col = colors[groupsf], pch=20, main = paste("PCA (prcomp), ", dataname))

legend("topright", inset=c(0,0), legend=levels(groupsf), col=colors,pch=20, ncol=1, cex=0.5)


#### PCA with function PCA() from package "FactoMineR", uses eig -----------------------

library(FactoMineR)

# takes in a dataframe X with n rows (individuals) and p columns (numeric variables)
X <- as.data.frame(X)
pca2 <- PCA(X)
coord <- pca2$ind$coord

## plot 2d
plot(pca2, choix="ind")
plot(x=coord[,1], y=coord[,2], habillage="ind", col.hab=colors[groupsf], label="none", legend=levels(groupsf), title=paste("2D PCA (FactoMineR), ", dataname))
legend("topright", inset = c(0,0), legend = levels(groupsf), pch = 20, col = colors, ncol=3, cex=0.6) # add legend

## plot 3d 
plot3d(x=coord[,1], y=coord[,2], z=coord[,3], xlab="PC1", ylab="PC2", zlab="PC3", col=colors[groupsf], main = paste("3D PCA (FactoMineR), ", dataname))
par3d(windowRect=c(0,0,1000,1000))
legend3d("topright", legend = levels(groupsf), inset=c(0.1,0.1), pch = 20, col = colors, ncol=1, cex=1) # add legend

