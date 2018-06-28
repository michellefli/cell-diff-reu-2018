# Michelle Li
# pca

# last updated: June 27, 2018
# tried 2 PCA packages, working with Nestorowa and Paul data
# added handwritten version


#### PREPARING DATA ----------------------------------------------

## read in data (Nestorowa, non-.RData file)
coordinates_gene_counts_flow_cytometry <- read.delim("~/Downloads/coordinates_gene_counts_flow_cytometry.txt", row.names=1)
rawdata <- coordinates_gene_counts_flow_cytometry
cleandata <- na.omit(rawdata) # remove rows with NA
data <- cleandata[,14:ncol(cleandata)] # choose relevant flow cytometry data

## load data 
load("/Users/Michelle/Downloads/data_Paul_Cell2015/Paul_Cell_MARSseq_GSE72857.RData")

dataname <- "Paul Data" # enter title of data

## MAKE SURE you have #rows>#columns
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
C <- cov(M1centered) # gives you Gram matrix (smaller dim than normal covariance matrix)


#### EIG/SVD ------------------------------------------------

## eig
technique <- "eig"
eig <- eigen(C) 
V <- eig$vectors # V=columns are the eig vectors of Gram matrix
V <- -V # by default, eig vectors in R point in the negative direction. Multiply by -1 to use the positive-pointing vector
U1 <- M1centered %*% V # U: columns are the eig vectors of original covariance matrix
U <- scale(U1, center=FALSE, scale=TRUE) # normalize U

# eig: U (eig vectors of original covariance matrix) = M1centered * V (eig vectors of fake covariance matrix (Gram Matrix))
# U from eig is a scale of the U from svd
# unclear on the exact relationship between U from eig and U from svd

## svd
technique <- "svd"
svd <- svd(M1centered)
U <- svd$u

## Calculate PC scores
scores <- t(U) %*% M1centered # rows of scores are the principal components


#### PLOTTING ----------------------------------------------

par(mar = c(5,5,5,5), xpd = "TRUE") # add extra space to margins for legend 

## plot scores 2d
plot(x=scores[1,], y=scores[2,], xlab="PC1 scores", ylab="PC2 scores", col = colors[groupsf], pch=20, main = paste("PCA Scores (Handwritten), ", dataname, ", ", technique))
legend("bottomright", inset = c(0,0), legend = levels(groupsf), pch = 20, col = colors, ncol=3, cex=0.4)

## plot scores 3d
plot3d(x=scores[1,], y=scores[2,], z=scores[3,], xlab="PC1 scores", ylab="PC2 scores", zlab="PC3 scores", col = colors[groupsf], pch=20, main = paste("PCA Scores (Handwritten), ", dataname, ", ", technique))
legend3d("bottomright", legend = levels(groupsf), inset=c(0.2,0.2), pch = 20, col = colors, ncol=3, cex=1)

## plot PCs 2d (eigenvectors of origianl covariance matrix)
plot(x=U[1,], y=U[2,], xlab="PC1", ylab="PC2", col = colors[groupsf], pch=20, main = paste("PCA (Handwritten), ", dataname, ", ", technique))



################### PCA Packages #######################

#### PCA with function prcomp, uses svd -----------------------

## do pca
pca1 <- prcomp(X, scale = TRUE)
# pca1t <- prcomp(t(datamat), scale=TRUE)
loadings <- pca1$rotation # loadings
scores <- pca1$x # PCs (scores)

## plot
par(mar = c(5,5,5,5), xpd = "TRUE") # add extra space to margins for legend 
plot(scores, col = colors[groupsf], pch=20, main = paste("PCA (prcomp), ", dataname))
legend("topright", inset=c(0,0), legend=levels(groupsf), col=colors,pch=20, ncol=1, cex=0.5)


#### PCA with function PCA() from package "FactoMineR", uses eig -----------------------

library(FactoMineR)
pca2 <- PCA(X)

## plot 2d
plot(pca2, habillage="ind", col.hab=colors[groupsf], label="none", legend=levels(groupsf), title=paste("2D PCA (FactoMineR), ", dataname))
legend("topright", inset = c(0,0), legend = levels(groupsf), pch = 20, col = colors, ncol=3, cex=0.6) # add legend

## plot 3d 
library(rgl)
coord <- pca2$ind$coord
plot3d(x=coord[,1], y=coord[,2], z=coord[,3], xlab="PC1", ylab="PC2", zlab="PC3", col=colors[groupsf], main = paste("3D PCA (FactoMineR), ", dataname))
legend3d("topright", legend = levels(groupsf), inset=c(0.1,0.1), pch = 20, col = colors, ncol=1, cex=1) # add legend

