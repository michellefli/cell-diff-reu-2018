# Michelle Li
# pca

# last updated: July 20, 2018
# tried 2 PCA packages, working with Nestorowa and Paul data
# added handwritten version

library(rgl)

#### PREPARING DATA ----------------------------------------------

## NOTES ON OUR DATA:
## Nestorowa data: 1645 cells x 4290 genes
## Paul data: 8716 genes x 2730 cells (need transpose for PCA)
## Grover data: 46175 genes x 135 cells (need transpose for PCA)
## Rockne data: 12674 genes x 101 cells (need transpose for PCA)

## PCA DIMENSIONS:
## cells = rows, genes = columns 

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


################### Handwritten PCA ##########################

## refer to above to choose whether you need to transpose your data or not
# data <- t(data) # may need to take transpose

## center data
M1centered <- scale(data, center=TRUE, scale=FALSE) # subtracts column means from each column


#### EIG/SVD ------------------------------------------------

## NOTE: Eig has not been looked at in a long time

## testing covariance
# C <- cov(M1centered) # returns ncol x ncol
# eigC <- eigen(C)
# eigvectors <- eigC$vectors

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
svd <- svd(M1centered, nu=3, nv=3) # only first 3 left/right singular vectors for efficiency
U <- svd$u # dim: nrowx3
V <- svd$v # dim: ncolx3
## Note: U = data %*% V %*% S^-1

## Calculate PC scores
# weights <- t(M1centered) %*% U # dim: ncolx3
scores <- M1centered %*% V # dim: nrowx3


#### PLOTTING -------------------------------------------------

# par(mar = c(5.1,4.1,4.1,8.1), xpd = "TRUE") # use if you need extra space in margins for legend 

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
       ncol=3, cex=0.4, pt.cex=1)

## plot scores 3d 
plot3d(x=scores[,1], y=scores[,2], z=scores[,3], 
       xlab="PC1 scores", ylab="PC2 scores", zlab="PC3 scores", 
       col = colors[groupsf],
       pch=20, 
       main = paste("PCA Scores (Handwritten),", dataname))
legend3d("bottomright", inset=c(0.2,0.2),
         legend = levels(groupsf),
         pch = 20,
         col = colors,
         ncol=3, cex=1)

## plot PCs 2d (eigenvectors of covariance matrix)
plot(x=V[,1], y=V[,2], 
     xlab="PC1", ylab="PC2", 
     col = colors[groupsf],
     pch=20, 
     main = paste("Principal Components (Handwritten PCA),", dataname))
legend("bottomleft", inset = c(0,0), 
       legend = levels(groupsf), 
       pch = 20, 
       col = colors, 
       ncol=3, cex=0.4, pt.cex=1)

## plot PCs 3d
plot3d(x=V[,1], y=V[,2], z=V[,3], 
       xlab="PC1", ylab="PC2", zlab="PC3", 
       col = colors[groupsf],
       pch=20, 
       main = paste("Principal Components (Handwritten PCA),", dataname))
legend3d("bottomright", inset=c(0.2,0.2),
         legend = levels(groupsf),
         pch = 20, 
         col = colors, 
         ncol=3, cex=1)


################### PCA Packages #######################

#### PCA with function prcomp, uses svd -----------------------

## run pca using prcomp
pca1 <- prcomp(data, center=TRUE, scale = TRUE)
rotation1 <- pca1$rotation # columns are the eigenvectors ncol x ncol (# cells)
x1 <- pca1$x # x = rotated data (centered and scaled if requested) %*% rotation matrix
## rotation1 is equivalent to V, x1 is equivalent to scores

## plot
par(mar = c(5,5,5,5), xpd = "TRUE") # add extra space to margins for legend 

## plot PCs (eigenvectors) 
## 2d
plot(x=rotation1[,1], y=rotation1[,2], 
     xlab="PC1", ylab="PC2",
     col = colors[groupsf],
     pch=20, 
     main = paste("Principal Components (prcomp), ", dataname))
legend("bottomright", inset = c(0,0),
       legend = levels(groupsf),
       pch = 20,
       col = colors,
       ncol=1, cex=0.4, pt.cex=1)

## 3d
plot3d(x=rotation1[,1], y=rotation1[,2], z=rotation1[,3],
     xlab="PC1", ylab="PC2", zlab="PC3",
     col = colors[groupsf],
     pch=20, 
     main = paste("Principal Components (prcomp), ", dataname))
legend3d("bottomright", inset = c(0,0),
       legend = levels(groupfs),
       pch = 20,
       col = colors,
       ncol=1, cex=1)

## plot scores
## 2d
plot(x=x1[,1], y=x1[,2],
     xlab="PC1 scores", ylab="PC2 scores",
     col = colors[groupsf],
     pch=20, 
     main = paste("PCA (prcomp),", dataname))
legend("bottomright", inset=c(0.2,0.2),
         legend = levels(groupsf),
         pch = 20, 
         col = colors, 
         ncol=3, cex=0.4, pt.cex=1)

## 3d
plot3d(x=x1[,1], y=x1[,2], z=x1[,3],
       xlab="PC1 scores", ylab="PC2 scores", zlab="PC3 scores",
       col = colors[groupsf],
       pch=20, 
       main = paste("PCA (prcomp),", dataname))
legend3d("bottomright", inset=c(0.2,0.2),
         legend = levels(groupsf),
         pch = 20, 
         col = colors, 
         ncol=3, cex=1)


#### PCA with function PCA() from package "FactoMineR", uses eig -----------------------

library(FactoMineR)

# takes in a dataframe X with n rows (individuals) and p columns (numeric variables)
data <- as.data.frame(data)
pca2 <- PCA(data) # n rows (individuals), p columns (numeric variables)
coord <- pca2$ind$coord

## plot 2d
plot(x=coord[,1], y=coord[,2], 
     xlab = "PC1", ylab = "PC2", 
     col=colors[groupsf],
     pch = 20,
     main=paste("2D PCA (FactoMineR),", dataname))
legend("topright", inset = c(0,0), 
       legend = levels(groupsf), 
       col = colors,
       pch = 20, 
       ncol=3, cex=0.6) # add legend

## plot 3d 
plot3d(x=coord[,1], y=coord[,2], z=coord[,3],
       xlab="PC1", ylab="PC2", zlab="PC3", 
       col=colors[groupsf], 
       main = paste("3D PCA (FactoMineR), ", dataname))
legend3d("topright", inset=c(0.1,0.1), 
         legend = levels(groupsf), 
         col = colors,
         pch = 20,
         ncol=1, cex=1)



