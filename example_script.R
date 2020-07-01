# Example Script

install.packages(c("ggplot2","RANN","igraph"))

# to install the "SyNet" the "tkrplot" need to be installed as well.
# if have have trouble doing it use terminal in Ubuntu systems using this command:
# sudo apt-get install tk-dev
install.packages(c("SyNet"))

# load required libraries

library(RANN)
library(igraph)
library(SyNet)


# to avoid functions get masked
mst <- SyNet::mst
clusters <- igraph::clusters


# source the file including all required functions for i-DivClu
source("i-divclu.R")


# load the toy 2 dimensional dataset
dd <- read.table("2d_data/toy_non_linear_new2018.dat")
dd_class <- read.table("2d_data/toy_non_linear_new2018_class.dat")


#get the true number of clusters
cl_n <- length(table(dd_class))

# Run both versions of i-DivClu and plot the data with respect to the retrieved cluster labels

# ------------------------------------
# iDivClu-D
res<-idivclu(dd, k=cl_n,SP="DEN")
plot(dd$V1, dd$V2, col = as.factor(res$Clusterids))

# ------------------------------------
# iDivClu-M
res<-idivclu(dd, k=cl_n,SP="MAR")
plot(dd$V1, dd$V2, col = as.factor(res$Clusterids))
# ------------------------------------


# Optional parameters for the algorithm:
# k the number of nearest neighbors used when building the kgraph
# adjust the multiplier of the bandwidth parameter used for the calculation of the kernel density estimation
# interpolation_points number of interpolation points upon which the calculation of the kernel density estimation takes place
# min_outlier minimum number of samples allowed to constitute a cluster
 
# Outputs:
# Clusterids vector of assigned cluster ids
# leafs information regarding the tree structure
# max_depth maximum depth of the resulting tree
