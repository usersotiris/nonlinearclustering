
# Nonlinear Dimensionality Reduction for Clustering

## Introduction

Clusters defined in low dimensional manifolds can have highly nonlinear
structure, which can cause linear dimensionality reduction methods to fail. We introduce an approach to divisive hierarchical clustering that is capable of identifying clusters in nonlinear manifolds. This approach uses the isometric mapping (Isomap) to recursively embed (subsets of) the data in one dimension, and then performs a binary partition designed to avoid the splitting of clusters.

The proposed algorithm uses the acronym "i-DivClu" (Isometric mapping for Divisive Clustering).
It is based on the idea that a suitably defined one-dimensional
representation is sufficient for constructing cluster boundaries that split the data
without breaking any of the clusters. Repeating the procedure recursively provides a
theoretically justified and efficient non-linear clustering technique.

We provide two variations of this methodology ,"i-DivClu-M" for maximum margin clustering and "i-DivClu-D"
for density based clustering.

![](https://github.com/usersotiris/nonlinearclustering/blob/master/output-2.png)
*Two steps of the hierarchical procedure using "i-DivClu-D"*

## Reference
If you use this code please reference the corresponding recently published article. 

Sotiris Tasoulis, Nicos G. Pavlidis, Teemu Roos, Nonlinear dimensionality reduction for clustering, Pattern Recognition, Volume 107, 2020, 107508, ISSN 0031-3203,
https://doi.org/10.1016/j.patcog.2020.107508.

[download bibtex](https://github.com/usersotiris/nonlinearclustering/blob/master/bibtex.txt)

## License
This project is licensed under the BSD-3-Clause License - see the LICENSE file for details.

## Example script
```r
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

```


