# Script for replicating the paper results

install.packages(c("ggplot2","RANN","igraph","MASS","e1071"))
install.packages(c("irlba","kernlab","PPCI","pdfCluster","densityClust","HDclassif","ClusterR"))
# to install the "SyNet" you will need to install "tkrplot"
# if have have trouble doing it use terminal in Ubuntu systems using this command:
# sudo apt-get install tk-dev
install.packages(c("SyNet"))

# load required libraries
# keep in mind that these libraries are reiquired to replicated the results
# presented in the paper, not to run the "i-divclu" algorithm


library(ggplot2)
library(RANN)
library(igraph)
library(SyNet)
library(MASS)
library(e1071)
library(irlba)
library(kernlab)
library(PPCI)
library(pdfCluster)
library(densityClust)
library(HDclassif)
library(ClusterR)


# to avoide functions get masked
mst <- SyNet::mst
clusters <- igraph::clusters

dataname <- list.files("./RdsData/")

dname <- dataname[seq(2, length(dataname), by = 2)]
cname <- dataname[seq(1, length(dataname), by = 2)]

# # for the bio data=====================
# =====================================

# -----------------------------------------
# prepare the data frame to add the results
datnum <- rep(0,length(dname))
# comres_frame <- data.frame(time = datnum, adj = datnum,pur = datnum, vme = datnum,  nmi = datnum, adj2 = datnum, pur2 = datnum ,nmi2 = datnum, clustn = datnum)
comres_frame <- data.frame(time = datnum, adj2 = datnum, pur2 = datnum ,nmi2 = datnum, clustn = datnum)
row.names(comres_frame) <- dname
# ---
idivclu_frame <- comres_frame
idivclu_frame2 <- comres_frame
spect_frame <- comres_frame
kkmean_frame <- comres_frame
agglo_frame <- comres_frame
agglo_frame2 <- comres_frame
dencl_frame <- comres_frame
hddc_frame <- comres_frame
ncutdc_frame <- comres_frame
# -----------------------------------------

# =====================
# --
for (ll in 1:length(dname)){
    # ll <- 10
  
    ddat <- paste("RdsData/",dname[ll], sep = "")
    cdat <- paste("RdsData/",cname[ll], sep = "")
    # -------
    dd <- readRDS(file = ddat)
    dd_class <- as.numeric(readRDS(file = cdat))
    
    # how many times to run the partitioning algorithms
    partimes <- 10
    parver <- rep(NA,partimes)

    # the true number of clusters given as input to all algorithms
    cl_n <- length(unique(dd_class))
    # ------------------------------------
    
    source("i-divclu.R")
    # iDivClu-D
    ptm <- proc.time()
    res<-idivclu(dd, k=cl_n, knei = 5, adjust=0.5,SP="DEN")
    idivclu_frame$time[ll] <- (proc.time() - ptm)[3]
    idivclu_frame$adj2[ll] <- external_validation(as.numeric(factor(res$Clusterids)), dd_class, method = "adjusted_rand_index", summary_stats = F)
    idivclu_frame$pur2[ll] <- external_validation(as.numeric(factor(res$Clusterids)), dd_class, method = "purity", summary_stats = F)
    idivclu_frame$nmi2[ll] <- external_validation(as.numeric(factor(res$Clusterids)), dd_class, method = "nmi", summary_stats = F)
    idivclu_frame$clustn[ll] <- length(table(res$Clusterids))
    # ------------------------------------
    # iDivClu-M
    ptm <- proc.time()
    res<-idivclu(dd, k=cl_n, knei = 5, adjust=0.5,SP="MAR")
    idivclu_frame2$time[ll] <- (proc.time() - ptm)[3]
    idivclu_frame2$adj2[ll] <- external_validation(as.numeric(factor(res$Clusterids)), dd_class, method = "adjusted_rand_index", summary_stats = F)
    idivclu_frame2$pur2[ll] <- external_validation(as.numeric(factor(res$Clusterids)), dd_class, method = "purity", summary_stats = F)
    idivclu_frame2$nmi2[ll] <- external_validation(as.numeric(factor(res$Clusterids)), dd_class, method = "nmi", summary_stats = F)
    idivclu_frame2$clustn[ll] <- length(table(res$Clusterids))
    # ------------------------------------
    # kernel MDDC
    ptm <- proc.time()
    x2 <- kernlab::kpca(as.matrix(dd), kernel = "rbfdot")@rotated
    sol2 <- ncutdc(as.matrix(x2),cl_n)
    ncutdc_frame$time[ll] <- (proc.time() - ptm)[3]
    ncutdc_frame$adj2[ll] <- external_validation(as.numeric(factor(sol2$cluster)), dd_class, method = "adjusted_rand_index", summary_stats = F)
    ncutdc_frame$pur2[ll] <- external_validation(as.numeric(factor(sol2$cluster)), dd_class, method = "purity", summary_stats = F)
    ncutdc_frame$nmi2[ll] <- external_validation(as.numeric(factor(sol2$cluster)), dd_class, method = "nmi", summary_stats = F)
    ncutdc_frame$clustn[ll] <- length(table(sol2$cluster))
    # ------------------------------------
    
    # spectral clustering
    sp_frame <- data.frame(time = parver, adj2 = parver, pur2 = parver ,nmi2 = parver, clustn = parver)
    for (i in 1:partimes){
        ptm <- proc.time()
        sc <- specc(as.matrix(dd), centers=cl_n)
        sp_frame$time[i] <- (proc.time() - ptm)[3]
        sp_frame$adj2[i] <- external_validation(sc, dd_class, method = "adjusted_rand_index", summary_stats = F)
        sp_frame$pur2[i] <- external_validation(sc, dd_class, method = "purity", summary_stats = F)
        sp_frame$nmi2[i] <- external_validation(sc, dd_class, method = "nmi", summary_stats = F)
        sp_frame$clustn[i] <- length(table(sc))
    }
    spect_frame[ll,] <- colMeans(sp_frame, na.rm = TRUE)
    
    # ------------------------------------
    # kernel kmeans
    km_frame <- data.frame(time = parver, adj2 = parver, pur2 = parver ,nmi2 = parver, clustn = parver)
    for (i in 1:partimes){
        ptm <- proc.time()
        
        er <- tryCatch({
        km <- kkmeans(as.matrix(dd), centers=cl_n)
        }, error=function(e){return(1)})
        
        length(er)
        if (length(er) != 1){
          km_frame$time[i] <- (proc.time() - ptm)[3]
          km_frame$adj2[i] <- external_validation(km, dd_class, method = "adjusted_rand_index", summary_stats = F)
          km_frame$pur2[i] <- external_validation(km, dd_class, method = "purity", summary_stats = F)
          km_frame$nmi2[i] <- external_validation(km, dd_class, method = "nmi", summary_stats = F)
          km_frame$clustn[i] <- length(table(km))
        }
    }
    kkmean_frame[ll,] <- colMeans(km_frame, na.rm = TRUE)
    # ------------------------------------

#   agglomerative clustering complete linkage
    ptm <- proc.time()
    distmat <- dist(dd)
    fit <- hclust(distmat, method="single")
    clusterCut <- cutree(fit,cl_n )
    agglo_frame$time[ll] <- (proc.time() - ptm)[3]
    agglo_frame$adj2[ll] <- external_validation(clusterCut, dd_class, method = "adjusted_rand_index", summary_stats = F)
    agglo_frame$pur2[ll] <- external_validation(clusterCut, dd_class, method = "purity", summary_stats = F)
    agglo_frame$nmi2[ll] <- external_validation(clusterCut, dd_class, method = "nmi", summary_stats = F)
    agglo_frame$clustn[ll] <- length(table(clusterCut))
    # ------------------------------------------------------
#   agglomerative clustering average linkage
    ptm <- proc.time()
    distmat <- dist(dd)
    fit <- hclust(distmat, method="average")
    clusterCut <- cutree(fit,cl_n )
    agglo_frame2$time[ll] <- (proc.time() - ptm)[3]
    agglo_frame2$adj2[ll] <- external_validation(clusterCut, dd_class, method = "adjusted_rand_index", summary_stats = F)
    agglo_frame2$pur2[ll] <- external_validation(clusterCut, dd_class, method = "purity", summary_stats = F)
    agglo_frame2$nmi2[ll] <- external_validation(clusterCut, dd_class, method = "nmi", summary_stats = F)
    agglo_frame2$clustn[ll] <- length(table(clusterCut))
    
    
    # ------------------------------------------------------
    # desnityclust
    dn_frame <- data.frame(time = parver, adj2 = parver, pur2 = parver ,nmi2 = parver, clustn = parver)
    for (i in 1:partimes){
      ptm <- proc.time()
      dd_dist <- dist(dd)
      dd_Clust <- densityClust(dd_dist, gaussian=TRUE)
      ff_Clust <- findClusters(dd_Clust, rho=mean(dd_Clust$rho), delta=mean(dd_Clust$delta))
      
      dn_frame$time[i] <- (proc.time() - ptm)[3]
      dn_frame$adj2[i] <- external_validation(ff_Clust$clusters, dd_class, method = "adjusted_rand_index", summary_stats = F)
      dn_frame$pur2[i] <- external_validation(ff_Clust$clusters, dd_class, method = "purity", summary_stats = F)
      dn_frame$nmi2[i] <- external_validation(ff_Clust$clusters, dd_class, method = "nmi", summary_stats = F)
      dn_frame$clustn[i] <- length(table(ff_Clust$clusters))
    }
    dencl_frame[ll,] <- colMeans(dn_frame, na.rm = TRUE)
      
    # ------------------------------------------------------
    # HDDC
    hd_frame <- data.frame(time = parver, adj2 = parver, pur2 = parver ,nmi2 = parver, clustn = parver)
    for (i in 1:partimes){
      ptm <- proc.time()
      prms1 <- hddc(dd, K=cl_n)
      
      hd_frame$time[i] <- (proc.time() - ptm)[3]
      hd_frame$adj2[i] <- external_validation(prms1$class, dd_class, method = "adjusted_rand_index", summary_stats = F)
      hd_frame$pur2[i] <- external_validation(prms1$class, dd_class, method = "purity", summary_stats = F)
      hd_frame$nmi2[i] <- external_validation(prms1$class, dd_class, method = "nmi", summary_stats = F)
      hd_frame$clustn[i] <- length(table(prms1$class))
    }
    hddc_frame[ll,] <- colMeans(hd_frame, na.rm = TRUE)
    # ------------------------------------------------------
    
}


save(idivclu_frame,idivclu_frame2, spect_frame,kkmean_frame, agglo_frame,agglo_frame2, dencl_frame,hddc_frame,ncutdc_frame, file = "example_results.RData")

