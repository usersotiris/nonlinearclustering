
install.packages(c("ggplot2","RANN","igraph","MASS","e1071"))
install.packages(c("irlba","kernlab","PPCI"))
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

    
# to avoide mst function get masked by SyNet
mst <- SyNet::mst


# source the file that includes all function required for executing "i-divclu"
source("i-divclu.R")
# ------------------------------------------

dname <- list.files("./Rdata/")

# -----------------------------------------
# prepare the data frame to add the results
datnum <- rep(0,length(dname))
comres_frame <- data.frame(time = datnum, pur = datnum, vme = datnum)
row.names(comres_frame) <- dname

idivclu_frame <- comres_frame
spect_frame <- comres_frame
kkmean_frame <- comres_frame
# -----------------------------------------

for (ll in 1:length(dname)){

    ddat <- paste("Rdata/",dname[ll], sep = "")
    load(ddat)
    # ------------------------------------------
    # define how many times to rerun the partitioning algorithms (re-initilization)
    partimes <- 20
    parver <- rep(NA,partimes)

    # the true number of clusters given as input to all algorithms
    cl_n <- length(unique(dd_class))
    # ------------------------------------
    # iDivClu
    ptm <- proc.time()
    res<-idivclu(dd, k=cl_n, knei = 5, adjust=0.5)
    idivclu_frame$time[ll] <- (proc.time() - ptm)[3]
    idivclu_frame$pur[ll] <- cluster_performance(res$Clusterids, dd_class)[2]
    idivclu_frame$vme[ll] <- cluster_performance(res$Clusterids, dd_class)[3]
    # ------------------------------------
    # spectral clustering
    sp_frame <- data.frame(time = parver, pur = parver, vme = parver)
    for (i in 1:partimes){
        ptm <- proc.time()
        sc <- specc(as.matrix(dd), centers=cl_n)
        sp_frame$time[i] <- (proc.time() - ptm)[3]
        sp_frame$pur[i] <- cluster_performance(sc,factor(dd_class))[2]
        sp_frame$vme[i] <- cluster_performance(sc,factor(dd_class))[3]
    }
    spect_frame[ll,] <- colMeans(sp_frame)
    # ------------------------------------
    # kernel kmeans
    km_frame <- data.frame(time = parver, pur = parver, vme = parver)
    for (i in 1:partimes){
        ptm <- proc.time()
        
        er <- tryCatch({
        km <- kkmeans(as.matrix(dd), centers=cl_n)
        }, error=function(e){return(1)})
        
        length(er)
        if (length(er) != 1){
            km_frame$time[i] <- (proc.time() - ptm)[3]
            km_frame$pur[i] <- cluster_performance(km,factor(dd_class))[2]
            km_frame$vme[i] <- cluster_performance(km,factor(dd_class))[3]
        }
    }
    kkmean_frame[ll,] <- colMeans(km_frame, na.rm = TRUE)
    # ------------------------------------
 
}

save(idivclu_frame, spect_frame,kkmean_frame, file = "real_results.RData")

# load("real_results.RData")
# =======================================================================

# results reports the performance of each algorithm with respect to purity and
# V-measure, as well as average computational time in seconds over 20 replications
# of each algorithm on each dataset.


