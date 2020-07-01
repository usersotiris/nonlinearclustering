

cal_dist_mat <- function(dd,kapa){
  res<- nn2(dd, query= dd, k=kapa, treetype = c("kd"), searchtype = c("standard"), eps = 0)
  em_mat <- matrix(0L, nrow = nrow(dd), ncol = nrow(dd))
  
  for (i in 1:nrow(dd)){
    em_mat[i,res$nn.idx[i,]] <- res$nn.dists[i,]
  }
  
  net3<-graph.adjacency(em_mat, mode=c("undirected"), weighted=TRUE, diag=FALSE)
  
  un_com <- clusters(net3)
  
  if (un_com$no>1){
    cpoint<-vector()
    for (i in 1: un_com$no){
      q_data <- dd[which(un_com$membership==i),]
      central <- colMeans(q_data)
      central <- matrix(central, nrow=1,ncol = dim(q_data)[2])
      kk <- nn2(q_data, query = central, k = 1, treetype = c("kd"), searchtype = c("standard"), eps = 0)
      cpoint[i] <- as.numeric(row.names(q_data)[kk$nn.idx])
    }  
    
    c_distmat <- as.matrix(dist(dd[cpoint,], diag = TRUE, upper = TRUE))
    aux <- mst(c_distmat)
    for (j in 1: length(aux$xy0)){
      em_mat[cpoint[aux$xy0[j]],cpoint[aux$xy1[j]]] <- c_distmat[aux$xy0[j],aux$xy1[j]] 
      em_mat[cpoint[aux$xy1[j]],cpoint[aux$xy0[j]]] <- c_distmat[aux$xy1[j],aux$xy0[j]] 
    }
    # ----------------
    net4<-graph.adjacency(em_mat, mode=c("undirected"), weighted=TRUE, diag=FALSE)  
    cal_dist_mat <- distances(net4)
  }
  else{
    cal_dist_mat <- distances(net3)
  }
}

# =============================================================================
center<-function(x){ x-mean(x) }

# =============================================================================
GetProjection<-function(x){
  print( sprintf( "r%d: c:%d", nrow(x) ,ncol(x) ) )
  return ( cmdscale(x, k = 1) );
}


# =============================================================================
SplitCriterion<-function(x,p,min_outlier){ 
	if(length(x$x)<3) {return("NoCluster") }
	cmins<-c()
	for(i in c(2:(length(x$x)-2)) ) {
		if( x$y[i-1]>= x$y[i]  & x$y[i]<=x$y[i+1] ) { cmins<-c(cmins,i) }
	}
	if(length(cmins)==0 ) { return("NoCluster") }
	no_points_split<-c()
	for( cm in cmins) { 
		no_points_split<-c(no_points_split,  min ( sum( p<x$x[cm]  ) , sum( p>=x$x[cm]  ) ));
	}
	cmins<-cmins[ no_points_split > min_outlier];
	if(length(cmins)==0 ) { return("NoCluster") } 

	selected_min<- which.min(x$y[cmins] )

	return( list( x=  (x$x[cmins])[ selected_min ] , y= min( (x$y[cmins])[selected_min])  )  );

}

# =============================================================================
ProcessClust<-function(x,cids,adjust,min_outlier,interpolation_points,SP){

	if( nrow(x) == 1  || is.vector(x)) {
		projections<-c(0);
	}
	print( sprintf(" dim: %d" ,dim(x[cids,]) ) ) ;
	projections<-GetProjection(x[cids,cids]);
	sr_projections <- projections[order(projections)]
	d<-c();
  if (SP == "DEN"){
	    if(length(projections) > 2) {
          d<-density(projections, bw="SJ",kernel = "rectangular" ,adjust=adjust, n=min(length(projections), interpolation_points))
          ord_mid <- sr_projections[-length(sr_projections)] + diff(sr_projections)/2
          newd <- approx(d$x,d$y,xout=ord_mid)
          d$x <- newd$x
          d$y <- newd$y
      } else {
        d<-density(projections, bw="SJ",kernel = "rectangular" ,adjust=adjust, n=min(length(projections), interpolation_points) )
      }
	    d$SplitPoint<-SplitCriterion(d,projections,min_outlier);
  }
	 # ------------------------------------------------------------------------------------------- 
	else if (SP == "MAR") {
	     if (length(sr_projections) > 2*min_outlier){
	       srange <- (min_outlier+1):(length(sr_projections)-min_outlier)
	       print(srange)
	       ss <- which.max(diff(sr_projections[srange]))
	       d$SplitPoint$x <- (sr_projections[srange][ss]+sr_projections[srange][ss+1])/2
	       d$SplitPoint$y <- 1/max(diff(sr_projections[srange]))
	       print("here")
	     } else {
	       d$SplitPoint[1] = "NoCluster"
	     }
	}
	  # ------------------------------------------------------------------------------------------- 
	 
	d$ids<-cids;
	if( d$SplitPoint[1]!="NoCluster") { 
		if( ( !any(projections< d$SplitPoint$x ) ) || (!any(projections>= d$SplitPoint$x) ) ) {
			d$SplitPoint<-"NoCluster";
		} else {
			d$kid1.ids<- d$ids[ projections< d$SplitPoint$x ] ;
			d$kid2.ids<- d$ids[ projections>= d$SplitPoint$x ] ;
		}
	}

	print( sprintf( " cids:%d left: %d right: %d , SplitPoint:%s ",   length(cids), length( d$kid1.ids) , length(d$kid2.ids) , d$SplitPoint[1] ) )
	return(d);
}

# =============================================================================

idivclu<-function(dd,k,adjust=0.5, min_outlier=round(max((dim(dd)[1]/k)/3,5)), interpolation_points = (dim(dd)[1]/2), knei=5,SP="DEN"){ 
	print(sprintf("min_outlier:%f",min_outlier))
  cdata <- cal_dist_mat(dd,knei)
  # ------------------------
  
	ids<-1:nrow(cdata);

	cres<-vector("list" , 1)
	max_depth<-1;
	cres[[1]]<-list(node = ProcessClust(cdata,ids,adjust,min_outlier,interpolation_points,SP) );

	cres[[1]]$node$depth<-0;
	cres[[1]]$node$NoNode<-1;
	cres[[1]]$node$res.id<-1;

	SplitPointsVals<-data.frame( SplitPoint=cres[[1]]$node$SplitPoint$y , id=1, NoPoints=length(ids));
	noClusters<-1;
	leafs<- data.frame ( id=1 , NoPoints=length(ids) );
	
	while( noClusters<k &  nrow(SplitPointsVals) > 0 ) {
		Cluster.To.Split<- SplitPointsVals$id[ which.min( SplitPointsVals$SplitPoint )  ][1] ;

		kid1<-vector("list" , 1); kid2<-vector("list" , 1)

		kid1[[1]]<- list(node = ProcessClust(cdata,  cres[[ Cluster.To.Split ]]$node$kid1.ids  ,adjust, min_outlier, interpolation_points,SP) );
		kid2[[1]]<- list(node = ProcessClust(cdata,  cres[[ Cluster.To.Split ]]$node$kid2.ids  ,adjust, min_outlier, interpolation_points,SP) );

        kid1[[1]]$node$depth <- cres[[ Cluster.To.Split ]]$node$depth+1
        kid2[[1]]$node$depth <- cres[[ Cluster.To.Split ]]$node$depth+1 

		kid1[[1]]$node$NoNode<-2*cres[[ Cluster.To.Split ]]$node$NoNode;
		kid2[[1]]$node$NoNode<-2*cres[[ Cluster.To.Split ]]$node$NoNode+1;

		kid1[[1]]$node$Splitted<-0;
		kid2[[1]]$node$Splitted<-0;
	
		kid1[[1]]$node$res.id<- length(cres)+1 ;
		kid2[[1]]$node$res.id<- length(cres)+2
		
		cres<- rbind(cres,kid1,kid2);
		
		if(  kid1[[1]]$node$SplitPoint[1] != "NoCluster" ) {
			SplitPointsVals<- rbind( SplitPointsVals ,  cbind( SplitPoint=kid1[[1]]$node$SplitPoint$y , id=length(cres)-1  ,  NoPoints=length( kid1[[1]]$node$ids) ) );
		}
		if(  kid2[[1]]$node$SplitPoint[1] != "NoCluster" ) {
			SplitPointsVals<- rbind( SplitPointsVals ,  cbind( SplitPoint=kid2[[1]]$node$SplitPoint$y , id=length(cres) , NoPoints=length( kid2[[1]]$node$ids) ) );
		}
	
		SplitPointsVals<-SplitPointsVals[ SplitPointsVals$id!=Cluster.To.Split,];

		cres[[ Cluster.To.Split ]]$node$kid1.id<- length(cres)-1
		cres[[ Cluster.To.Split ]]$node$kid2.id<- length(cres)
		cres[[ Cluster.To.Split ]]$node$Splitted<- 1 
		cleafs<- data.frame( id= c( length(cres) -1, length(cres) ) , NoPoints=c( length(kid1[[1]]$node$ids), length(kid2[[1]]$node$ids)  ) ) ;
		leafs<-rbind(leafs, cleafs );
		leafs<- leafs[ leafs$id!=Cluster.To.Split, ];

		noClusters<-noClusters+1;

		max_depth<- max( kid2[[1]]$node$depth , max_depth);
	}
	Clusters<-rep("Cluster" ,length(ids) )
	leafs$lab<-rep("",length(leafs$id))
	for( i in seq_along(leafs$id) ) { 
		Clusters[ cres[[ leafs$id[i] ]]$node$ids ]<- sprintf("Cluster%d",i ) 
		leafs$lab[ i]    <- sprintf("Cluster%d",i )

	}

	return( list(  Clusterids=Clusters,  leafs=leafs, max_depth=max_depth) ) ;
}
