Wvalid <- function(x, kmax, kmin=2, method='kmeans',
                    corr='pearson', nstart=100, NCstart = TRUE) {
  if(missing(x))
    stop("Missing input argument. A numeric data frame or matrix is required")
  if(missing(kmax))
    stop("Missing input argument. A maximum number of clusters is required")
  if(!is.numeric(kmax))
    stop("Argument 'kmax' must be numeric")
  if(kmax > nrow(x))
    stop("The maximum number of clusters for consideration should be less than or equal to the number of data points in dataset.")
  if(!any(method  == c("kmeans", "hclust_complete", "hclust_average", "hclust_single")))
    stop("Argument 'method' should be one of 'kmeans', 'hclust_complete', 'hclust_average' or 'hclust_single'")
  if(!any(corr == c("pearson","kendall","spearman")))
    stop("Argument 'corr' should be one of 'pearson', 'kendall', 'spearman'")
  if(method == "kmeans"){
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
  }
  if(!is.logical(NCstart))
    stop("Argument 'NCstart' must be logical")
  dm=dim(x)
  d = as.vector(dist(x))
  crr = rep(0,kmax-kmin+2)
  if (NCstart) {
    dtom = sqrt(rowSums((x-colMeans(x))^2))
    crr[1] = sd(dtom)/(max(dtom)-min(dtom))
  }
  if(startsWith(method,"hclust_")){ # Check  hierarchical method
    H.model = hclust(dist(x),method = sub("hclust_", "", method))
  }
  if (kmin == 2){
    lb = 2
  } else {
    lb = kmin-1
  }
  for (k in lb:(kmax+1)){
    xnew = matrix(0,dm[1],dm[2])
    centroid = matrix(0,k,dm[2])
    if (method == 'kmeans'){
      K.model=kmeans(x,k,nstart =nstart)
      cluss = K.model$cluster
      xnew = K.model$centers[cluss,]
    }
    else if(startsWith(method,"hclust_")){
      cluss = cutree(H.model,k)
      for (j in 1:k){
        if (is.null(nrow(x[cluss==j,])) | sum(nrow(x[cluss==j,]))==1){
          centroid[j,] = as.numeric(x[cluss==j,])
        } else {
          centroid[j,] = colMeans(x[cluss==j,])
        }
      }
      xnew = centroid[cluss,]
    }
    if(!all(seq(k) %in% unique(cluss)))
      warning("Some clusters are empty.")
    d2 = as.vector(dist(xnew))
    crr[k-kmin+2]= cor(d,d2,method=corr)
  }
  K = length(crr)
  NWI = ((crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)]))/pmax(0,(crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)]))
  NWI2 = (crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)])-(crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)])
  NWI3 = NWI
  if (max(NWI)<Inf){
    if (min(NWI) == -Inf){
      NWI3[NWI==-Inf] = min(NWI[is.finite(NWI)])
    }
  }
  if (max(NWI)==Inf){
    NWI3[NWI==Inf] = max(NWI[is.finite(NWI)])+NWI2[NWI==Inf]
    NWI3[NWI<Inf] = NWI[NWI<Inf] + NWI2[NWI<Inf]
    if (min(NWI) == -Inf){
      NWI3[NWI==-Inf] = min(NWI[is.finite(NWI)])+NWI2[NWI==-Inf]
    }
  }
  NWI = data.frame(cbind("k"= kmin:kmax, "NCI1" = NWI))
  NWI2 = data.frame(cbind("k"=kmin:kmax,"NCI2"=NWI2))
  NWI3 = data.frame(cbind("k"=kmin:kmax,"NCI"=NWI3))
  crr = data.frame(cbind("k" = (kmin-1):(kmax+1), "NC" = crr))
  my_list <- list("NC" = crr, "NCI" = NWI3, "NCI1" = NWI, "NCI2" = NWI2)
  return(my_list)
}
