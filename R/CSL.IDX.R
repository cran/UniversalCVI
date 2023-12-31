CSL.IDX <- function(x, kmax, kmin=2,
                      method = 'kmeans',
                      nstart = 100){
  if(missing(x))
    stop("Missing input argument. A numeric data frame or matrix is required")
  if(missing(kmax))
    stop("Missing input argument. A maximum number of clusters is required")
  if(!is.numeric(kmax))
    stop("Argument 'kmax' must be numeric")
  if(kmax > nrow(x))
    stop("The maximum number of clusters for consideration should be less than or equal to the number of data points in dataset.")
  if(!is.numeric(kmin))
    stop("Argument 'kmin' must be numeric")
  if(kmin <=1)
    warning("The minimum number of clusters for consideration should be more than 1",immediate. = TRUE)
  if(!any(method  == c("kmeans","hclust_complete","hclust_average","hclust_single")))
    stop("Argument 'method' should be one of 'kmeans', 'hclust_complete', 'hclust_average', 'hclust_single'")
  if(method == "kmeans"){
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
  }
  if(startsWith(method,"hclust_")){
    H.model = hclust(dist(x),method = sub("hclust_", "", method))
  }
  dm = dim(x)
  csl = rep(0,kmax-kmin+1)
  for(k in kmin:kmax){
    centroid = matrix(0,k,dm[2])
    if(method == "kmeans"){
      K.model = kmeans(x,k,nstart =nstart)
      cluss = K.model$cluster
      centroid = K.model$centers
      for(j in 1:k){
        dd = as.matrix(dist(x[cluss==j,]))
        csl[k-kmin+1] = csl[k-kmin+1] + sum(apply(dd,2,max))/sum(cluss==j)
      }
    } else if(startsWith(method,"hclust_")){
      cluss = cutree(H.model,k)
      for (j in 1:k){
        if (is.null(nrow(x[cluss==j,])) | sum(nrow(x[cluss==j,]))==1){
          centroid[j,] = as.numeric(x[cluss==j,])
        } else {
          centroid[j,] = colMeans(x[cluss==j,])
        }
        dd = as.matrix(dist(x[cluss==j,]))
        csl[k-kmin+1] = csl[k-kmin+1] + sum(apply(dd,2,max))/sum(cluss==j)
      }
    }
    if(!all(seq(k) %in% unique(cluss)))
      warning("Some clusters are empty.")
    dd2 = as.matrix(dist(centroid))
    dd2 = matrix(dd2[dd2 > 0],k-1,k)
    csl[k-kmin+1] = csl[k-kmin+1]/sum(apply(dd2,2,min))
  }
  CSL.data = data.frame(cbind("k"=kmin:kmax,"CSL"=csl))
  return(CSL.data)
}
