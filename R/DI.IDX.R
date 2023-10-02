DI.IDX <- function(x, kmax, kmin=2,
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
  di = vector()
  for(k in kmin:kmax){
    if(method == "kmeans"){
      K.model = kmeans(x,k,nstart =nstart)
      cluss = K.model$cluster
    } else if(startsWith(method,"hclust_")){
      cluss = cutree(H.model,k)
    } # End check algorithm
    if(!all(seq(k) %in% unique(cluss)))
      warning("Some clusters are empty.")
    size = table(cluss)
    dunn.dem = 0
    dunn.num = 1e10
    for (i in 1:(k-1)){
      dunn.dem = max(dunn.dem,max(dist(x[cluss==i,])))
      for(j in (i+1):k){
        dunn.num = min(dunn.num,min(as.matrix(dist(rbind(x[cluss==i,],x[cluss==j,])))[(size[i]+size[j]):(size[i]+1),1:size[i]]))
      }
    }
    dunn.dem = max(dunn.dem,max(dist(x[cluss==k,])))
    di[k-kmin+1] = dunn.num / dunn.dem
  }
  DI.data = data.frame(cbind("k"=kmin:kmax,"DI"=di))
  return(DI.data)
}
