CH.IDX <- function(x, kmax, kmin=2,
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
  if (!any(method == c("kmeans", "hclust_ward.D", "hclust_ward.D2", "hclust_complete", "hclust_average", 
                       "hclust_single"))) 
    stop("Argument 'method' should be one of 'kmeans', 'hclust_ward.D', 'hclust_ward.D2', 'hclust_complete', 'hclust_average', 'hclust_single'")
  if(method == "kmeans"){
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
  }
  if(startsWith(method,"hclust_")){
    H.model = hclust(dist(x),method = sub("hclust_", "", method))
  }
  dm = dim(x)
  ch = vector()
  for(k in kmin:kmax){
    xnew = matrix(0,dm[1],dm[2])
    centroid = matrix(0,k,dm[2])
    if(method == "kmeans"){
      K.model = kmeans(x,k,nstart =nstart)
      cluss = K.model$cluster
      centroid = K.model$centers
      xnew = centroid[cluss,]
    } else if(startsWith(method,"hclust_")){
      cluss = cutree(H.model,k)
      for (j in 1:k){
        if (is.null(nrow(x[cluss==j,])) | sum(nrow(x[cluss==j,]))==1){
          centroid[j,] = as.numeric(x[cluss==j,])
        } else {
          centroid[j,] = colMeans(x[cluss==j,])
        }
      }
      xnew = centroid[cluss,]
    } # End check algorithm
    if(!all(seq(k) %in% unique(cluss)))
      warning("Some clusters are empty.")
    d.cen = rowSums((x - xnew)^2)
    num = sapply(1:k, function(i) {
      ck = x[cluss == i, ]
      cen.k = centroid[i, ]
      n.clust_k = nrow(ck)
      n.clust_k * sum((cen.k - colMeans(x))^2)
    })
    num = as.numeric(num)
    dem = sapply(1:k, function(i) sum(d.cen[cluss == i]))
    ch[k-kmin+1] = ((nrow(x) - k) / (k-1)) * (sum(num[!is.na(num)]) / sum(dem))
  }
  CH.data = data.frame(cbind("k"=kmin:kmax,"CH"=ch))
  return(CH.data)
}
