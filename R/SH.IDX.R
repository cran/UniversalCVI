SH.IDX = function (x, kmax, kmin = 2, method = "kmeans", nstart = 100)
{
  if (missing(x))
    stop("Missing input argument. A numeric data frame or matrix is required")
  if (missing(kmax))
    stop("Missing input argument. A maximum number of clusters is required")
  if (!is.numeric(kmax))
    stop("Argument 'kmax' must be numeric")
  if (kmax > nrow(x))
    stop("The maximum number of clusters for consideration should be less than or equal to the number of data points in dataset.")
  if (!is.numeric(kmin))
    stop("Argument 'kmin' must be numeric")
  if (kmin <= 1)
    warning("The minimum number of clusters for consideration should be more than 1",
            immediate. = TRUE)
  if (!any(method == c("kmeans", "hclust_ward.D", "hclust_ward.D2", "hclust_complete", "hclust_average", 
                       "hclust_single"))) 
    stop("Argument 'method' should be one of 'kmeans', 'hclust_ward.D', 'hclust_ward.D2', 'hclust_complete', 'hclust_average', 'hclust_single'")
  if (method == "kmeans") {
    if (!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
  }
  if (startsWith(method, "hclust_")) {
    H.model = hclust(dist(x), method = sub("hclust_", "",
                                           method))
  }
  sc = vector()
  n = nrow(x)
  for (k in kmin:kmax) {
    if (method == "kmeans") {
      K.model = kmeans(x, k, nstart = nstart)
      cluss = K.model$cluster
    }
    else if (startsWith(method, "hclust_")) {
      cluss = cutree(H.model, k)
    }
    if (!all(seq(k) %in% unique(cluss)))
      warning("Some clusters are empty.")
    size = table(cluss)
    sss = 1:n
    s=1
    for (i in 1:k) {
      mm = as.matrix(dist(x[cluss == i,]))
      ll = as.numeric(size[i])
      if (ll > 1){
        for (u in 1:ll){
          b = Inf
          a = mean(mm[u,-u])
          for (ii in setdiff(1:k,i)){
            b = min(b,mean(as.matrix(dist(rbind(x[cluss == i,][u,], x[cluss == ii,])))[1,-1]))
          }
          sss[s] = (b-a)/max(a,b)
          s= s+1
        }
      }else{
        sss[s] = 0
        s=s+1
      }
    }
    sc[k - kmin + 1] = mean(sss)
  }
  SH.data = data.frame(cbind(k = kmin:kmax, SH = sc))
  return(SH.data)
}
