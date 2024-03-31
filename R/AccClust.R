AccClust <- function(x, label.names = "label",algorithm = "FCM", fzm = 2, scale = TRUE,
                     nstart = 100,iter = 100){
  if(missing(x))
    stop("Missing input argument. A numeric data frame or matrix is required")
  if(!is.character(label.names))
    stop("Argument 'label.names' must be a character string indicating the true label column name")
  if(!(label.names %in% colnames(x)))
    stop(paste("Column", label.names , "is not in the data frame or matrix"))
  if(!any(algorithm  %in% c("FCM", "EM", "Kmeans")))
    stop("Argument 'algorithm' should be 'FCM', 'EM' and 'Kmeans'")
  if(!is.logical(scale))
    stop("Argument 'scale' must be a logical")
  if(any(algorithm  %in% "FCM")){
    if(fzm <= 1)
      stop("Argument 'fcm' should be the number greater than 1",call. = FALSE)
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
    if(!is.numeric(iter))
      stop("Argument 'iter' must be numeric")
  }
  wd = 100000
  x = data.frame(x)
  x[,label.names] = as.numeric(as.factor(x[,label.names])) #new order label
  k = length(unique(x[,label.names])) # True number of cluster
  tot.alg = length(algorithm) # Number of algorithms
  u = order(table(x[,label.names]),decreasing = T) # Order labels
  result = data.frame() # For recording result

  # scale option : TRUE / FALSE
  if(scale){
    alg.data = scale(x[,!(names(x)) == label.names])
  }else{
    alg.data = x[,!(names(x)) == label.names]
  }

  # Algorithm selected : Kmeans, EM and FCM
  for (n in 1:tot.alg) {
    result[n,"ALGORITHM"] = algorithm[n]

    if(algorithm[n] == "Kmeans"){
      alg.result = kmeans(alg.data,k,nstart = nstart)
      assign(paste0(algorithm[n],".cluster"),alg.result$cluster)
    }else if(algorithm[n] == "EM"){
      alg.result = Mclust(alg.data,G=k,verbose = FALSE)
      assign(paste0(algorithm[n],".cluster"),alg.result$classification)
    }else if(algorithm[n] == "FCM"){

      for (nr in 1:nstart){
        alg.result = cmeans(alg.data,k, iter.max = iter, dist = "euclidean", m = fzm)
        assign("wd.error",alg.result$withinerror)
        if (wd.error < wd){
          wd = wd.error
          assign(paste0(algorithm[n],".cluster"),alg.result$cluster)
        }
      }
    }
    cluster = get(paste0(algorithm[n],".cluster"))+100
    cluster2 = cluster
    wold = NULL
    for(i in 1:k){
      ttt = sort(table(cluster[x[,label.names]==u[i]]),decreasing = T)
      for (r in 1:length(ttt)){
        w = as.numeric(names(ttt[r]))
        if (!any(wold==w)){
          break
        }
      }
      cluster2[cluster2==w]=u[i]
      wold = c(wold,w)
    }
    result[n,"ACC"] = mean(cluster2==x[,label.names])
  }
  return(result)
}






