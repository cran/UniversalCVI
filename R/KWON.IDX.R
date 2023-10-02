KWON.IDX <- function(x, cmax, cmin = 2, method = "FCM",
                     fzm = 2, nstart = 20, iter = 100){
  if(missing(x))
    stop("Missing input argument. A numeric data frame or matrix is required")
  if(missing(cmax))
    stop("Missing input argument. A maximum number of clusters  is required")
  if(!is.numeric(cmax))
    stop("Argument 'cmax' must be numeric")
  if(cmax > nrow(x))
    stop("The maximum number of clusters for consideration should be less than or equal to the number of data points in dataset.")
  if(!is.numeric(cmin))
    stop("Argument 'cmin' must be numeric")
  if(cmin <=1)
    warning("The minimum number of clusters for consideration should be more than 1",immediate. = TRUE)
  if(!any(method  == c("FCM","EM")))
    stop("Argument 'method' should be one of 'FCM','EM' ")
  if(method == "FCM"){
    if(fzm <= 1)
      stop("Argument 'fcm' should be the number greater than 1",call. = FALSE)
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
    if(!is.numeric(iter))
      stop("Argument 'iter' must be numeric")
  }
  # Defined vector
  kwon = vector()
  # start k loop

  for(k in cmin:cmax){
    if(method == "EM"){ # EM Algorithm
      EM.model <- Mclust(x,G=k,verbose=FALSE)
      assign("m",EM.model$z)
      assign("c",t(EM.model$parameters$mean))
    }else if(method == "FCM"){ # FCM Algorithm
      wd = Inf
      # cm.out = list()
      for (nr in 1:nstart){
        FCM.model = cmeans(x,k,iter,verbose=FALSE,method="cmeans",m=fzm)
        if (FCM.model$withinerror < wd){
          wd = FCM.model$withinerror
          FCM.model2 =FCM.model
        }
      }
      assign("m",FCM.model2$membership)
      assign("c",FCM.model2$centers)
    }
    d1= vector()
    d2 = rowSums((c-matrix(colMeans(x),k,ncol(x),byrow=T))^2) #NW
    d3 = vector()
    n = nrow(x)
    for (j in 1:k){
      center = matrix(c[j,],n,ncol(x),byrow = T) #NW
      d1[j] = (m[,j])^2%*%rowSums((x-center)^2)
    }

    s=1
    for(i in 1:(k-1)){
      for(j in (i+1):k){
        d3[s]=sum((c[i,]-c[j,])^2)
        s=s+1
      }
    }
    kwon[k-cmin+1] = (sum(d1)+ mean(d2))/min(d3)
  }
  KWON.data = data.frame(cbind("c"=cmin:cmax,"KWON"=kwon))
  return(KWON.data)
}
