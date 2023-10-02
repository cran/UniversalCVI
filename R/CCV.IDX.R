CCV.IDX <- function(x, cmax, cmin = 2, indexlist = "all", method = 'FCM', fzm = 2,
                    iter = 100, nstart = 20){
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
  if(!any(indexlist %in% c("all","CCVP", "CCVS")))
    stop("Argument 'indexlist' is not in 'all', 'CCVP', 'CCVS'")
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
  ccvp = vector()
  ccvs = vector()
  distance =dist(x,diag = TRUE,upper= TRUE)
  # FOR CCVP CCVS
  distc = as.vector(as.matrix(distance))
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
    uut = m%*%t(m)
    vnew = as.vector(1-(uut/max(uut)))

    if(sum(indexlist %in% c("all","CCVP"))>=1){
      ccvp[k-cmin+1] = cor(distc-mean(distc),vnew-mean(vnew),method = "pearson") #NW
    }
    if(sum(indexlist %in% c("all","CCVS"))>=1){
      ccvs[k-cmin+1] = cor(distc,vnew,method = "spearman") #NW
    }
  } # END CCVP CCVS index
  CCVP = data.frame(cbind("c"=cmin:cmax,"CCVP"=ccvp))
  CCVS = data.frame(cbind("c"=cmin:cmax,"CCVS"=ccvs))
  CCV.list = list("CCVP"= CCVP, "CCVS" = CCVS)
  if (sum(indexlist %in% "all")>=1){
    return(CCV.list)
  } else {
    return(CCV.list[indexlist])
  }
}
