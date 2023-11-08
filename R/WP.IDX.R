WP.IDX <- function(x, cmax, cmin = 2, corr = 'pearson',
                   method = 'FCM', fzm = 2,
                   gamma = (fzm^2*7)/4,
                   sampling = 1,
                   iter = 100,
                   nstart = 20,
                   NCstart = TRUE){
  if(missing(x))
    stop("Missing input argument. A numeric data frame or matrix is required")
  if(missing(cmax))
    stop("Missing input argument. A maximum number of clusters is required")
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
  if(!any(corr == c("pearson","kendall","spearman")))
    stop("Argument 'corr' should be one of 'pearson', 'kendall', 'spearman'")
  if(!is.logical(NCstart))
    stop("Argument 'NCstart' must be logical")
  if(!is.numeric(gamma))
    stop("Argument 'gamma' must be numeric or leave blank for default value")
  if(method == "FCM"){
    if(fzm <= 1)
      stop("Argument 'fcm' should be the number greater than 1",call. = FALSE)
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
    if(!is.numeric(iter))
      stop("Argument 'iter' must be numeric")
  }
  if(!is.numeric(sampling))
    stop("Argument 'sampling' must be numeric")
  if(!(sampling > 0 & sampling <= 1))
    stop("'sampling' must be greater than 0 and less than or equal to 1")
  if(sampling == 1){
    x = x
  }else {
    sample = sample(1:(nrow(x)),ceiling(nrow(x)*sampling),replace = FALSE)
    x = x[sample,]
  }
  crr = vector()
  WPI = vector()
  WPCI2 = vector()
  WPCI3 = vector()
  # Distance part
  distance = dist(x,diag = TRUE,upper= TRUE)
  # FOR WP idx (single distance)
  distx = as.vector(distance)
  # Algorithm method
  if(method == "EM"){
    if(NCstart & (cmin<=2)){
      dtom = sqrt(rowSums((x-colMeans(x))^2))
      crr[1] = sd(dtom)/(max(dtom)-min(dtom))
    } else{
      EM.model <- Mclust(x,G=cmin-1,verbose = FALSE)
      xnew = ((EM.model$z^gamma)/rowSums(EM.model$z^gamma))%*%t(EM.model$parameters$mean)
      crr[1]= cor(distx,as.vector(dist(xnew)),method=corr)
    }
    EM.model <- Mclust(x,G = cmax+1,verbose = FALSE)
    xnew = ((EM.model$z^gamma)/rowSums(EM.model$z^gamma))%*%t(EM.model$parameters$mean)
    crr[cmax-cmin+3]= cor(distx,as.vector(dist(xnew)),method=corr)
  }else if(method == "FCM"){
    if(NCstart & cmin<=2){
      dtom = sqrt(rowSums((x-colMeans(x))^2))
      crr[1] = sd(dtom)/(max(dtom)-min(dtom))
    }else{
      wd = Inf
      for (nr in 1:nstart){
        minFCM.model = cmeans(x,cmax+1,iter,verbose=FALSE,method="cmeans",m=fzm)
        if (minFCM.model$withinerror < wd){
          wd = minFCM.model$withinerror
          minFCM.model2 =minFCM.model
        }
      }
      xnew = ((minFCM.model2$membership^gamma)/rowSums(minFCM.model2$membership^gamma))%*%minFCM.model2$center
      crr[1]= cor(distx,as.vector(dist(xnew)),method=corr)
    }
    # crr cmax + 1
    wd = Inf
    for (nr in 1:nstart){
      maxFCM.model = cmeans(x,cmax+1,iter,verbose=FALSE,method="cmeans",m=fzm)
      if (maxFCM.model$withinerror < wd){
        wd = maxFCM.model$withinerror
        maxFCM.model2 = maxFCM.model
      }
    }
    xnew = ((maxFCM.model2$membership^gamma)/rowSums(maxFCM.model2$membership^gamma))%*%maxFCM.model2$center
    crr[cmax-cmin+3]= cor(distx,as.vector(dist(xnew)),method=corr)
  } # END if first process defined

  # start k loop
  for(k in cmin:cmax){
    if(method == "EM"){ # EM Algorithm
      EM.model <- Mclust(x,G=k,verbose = FALSE)
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

    xnew = ((m^gamma)/rowSums(m^gamma))%*%c
    crr[k-cmin+2]= cor(distx,as.vector(dist(xnew)),method=corr)
  }

  K = length(crr)
  WPI = ((crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)]))/pmax(0,(crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)]))
  WPCI2 = (crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)])-(crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)])
  WPCI3 = WPI
  if(sum(is.finite(WPI))==0){
    WPCI3[WPI==Inf] = WPCI2[WPI==Inf]
    WPCI3[WPI==-Inf] = pmin(0,WPCI2[WPI==-Inf])
  }else{
    if (max(WPI)<Inf){
      if (min(WPI) == -Inf){
        WPCI3[WPI==-Inf] = min(WPI[is.finite(WPI)])
      }
    }
    if (max(WPI)==Inf){
      WPCI3[WPI==Inf] = max(WPI[is.finite(WPI)])+WPCI2[WPI==Inf]
      WPCI3[WPI<Inf] = WPI[WPI<Inf] + WPCI2[WPI<Inf]  #added
      if (min(WPI) == -Inf){
        WPCI3[WPI==-Inf] = min(WPI[is.finite(WPI)])+WPCI2[WPI==-Inf]
      }
    }
  }
  # Data frame for result
  c = cmin:cmax
  WPI = data.frame(cbind("c"= c, "WPI1" = WPI))
  WPCI2 = data.frame(cbind("c"=c,"WPCI2"=WPCI2))
  WPCI3 = data.frame(cbind("c"=c,"WPI"=WPCI3))
  crr = data.frame(cbind("c" = (cmin-1):(cmax+1), "NC" = crr))

  WP.list = list("WPC" = crr, "WP" = WPCI3, "WPCI1" = WPI ,"WPCI2" = WPCI2)
  return(WP.list)
}

