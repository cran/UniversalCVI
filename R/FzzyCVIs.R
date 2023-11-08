FzzyCVIs <- function(x, cmax, cmin = 2, indexlist = 'all', corr = 'pearson',
                     method = 'FCM', fzm = 2,
                     gamma = (fzm^2*7)/4,
                     sampling = 1,
                     iter = 100,
                     nstart = 20,
                     NCstart = TRUE){
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
  if(!any(indexlist %in% c("all", "WPC", "WP", "WPCI1", "WPCI2", "XB", "KWON", "KWON2", "TANG", "HF", "WL", "PBM", "KPBM", "CCVP", "CCVS", "GC1", "GC2", "GC3", "GC4")))
    stop("Argument 'indexlist' is not in the possible value")
  if(!any(method  == c("FCM","EM")))
    stop("Argument 'method' should be one of 'FCM', 'EM' ")
  if(method == "FCM"){
    if(fzm <= 1)
      stop("Argument 'fcm' should be the number greater than 1",call. = FALSE)
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
    if(!is.numeric(iter))
      stop("Argument 'iter' must be numeric")
  }
  if(any(indexlist %in% c("all","WPC","WP","WPCI1","WPCI2"))){ # For WP index check
    if(!any(corr == c("pearson","kendall","spearman")))
      stop("Argument 'corr' should be one of 'pearson', 'kendall', 'spearman'")
    if(!is.logical(NCstart))
      stop("Argument 'NCstart' must be logical")
    if(!is.numeric(gamma))
      stop("Argument 'gamma' must be numeric or leave blank for default value")
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
  # Defined vector
  xb = vector()
  kwon = vector()
  kwon2 = vector()
  tang = vector()
  hf = vector()
  wl = vector()
  pbm = vector()
  kpbm = vector()
  ccvp = vector()
  ccvs = vector()
  crr = vector()
  WPI = vector()
  WPCI2 = vector()
  WPCI3 = vector()
  gc1 = vector()
  gc2 = vector()
  gc3 = vector()
  gc4 = vector()
  # Algorithm method
  if(sum(indexlist %in% c("all","WPC","WP","WPCI1","WPCI2","CCVP","CCVS","GC1","GC2","GC3","GC4"))>=1){
    distance =dist(x,diag = TRUE,upper= TRUE)
    # FOR WP idx (single distance)
    distx = as.vector(distance)
    # FOR CCVP CCVS
    distc = as.vector(as.matrix(distance))
    # GC distance matrix defined
    if(sum(indexlist %in% c("all","GC1","GC2","GC3","GC4"))>=1){ # FOR GC IDX
      dd = as.matrix(distance)
      dd[lower.tri(dd)] = 0
      diag(dd) = 0
      Check_GC1 = sum(indexlist %in% c("all","GC1"))>=1
      Check_GC2 = sum(indexlist %in% c("all","GC2"))>=1
      Check_GC3 = sum(indexlist %in% c("all","GC3"))>=1
      Check_GC4 = sum(indexlist %in% c("all","GC4"))>=1
    }
    # IF EM Algorithm
    if(method == "EM"){
      if(sum(indexlist %in% c("all","WPC","WP","WPCI1","WPCI2"))>=1){
        if(cmin<=2){
          dtom = sqrt(rowSums((x-colMeans(x))^2))
          if(NCstart){
            crr[1] = sd(dtom)/(max(dtom)-min(dtom))
          }else{
            crr[1] = 0
          }
        }else{
          EM.model <- Mclust(x,G=cmin-1,verbose = FALSE)
          xnew = ((EM.model$z^gamma)/rowSums(EM.model$z^gamma))%*%t(EM.model$parameters$mean)
          crr[1]= cor(distx,as.vector(dist(xnew)),method=corr)
        }
        EM.model <- Mclust(x,G = cmax+1,verbose = FALSE)
        xnew = ((EM.model$z^gamma)/rowSums(EM.model$z^gamma))%*%t(EM.model$parameters$mean)
        crr[cmax-cmin+3]= cor(distx,as.vector(dist(xnew)),method=corr)
      }
    }else if(method == "FCM"){# IF FCM Algorithm
      if (sum(indexlist %in% c("all","WPC","WP","WPCI1","WPCI2"))>=1){
        if(cmin<=2){
          dtom = sqrt(rowSums((x-colMeans(x))^2))
          if(NCstart){
            crr[1] = sd(dtom)/(max(dtom)-min(dtom))
          }else{
            crr[1] = 0
          }
        }else{
          wd = Inf
          for (nr in 1:nstart){
            minFCM.model = cmeans(x,cmax+1,iter,verbose=FALSE,method="cmeans",m=fzm)
            if (minFCM.model$withinerror < wd){
              wd = minFCM.model$withinerror
              minFCM.model2 =minFCM.model
            }
          }
          xnew = ((minFCM.model2$membership^gamma)/rowSums(minFCM.model2$membership^gamma)) %*% minFCM.model2$center
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
      }
    }
  } # END if first process defined
  # start k loop
  for(k in cmin:cmax){
    if(method == "EM"){ # EM Algorithm
      EM.model <- Mclust(x,G=k,verbose = FALSE)
      m <-  EM.model$z # m: membership
      c <-  t(EM.model$parameters$mean) # c: centroid
    }else if(method == "FCM"){ # FCM Algorithm
      wd = Inf
      for (nr in 1:nstart){
        FCM.model = cmeans(x,k,iter,verbose=FALSE,method="cmeans",m=fzm)
        if (FCM.model$withinerror < wd){
          wd = FCM.model$withinerror
          FCM.model2 =FCM.model
        }
      }
      m <- FCM.model2$membership
      c <- FCM.model2$centers
    }
    if(sum(indexlist %in% c("all","KWON","KWON2","TANG","HF","XB","WL","PBM","KPBM"))>=1){
      # Defined vector
      d1 = vector()
      d3 = vector()
      d4 = vector()
      d5 = vector()
      d6 = vector()
      d8 = vector()
      d9 = vector()
      d2 = rowSums((c-matrix(colMeans(x),k,ncol(x),byrow=T))^2) #NW
      n = nrow(x)
      d7 = sqrt(rowSums((x-matrix(colMeans(x),n,ncol(x),byrow=T))^2)) #NW
      w1 = (n-k+1)/n
      w2 = (k/(k-1))^sqrt(2)
      w3 = (n*k)/(n-k+1)^2
      for (j in 1:k){
        center = matrix(c[j,],n,ncol(x),byrow = T) #NW
        d1[j] = (m[,j])^2%*%rowSums((x-center)^2)
        d4[j] = ((m[,j])^(2^sqrt(fzm/2))%*%rowSums((x-center)^2))
        d5[j] = (m[,j])^fzm%*%rowSums((x-center)^2)
        d6[j] = sum((m[,j]))
        d8[j] = (m[,j])%*%sqrt(rowSums((x-center)^2))  #NW
      }

      s=1
      for(i in 1:(k-1)){
        for(j in (i+1):k){
          d3[s]=sum((c[i,]-c[j,])^2)
          d9[s]= sqrt(d3[s]) #NW
          s=s+1
        }
      }

      adp = (1/(k*(k-1)))*sum(d3)
      kwon[k-cmin+1] = (sum(d1)+ mean(d2))/min(d3) #NW
      kwon2[k-cmin+1] = w1*((w2*sum(d4)) + (sum(d2)/max(d2)) + w3) / (min(d3) + 1/k + 1/(k^fzm - 1))
      tang[k-cmin+1] = (sum(d1) + 2*adp)/ (min(d3) + 1/k) #NW
      hf[k-cmin+1] =  (sum(d5) + 2*adp)/ ((n/2*k)* (min(d3) + median(d3))) #NW
      xb[k-cmin+1] = (sum(d1))/(n*min(d3))
      wl[k-cmin+1] = sum(d1/d6) / (min(d3) + median(d3))
      pbm[k-cmin+1] = ((1/k)*(sum(d7)*max(d9)/sum(d8)))^2 #NW
      kpbm[k-cmin+1] = ((1/k)*(max(d9)/sum(d8)))^2 #NW
    } #end loop if for indexes based on compactness and seperated

    if (sum(indexlist %in% c("all","WPC","WP","WPCI1","WPCI2"))>=1){
      xnew = ((m^gamma)/rowSums(m^gamma))%*%c
      crr[k-cmin+2]= cor(distx,as.vector(dist(xnew)),method=corr)
    } # END WP index

    if (sum(indexlist %in% c("all","CCVP","CCVS"))>=1){
      uut = m%*%t(m)
      vnew = as.vector(1-(uut/max(uut)))
      ccvp[k-cmin+1] = cor(distc-mean(distc),vnew-mean(vnew),method = "pearson") #NW
      ccvs[k-cmin+1] = cor(distc,vnew,method = "spearman") #NW
    } # END CCVP CCVS index

    if (sum(indexlist %in% c("all","GC1","GC2","GC3","GC4"))>=1){
      n = nrow(x)
      mt = t(m)
      NI = colSums(m)
      nws = floor(sum(NI*(NI-1)/2))
      # GC1: sum product
      if(Check_GC1){
        PD1 = m%*%mt
        G1 = sum(PD1*dd)
        Gmax1 = sum(sort(dd,decreasing = T)[1:nws]*sort(PD1[lower.tri(PD1)],decreasing = T)[1:nws])
        Gmin1 = sum(sort(dd,decreasing = T)[1:nws]*sort(PD1[lower.tri(PD1)],decreasing = F)[1:nws])
        gc1[k-cmin+1] = (G1-Gmin1)/(Gmax1-Gmin1)
      }

      if(sum(indexlist %in% c("all","GC2","GC3","GC4"))>=1){

        PD2 = matrix(rep(0,n^2),n,n)
        PD3 = matrix(rep(0,n^2),n,n)
        PD4 = matrix(rep(0,n^2),n,n)

        # check
        for (s in 1:n){
          if(Check_GC2){
            #sum-min
            PD2[s,] =  colSums(pmin(mt,m[s,]))
          }
          if(Check_GC3){
            #max-product
            PD3[s,] = apply(m[s,]*mt,2,max)
          }
          if(Check_GC4){
            #max-min
            PD4[s,] = apply(pmin(mt,m[s,]),2,max)
          }
        }# end loop s

        if(Check_GC2){
          G2 = sum(PD2*dd)
          Gmax2 = sum(sort(dd,decreasing = T)[1:nws]*sort(PD2[lower.tri(PD2)],decreasing = T)[1:nws])
          Gmin2 = sum(sort(dd,decreasing = T)[1:nws]*sort(PD2[lower.tri(PD2)],decreasing = F)[1:nws])
          gc2[k-cmin+1] = (G2-Gmin2)/(Gmax2-Gmin2)
        }
        if(Check_GC3){
          G3 = sum(PD3*dd)
          Gmax3 = sum(sort(dd,decreasing = T)[1:nws]*sort(PD3[lower.tri(PD3)],decreasing = T)[1:nws])
          Gmin3 = sum(sort(dd,decreasing = T)[1:nws]*sort(PD3[lower.tri(PD3)],decreasing = F)[1:nws])
          gc3[k-cmin+1] = (G3-Gmin3)/(Gmax3-Gmin3)
        }
        if(Check_GC4){
          G4 = sum(PD4*dd)
          Gmax4 = sum(sort(dd,decreasing = T)[1:nws]*sort(PD4[lower.tri(PD4)],decreasing = T)[1:nws])
          Gmin4 = sum(sort(dd,decreasing = T)[1:nws]*sort(PD4[lower.tri(PD4)],decreasing = F)[1:nws])
          gc4[k-cmin+1] = (G4-Gmin4)/(Gmax4-Gmin4)
        }
      }
    } # END GC index
  }#end k loop
  # Create data frame for IDX result
  XB = data.frame(cbind("c"=cmin:cmax,"XB"=xb))
  KWON = data.frame(cbind("c"=cmin:cmax,"KWON"=kwon))
  KWON2 = data.frame(cbind("c"=cmin:cmax,"KWON2"=kwon2))
  TANG = data.frame(cbind("c"=cmin:cmax,"TANG"=tang))
  HF = data.frame(cbind("c"=cmin:cmax,"HF"=hf))
  WL = data.frame(cbind("c"=cmin:cmax,"WL"=wl))
  PBM = data.frame(cbind("c"=cmin:cmax,"PBM"=pbm))
  KPBM = data.frame(cbind("c"=cmin:cmax,"KPBM"=kpbm))
  CCVP = data.frame(cbind("c"=cmin:cmax,"CCVP"=ccvp))
  CCVS = data.frame(cbind("c"=cmin:cmax,"CCVS"=ccvs))
  GC1 = data.frame(cbind("c"=cmin:cmax,"GC1"=gc1))
  GC2 = data.frame(cbind("c"=cmin:cmax,"GC2"=gc2))
  GC3 = data.frame(cbind("c"=cmin:cmax,"GC3"=gc3))
  GC4 = data.frame(cbind("c"=cmin:cmax,"GC4"=gc4))
  if (sum(indexlist %in% c("all","WPC","WP","WPCI1","WPCI2"))>=1){
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
    WPI = data.frame(cbind("c"= cmin:cmax, "WPI1" = WPI))
    WPCI2 = data.frame(cbind("c"=cmin:cmax,"WPCI2"=WPCI2))
    WPCI3 = data.frame(cbind("c"=cmin:cmax,"WPI"=WPCI3))
    crr = data.frame(cbind("c" = (cmin-1):(cmax+1), "WPC" = crr))
  }else{
    crr = 0
    WPI = 0
    WPCI2 = 0
    WPCI3 = 0
  }
  IDX.list = list("WPC" = crr,"WP"= WPCI3,"WPCI1" = WPI, "WPCI2" = WPCI2,"XB" = XB, "KWON" = KWON, "KWON2" = KWON2,"TANG" = TANG,"HF"=HF,"WL" = WL,"PBM" =PBM,
                  "KPBM" = KPBM, "CCVP"= CCVP, "CCVS" = CCVS, "GC1"= GC1,"GC2"= GC2,"GC3"= GC3,"GC4"= GC4)
  if (sum(indexlist == "all")==1){
    return(IDX.list)
  } else {
    return(IDX.list[indexlist])
  }
}
