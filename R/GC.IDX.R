GC.IDX <- function(x, cmax, cmin = 2, indexlist = "all", method = 'FCM', fzm = 2,
                   iter = 100, nstart = 20){

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
  if(!any(indexlist %in% c("all", "GC1", "GC2", "GC3", "GC4")))
    stop("Argument 'indexlist' is not in 'all', 'GC1', 'GC2', 'GC3', 'GC4'")
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
  gc1 = vector()
  gc2 = vector()
  gc3 = vector()
  gc4 = vector()
  distance =dist(x,diag = TRUE,upper= TRUE)
  dd = as.matrix(distance)
  dd[lower.tri(dd)] = 0
  diag(dd) = 0
  # For checking IDX
  Check_GC1 = sum(indexlist %in% c("all","GC1"))>=1
  Check_GC2 = sum(indexlist %in% c("all","GC2"))>=1
  Check_GC3 = sum(indexlist %in% c("all","GC3"))>=1
  Check_GC4 = sum(indexlist %in% c("all","GC4"))>=1

  # Start k loop
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
    # Defined variable
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
      for (s in 1:n){
        #GC2: sum-min
        if(Check_GC2){
          PD2[s,] =  colSums(pmin(mt,m[s,]))
        }
        # GC3: max-product
        if(Check_GC3){
          PD3[s,] = apply(m[s,]*mt,2,max)
        }
        # GC4: max-min
        if(Check_GC4){
          PD4[s,] = apply(pmin(mt,m[s,]),2,max)
        }
      } # end loop s
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
    }# END GC index
  }#end k loop
  # Create data frame for IDX result
  GC1 = data.frame(cbind("c"=cmin:cmax,"GC1"=gc1))
  GC2 = data.frame(cbind("c"=cmin:cmax,"GC2"=gc2))
  GC3 = data.frame(cbind("c"=cmin:cmax,"GC3"=gc3))
  GC4 = data.frame(cbind("c"=cmin:cmax,"GC4"=gc4))
  GC.list = list("GC1"= GC1,"GC2"= GC2,"GC3"= GC3,"GC4"= GC4)
  if (sum(indexlist == "all")==1){
    return(GC.list)
  } else {
    return(GC.list[indexlist])
  }
}
