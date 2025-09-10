Hvalid <-function(x, kmax, kmin=2,
                  indexlist = 'all',
                  method = 'kmeans',
                  p = 2, q = 2, # For DB ans DBs index
                  corr = 'pearson',
                  nstart = 100,
                  sampling = 1,
                  NCstart = TRUE){
  if(missing(x))
    stop("Missing input argument. A numeric data frame or matrix is required")
  if(missing(kmax))
    stop("Missing input argument. A maximum number of clusters  is required")
  if(!is.numeric(kmax))
    stop("Argument 'kmax' must be numeric")
  if(kmax > nrow(x))
    stop("The maximum number of clusters for consideration should be less than or equal to the number of data points in dataset.")
  if(!is.numeric(kmin))
    stop("Argument 'kmin' must be numeric")
  if(kmin <=1)
    warning("The minimum number of clusters for consideration should be more than 1",immediate. = TRUE)
  if(!any(indexlist %in% c("all", "PB","CH","CSL","DB","DBs","SF","DI","NC","NCI","NCI1","NCI2","NCI3","STR","PBM")))
    stop("Argument 'indexlist' is not in the possible value")
  if (!any(method == c("kmeans", "hclust_ward.D", "hclust_ward.D2", "hclust_complete", "hclust_average", 
                       "hclust_single"))) 
    stop("Argument 'method' should be one of 'kmeans', 'hclust_ward.D', 'hclust_ward.D2', 'hclust_complete', 'hclust_average', 'hclust_single'")
  if(method == "kmeans"){
    if(!is.numeric(nstart))
      stop("Argument 'nstart' must be numeric")
  }
  if(any(indexlist %in% c("all","DB","DBs"))){
    if(!is.numeric(p))
      stop("Argument 'p' must be numeric")
    if(!is.numeric(q))
      stop("Argument 'q' must be numeric")
  }
  if(any(indexlist %in% c("all","NC", "NCI", "NCI1", "NCI2", "PB","STR"))){
    if(!any(corr == c("pearson","kendall","spearman")))
      stop("Argument 'corr' should be one of 'pearson', 'kendall', 'spearman'")
    if(!is.logical(NCstart))
      stop("Argument 'NCstart' must be logical")
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
    dis = dist(x) #Distance of data
    d = as.vector(dis)
    dtom = sqrt(rowSums((x-colMeans(x))^2))
    E0 = sum(dtom)
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
  dm=dim(x) # Dimension of data
  # Defined vector
  pb = vector()
  ch = vector()
  db = vector()
  dbs = vector()
  sf = vector()
  di = vector()
  crr= vector()
  str = rep(0,kmax-kmin+1)
  pbm = rep(0,kmax-kmin+1)
  EK = rep(0,kmax-kmin+2)
  DK = rep(0,kmax-kmin+3)
  md = rep(0,kmax-kmin+3)
  csl = rep(0,kmax-kmin+1) #except
  NCI1 = vector()
  NCI2 = vector()
  NCI3 = vector()
  # Check  hierarchical method
  if(startsWith(method,"hclust_")){
    H.model = hclust(dist(x),method = sub("hclust_", "", method))
  }
  # Check for "NC", "NCI", "NCI1", "NCI2","STR","PBM" index
  if(any(indexlist %in% c("all","NC", "NCI", "NCI1", "NCI2", "STR", "PBM"))){
    if(kmin<=2){
      if(any(indexlist %in% c("all","NC", "NCI", "NCI1", "NCI2"))){
        if(NCstart){
          crr[1] = sd(dtom)/(max(dtom)-min(dtom))
        }else{
          crr[1] = 0
        }
      }
      if(any(indexlist %in% c("all", "STR", "PBM"))){
        EK[1] = E0
      }
    }else{
      if(method == "kmeans"){
        K.model = kmeans(x,kmin-1,nstart = nstart)
        cluss = K.model$cluster
        xnew = K.model$centers[cluss,]
        if(any(indexlist %in% c("all","NC", "NCI", "NCI1", "NCI2"))){
          crr[1]= cor(d,as.vector(dist(xnew)),method=corr)
        }
        if(any(indexlist %in% c("all","STR", "PBM"))){ #For str and pbm index
          EK[1] = sum(sqrt(rowSums((x - xnew)^2)))
        }
      }else if(startsWith(method,"hclust_")){
        cluss = cutree(H.model,kmin-1)
        centroid = matrix(0,(kmin-1),dm[2])
        for (j in 1:(kmin-1)) {
          if (is.null(nrow(x[cluss==j,])) | sum(nrow(x[cluss==j,]))==1){
            centroid[j,] = as.numeric(x[cluss==j,])
          } else {
            centroid[j,] = colMeans(x[cluss==j,])
          }
        }
        xnew = centroid[cluss,]
        if(any(indexlist %in% c("all","NC", "NCI", "NCI1", "NCI2"))){
          crr[1]= cor(d,as.vector(dist(xnew)),method=corr)
        }
        if(any(indexlist %in% c("all","STR", "PBM"))){ #For str and pbm index
          EK[1] = sum(sqrt(rowSums((x - xnew)^2)))
        }
      }
    }
    if(method == "kmeans"){
      K.model <- kmeans(x,kmax+1,nstart = nstart)
      xnew = K.model$centers[K.model$cluster,]
      if(any(indexlist %in% c("all","NC", "NCI", "NCI1", "NCI2"))){
        crr[kmax-kmin+3]= cor(d,as.vector(dist(xnew)),method=corr)
      }
      if(any(indexlist %in% c("all","STR", "PBM"))){
        ddd = dist(K.model$centers)
        md[kmax-kmin+3] = max(ddd)
        DK[kmax-kmin+3] = max(ddd)/min(ddd)
      }
    }else if(startsWith(method,"hclust_")){
      cluss = cutree(H.model,kmax+1)
      xnew = matrix(0,dm[1],dm[2])
      centroid = matrix(0,(kmax+1),dm[2])
      for (j in 1:(kmax+1)) {
        if (is.null(nrow(x[cluss==j,])) | sum(nrow(x[cluss==j,]))==1){
          centroid[j,] = as.numeric(x[cluss==j,])
        } else {
          centroid[j,] = colMeans(x[cluss==j,])
        }
      }
      xnew = centroid[cluss,]
      if(any(indexlist %in% c("all","NC", "NCI", "NCI1", "NCI2"))){
        crr[kmax-kmin+3]= cor(d,as.vector(dist(xnew)),method=corr)
      }
      if(any(indexlist %in% c("all","STR", "PBM"))){
        ddd = dist(centroid)
        md[kmax-kmin+3] = max(ddd)
        DK[kmax-kmin+3] = max(ddd)/min(ddd)
      }
    }
  } #End check condition for NCvalid
  # Start run loop kmin:kmax
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
    } # End check algorithm (cluss, xnew and centroid defined)
    if(!all(seq(k) %in% unique(cluss)))
      warning("Some clusters are empty.")
    #start index
    if(any(indexlist %in% c("all","CSL"))){
      for (j in 1:k){
        dd = as.matrix(dist(x[cluss==j,]))
        csl[k-kmin+1] = csl[k-kmin+1] + sum(apply(dd,2,max))/sum(cluss==j)
      } # end 1:k
      dd2 = as.matrix(dist(centroid))
      dd2 = matrix(dd2[dd2 > 0],k-1,k)
      csl[k-kmin+1] = csl[k-kmin+1]/sum(apply(dd2,2,min))
    }
    if(any(indexlist %in% c("all", "CH"))){
      d.cen = rowSums((x - xnew)^2)
      num = sapply(1:k, function(i) {
        ck = x[cluss == i, ]
        cen.k = centroid[i, ]
        n.clust_k = nrow(ck)
        n.clust_k * sum((cen.k - colMeans(x))^2)
      })
      dem = sapply(1:k, function(i) sum(d.cen[cluss == i]))
      ch[k-kmin+1] = ((nrow(x) - k) / (k-1)) * (sum(num) / sum(dem))
    }
    if(any(indexlist %in% c("all","DB","DBs","SF"))){
      # Si,q
      S = vector()
      sizecluss  = as.vector(table(cluss))
      for(i in 1:k){
        C = sizecluss[i]
        if(C>1){
          cenI = xnew[cluss ==i,]
          S[i] = (sum(sqrt(rowSums((x[cluss==i,]-cenI)^2))^q)/C)^(1/q)
        }else{
          S[i] = 0
        }
      }
      m = as.matrix(dist(centroid,method="minkowski",p=p))
      R = matrix(0,k,k)
      r = vector()
      rs = vector()
      wcdd = vector()
      for (i in 1:k){
        C = sizecluss[i]
        r[i] = max((S[i] + S[-i])/m[i,][m[i,]!=0])
        rs[i] = max(S[i] + S[-i])/min(m[i,][m[i,]!=0])
        # for SF index
        wcdd[i] = sum(dist(rbind(centroid[i,],x[cluss==i,]))[1:C])/C
      }
      db[k-kmin+1] = mean(r)
      dbs[k-kmin+1] = mean(rs)
      # SF index
      bcd = sum(dist(rbind(colMeans(x),centroid))[1:k]*sizecluss)/(dim(x)[1]*k)
      wcd = sum(wcdd)
      sf[k-kmin+1] = 1-1/exp(bcd+wcd)
    }
    if(any(indexlist %in% c("all","STR", "PBM"))){
      EK[k-kmin+2] = sum(sqrt(rowSums((x - xnew)^2)))
      ddd = dist(centroid)
      md[k-kmin+2] = max(ddd)
      DK[k-kmin+2] = max(ddd)/min(ddd)
    }
    if(any(indexlist %in% c("all","DI"))){
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
    if(any(indexlist %in% c("all","PB"))){
      d3 = as.vector(dist(xnew))
      d3[d3>0] = 1
      pb[k-kmin+1] = cor(d,d3,method=corr)
    }
    if(any(indexlist %in% c("all","NC", "NCI", "NCI1", "NCI2"))){
      d2 = as.vector(dist(xnew))
      crr[k-kmin+2]= cor(d,d2,method=corr)
    }
  }# End kloop
  if(any(indexlist %in% c("all","STR", "PBM"))){
    EKK = E0/EK
    str = (EKK[2:length(EKK)]-EKK[1:(length(EKK)-1)])*(DK[3:(length(DK))]-DK[2:(length(DK)-1)])
    pbm = EKK[2:(length(EKK))]*md[2:(length(EKK))]/(kmin:kmax)
  }

  # Create data frame for IDX result
  PB = data.frame(cbind("k"=kmin:kmax,"PB"=pb))
  CSL = data.frame(cbind("k"=kmin:kmax,"CSL"=csl))
  CH = data.frame(cbind("k"=kmin:kmax,"CH"=ch))
  DB = data.frame(cbind("k"=kmin:kmax,"DB"=db))
  DBs = data.frame(cbind("k"=kmin:kmax,"DBs"=dbs))
  SF = data.frame(cbind("k"=kmin:kmax,"SF"=sf))
  DI = data.frame(cbind("k"=kmin:kmax,"DI"=di))
  STR = data.frame(cbind("k"= kmin:kmax, "STR" = str))
  PBM = data.frame(cbind("k"= kmin:kmax, "PBM" = pbm))
  # End defined data frame
  if(any(indexlist %in% c("all","NC", "NCI", "NCI1", "NCI2"))){
    K = length(crr)
    NWI = ((crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)]))/pmax(0,(crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)]))
    NWI2 = (crr[2:(K-1)]-crr[1:(K-2)])/(1-crr[1:(K-2)])-(crr[3:K]-crr[2:(K-1)])/(1-crr[2:(K-1)])
    NWI3 = NWI
    if (max(NWI)<Inf){
      if (min(NWI) == -Inf){
        NWI3[NWI==-Inf] = min(NWI[is.finite(NWI)])
      }
    }
    if (max(NWI)==Inf){
      NWI3[NWI==Inf] = max(NWI[is.finite(NWI)])+NWI2[NWI==Inf]
      NWI3[NWI<Inf] = NWI[NWI<Inf] + NWI2[NWI<Inf]
      if (min(NWI) == -Inf){
        NWI3[NWI==-Inf] = min(NWI[is.finite(NWI)])+NWI2[NWI==-Inf]
      }
    }
    NWI = data.frame(cbind("k"= kmin:kmax, "NCI1" = NWI))
    NWI2 = data.frame(cbind("k"=kmin:kmax,"NCI2"=NWI2))
    NWI3 = data.frame(cbind("k"=kmin:kmax,"NCI"=NWI3))
    crr = data.frame(cbind("k" = (kmin-1):(kmax+1), "NC" = crr))
  }else{
    crr = 0
    NWI = 0
    NWI2 = 0
    NWI3 = 0
  }
  my_list <- list("NC"=crr, "NCI" = NWI3, "NCI1" = NWI, "NCI2" = NWI2, "CSL"=CSL, "STR"= STR, "CH"=CH,"DB"=DB, "DBs"=DBs, "PBM" = PBM, "DI"=DI, "PB"=PB, "SF"=SF)
  if (sum(indexlist == "all")==1){
    return(my_list)
  } else {
    return(my_list[indexlist])
  }
}
