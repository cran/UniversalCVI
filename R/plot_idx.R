plot_idx <- function(idxresult, selected.idx = NULL){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if(!any(class(idxresult) %in% c("list","data.frame")))
    stop("Argument 'idxresult' must be list or data frame.")
  if(sum(!(names(idxresult) %in% c("c","WPC","WP","WPCI1","WPCI2","XB","KWON","KWON2","TANG","HF","WL","PBM","KPBM","CCVP","CCVS","GC1","GC2","GC3","GC4",
                                   "NC","NCI","NCI1","NCI2","CSL","CH","DB","DBs","DI","PB","SF","STR"))) >0)
    stop("Bad input data, 'idxresult' is not result from our package function.")
  if(is.data.frame(idxresult)){
    par(mar = c(4, 4, 0.5, 0.5))
    par(mfrow=c(1,1))
    IDX.name = names(idxresult)[2]
    if(any(IDX.name %in% c("WPC","NC"))){
      plot(data.frame(idxresult),ylab = IDX.name, xlab = "c",type='b')
    }
    # Plot the indexes that the largest value indicates a valid optimal partition.
    else if(any(IDX.name %in% c("WP","WPCI1","WPCI2","PBM","KPBM","CCVP","CCVS","CH","DI","PB","NCI","NCI1","NCI2","STR"))){
      plot(data.frame(idxresult),ylab = IDX.name, xlab = "c",type='b')
      points(data.frame(idxresult)[which.max(data.frame(idxresult)[,2]),1],max(data.frame(idxresult)[,2]),col='red',pch=20)
    }else { # Plot the indexes that the smallest value indicates a valid optimal partition.
      plot(data.frame(idxresult),ylab = IDX.name, xlab = "c",type='b')
      points(data.frame(idxresult)[which.min(data.frame(idxresult)[,2]),1],min(data.frame(idxresult)[,2]),col='red',pch=20)
    }
  }else{
    if(!any(class(selected.idx) %in% c("integer","numeric","NULL")))
      stop("Argument 'selected.idx' must be in one of the classes: 'integer', 'numeric' or 'NULL'.")
    if(length(selected.idx) > length(idxresult) | sum(selected.idx > length(idxresult)) > 0)
      stop("The largest number or length in argument 'selected.idx' must be less than or equal to the number of indices in idxresult")
    if(is.null(selected.idx)){
      idxresult = idxresult
    } else {
      idxresult = idxresult[names(idxresult)[selected.idx]]
    }
    n.idx = length(idxresult)
    name.idx = names(idxresult)
    if(n.idx %in% c(13,18)){
      name.idx = name.idx[!name.idx %in% c('WPC', 'WPCI1', 'WPCI2','NC','NCI1','NCI2')][1:8]
      n.idx = length(name.idx)
    }else{
      n.idx = min(8,n.idx)
      name.idx = name.idx[1:n.idx]
    }
    par(mar = c(4, 4, 0.5, 0.5))
    if(n.idx <= 3){
      par(mfrow=c(1,n.idx))
    }else {
      par(mfrow=c(2,ceiling(n.idx/2)))
    }
    for(np in 1:n.idx){
      IDX.name = name.idx[np] #name of index
      IDX = idxresult[[IDX.name]]
      if(any(IDX.name %in% c("WPC","NC"))){
        plot(data.frame(IDX),ylab = IDX.name, xlab = "c",type='b')
      }
      # Plot the indexes that the largest value indicates a valid optimal partition.
      else if(any(IDX.name %in% c("WP","WPCI1","WPCI2","PBM","KPBM","CCVP","CCVS","CH","DI","PB","NCI","NCI1","NCI2","STR"))){
        plot(data.frame(IDX),ylab = IDX.name, xlab = "c",type='b')
        points(data.frame(IDX)[which.max(data.frame(IDX)[,2]),1],max(data.frame(IDX)[,2]),col='red',pch=20)
      }else { # Plot the indexes that the smallest value indicates a valid optimal partition.
        plot(data.frame(IDX),ylab = IDX.name, xlab = "c",type='b')
        points(data.frame(IDX)[which.min(data.frame(IDX)[,2]),1],min(data.frame(IDX)[,2]),col='red',pch=20)
      }
    }
  }
}

