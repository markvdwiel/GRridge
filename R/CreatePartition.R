CreatePartition <- function(vec,varnamesdata=NULL,
                            subset=NULL,grsize=NULL,
                            decreasing=TRUE,uniform=FALSE,
                            ngroup=10,mingr=25){
    if(is.factor(vec)){
      firstcl <- lapply(as.character(levels(vec)),function(xg) which(vec==xg))
      names(firstcl) <- levels(vec)
      } else {
        if(is.numeric(vec)){
          if(uniform){
            if(is.null(grsize)){
              grsize <- floor(length(vec)/ngroup)
              print(paste("Group size set to:",grsize))
              } else {
                print(paste("Group size",grsize))
                }
            if(decreasing) {
              print("Sorting vec in decreasing order, assuming small values are LESS relevant")
                orderp2 <- order(vec,decreasing=TRUE) 
                lroep <- length(vec)
                } else {
                print("Sorting vec in increasing order, assuming small values are MORE relevant")
                orderp2 <- order(vec,decreasing=FALSE) 
                lroep <- length(vec)   
                }
            
            ngr <- floor(lroep/grsize)
            firstcl <- lapply(1:ngr,function(xg) {
            if(xg < ngr) els <- orderp2[(1+(xg-1)*grsize):(xg*grsize)] else 
            els <- orderp2[(1+(xg-1)*grsize):lroep]
            return(els)
            }
            )
            names(firstcl) <- sapply(1:length(firstcl),function(i) paste("group",i,sep=""))
    } else {
    if(decreasing) {
                print("Sorting vec in decreasing order, assuming small values are LESS relevant")
                orderp2 <- order(vec,decreasing=TRUE) 
                lroep <- length(vec)
                } else {
                print("Sorting vec in increasing order, assuming small values are MORE relevant")
                orderp2 <- order(vec,decreasing=FALSE) 
                lroep <- length(vec)   
                }
    p <- length(vec) 
    if(ngroup*mingr >= p) {
    print("ERROR: Number of groups (ngroup) times minimal group size (mingr) is larger 
          than number of variables. Please use uniform = TRUE or decrease either ngroup or mingr.")
    return(NULL)  
    }
    povermin <- p/mingr
    parint <-povermin^{1/ngroup}
      
    lefts <- povermin+1
    gfun2 <- function(x){1-x^(ngroup+1) - lefts*(1-x)}
    root <- uniroot(f=gfun2, lower=1.000001,upper=parint)$root

    grs <- sapply(1:ngroup,function(i) if(i==1) floor(mingr*root^i) else round(mingr*root^i)) 
    sm <- sum(grs)
    grs[ngroup] <- grs[ngroup] -(sm-p)
    print("Summary of group sizes:")
    print(summary(grs))
    cumul <- cumsum(c(0,grs))
    firstcl <- lapply(1:ngroup,function(xg) {
            els <- orderp2[(cumul[xg]+1):cumul[xg+1]]
            return(els)
            }
            )
    names(firstcl) <- sapply(1:length(firstcl),function(i) paste("group",i,sep=""))
    }
} else {  #assume character
    if(!is.character(vec)){
    print("Argument vec is not correctly specified")
    return(NULL)
    } else {
        if(is.null(varnamesdata)){
        print("Please specify a character vector for varnamesdata")
        return(NULL)
        } else {
        whin <- match(vec,varnamesdata)
        whin <- unique(whin[!is.na(whin)])
        firstcl <- list(VarIn=whin,VarOut=(1:length(varnamesdata))[-whin])
        }
    }
}
}
if(!is.character(vec) & !is.null(subset)){ #remapping 
if(is.null(varnamesdata)){
print("ERROR: varnamesdata required for subsetting")
return(NULL)
}
if(length(vec) != length(subset)){
print("ERROR: Length of vec does not match length of subset.")
return(NULL)
} else {
matchss <- match(subset,varnamesdata)
firstcl <- lapply(firstcl,function(vector) matchss[vector])
}
}
print("Summary of group sizes:")
print(unlist(lapply(firstcl,length)))
return(firstcl)
}
