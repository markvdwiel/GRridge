CreatePartition <-
function(vec,varNames=NULL,grsize=NULL,decreasing=TRUE,uniform=FALSE,ngroup=100,mingr=10){
#vec <- pvalFarkas;mingr=10;ngroup=100;decreasing=FALSE;uniform=FALSE
 #vec <- pvalonly;mingr=10;ngroup=3;decreasing=FALSE;uniform=FALSE
if(is.factor(vec)){
firstcl <- lapply(as.character(levels(vec)),function(xg) which(vec==xg))
names(firstcl) <- levels(vec)
} else {
if(is.numeric(vec)){
    if(uniform){
        if(is.null(grsize)){
        print("Please specify grsize for numeric input of vec")
        return(NULL)
        } else {
            if(decreasing) {
                print("Sorting vec in decreasing order, assuming small values are LESS relevant")
                orderp2 <- order(vec,decreasing=T) 
                lroep <- length(vec)
                } else {
                print("Sorting vec in increasing order, assuming small values are MORE relevant")
                orderp2 <- order(vec,decreasing=F) 
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
        }
    } else {
    if(decreasing) {
                print("Sorting vec in decreasing order, assuming small values are LESS relevant")
                orderp2 <- order(vec,decreasing=T) 
                lroep <- length(vec)
                } else {
                print("Sorting vec in increasing order, assuming small values are MORE relevant")
                orderp2 <- order(vec,decreasing=F) 
                lroep <- length(vec)   
                }
    p <- length(vec) 
    if(ngroup*mingr >= p) {
    print("ERROR: Number of groups (ngroup) times minimal group size (mingr) is larger than number of variables. Please use uniform = TRUE or decrease either ngroup or mingr.")
    return(NULL)  
    }
    povermin <- p/mingr
    parint <-povermin^{1/ngroup}
      
    lefts <- povermin+1
    gfun2 <- function(x){1-x^(ngroup+1) - lefts*(1-x)}
    root <- uniroot(f=gfun2, lower=1.000001,upper=parint)$root
    
    
    
    grs <- sapply(1:ngroup,function(i) if(i==1) floor(mingr*root^i) else round(mingr*root^i)) #garantees the smallest size equals mingr
    sm <- sum(grs)
    grs[ngroup] <- grs[ngroup] -(sm-p)
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
        if(is.null(varNames)){
        print("Please specify a character vector for varNames")
        return(NULL)
        } else {
  #      print("Please verify that length(varNames) equals the number of variables (features; usually rows) in your data")
  #vec=genesroepman;varNames=genes
        whin <- match(vec,varNames)
        whin <- unique(whin[!is.na(whin)])
        firstcl <- list(VarIn=whin,VarOut=(1:length(varNames))[-whin])
        }
    }
}
}
return(firstcl)
}
