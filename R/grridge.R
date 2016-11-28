grridge <- function(highdimdata, response, partitions, unpenal = ~1, 
                    offset=NULL, method="exactstable",
                    niter=10, monotone=NULL, optl=NULL, innfold=NULL, 
                    fixedfoldsinn=TRUE, selectionForward=FALSE,
                    maxsel=100,selectionEN=FALSE,stepsel=1,cvlmarg=1,
                    savepredobj="all", dataunpen=NULL, ord = 1:length(partitions),
                    comparelasso=FALSE,optllasso=NULL,cvllasso=TRUE,
                    compareEN=FALSE,compareunpenal=FALSE,trace=FALSE,modus=1){#

    if(method=="adaptridge" | method== "exact") niter <- 1
  
  if(class(partitions[[1]]) =="integer"){
    partitions=list(group=partitions)
  }
  nclass <- length(partitions)
  
  if(is.null(monotone)) monotone <- rep(FALSE,nclass)
  if(length(monotone) != length(partitions)) {
    print(paste("ERROR: length 'monotone' unequal to length 'partitions' "))
    return(NULL) 
  }
  partitions <- partitions[ord]
  monotone <- monotone[ord]
  nr <- nrow(highdimdata)
  
  for(ncl in 1:nclass){
    indexset <- unlist(partitions[[ncl]])
    if(length(indexset) < nr){
      print(paste("Warning: partition",ncl,"does not contain all row indices of the data"))
    }
    if(max(indexset) > nr | min(indexset)<1){
    print(paste("ERROR: partition",ncl,"contains an invalid index, e.g. larger than number of data rows"))
    return(NULL)
    }
  }
  
  overlap <- c()
  Wmat <- c()
  nfeattot <- c()
  for(ncl in 1:nclass){
    indexset <- unlist(partitions[[ncl]])
    nfeatcl <- length(unique(indexset))
    nfeattot <- c(nfeattot,nfeatcl)
    if(length(indexset) > nfeatcl){
      print(paste("Grouping",ncl,"contains overlapping groups"))
      overlap <- c(overlap,TRUE)
      whgroup <- partitions[[ncl]]
      nover <- rep(0,nr)
      for(k in 1:length(whgroup)){
        wh <- whgroup[[k]]  
        nover[wh] <- nover[wh] + 1
      }
      Wmat <- cbind(Wmat,sqrt(1/nover))
    } else {
      print(paste("Grouping",ncl,"contains mutually exclusive groups"))
      overlap <- c(overlap,FALSE)  
      Wmat <- cbind(Wmat,rep(1,nr))
    }
  }
  
  arguments <- list(partitions=partitions,unpenal=unpenal, offset=offset, method=method, 
                    niter=niter, monotone=monotone, optl=optl, innfold=innfold,
                    fixedfoldsinn=fixedfoldsinn,selectionForward=selectionForward,
                    selectionEN=selectionEN,maxsel=maxsel,stepsel=stepsel, 
                    cvlmarg=cvlmarg, dataunpen=dataunpen,savepredobj=savepredobj, ord=ord, 
                    comparelasso=comparelasso, optllasso=optllasso, compareEN=compareEN, 
                    compareunpenal=compareunpenal, modus=modus)
  
  
  if(nr > 10000 & is.null(innfold)) print("NOTE: consider setting innfold=10 to save computing time")
  
  nmp0 <- names(partitions)
  if(is.null(nmp0)) nmp0 <- sapply(1:length(partitions),function(i) paste("Grouping",i))  
  nmp0 <- sapply(1:length(partitions),function(i) {
    if(nmp0[i]=="") return(paste("Grouping",i)) else return(nmp0[i])})
  
  nmp <- c("NoGroups","GroupRegul")
  nmpweight <- nmp #to be used later
  if(selectionForward) nmp <- c(nmp,"ForwSel") 
  if(comparelasso) nmp <- c(nmp,"lasso") 
  if(compareEN) nmp <- c(nmp,"EN")
  if(compareunpenal) nmp <- c(nmp,"modelunpen")
  
  if(class(response) =="factor") {
    nlevel <- length(levels(response))
    if(nlevel != 2){
      print("Response is not binary, so not suitable for two-class classification.")
      return(NULL)
    } else {
      model = "logistic"
      print("Binary response, executing logistic ridge regression")
      lev <- levels(response)
      print(paste("Predicting probability on factor level",lev[2]))
    }} else {
      if(class(response) == "numeric" | class(response)=="integer"){
        valresp <- sort(unique(response)) 
        if(length(valresp)==2 & valresp[1]==0 & valresp[2]==1) {
          model = "logistic"
          print("Binary response, executing logistic ridge regression")
        } else {
          model = "linear"
          print("Numeric continuous response, executing linear ridge regression")
        }
      } else {
        if(class(response) == "Surv"){
          model="survival"
          print("Survival response, executing cox ridge regression")
        } else {
          print("Non-valid response. Should be binary, numeric or survival.")
          return(NULL)
        }
      }
    }
  
  if((unpenal != ~0) & (unpenal != ~1)) {
    if(is.null(dataunpen)) {print("If unpenal contains variables, data of 
                                  the unpenalized variables should be specified in the data slot!")
      return(NULL)
    } 
  }
  nsam <- ncol(highdimdata)
  if(!is.null(offset)){
    noffs <- length(offset)
    offsets <-"c("
    if(noffs==1) {for(i in 1:(nsam-1)) offsets <- paste(offsets,offset,",",sep="")} else {
      for(i in 1:(nsam-1)) offsets <- paste(offsets,offset[i],",",sep="")}
    if(noffs==1) offsets <- paste(offsets,offset,")",sep="")  else offsets <- paste(offsets,offset[nsam],")",sep="")
    if((unpenal != ~0) & (unpenal != ~1)){
      unpenal <- formula(paste(deparse(unpenal),"+ offset(",offsets,")",sep=""))
    } else {
      unpenal <- formula(paste("~","offset(",offsets,")",sep=""))
    }
  }
  
  
  if(is.null(dataunpen)) datapred <- data.frame(fake=rep(NA,ncol(highdimdata))) else datapred <- dataunpen
  nopen <- unpenal
  
  if(is.null(innfold)) foldinit <- nsam else foldinit <- innfold
  pmt0<- proc.time()
  optl0 <- optl
  if(is.null(optl)){
    print("Finding lambda for initial ridge regression")
    if(fixedfoldsinn) set.seed(346477)
    opt <- optL2(response, penalized = t(highdimdata),fold=foldinit,unpenalized=nopen,data=datapred,trace=trace)
    time1 <- proc.time()-pmt0
    print(opt$cv)
    print(paste("Computation time for cross-validating main penalty parameter:",time1[3]))
    optl <- opt$lambda
    print(paste("lambda2",optl))
    arguments$optl <- optl
  }
  pmt <- proc.time()
  nsam <- ncol(highdimdata)
  nfeat <- nrow(highdimdata)
  XM0 <- t(highdimdata)
  response0 <- response
  
  if(!selectionForward) cvlnstot <- rep(0,(nclass+1)) else cvlnstot <- rep(0,(nclass+2))
  allpreds <- c()
  whsam <- 1:nsam
  responsemin <- response0
  pen0 <- penalized(responsemin, penalized = XM0, lambda2 = optl, unpenalized=nopen,data=cbind(XM0,datapred))
  nmunpen <- names(pen0@unpenalized)
  if(is.element("(Intercept)",nmunpen)) addintercept <- TRUE else addintercept <- FALSE 
  if(is.null(innfold)) {nf <- nrow(XM0)} else {
    if(!is.null(optl0)) {
      nf <- innfold
      if(fixedfoldsinn) set.seed(346477)
    } else {nf <- opt$fold}
  }
  opt2 <- cvl(responsemin, penalized = XM0,fold=nf, lambda2 = optl,unpenalized=nopen,data=datapred,trace=trace)
  nf <- opt2$fold 
  cvln0 <- opt2$cvl
  
  cvlnprev <- cvln0
  penprev <- pen0
  pen <- pen0
  print(cvln0)
  XMw0 <- XM0
  XMw0prev <- XM0
  converged <- FALSE
  conv <- rep(FALSE,nclass)
  controlbound1 <-1000
  controlbound2 <- controlbound3 <- 10
  
  almvecall <- rep(1,nfeat)
  lambdas <- lapply(partitions, function(cla) {
    ngroup <- length(cla)
    return(rep(1,ngroup))
  })
  lmvec <- lmvecprev <- array(1,nfeat)
  i <- 1
  
  while(!converged & i <= niter){
    cl <- 1
    if(method=="adaptridge") cl <- nclass
    while(cl <= nclass){
      convcl <- conv[cl]
      
      if(!convcl){
        whgr <- partitions[[cl]]
        lenggr <- unlist(lapply(whgr,length))
        ngroup1 <- length(whgr)
        names(lambdas[[cl]]) <- names(whgr)
        coeff <- penprev@penalized
        if(model == "survival"){
          preds <- predict(penprev,XMw0,data=datapred)
        } else {
          preds <- predict(penprev,XMw0,data=datapred)[1:nsam]   
        }  
        coeffsq <- coeff^2
        
        if(model=="logistic") {
          Wi <- sqrt(preds*(1-preds)) 
          constlam <- 2
        }
        if(model == "linear"){
          Wi <- rep(1,length(preds))
          constlam <- 1
        }
        if(model == "survival"){
          resptime <- response[,1]
          predsnew <- -log(sapply(1:nsam,function(k) survival(preds,time=resptime[k])[k]))
          Wi <- sqrt(predsnew)
          constlam <- 2
        }
        
        if(!is.null(dataunpen)) {
          mm <- model.matrix(nopen,dataunpen) 
          XMW <- t(t(cbind(XMw0,10^5*mm)) %*% diag(Wi)) 
        } else {
          if(addintercept) XMW <- t(t(cbind(XMw0,rep(10^5,nsam))) %*% diag(Wi)) else XMW <- t(t(XMw0) %*% diag(Wi))
        }
        SVD <- svd(XMW)
        leftmat <- SVD$v %*% diag(1/((SVD$d)^2+constlam*optl)) %*% diag(SVD$d) %*% t(SVD$u)
        
        if(model=="linear"){
          Hatm <- XMW %*% leftmat 
          df <- nsam - sum(diag(2*Hatm - Hatm %*% t(Hatm)))
          VarRes <- sum((response - preds)^2)/df
          print(paste("Sigma^2 estimate:",VarRes))
          vars3 <- VarRes*rowSums(leftmat^2)
        } else {
          vars3 <- rowSums(leftmat^2)
        }
        
        which0 <- which(vars3==0)
        vars3[which0] <- 10^{-30}
        
        if(model=="linear"){
          mycoeff2svd <- (leftmat %*% response)^2
        } 
        if(model=="logistic"){
          if(is.factor(response)) respnum <- as.numeric(response)-1 else respnum <- response
          z <- matrix(log(preds/(1-preds))+(respnum - preds)/(preds*(1-preds)),ncol=1) 
          if(modus == 1) mycoeff2svd <- coeffsq
          if(modus == 2) mycoeff2svd <- (leftmat %*% z)^2 
        }
        if(model=="survival"){
          mycoeff2svd <- coeffsq 
        } 
        cii2 <- (rowSums(leftmat * t(XMW)))^2
        
        
        leftmat <- leftmat/sqrt(vars3)
        lowerleft <- 10^(-30)
        lefts2 <- function(group){
          ind <- whgr[[group]]
          ngr <- lenggr[group]
          coefftau2 <- sum(sapply(mycoeff2svd[ind]/vars3[ind],function(x) max(x,1)))-length(ind)
          return(max(lowerleft,coefftau2/ngr))
        } 
        
        leftside <- sapply(1:length(whgr),lefts2)
        ellarg0 <- length(leftside[leftside>lowerleft])/length(leftside)
        
         if(ellarg0 <=.5){ 
           print(paste("Partition",nmp0[cl],"NOT ITERATED")) 
           conv[cl] <- TRUE
           cvln1 <- cvlnprev
           XMw0 <- XMw0prev
           pen <- penprev
         } else {
        lefts2ran <- function(group,randomind){
          ind <- whgr[[group]]
          ngr <- lenggr[group]
          coefftau2 <- sum(sapply(mycoeff2svd[1:nfeat][randomind[ind]]/vars3[1:nfeat][randomind[ind]],
                                  function(x) max(x,1)))-length(randomind[ind])
          return(max(lowerleft,coefftau2/ngr))
        } 
        randomiz <- function(fakex){
          randomind2 <- sample(1:nfeat)
          leftsideran <- sapply(1:length(whgr),lefts2ran,randomind=randomind2)
          return(leftsideran)}
        nlefts <- 100
        leftsran <-  sapply(1:nlefts,randomiz)
        means <- apply(leftsran,1,mean)
        leftsrancen <- t(t(leftsran)-means)
        relerror <- sum(abs(leftsrancen))/(nlefts*sum(abs(means)))
        if(cl==1 & i==1) print(cvln0)
        print(paste("Relative error:",relerror))  
        if(relerror>= 0.1) print("WARNING: large relative error (>=0.1). Consider using larger groups of variable.")
        nadd <- ncol(XMW) - ncol(XMw0)
        
        rightmat =  t(t(XMW) * c(Wmat[,cl],rep(1,nadd)))
        
        rightmats <- lapply(1:length(whgr),function(j){
          rightj2 <- rightmat[,whgr[[j]]]
          rcp <- rightj2 %*% t(rightj2)
          return(rcp)
        })
        
        coefmatfast <- t(apply(matrix(1:length(whgr),nrow=length(whgr)),1,  
                               function(i){
                                 lefti2  <- leftmat[whgr[[i]],]
                                 lcp <- t(lefti2) %*% lefti2
                                 ckls <- sapply(1:length(whgr),function(j){
                                   rcp <- rightmats[[j]]
                                   return(sum(lcp*rcp))
                                 })
                                 return(ckls)
                               }))
        coefmatfast <- coefmatfast/lenggr   
        
        CNfun <- function(lam,cfmmat=coefmatfast){
          ng <- nrow(cfmmat)
          dmax <- max(diag(cfmmat))
          cfmlam <- (1-lam)*cfmmat + lam*diag(dmax,nrow=ng)
          eigenvals <- eigen(cfmlam,only.values=TRUE)$values
          CN <- eigenvals[1]/eigenvals[ng]
          return(Re(CN)) 
        }
        
        lams <- seq(0,1,by=0.005)
        CNsRan <- sapply(lams, CNfun,cfmmat=coefmatfast)
        CNsRanre <- CNsRan*relerror
        if(relerror<=0.1){
          lam <- lams[which(CNsRanre<=0.1)[1]]
        } else lam <- 1
        print(paste("Shrink Factor coefficient matrix",lam))
        cfmmat <- coefmatfast;
        ng <- nrow(cfmmat)
        dmax <- max(diag(cfmmat))
        cfmlam <- (1-lam)*cfmmat + lam*diag(dmax,nrow=ng)

        if(method=="exactstable"){
          soltau = solve(sum(cfmlam),sum(leftside))
          sol = solve(cfmlam,leftside)
          low <- soltau/controlbound1;up = soltau*controlbound1
          parinint <- sapply(sol,function(x) min(max(low,x),up))
          minopt <- optim(par = parinint,fn = function(pars=c(parinint)) sum(leftside - cfmlam %*% pars)^2,
                          method="L-BFGS-B",
                          lower=rep(low,ngroup1),upper=rep(up,ngroup1))
          tausqest0 <- minopt$par
        }
        
        if(method=="exact"){
          soltau = solve(sum(coefmatfast),sum(leftside))
          sol = solve(coefmatfast,leftside)
          low <- soltau/controlbound2;up = soltau*controlbound2
          parinint <- sapply(sol,function(x) min(max(low,x),up))
          minopt <- optim(par = parinint,fn = function(pars=c(parinint)) sum(leftside- cfmlam %*% pars)^2,
                          method="L-BFGS-B",
                          lower=rep(low,ngroup1),upper=rep(up,ngroup1))
          tausqest0 <- minopt$par
        }
        
        if(method=="stable"){
          soltau = solve(sum(coefmatfast),sum(leftside))
          solhyb <- sapply(1:ngroup1,function(i){
            leftsidei <- leftside[i] 
            rightsidei <- c(coefmatfast[i,i], sum(coefmatfast[i,-i]*soltau))
            soli <- (leftsidei-rightsidei[2])/rightsidei[1]
            return(max(min(soltau*controlbound3,soli),soltau/controlbound3))
          })
        }
        
        if(method=="simple"){
          solsim = leftside
          tausqest <- solsim
          print("simple")
        }
        
        if(method=="exact") {
          tausqest <- tausqest0
          print("exact")
        }
        
        if(method=="exactstable") {
          tausqest <- tausqest0
          print("exactstable")
        }
        if(method=="stable") {
          print("stable")
          tausqest <- solhyb
        } 
        
        if(method=="adaptridge") {
          print("adaptive ridge")
        } 
        
        if(method=="stable" | method=="exact" | method=="exactstable" | method=="simple"){
          
          lambdanoncal <- 1/tausqest
          
          if(monotone[cl]){
            weigh = unlist(lapply(whgr,length))
            lambdamultnoncal <- pava(lambdanoncal,w=weigh)
          } else lambdamultnoncal <- lambdanoncal
          
          tausqest<-1/lambdamultnoncal
          nfeatcl <- nfeattot[cl]
          overl <- overlap[cl]
          if(!overl){
            con3 <- sum(sapply(1:length(whgr),function(gr){return(length(whgr[[gr]])*tausqest[gr])}))
            tausqestcal<-  nfeatcl/con3*tausqest
            lambdamult <- 1/tausqestcal
            print(lambdamult)
            
            for(k in 1:length(whgr)){
              wh <- whgr[[k]]
              XMw0[,wh] <- XMw0[,wh]/sqrt(lambdamult[k])
            }
          } else {   
            tauk <- rep(0,nfeat)
            Wsq <- (Wmat[,cl])^2
            for(k in 1:length(whgr)){
              wh <- whgr[[k]]  
              tauk[wh] <- tauk[wh] + tausqest[k]
            }
            tauk <- tauk*Wsq
            whna <- which(is.na(tauk))
            if(length(whna)>0) con3 <- sum(tauk[-whna]) else  con3 <- sum(tauk) 
            
            tausqestcal0<-  nfeatcl/con3*tausqest
            lambdamult <- 1/tausqestcal0
            print(lambdamult)
            
            tausqestcal<-  (nfeatcl/con3)*tauk
            lambdamultperk <- 1/tausqestcal
            lambdamultperk[whna]<- 1
            XMw0 <- t(t(XMw0)/sqrt(lambdamultperk))
          }
        } else {
          tausqest <- coeffsq
          con3 <- sum(tausqest)
          tausqestcal<-  nfeat/con3*tausqest
          lambdamult <- 1/tausqestcal
          XMw0 <- t(t(XMw0)/sqrt(lambdamult))
        }
        
        opt2w <- cvl(responsemin, penalized = XMw0,fold=nf,lambda2=optl,
                     unpenalized=nopen,data=datapred, trace=trace)
        cvln1 <- opt2w$cvl
        print(cvln1)
        
        if((cvln1 - cvlnprev)/abs(cvlnprev) > 1/100  |  ((cvln1 - cvlnprev)/abs(cvlnprev) >= 0 & i==1)){ 
          pen <- penalized(responsemin, penalized = XMw0, 
                           lambda2 = optl,unpenalized=nopen,data=datapred)
          if(niter>1){
            if(!overl){
              for(group in 1:ngroup1){
                ind <- whgr[[group]]
                lmvec[ind] <- lmvec[ind]*lambdamult[group]
              }} else {
                lmvec <- lmvec * lambdamultperk
              }
          }
          lambdas[[cl]] <- lambdas[[cl]]*lambdamult
          cvlnprev <- cvln1
          penprev <- pen
          XMw0prev <- XMw0
          print(paste("Partition",nmp0[cl],"improved results"))
        } else {
          if(niter>1) print(paste("Partition",nmp0[cl],"CONVERGED after",i,"iterations")) 
          else print(paste("Partition",nmp0[cl],"did not improve results"))
          
          conv[cl] <- TRUE
          cvln1 <- cvlnprev
          XMw0 <- XMw0prev
          pen <- penprev
        }
         } 
      } 
      cl <- cl+1
    }
    if(sum(conv)==nclass){
      converged <- TRUE
      if(niter>1) print(paste("All partitions CONVERGED after",i,"iterations"))
    }
    i <- i+1
  }

if(niter==0) {pen <- pen0;XMw0<-XM0;cvln1<-cvln0;soltau <- NULL}
if(model=="survival"){
pred0 <- predict(pen0,XM0,unpenalized=nopen,data=datapred)
predw <- predict(pen,XMw0,unpenalized=nopen,data=datapred)
} else {
pred0 <- predict(pen0,XM0,unpenalized=nopen,data=datapred)[1:nsam]  
predw <- predict(pen,XMw0,unpenalized=nopen,data=datapred)[1:nsam]
}
predshere <- cbind(pred0,predw)
cvlnssam <- c(cvln0,cvln1)
lmvecall <- lmvec
almvecall <- cbind(almvecall,lmvecall)
predobj <- c(pen0,pen)

allpreds <- predshere
whichsel <- NULL
betassel <- NULL

npr <- length(predobj)
pred2 <-predobj[[npr]]
lambs <- lmvecall
oldbeta <- pred2@penalized
newbeta <- oldbeta/sqrt(lambs)
time2 <- proc.time()-pmt
print(paste("Computation time for adaptive weigthing:",time2[3]))

if(selectionForward){
print("Start posthoc variable selection")
pmt <- proc.time()
ord <- order(abs(newbeta),decreasing=TRUE)
sequ <- seq(0,maxsel,by=stepsel)
cvlsels <- c()
for(nsel in sequ){
whsel <- ord[1:nsel]
datwsel <- XMw0[,whsel]

pensel <- penalized(responsemin,datwsel,lambda2=optl, unpenalized=nopen,data=datapred,trace=FALSE)
optsel <- cvl(responsemin,datwsel,fold=nf,lambda2=optl,unpenalized=nopen,data=datapred, trace=trace)
cvlsel <- optsel$cvl
cvlsels <- c(cvlsels,cvlsel)
}
whbest <- which.max(cvlsels)
cvmax <- cvlsels[whbest]
nsel2 <- sequ[(which((cvlsels - cvmax) >= cvlmarg/100*(cvmax)))[1]]
whsel2 <- ord[1:nsel2]
whichsel <- whsel2
betassel <- newbeta[whsel2]
datwsel2 <- XMw0[,whsel2]
pensel2 <- penalized(responsemin,datwsel2,lambda2=optl,unpenalized=nopen,data=datapred)
optsel2 <- cvl(responsemin,datwsel2,fold=nf,lambda2=optl,unpenalized=nopen,data=datapred, trace=trace)
cvlsel2 <- optsel2$cvl
if(model=="survival"){
  predsel <- predict(pensel2,penalized=XMw0[,whsel2,drop=FALSE],unpenalized=nopen,data=datapred)
} else {
  predsel <- predict(pensel2,penalized=XMw0[,whsel2,drop=FALSE],unpenalized=nopen,data=datapred)[1:nsam]
}

allpreds <- cbind(allpreds,predsel)
predobj <- c(predobj,pensel2)
cvlnssam <- c(cvlnssam,cvlsel2)
print(paste("Number of selected markers by forward selection:",nsel2))
time3 <- proc.time()-pmt
print(paste("Computation time for feature selection by forward selection:",time3[3]))
}
if(selectionForward) cvlssel <- data.frame(nsel=sequ,cvl=cvlsels) else cvlssel <- NULL

cvlnstot <- cvlnssam
reslasso <- NULL
if(comparelasso){
print("Starting lasso")
    if(is.null(optllasso)){
    print("Finding lambda for lasso regression")
    opt <- optL1(response, penalized = t(highdimdata),fold=nf,unpenalized=nopen,data=datapred,trace=trace)
    print(opt$cv)
    optllasso <- opt$lambda
    print(paste("lambda1",optllasso))
    arguments$optllasso <- optllasso
    cvllasso <- opt$cv
    } else {
    cvliklasso <- if(cvllasso) try(cvl(response,penalized = t(highdimdata), lambda1 = optllasso,
                                       fold=nf,unpenalized=nopen,data=datapred, trace=trace))
    if(class(cvliklasso) == "try-error" | !cvllasso) cvllasso <- NA else cvllasso <- cvliklasso$cvl 
        }
cvlnstot <- c(cvlnstot,cvllasso)
penlasso <- penalized(response, penalized = t(highdimdata), lambda1 = optllasso, 
                      unpenalized=nopen,data=cbind(XM0,datapred))
whichlasso <- which(penlasso@penalized != 0)
betaslasso <- penlasso@penalized[whichlasso]
predobj <- c(predobj,penlasso)
reslasso <- list(cvllasso=cvllasso,whichlasso=whichlasso,betaslasso=betaslasso)
}

resEN <- NULL 
if(selectionEN){
  
  if(selectionForward) maxsel <- length(whsel2)
  
  fsel <- function(lam1,maxselec=maxsel,lam2){
    if(lam1==0) return(nfeat-maxselec) else {
      penselEN <- penalized(responsemin,XMw0,lambda1=lam1,lambda2=lam2, 
                            unpenalized=nopen,data=datapred,trace=FALSE,maxiter=100)
      coef <- penselEN@penalized
      return(length(coef[coef!=0])-maxselec)
    }
  }
  lam1 <- uniroot(fsel,interval=c(0,optl*10),maxiter=50,lam2=optl)$root
   penselEN0 <- penalized(responsemin,XMw0,lambda1=lam1,lambda2=optl, unpenalized=nopen,data=datapred,
                         trace=FALSE,maxiter=100)
   coefEN0 <- penselEN0@penalized
   whichEN <- which(coefEN0 != 0)
   penselEN <- penalized(responsemin,XMw0[,whichEN,drop=FALSE],lambda2=optl, unpenalized=nopen,data=datapred,
                         trace=FALSE,maxiter=100)
  coefEN <- penselEN@penalized
  predobj <- c(predobj,penselEN)
  resEN <- list(whichEN=whichEN,betasEN=coefEN)
}


if(compareunpenal){
if(model=="survival"){
  print("Starting unpenalized Cox-model")
  bogus <- matrix(rnorm(nsam),ncol=1)
  print(dim(datapred))
  penlambdas0 <- penalized(response,penalized = bogus,unpenalized = nopen,lambda1=0,lambda2=10^8, data=datapred) 
  predobj <- c(predobj,penlambdas0)
} else {
if(model == "logistic") famglm <- "binomial" else famglm <- "gaussian"
print("Starting unpenalized glm")
form <- formula(paste("response","~",as.character(unpenal)[2]))
modelglm <- glm(form,family=famglm,data=dataunpen)
predobj <- c(predobj,list(modelglm))
}
}

printlam <- function(lambs) {if(length(lambs)<=10) return(lambs) else return(summary(lambs))}

suml <- lapply(lambdas,printlam)
print("Final lambda multipliers (summary):")
print(suml)
print(paste("CVLs",cvlnstot))
timetot <- proc.time()-pmt0
print(paste("Total computation time:",timetot[3]))
names(predobj) <- nmp
colnames(almvecall) <- nmpweight
if(savepredobj=="last") {
predobj <- predobj[length(predobj)]
almvecall <- matrix(lmvecall,ncol=1)
}
if(savepredobj=="none") predobj <- NULL
return(list(true=response,cvls = cvlnstot,lambdamults = lambdas, optl=optl, lambdamultvec = almvecall, 
            predobj=predobj,betas=newbeta, whichsel = whichsel,cvlssel = cvlssel,reslasso=reslasso, 
            resEN = resEN, model=model, arguments=arguments,allpreds=allpreds))
}
