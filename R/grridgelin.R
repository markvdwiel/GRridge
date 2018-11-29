.grridgelin <- function(highdimdata, response, partitions, unpenal = ~1, 
                    offset=NULL, method="exactstable",
                    niter=10, monotone=NULL, optl=NULL, innfold=NULL, 
                    fixedfoldsinn=TRUE, maxsel=c(25,100),selectionEN=FALSE,cvlmarg=1,
                    savepredobj="all", dataunpen=NULL, ord = 1:length(partitions),
                    comparelasso=FALSE,optllasso=NULL,cvllasso=TRUE,
                    compareunpenal=FALSE,trace=FALSE,modus=1,EBlambda=FALSE,standardizeX = TRUE){#
  # highdimdata=simdata; response=Y; partitions=part5; unpenal = ~1; innfold=10
  # offset=NULL; method="exactstable";
  # niter=10; monotone=NULL; optl=NULL; innfold=NULL;
  # fixedfoldsinn=TRUE; maxsel=c(25,100);selectionEN=FALSE;cvlmarg=1;
  # savepredobj="all"; dataunpen=NULL; ord = 1:length(partitions);
  # comparelasso=FALSE;optllasso=NULL;cvllasso=TRUE;
  # compareunpenal=FALSE;trace=FALSE;modus=1;
  # EBlambda=FALSE;standardizeX = TRUE
  #grl <- .grridgelin(highdimdata=simdata, response=Y, partitions=part5,maxsel=25, unpenal = ~1, innfold=10,selectionEN=TRUE)
  print("Using GLMnet for fitting")

  if(standardizeX) {
    print("Covariates are standardized")
    sds <- apply(highdimdata,1,sd)
    sds2 <- sapply(sds,function(x) max(x,10^{-5}))
    highdimdata <- (highdimdata-apply(highdimdata,1,mean))/sds2
  }
  
  nsam <- ncol(highdimdata)
  
  
  if(method=="adaptridge" | method== "exact") niter <- 1
  
  if(class(partitions[[1]]) =="integer"){
    partitions=list(group=partitions)
    ord <- 1
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
    #ncl <-2
    indexset <- unlist(partitions[[ncl]])
    if(length(indexset) < nr){
      print(paste("Warning: partition",ncl,"does not contain all row indices of the data"))
    }
    if(max(indexset) > nr | min(indexset)<1 | is.na(max(indexset))){
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
                    fixedfoldsinn=fixedfoldsinn,
                    selectionEN=selectionEN,maxsel=maxsel,
                    cvlmarg=cvlmarg, dataunpen=dataunpen,savepredobj=savepredobj, ord=ord, 
                    comparelasso=comparelasso, optllasso=optllasso, 
                    compareunpenal=compareunpenal, modus=modus, EBlambda=EBlambda,standardizeX=standardizeX)
  
  
  if(nr > 10000 & is.null(innfold)) print("NOTE: consider setting innfold=10 to save computing time")
  
  nmp0 <- names(partitions)
  if(is.null(nmp0)) nmp0 <- sapply(1:length(partitions),function(i) paste("Grouping",i))  
  nmp0 <- sapply(1:length(partitions),function(i) {
    if(nmp0[i]=="") return(paste("Grouping",i)) else return(nmp0[i])})
  
  nmp <- c("NoGroups","GroupRegul")
  nmpweight <- nmp #to be used later
  
  if(comparelasso) nmp <- c(nmp,"lasso") 
  #new 29/11
  if(selectionEN) nmp <- c(nmp,paste("EN",maxsel,sep=""))
  if(compareunpenal) nmp <- c(nmp,"modelunpen")
  
  if(class(response) =="factor") {
    nlevel <- length(levels(response))
    if(nlevel != 2){
      print("Response is not binary, so not suitable for two-class classification.")
      return(NULL)
    } else {
      model = "logistic"
      fam="binomial"
      print("Binary response, executing logistic ridge regression")
      lev <- levels(response)
      print(paste("Predicting probability on factor level",lev[2]))
    }} else {
      if(class(response) == "numeric" | class(response)=="integer"){
        valresp <- sort(unique(response)) 
        if(length(valresp)==2 & valresp[1]==0 & valresp[2]==1) {
          model = "logistic"
          fam="binomial"
          print("Binary response, executing logistic ridge regression")
        } else {
          model = "linear"
          fam="gaussian"
          print("Numeric continuous response, executing linear ridge regression")
          # print("RESPONSE IS STANDARDIZED")
          # response <- (response-mean(response))/sd(response)
        }
      } else {
        if(class(response) == "Surv"){
          model="survival"
          fam="cox"
#          print("Survival response, executing cox ridge regression")
          print("Survival response is not yet included in GRridge2. It will be soon.")
          return(NULL)
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
  
  if(is.null(offset)) offset <- rep(0,nsam)
  
  interc <- TRUE
  if(unpenal == ~0) interc <- FALSE
  
  if((is.null(dataunpen)) | (unpenal == ~0) | (unpenal == ~1)) {
    X0 <- t(highdimdata)
    pf <- rep(1,nr)
    mm <- NULL
    nunpen <- 0
  } else  {
    mm <- model.matrix(unpenal,dataunpen)
    if(prod(mm[,1]==rep(1,nsam))==1) {
      interc <- TRUE
      mm <- mm[,-1,drop=FALSE]
    } else {
      interc <- FALSE
    }
    nunpen <- ncol(mm)
    pf <- c(rep(1,nr),rep(0,nunpen))
  }
  X0 <- cbind(t(highdimdata),mm) 
  
  #EB estimation lambda, linear model
  mlestlin <- function(Y,X){
    maxv <- var(Y)
    sim2 = function(ts){
      #ts <- c(log(0.01),log(maxv));X <- highdimdata
      tausq<-ts[1];sigmasq<-ts[2]
      n<- nrow(X)
      varY = X %*% t(X) * exp(tausq) + diag(rep(1,n))*exp(sigmasq)
      mlk <- -dmvnorm(Y,mean=rep(0,n),sigma=varY,log=TRUE)
      return(mlk)
    }
    op <- optim(c(log(0.01),log(maxv)),sim2)
    mlests <- exp(op$par)
    return(mlests)
  }
  
  if(is.null(innfold)) foldinit <- nsam else foldinit <- innfold
  #EBlambda <- TRUE;optl<-NULL
  pmt0<- proc.time()
  optl0 <- optl
  if(is.null(optl)){
    if(EBlambda){ #TO BE UPDATED
      if(model != "linear"){
      print("EB estimation of lambda currently only available for linear models")
      print("Applying CV instead")
      EBlambda <- FALSE
      } else {
      ts <- mlestlin(response,t(highdimdata))
      tausqest <- ts[1]
      sigmasq <- ts[2]
      #division by N/2 makes it suitable for glmnet
      optl <- (sigmasq/tausqest)/(nsam/2)
      print(sigmasq)
      print(tausqest)
      print(paste("lambda2",optl*(nsam/2)))
      time1 <- proc.time()-pmt0
      print(paste("Computation time EB estimation penalty parameter:",time1[3]))
      arguments$optl <- optl
      }
    } 
   opt <- NULL 
    if(!EBlambda){ #CV
      print("Finding lambda for initial ridge regression")
      if(fixedfoldsinn) set.seed(346477)
      #opt <- optL2(response, penalized = t(highdimdata),fold=foldinit,unpenalized=unpenal,data=dataunpen,trace=trace,standardize=TRUE)
      opt <- cv.glmnet(x=X0,y=response,offset=offset,nfolds=foldinit,penalty.factor=pf,alpha=0,standardize=FALSE,family=fam,
                       intercept=interc,keep=TRUE) 
      time1 <- proc.time()-pmt0
      print(paste("Computation time for cross-validating main penalty parameter:",time1[3]))
      optl <- opt$lambda.min
      if(optl == Inf) optl <- 10^10
      print(paste("lambda2 (multiplied by N/2):",optl*nsam/2))
      arguments$optl <- optl
    }
  }
  
  pmt <- proc.time()
  nsam <- ncol(highdimdata)
  nfeat <- nrow(highdimdata)
  XM0 <- X0
  response0 <- response
  
  cvlnstot <- rep(0,(nclass+1)) 
  allpreds <- c()
  whsam <- 1:nsam
  responsemin <- response0
  
  pen0 <- glmnet(x=X0,y=response,offset=offset,nlambda=1,lambda=optl,penalty.factor=pf,alpha=0,
                 family=fam,intercept=interc)
  
  addintercept <- interc
  
  if(is.null(innfold)) {nf <- nrow(XM0);myfold<-NULL} else {
    if(!is.null(optl0) | EBlambda) {  #NEW 7/2017
      nf <- innfold
      myfold <- NULL
      if(fixedfoldsinn) set.seed(346477)
    } else {
      nf <- innfold
      myfold <- opt$foldid
    }
  }
  
  if(is.null(myfold)) opt2 <- cv.glmnet(x=X0,y=response,offset=offset,nfolds=nf,lambda=c(optl,optl/2),
                                        penalty.factor=pf,alpha=0,family=fam,intercept=interc,keep=TRUE) else   
                                          opt2 <- cv.glmnet(x=X0,y=response,offset=offset,nfolds=nf,lambda=c(optl,optl/2),
                                                            foldid=myfold,penalty.factor=pf,alpha=0,family=fam,intercept=interc,keep=TRUE) 
  loss <- opt2$cvm[1]
  myfold <- opt2$foldid
  cvln0 <- loss
  
  #needed because glmnet lambda equals real lambda/(N/2)
  myoptl <- optl*nsam/2
  
  #betaexact <- solve(t(X0) %*% X0 + myoptl*diag(rep(1,nr))) %*% t(X0) %*% matrix(response,nrow=nsam)

#################### BEGIN LOOP ################   
  
  cvlnprev <- cvln0
  penprev <- pen0
  pen <- pen0
  print(cvln0)
  #contains highdim data only
  XMw0 <- t(highdimdata)
  XMw0prev <- t(highdimdata)
  Xglmnet <- cbind(XMw0,mm)
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
        coeff <- penprev$beta[1:nr]
        #preds <- as.numeric(predict(penprev,Xglmnet,s=c(optl),type="response",offset=offset))
        
        #changed 31/10/2018
        preds <- try(as.numeric(predict(penprev,Xglmnet,s=c(optl),type="response",offset=offset)),silent=T)
        if(class(preds) == "try-error") preds <- as.numeric(predict(penprev,Xglmnet,s=c(optl),type="response",newoffset=offset))
        
        
      coeffsq <- coeff^2
      
      if(model=="logistic") {
        Wi <- sqrt(preds*(1-preds)) 
        constlam <- 2
      }
      if(model == "linear"){
        Wi <- rep(1,length(preds))
        constlam <- 1
      }
      if(model == "survival"){  #TODO
        resptime <- response[,1]
        predsnew <- -log(sapply(1:nsam,function(k) survival(preds,time=resptime[k])[k]))
        Wi <- sqrt(predsnew)
        constlam <- 2
      }
      
      if(!((is.null(dataunpen)) | (unpenal == ~0) | (unpenal == ~1))) {
        #mm <- model.matrix(unpenal,dataunpen) 
        XMW <- t(t(cbind(XMw0,10^5*mm)) %*% diag(Wi)) 
      } else {
        if(addintercept) XMW <- t(t(cbind(XMw0,rep(10^5,nsam))) %*% diag(Wi)) else XMW <- t(t(XMw0) %*% diag(Wi))
      }
      pmt <- proc.time()
      SVD <- svd(XMW)
      proc.time()-pmt
      leftmat <- SVD$v %*% diag(1/((SVD$d)^2+constlam*myoptl)) %*% diag(SVD$d) %*% t(SVD$u)
      
      if(model=="linear"){
        #depreciated because bad estimator; 
        # Hatm <- XMW %*% leftmat
        # df <- nsam - sum(diag(2*Hatm - Hatm %*% t(Hatm)))
        # VarRes <- sum((response - preds)^2)/df
        #NEW: 10-7-2017
        VarRes <-  mlestlin(response,XMW)[2]
        print(paste("Sigma^2 estimate:",VarRes))
        #print(paste("Sigma^2 estimate:",sigmasq))
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
        # 
        # opt2w <- cvl(responsemin, penalized = XMw0,fold=nf,lambda2=optl,
        #              unpenalized=unpenal,data=datapred, trace=trace)
        
        #vervang door data.frame(XMw0,mm)
        opt2w <- cv.glmnet(x=cbind(XMw0,mm),y=responsemin,offset=offset,nfolds=nf,lambda=c(optl,optl/2),foldid=myfold,
                           penalty.factor=pf,alpha=0,family=fam,intercept=interc,standardize=FALSE) 
        loss <- opt2w$cvm[1]
        cvln1 <- loss
        print(cvln1)
        
        if((cvln1 - cvlnprev)/abs(cvlnprev) < -1/100  |  ((cvln1 - cvlnprev)/abs(cvlnprev) <= 0 & i==1)){ 
          
          pen <- glmnet(x=cbind(XMw0,mm),y=response,offset=offset,nlambda=1,lambda=optl,penalty.factor=pf,alpha=0,family=fam,
                        intercept=interc,standardize=FALSE)
          
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
          conv[cl] <- TRUE
          cvln1 <- cvlnprev
          XMw0 <- XMw0prev
          pen <- penprev
          if(niter>1) print(paste("Partition",nmp0[cl],"CONVERGED after",i,"iterations")) else print(paste("Partition",nmp0[cl],"did not improve results"))
        }
      }
      } #end if(!convcl)
      cl <- cl+1
    }
    if(sum(conv)==nclass){
      converged <- TRUE
      if(niter>1) print(paste("All partitions CONVERGED after",i,"iterations"))
    }
    i <- i+1
  }
  
#################### END LOOP ################  
  
Xglmnet <- cbind(XMw0,mm)

  #changed 31/10  
if(niter==0) {pen <- pen0;XMw0<-XM0;cvln1<-cvln0;soltau <- NULL}
pred0 <- try(as.numeric(predict(pen0,XM0,s=c(optl),offset=offset,type="response")),silent=T)
predw <- try(as.numeric(predict(pen,Xglmnet,s=c(optl),offset=offset,type="response")),silent=T)
             
if(class(pred0) == "try-error"){
pred0 <- predict(pen0,XM0,s=c(optl),newoffset=offset,type="response") 
predw <- predict(pen,Xglmnet,s=c(optl),newoffset=offset,type="response") 
}


predshere <- cbind(pred0,predw)
cvlnssam <- c(cveridge=cvln0,cvegrridge=cvln1)
lmvecall <- lmvec
almvecall <- cbind(almvecall,lmvecall)
predobj <- list(pen0)
predobj <- c(predobj,list(pen))

allpreds <- predshere
whichsel <- NULL
betassel <- NULL

npr <- length(predobj)
pred2 <-predobj[[npr]]
lambs <- lmvecall

#only the highdim variables
oldbeta <- pred2$beta[1:nr]
newbeta <- oldbeta/sqrt(lambs)
time2 <- proc.time()-pmt
print(paste("Computation time for adaptive weigthing:",time2[3]))


cvlnstot <- cvlnssam
reslasso <- NULL
if(comparelasso){
print("Starting lasso")
    if(is.null(optllasso)){
    print("Finding lambda for lasso regression")
    if(fixedfoldsinn) set.seed(346477)
    #alpha=1 implies lasso  
    opt <- cv.glmnet(x=X0,y=response,offset=offset,nfolds=nf,foldid=myfold,penalty.factor=pf,alpha=1,standardize=FALSE,family=fam,
                     intercept=interc) 
    optllasso <- opt$lambda.min
    whmin <- which(opt$lambda==optllasso)
    print(paste("lambda1 (multiplied by N):",optllasso*nsam))
    arguments$optllasso <- optllasso
    cvliklasso <- opt$cvm[whmin]
    } else {
    cvliklasso0 <- if(cvllasso) try(cv.glmnet(x=X0,y=response,offset=offset,nfolds=nf,lambda=c(optllasso,optllasso/2),
                                              foldid=myfold, penalty.factor=pf,alpha=1,family=fam,intercept=interc,
                                              standardize=FALSE) )
    if(class(cvliklasso0) == "try-error" | !cvllasso) cvliklasso <- NA else cvliklasso <- cvliklasso0$cvm[1]
        }
cvlnstot <- c(cvlnstot,cvelasso=cvliklasso)
penlasso <- glmnet(x=X0,y=response,offset=offset,nlambda=1,lambda=optllasso,penalty.factor=pf,alpha=1,family=fam,
                           intercept=interc,standardize=FALSE)
betaspenalizedlasso <- penlasso$beta[1:nr]  
whichlasso <- which(betaspenalizedlasso != 0)
betaslasso <- betaspenalizedlasso[whichlasso]
predobj <- c(predobj,list(penlasso))
reslasso <- list(cvllasso=cvliklasso,whichlasso=whichlasso,betaslasso=betaslasso)
print(paste("lasso uses",length(whichlasso),"penalized variables"))
}

resEN <- list() 

#one cannot select more than the nr of variables
#new 29/11
if(selectionEN){
    print("Variable selection by elastic net started...")  
    for(maxsel0 in maxsel){
      maxsel2 <- min(maxsel0,nr)
      print(paste("Maximum nr of variables",maxsel2))
      
      fsel <- function(lam1,maxselec=maxsel2,lam2){
        if(lam1==0) return(nfeat-maxselec) else {
          alp = lam1/(lam1+lam2)
          lam <- lam1+lam2
          penselEN <- glmnet(x=cbind(XMw0,mm),y=response,offset=offset,nlambda=1,lambda=lam,
                             penalty.factor=pf,alpha=alp,family=fam, intercept=interc,standardize=FALSE)
          coef <- penselEN$beta[1:nr]
          return(length(coef[coef!=0])-maxselec)
        }
      }
      lam1 <- uniroot(fsel,interval=c(0,optl),maxiter=50,lam2=optl)$root
      
      cvEN <- cv.glmnet(x=cbind(XMw0,mm),y=response,offset=offset,nfolds=nf,
                foldid=myfold, penalty.factor=pf,alpha=lam1/(lam1+optl),family=fam,intercept=interc,
                standardize=FALSE)
      
      whmin <- which(cvEN$lambda==cvEN$lambda.min)
      #glmnet output includes number of non-zero unpen covariates
      nsel <-  cvEN$nzero[whmin]-nunpen
      
      #if model with fewer than maxsel vars is better, use that
      if(nsel > maxsel2) lambdaEN <- lam1+optl else lambdaEN <- cvEN$lambda.min
      penselEN0 <- glmnet(x=cbind(XMw0,mm),y=response,offset=offset,penalty.factor=pf,nlambda=1,lambda=lambdaEN,
                     alpha=lam1/(lam1+optl),family=fam, intercept=interc,standardize=FALSE)
              
      coefEN0 <- penselEN0$beta[1:nr]
      whichEN <- which(coefEN0 != 0)
      nselEN <- length(whichEN)
      print(paste("Model with EN selection uses",nselEN,"penalized variables"))
      
      #now re-fit GRridge model with selected variables only
      penselEN <- glmnet(x=cbind(XMw0[,whichEN,drop=FALSE],mm),y=response,offset=offset,
                         penalty.factor=pf,nlambda=1,lambda=optl,alpha=0,family=fam, intercept=interc,standardize=FALSE)
      coefEN <- penselEN$beta[1:nselEN]
      predobj <- c(predobj,list(penselEN))
      resEN <- c(resEN,list(list(whichEN=whichEN,betasEN=coefEN)))
  }
  names(resEN) <- paste("resEN",maxsel,sep="")
}


if(compareunpenal){
if(model=="survival"){
  print("Starting unpenalized Cox-model")
  bogus <- matrix(rnorm(nsam),ncol=1)
  print(dim(datapred))
  penlambdas0 <- penalized(response,penalized = bogus,unpenalized = unpenal,lambda1=0,lambda2=10^8, data=datapred) 
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
print("CV error (high-dim models only):")
print(cvlnstot)
timetot <- proc.time()-pmt0
print(paste("Total computation time:",timetot[3]))
names(predobj) <- nmp
colnames(almvecall) <- nmpweight
if(savepredobj=="last") {
predobj <- predobj[length(predobj)]
almvecall <- matrix(lmvecall,ncol=1)
}
if(savepredobj=="none") predobj <- NULL

#added mm (model matrix), to be used by grridgeCV
return(list(true=response,cvfit = -cvlnstot,lambdamults = lambdas, optl=optl, lambdamultvec = almvecall, 
            predobj=predobj,betas=newbeta, whichsel = whichsel,reslasso=reslasso, 
            resEN = resEN, model=model, arguments=arguments,allpreds=allpreds, mm=mm))
}
