grridge <-
function(highdimdata, response, partitions, unpenal = ~1, offset=NULL, method="exactstable", niter=10, monotone=NULL, optl=NULL, innfold=NULL, 
fixedfoldsinn=TRUE, selection=FALSE,maxsel=100,stepsel=1,cvlmarg=1,savepredobj="all", dataunpen=NULL, ord = 1:length(partitions),
comparelasso=FALSE,optllasso=NULL,cvllasso=T,compareEN=F,compareunpenal=FALSE,trace=FALSE,modus=1){#
#method: "stable", "exact","adaptridge"
#highdimdata<-betasfcen;response<- fl;unpenal=~0;offset=NULL;partitions<-partitionsVerlaat;monotone = c(TRUE,FALSE);
#selection=TRUE;savepredobj="all"; optl=3;innfold=NULL; selection=TRUE;maxsel=100;stepsel=1;cvlmarg=1
#dataunpen=NULL;niter<-10;ord <- 1:2
#partitions <- list(partitionsVerlaat[[2]],partitionsVerlaat[[1]]);monotone=c(FALSE,TRUE)
#highdimdata <- t(X);response <- y;calssif <- partitions; unpenal = ~0;printl=FALSE; niter=10; monotone=NULL; optl=NULL; innfold=NULL;selection=FALSE;maxsel=100;stepsel=1;savepredobj="last"; dataunpen=NULL
#highdimdata <- gesel;response <- nodestat;partitions <- clroepman; printl=T;printp=T;unpenal = ~0; niter=10; monotone=NULL; optl=NULL; innfold=10;selection=FALSE;maxsel=100;stepsel=1;savepredobj="last"; dataunpen=NULL
#response: factor with two levels (logistic ridge regression) or numeric (linear ridge regression, except when only 0-1 values appear). 
#savepredobj = "last": save only last pred obcject; "all": save all; "none": save none
#highdimdata <- betasfcen;response<-fl;monotone <- c(TRUE);niter <- 10;partitions <- partitionwina;innfold=NULL;
#unpenal=~0+hpv16;selection=TRUE;maxsel=100;stepsel=1;optl=37.32;dataunpen=hpv16dat
#addintercept=TRUE;addunpenal=TRUE
#... futher arguments for penalized, optL2 and cvl, such as unpenalized
#model is determined form the response
#cvlmarg: Margin for CVL in percentages when selecting
#printl: print update for lambda multiplier
#niter: setting it to zero renders results for ordinary ridge regression
#partitions<-partitionsVerlaat;monotone = c(TRUE,FALSE,FALSE);ord <- c(1,2,3);highdimdata<-datcenVerlaat
#highdimdata<-datcenFarkas;response<-respFarkas;partitions<- partitionFarkas;offset=NULL;monotone = NULL;
#selection=TRUE;savepredobj="all";optl=5.68;innfold=NULL; selection=TRUE;maxsel=100;stepsel=1;cvlmarg=1
#dataunpen=NULL;niter<-10;ord <- 1;unpenal = ~1;fixedfoldsinn=TRUE;method<-"adapt ridge";controlbound<-5
#highdimdata<-mirnormcen;response<-resp;partitions<- partkeep;offset=NULL;monotone = c(TRUE,TRUE);
#selection=FALSE;savepredobj="all";optl=192;innfold=NULL; maxsel=100;stepsel=1;cvlmarg=1
#dataunpen=NULL;niter<-10;ord <- 1:2;unpenal = ~1;fixedfoldsinn=TRUE;comparelasso=TRUE;optllasso=NULL
#highdimdata<-mirnormcen;response<-resp;partitions<- partkeep;offset=NULL;monotone = c(TRUE);
#selection=FALSE;savepredobj="all";optl=664;innfold=10; maxsel=100;stepsel=1;cvlmarg=1;compareunpenal=T;
#dataunpen=data.frame(therapy=therapy);niter<-10;ord <- 1;unpenal = ~1+therapy;fixedfoldsinn=TRUE;comparelasso=TRUE;optllasso=7.06;

#see grridgeExample for creation of data

#sizes = c(5000,2000,1300,30,50,300)
#selectuneq <- sapply(1:6,function(i){part <- partitionFarkas[[1]][[i]]; return(sample(part,size=sizes[i]))})
#wh <- unlist(selectuneq)
#datcenF <- datcenF
#CpGF<- CpGannFarkas
#partitionF <- list(cpg=CreatePartition(CpGF))
#
  
# data(dataFarkas);partitionFarkas <- list(cpg=CreatePartition(CpGannFarkas))
# data(dataFarkas);partitionFarkasr <- list(cpg=CreatePartition(sample(CpGannFarkas)))
# highdimdata<-datcenFarkas;response<-respFarkas;partitions<- partitionFarkas;offset=NULL;monotone = FALSE;
# selection=TRUE;savepredobj="all";innfold=10; selection=FALSE;maxsel=100;stepsel=1;cvlmarg=1;
# dataunpen=NULL;niter<-10;ord <- 1;unpenal = ~0;fixedfolds=TRUE;compareunpenal=F;comparelasso=F;optllasso=NULL;
# controlbound<-1000;optl=5.68;trace=FALSE;method="stable"; modus=1;
  
# library(GRridge); load("C:/VUData/Wina/MethInfFull/partitionsVerlaat.Rdata"); highdimdata<-datcenVerlaat; response<-respVerlaat;unpenal=~0;partitions <- partitionsVerlaat;
# monotone = c(TRUE,FALSE,FALSE,TRUE);selection=T;optl=548; savepredobj="all";modus=1;method="stable";niter=5;
# dataunpen=NULL;ord <- 1:length(partitions);fixedfoldsinn=TRUE;compareunpenal=F;comparelasso=T;optllasso=NULL;offset=NULL;
# controlbound<-1000;trace=FALSE;maxsel=25;stepsel=1;cvlmarg=1;savepredobj="all";innfold=10;

  
# setwd("C:\\VUData\\Maarten\\AnalysePlanProject2");load("GrridgePFSMaarten.Rdata");load("mirnorms.Rdata")
# highdimdata<-mirnormcenTS_PFS;response<-survPFS;outerfold=10;fixedfolds<-T;partitions=partkeep;ord <- 1:length(partitions);
# cvlmarg=1;offset = NULL; fixedfoldsinn=TRUE;dataunpen = datfr;
# unpenal = ~1 + adjth + thscheme + age + neoadj + pcrcdiff;selection=T;stepsel=2;monotone=c(TRUE,TRUE);niter=1;
# savepredobj = "all";innfold=10;method="exact";
# optl=NULL;maxsel=100;comparelasso=TRUE;optllasso=NULL;compareunpenal=TRUE;trace=FALSE;modus=1
#   
# load("C:\\ExternData\\NeuroBlastoma\\seqNB.Rdata")
# library(GSEABase)
# oncocursym <- getGmt("C:\\ExternData\\CensusGenes\\GSEA\\OncoCurated.all.v5.0.symbols.gmt")
# # library(GRridge); source('C:/Synchr/RPackages/GRridge/V21/GRridge/R/matchGeneSets.R')  
# sds <- apply(seqTrcen100,1,sd);genes <- seqann[,1];genesets <- oncocursym;
#   partitions=list(gs=matchGeneSets(genes,genesets,remain=T),sd=CreatePartition(sds,decreasing=FALSE,uniform=T,grsize=5000));
# highdimdata<-seqTrcen100;response<-clinOSTr100;outerfold=5;fixedfolds<-T;ord <- 1:length(partitions);
#   cvlmarg=1;offset = NULL; fixedfoldsinn=TRUE;dataunpen=NULL;
#   unpenal = ~1;selection=F;stepsel=2;monotone=c(FALSE,FALSE);niter=3;
#   savepredobj = "all";innfold=5;method="stable";
#    optl=7500; maxsel=100;comparelasso=F;optllasso=NULL;compareunpenal=TRUE;trace=TRUE;modus=1;
  
# LINEAR CASE 
# setwd("C:/VUData/Tonje/"); datcen <- read.table("methylValues.txt"); pvals <-  read.table("wValues_LowMeansRelevant.txt")[,1];resp <- read.table("Y.txt")[,1]
# firstPartition <- CreatePartition(pvals,grsize=10,uniform=T,decreasing=FALSE); partitionP <- list(wPvals=firstPartition)
# highdimdata<-datcen;response <- resp;partitions<- partitionP;monotone=c(TRUE);savepredobj="all";selection=FALSE;
# outerfold=5;fixedfolds<-T;ord <- 1:length(partitions);
#     cvlmarg=1;offset = NULL; fixedfoldsinn=TRUE;dataunpen=NULL;
#     unpenal = ~1;niter=3; savepredobj = "all";innfold=5;method="stable";stepsel <- 5;
#      optl=NULL; maxsel=100;comparelasso=F;optllasso=NULL;compareunpenal=TRUE;trace=TRUE;modus=1 
# 
  
# PATHWAYS 
# load("C:/VUData/Wurdinger/WurdingerObjects.RData"); data <- WurdingerObjects$dataSqrt; 
# datamc <- data - apply(data,1,mean);sds <- apply(data,1,sd);datamc <- datamc/sds;
#   sdpartnew <- CreatePartition(sds,uniform=TRUE,grsize=1000,decreasing=FALSE)
# resp <- WurdingerObjects$group; part <- WurdingerObjects$partitions;part2 <- part[c(1,5,6,7)] #use only the 'pathways'
# highdimdata<-datamc;response <- resp;partitions<- part2;monotone=c(TRUE,FALSE,FALSE,FALSE);savepredobj="all";selection=FALSE;
# outerfold=5;fixedfolds<-T;ord <- 1:length(partitions);
#     cvlmarg=1;offset = NULL; fixedfoldsinn=TRUE;dataunpen=NULL;
#     unpenal = ~1;niter=3; savepredobj = "all";innfold=5;method="exact";stepsel <- 5;
#      optl=780; maxsel=100;comparelasso=F;optllasso=NULL;compareunpenal=FALSE;trace=TRUE;modus=1 
#  
#  RANDOM groups 
#   load("C:/VUData/Wurdinger/WurdingerObjects.RData"); data <- WurdingerObjects$dataSqrt; 
#   datamc <- data - apply(data,1,mean);sds <- apply(data,1,sd);datamc <- datamc/sds;
#   resp <- WurdingerObjects$group; 
#  nfeat <- nrow(datamc);ngr <- 10; npgr <- ceiling(nfeat/ngr);sams <- sample(rep(1:ngr,npgr))[1:nfeat]
#  sams <- sample(rep(1:ngr,npgr))[1:nfeat]; partrand=list(rand=CreatePartition(factor(sams)))
#   highdimdata<-datamc;response <- resp;partitions<- partrand;monotone=c(FALSE);savepredobj="all";selection=FALSE;
#   outerfold=5;fixedfolds<-T;ord <- 1:length(partitions);
#       cvlmarg=1;offset = NULL; fixedfoldsinn=TRUE;dataunpen=NULL;
#       unpenal = ~1;niter=3; savepredobj = "all";innfold=5;method="stable";stepsel <- 5;
#        optl=780; maxsel=100;comparelasso=F;optllasso=NULL;compareunpenal=FALSE;trace=TRUE;modus=1 
#   
#   highdimdata<- datStd;response <- resp; partitions<-parSym; monotone = c(TRUE,FALSE);innfold=10;
#  selection=FALSE;savepredobj="all";method="exactstable";stepsel <- 5;optl=115;niter=1;unpenal = ~1;
# maxsel=100;comparelasso=T;optllasso=NULL;compareunpenal=FALSE;trace=TRUE;modus=1;
# fixedfolds<-T;ord <- 1:length(partitions); cvlmarg=1;offset = NULL; fixedfoldsinn=TRUE;dataunpen=NULL;
  
  #### START NEW2 ###

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
  
  #check whether each partition is correct; update 30/11/2015: is allowed now
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
  
  #### END NEW2 ###
  
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
      #     pmt <- proc.time()
      #     nover <- sapply(1:nr,function(x) sum(indexset==x)) #can possibly be done more efficiently
      #     proc.time()-pmt
      #     pmt <- proc.time()
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
  
  #NEW
  arguments <- list(partitions=partitions,unpenal=unpenal, offset=offset, method=method, niter=niter, monotone=monotone, optl=optl, innfold=innfold, fixedfoldsinn=fixedfoldsinn, 
                    selection=selection, maxsel=maxsel,stepsel=stepsel, cvlmarg=cvlmarg, dataunpen=dataunpen,savepredobj=savepredobj, ord=ord, 
                    comparelasso=comparelasso, optllasso=optllasso, compareEN=compareEN, compareunpenal=compareunpenal, modus=modus)
  
  
  if(nr > 10000 & is.null(innfold)) print("NOTE: consider setting innfold=10 to save computing time")
  
  nmp0 <- names(partitions)
  if(is.null(nmp0)) nmp0 <- sapply(1:length(partitions),function(i) paste("Grouping",i))  
  nmp0 <- sapply(1:length(partitions),function(i) {
    if(nmp0[i]=="") return(paste("Grouping",i)) else return(nmp0[i])})
  
  
  nmp <- c("NoGroups","GroupRegul")
  nmpweight <- nmp #to be used later
  if(selection) nmp <- c(nmp,"ForwSel") 
  if(comparelasso) nmp <- c(nmp,"lasso") #NEW
  if(compareEN) nmp <- c(nmp,"EN") #NEW
  if(compareunpenal) nmp <- c(nmp,"modelunpen") #NEW
  
  
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
    if(is.null(dataunpen)) {print("If unpenal contains variables, data of the unpenalized variables should be specified in the data slot!")
      return(NULL)
    } 
  }
  
  
  nsam <- ncol(highdimdata)
  #offset <- log(1/9);unpenal= ~0;nsam=10
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
  optl0 <- optl #optl0 for later use
  if(is.null(optl)){
    print("Finding lambda for initial ridge regression")
    if(fixedfoldsinn) set.seed(346477)
    opt <- optL2(response, penalized = t(highdimdata),fold=foldinit,unpenalized=nopen,data=datapred,trace=trace)
    #opt <- optL2(response, penalized = t(highdimdata),fold=foldinit,trace=T,unpenalized=nopen,data=datapred,...)
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
  
#NEW2  if(is.null(monotone)) monotone <- rep(FALSE,nclass)
  if(!selection) cvlnstot <- rep(0,(nclass+1)) else cvlnstot <- rep(0,(nclass+2))
  allpreds <- c()
  whsam <- 1:nsam
  
  
  #whsam <- c(11,12,31,35)
  
  responsemin <- response0
  pen0 <- penalized(responsemin, penalized = XM0, lambda2 = optl, unpenalized=nopen,data=cbind(XM0,datapred))
  nmunpen <- names(pen0@unpenalized)
  if(is.element("(Intercept)",nmunpen)) addintercept <- TRUE else addintercept <- FALSE #needed later for matrix operations
  #pen0 <- penalized(responsemin, penalized = XM, lambda2 = optl, unpenalized=nopen,data=datapred,...)
  if(is.null(innfold)) {nf <- nrow(XM0)} else {
    if(!is.null(optl0)) {
      nf <- innfold
      if(fixedfoldsinn) set.seed(346477)
    } else {nf <- opt$fold}
  }
  opt2 <- cvl(responsemin, penalized = XM0,fold=nf, lambda2 = optl,unpenalized=nopen,data=datapred,trace=trace)
  nf <- opt2$fold  #folds are fixed for all cvs hereafter
  #opt2 <- cvl(responsemin, penalized = XM,fold=nf, lambda2 = optl, unpenalized=nopen,data=datapred,...)
  cvln0 <- opt2$cvl
  
  cvlnprev <- cvln0
  penprev <- pen0
  pen <- pen0
  print(cvln0)
  XMw0 <- XM0
  XMw0prev <- XM0
  converged <- FALSE
  conv <- rep(FALSE,nclass)
  controlbound1 <-1000 #maximal relative deviation of the "exactstable" estimates
  controlbound2 <- controlbound3 <- 10 #maximal relative deviation of the "stable" and "exact" estimates
  
  almvecall <- rep(1,nfeat) #lambda multipliers for unweighted ridge 
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
        lenggr <- unlist(lapply(whgr,length))  #NEW
        ngroup1 <- length(whgr)
        names(lambdas[[cl]]) <- names(whgr)
        coeff <- penprev@penalized
        if(model == "survival"){
          preds <- predict(penprev,XMw0,data=datapred)
        } else {
          preds <- predict(penprev,XMw0,data=datapred)[1:nsam]  #needed because linear model returns matrix   
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
          mm <- model.matrix(nopen,dataunpen) #model.matrix automatically adds the intercept
          XMW <- t(t(cbind(XMw0,10^5*mm)) %*% diag(Wi)) #10^8: large number that 'undoes' the generic penalty 
        } else {
          if(addintercept) XMW <- t(t(cbind(XMw0,rep(10^5,nsam))) %*% diag(Wi)) else XMW <- t(t(XMw0) %*% diag(Wi))
        }
        
        SVD <- svd(XMW)
        
        leftmat <- SVD$v %*% diag(1/((SVD$d)^2+constlam*optl)) %*% diag(SVD$d) %*% t(SVD$u)
        
        
        #XMWn <- SVD$u %*% diag(SVD$d) %*% t(SVD$v) #checked for large matrix (80x18845: fine, equals XMW)
        
        if(model=="linear"){
          Hatm <- XMW %*% leftmat 
          df <- nsam - sum(diag(2*Hatm - Hatm %*% t(Hatm)))
          VarRes <- sum((response - preds)^2)/df
          print(paste("Sigma^2 estimate:",VarRes))
          vars3 <- VarRes*rowSums(leftmat^2)
        } else {
          vars3 <- rowSums(leftmat^2)
        }
        
        #to solve numerical problems
        which0 <- which(vars3==0)
        vars3[which0] <- 10^{-30}
        
        #XMW <- XMW  #checked; coincides with svd solution above
        #p1 <- t(XMW) %*% XMW + constlam*optl*diag(x=1,nrow=P)
        #inv1 <- solve(a=p1,b=diag(x=1,nrow=P))
        #vars1 <- inv1 %*% t(XMW) %*% XMW %*% inv1
        #vs <- VarRes*diag(vars1)
        
        #BELOW; does not apply to survival
        if(model=="linear"){
          mycoeff2svd <- (leftmat %*% response)^2 #for linear ridge the formula is exact
        } 
        if(model=="logistic"){
          if(is.factor(response)) respnum <- as.numeric(response)-1 else respnum <- response
          z <- matrix(log(preds/(1-preds))+(respnum - preds)/(preds*(1-preds)),ncol=1) 
          #
          #if(modus ==2) mycoeff2svd <- (leftmat %*% diag(Wi) %*% z)^2 
          if(modus == 1) mycoeff2svd <- coeffsq
          if(modus == 2) mycoeff2svd <- (leftmat %*% z)^2 
        }
        
        
        #NEW 11/1/2016
        if(model=="survival"){
          mycoeff2svd <- coeffsq #for linear ridge the formula is exact
        } 
        #
        #coefmatrix3 <- (vars2left %*% XMW)^2
        #coefmatrix2 <- (inv1 %*% t(XMW) %*% XMW)^2
        
        #mycoeff2 <- (inv1 %*% t(XMW) %*% z)^2 #checked: coincides with svd solution
        
        #cii <- diag(coefmatrix2)
        cii2 <- (rowSums(leftmat * t(XMW)))^2  #checked
        
        
        leftmat <- leftmat/sqrt(vars3)
        
        #lefts <- function(group){
        ###group <- 1
        #ind <- whgr[[group]]
        #coefftau2 <- sum(mycoeff2svd[ind]-vars3[ind])
        #return(coefftau2)
        #}
        #leftside <- sapply(1:length(whgr),lefts)
        
        #random:
        #whgr <- lapply(1:ngroup1,function(i){ni <- length(whgr[[i]]);return(sample(1:nr,ni))})
        
        
        
        #NEW 11/1/2016; also includes survival
        #10^(-30): just an arbitrary small number
        lowerleft <- 10^(-30)
        lefts2 <- function(group){
          ##group <- 3
          ind <- whgr[[group]]
          ngr <- lenggr[group]
          coefftau2 <- sum(sapply(mycoeff2svd[ind]/vars3[ind],function(x) max(x,1)))-length(ind)
          return(max(lowerleft,coefftau2/ngr))
        } 
        
        leftside <- sapply(1:length(whgr),lefts2)
        ellarg0 <- length(leftside[leftside>lowerleft])/length(leftside)
        
         if(ellarg0 <=.5){ #NO ITERATION WHEN MANY LEFT SIDES [MORE THAN HALF] EQUAL ZERO 
           print(paste("Partition",nmp0[cl],"NOT ITERATED")) 
           conv[cl] <- TRUE
           cvln1 <- cvlnprev
           XMw0 <- XMw0prev
           pen <- penprev
         } else {
        #use 1:nfeat, because last ones correspond to non-penalized
        lefts2ran <- function(group,randomind){
          ##group <- 1;randomind <- sample(1:nfeat)
          ind <- whgr[[group]]
          ngr <- lenggr[group]
          coefftau2 <- sum(sapply(mycoeff2svd[1:nfeat][randomind[ind]]/vars3[1:nfeat][randomind[ind]],function(x) max(x,1)))-
            length(randomind[ind])
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
        #relerror2 <- sum(sqrt(apply(leftsrancen^2,2,sum)))/(nlefts*sqrt(sum(means^2)))
        
        if(cl==1 & i==1) print(cvln0)
        print(paste("Relative error:",relerror))  
        if(relerror>= 0.1) print("WARNING: large relative error (>=0.1). Consider using larger groups of variable.")
        #lefts2 <- function(group){
        ###group <- 1
        #ind <- whgr[[group]]
        #coefftau2 <- sum(coeffsq[ind]/vars3[ind])-length(ind)
        #return(coefftau2)
        #}
        #leftside <- sapply(1:length(whgr),lefts2)
        #
        
        #leftmat of dimension p*n
        
        #note: XMW also contains additional covariates; therefore add weigts to Wmat to match dimensions; those weigths have no impact
        
        nadd <- ncol(XMW) - ncol(XMw0)
        
        rightmat =  t(t(XMW) * c(Wmat[,cl],rep(1,nadd)))
        
        #store the rightmat products (dim n*n) to avoid recomputing
        rightmats <- lapply(1:length(whgr),function(j){
          rightj2 <- rightmat[,whgr[[j]]]
          rcp <- rightj2 %*% t(rightj2)
          return(rcp)
        })
        
        
        coefmatfast <- t(apply(matrix(1:length(whgr),nrow=length(whgr)),1,  
                               function(i){
                                 #i<-1;j<-1
                                 lefti2  <- leftmat[whgr[[i]],]
                                 #lefti2 <- matrix(c(1,2,3,4,3,8),nrow=3)
                                 lcp <- t(lefti2) %*% lefti2
                                 #lcp <- leftcrosspr[upper.tri(leftcrosspr)]
                                 ckls <- sapply(1:length(whgr),function(j){
                                   #rightj2 <- rightmat[,whgr[[j]]]
                                   ##rightj2 <- matrix(c(1,0,2,4,4,1),nrow=2)
                                   #rcp <- rightj2 %*% t(rightj2)
                                   rcp <- rightmats[[j]]
                                   return(sum(lcp*rcp))
                                 })
                                 return(ckls)
                               }))
        coefmatfast <- coefmatfast/lenggr   
        
        #sum((lefti2 %*% rightj2)^2)
        
        #coefmatFast <- function(L,R){
        #lcp <- t(L) %*% L
        #rcp <- R %*% t(R)
        #return(sum(lcp*rcp))
        #}
        
        #Determines condition numbers (CN)   
        CNfun <- function(lam,cfmmat=coefmatfast){
          ng <- nrow(cfmmat)
          dmax <- max(diag(cfmmat))
          cfmlam <- (1-lam)*cfmmat + lam*diag(dmax,nrow=ng)
          eigenvals <- eigen(cfmlam,only.values=T)$values
          CN <- eigenvals[1]/eigenvals[ng]
          return(Re(CN)) 
        }
        
        #Determine shirnkage factor lam such that CN*relerror <- 0.1, hence relerror in solution <= 0.1
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
          #leftsidetemp <- leftsidefark/lengfark;leftsidetemp <- leftside
          soltau = solve(sum(cfmlam),sum(leftside))
          sol = solve(cfmlam,leftside)
          #sol=mygeninv(coefmatfast/leng) %*% (leftside/leng)
          low <- soltau/controlbound1;up = soltau*controlbound1
          parinint <- sapply(sol,function(x) min(max(low,x),up))
          minopt <- optim(par = parinint,fn = function(pars=c(parinint)) sum(leftside - cfmlam %*% pars)^2,
                          method="L-BFGS-B",
                          lower=rep(low,ngroup1),upper=rep(up,ngroup1))
          tausqest0 <- minopt$par
        }
        
        
        if(method=="exact"){
          #leftsidetemp <- leftsidefark/lengfark;leftsidetemp <- leftside
          soltau = solve(sum(coefmatfast),sum(leftside))
          sol = solve(coefmatfast,leftside)
          #sol=mygeninv(coefmatfast/leng) %*% (leftside/leng)
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
            #i <-3
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
          #tausqest <- solhybrec
          tausqest <- solhyb
        } 
        
        if(method=="adaptridge") {
          print("adaptive ridge")
          #tausqest <- solhybrec
        } 
        
        if(method=="stable" | method=="exact" | method=="exactstable" | method=="simple"){
          
          lambdanoncal <- 1/tausqest
          
          if(monotone[cl]){
            weigh = unlist(lapply(whgr,length))
            lambdamultnoncal <- pava(lambdanoncal,w=weigh)
          } else lambdamultnoncal <- lambdanoncal
          
          #lambdamultnoncal<-1/rev(cumsum(rev(log(1/lambdanoncal))))/(ngroup1:1)
          
          #logtausqest<-log(1/lambdamultnoncal)
          #con3 <- sum(sapply(1:length(whgr),function(gr){return(length(whgr[[gr]])*logtausqest[gr])}))/nfeat
          #loglambdamult<-  - logtausqest + con3
          #lambdamult <- exp(loglambdamult)
          
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
          } else {   #overlap
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
        
        opt2w <- cvl(responsemin, penalized = XMw0,fold=nf,lambda2=optl,unpenalized=nopen,data=datapred, trace=trace)
        #opt2w <- cvl(responsemin, penalized = XMw,fold=nf,lambda2=optl, unpenalized=nopen,data=datapred,...)
        #pred <- opt2w$predictions
        cvln1 <- opt2w$cvl
        print(cvln1)
        
        #### START NEW2 ###
        if((cvln1 - cvlnprev)/abs(cvlnprev) > 1/100  |  ((cvln1 - cvlnprev)/abs(cvlnprev) >= 0 & i==1)){ 
          #if((cvln1 - cvlnprev)/abs(cvlnprev) >= cvlmarg) { 
          pen <- penalized(responsemin, penalized = XMw0, lambda2 = optl,unpenalized=nopen,data=datapred)
          if(niter>1){
            if(!overl){
              for(group in 1:ngroup1){
                ind <- whgr[[group]]
                lmvec[ind] <- lmvec[ind]*lambdamult[group]
              }} else {  #overlap
                lmvec <- lmvec * lambdamultperk
              }
          }
          lambdas[[cl]] <- lambdas[[cl]]*lambdamult
          cvlnprev <- cvln1
          penprev <- pen
          XMw0prev <- XMw0
          print(paste("Partition",nmp0[cl],"improved results"))
        } else {
        #### END NEW2 ###  
          if(niter>1) print(paste("Partition",nmp0[cl],"CONVERGED after",i,"iterations")) 
          else print(paste("Partition",nmp0[cl],"did not improve results"))
          
          conv[cl] <- TRUE
          cvln1 <- cvlnprev
          XMw0 <- XMw0prev
          pen <- penprev
        }
         } #end if(ellarg <=.5)
      } #end if(!convcl)
      cl <- cl+1
    }
    if(sum(conv)==nclass){
      converged <- TRUE
      if(niter>1) print(paste("All partitions CONVERGED after",i,"iterations"))
    }
    i <- i+1
  }


#if(reCV){
#print("Finding NEW lambda")
#if(fixedfoldsinn) set.seed(346477)
#optnew <- optL2(response, penalized = XMw0,fold=foldinit,trace=T,unpenalized=nopen,data=datapred)
##opt <- optL2(response, penalized = t(highdimdata),fold=foldinit,trace=T,unpenalized=nopen,data=datapred,...)
#time1 <- proc.time()-pmt0
#print(optnew$cv)
#print(paste("Computation time for re-cross-validating main penalty parameter:",time1[3]))
#optlnew <- optnew$lambda
#pen <- penalized(responsemin, penalized = XMw0, lambda2 = optlnew, unpenalized=nopen,data=datapred)
#}

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


#load("C:/VUData/Wina/MethInfFull/workspaceVerlaat.Rdata") 

if(selection){
print("Start posthoc variable selection")
pmt <- proc.time()
ord <- order(abs(newbeta),decreasing=T)
sequ <- seq(0,maxsel,by=stepsel)
cvlsels <- c()
for(nsel in sequ){
whsel <- ord[1:nsel]
datwsel <- XMw0[,whsel]


pensel <- penalized(responsemin,datwsel,lambda2=optl, unpenalized=nopen,data=datapred,trace=FALSE)
#pensel <- penalized(responsemin,datwsel,lambda2=optl, unpenalized=nopen,data=datapred,...)
optsel <- cvl(responsemin,datwsel,fold=nf,lambda2=optl,unpenalized=nopen,data=datapred, trace=trace)
#optsel <- cvl(responsemin,datwsel,fold=nf,lambda2=optl, unpenalized=nopen,data=datapred,...)
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
#pensel2 <- penalized(responsemin,datwsel2,lambda2=optl, unpenalized=nopen,data=datapred,...)
#optsel2 <- cvl(responsemin,datwsel2,fold=nf,lambda2=optl, unpenalized=nopen,data=datapred,...)
cvlsel2 <- optsel2$cvl
#predsel <- predict(pensel2,penalized=XMw0[,whsel2,drop=FALSE], unpenalized=nopen)
#predsel <- predict(pensel2,penalized=XMw0[,whsel2,drop=FALSE], unpenalized=nopen)
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
if(selection) cvlssel <- data.frame(nsel=sequ,cvl=cvlsels) else cvlssel <- NULL

cvlnstot <- cvlnssam




#mat <- cbind(response0,allpreds) 
#if(printp){
#print("True:")
#print(response0)
#print(round(allpreds,3))
#}

#NEW
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
    cvliklasso <- if(cvllasso) try(cvl(response,penalized = t(highdimdata),lambda1 = optllasso,fold=nf,unpenalized=nopen,data=datapred, trace=trace))
    if(class(cvliklasso) == "try-error" | !cvllasso) cvllasso <- NA else cvllasso <- cvliklasso$cvl  #optsel <- cvl(responsemin,datwsel,fold=nf,lambda2=optl, unpenalized=nopen,data=datapred,...)
    }
cvlnstot <- c(cvlnstot,cvllasso)
penlasso <- penalized(response, penalized = t(highdimdata), lambda1 = optllasso, unpenalized=nopen,data=cbind(XM0,datapred))
whichlasso <- which(penlasso@penalized != 0)
betaslasso <- penlasso@penalized[whichlasso]
predobj <- c(predobj,penlasso)
reslasso <- list(cvllasso=cvllasso,whichlasso=whichlasso,betaslasso=betaslasso)
}

resEN <- NULL 
if(compareEN){
  
  if(selection) maxsel <- length(whsel2)
  
  fsel <- function(lam1,maxselec=maxsel,lam2){
    if(lam1==0) return(nfeat-maxselec) else {
      penselEN <- penalized(responsemin,XMw0,lambda1=lam1,lambda2=lam2, unpenalized=nopen,data=datapred,trace=FALSE,maxiter=100)
      coef <- penselEN@penalized
      return(length(coef[coef!=0])-maxselec)
    }
  }
  lam1 <- uniroot(fsel,interval=c(0,optl*10),maxiter=50,lam2=optl)$root
   penselEN0 <- penalized(responsemin,XMw0,lambda1=lam1,lambda2=optl, unpenalized=nopen,data=datapred,
                         trace=FALSE,maxiter=100)
   coefEN0 <- penselEN0@penalized
   whichEN <- which(coefEN0 != 0)
#   newl2 <- optL2(responsemin,XMw0,lambda1=lam1,fold=nf,unpenalized=nopen,data=datapred,trace=TRUE)
#   optlnew <- max(0.001,newl2$lambda)
#   print(optlnew)
#   lam1new <- uniroot(fsel,interval=c(0,optl/2),maxiter=50,lam2=optlnew)$root
#   print(lam1new)
  
   penselEN <- penalized(responsemin,XMw0[,whichEN,drop=FALSE],lambda2=optl, unpenalized=nopen,data=datapred,
                         trace=FALSE,maxiter=100)
#   penselEN <- penalized(responsemin,XMw0,lambda1=lam1new,lambda2=optlnew, unpenalized=nopen,data=datapred,
#                         trace=FALSE,maxiter=100)
  coefEN <- penselEN@penalized
#   coefEN <- penselEN@penalized
#   whichEN <- which(coefEN!=0)
  
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
#colnames(allpreds) <- nmp
#NEW
return(list(true=response,cvls = cvlnstot,lambdamults = lambdas, optl=optl, lambdamultvec = almvecall, 
            predobj=predobj,betas=newbeta, whichsel = whichsel,cvlssel = cvlssel,reslasso=reslasso, 
            resEN = resEN, model=model, arguments=arguments,allpreds=allpreds))
}
