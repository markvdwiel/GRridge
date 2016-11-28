grridgeCV <- function(grr,highdimdata,response,
                       outerfold=length(response),
                       fixedfolds=TRUE,recalibrate=FALSE){
model <- grr$model

if(model=="survival") allobstimes <- sort(response[,1])
balance <- FALSE
nsam <- ncol(highdimdata)
if(is.null(colnames(highdimdata))){
print("No sample names available. Creating names S1, ..., Sn.")
cnames <- sapply(1:nsam,function(i) paste("S",i,sep=""))
colnames(highdimdata) <- cnames
}
if(outerfold < 10 & model != "linear") {
  print("Stratified splits for events applied")
  balance <- TRUE
}
if(fixedfolds) set.seed(3534) else set.seed(NULL)
if(!balance){
rand<-sample(1:nsam)
grs1 <- floor(nsam/outerfold)
grs2 <- grs1+1
ngr1 <- outerfold*grs2 - nsam
folds <- lapply(1:outerfold,function(xg) {
if(xg <= ngr1) els <- rand[(1+(xg-1)*grs1):(xg*grs1)] else els <- rand[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
return(els)
}
)} else {
if(model=="logistic") if(class(response)=="factor") nev <- which((as.numeric(response)-1)==1) else nev <- which(response==1)  
if(model=="survival") nev <- which(response[,1]==1)    
nsamev <- length(nev) 
randev<-sample(nev)
grs1 <- floor(nsamev/outerfold)
grs2 <- grs1+1
ngr1 <- outerfold*grs2 - nsamev
foldsev <- lapply(1:outerfold,function(xg) {
  if(xg <= ngr1) els <- randev[(1+(xg-1)*grs1):(xg*grs1)] else els <- randev[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
  return(els)
}
)
nonev <- setdiff(1:nsam,nev)
nsamnonev <- length(nonev) 
randnonev<-sample(nonev)
grs1 <- floor(nsamnonev/outerfold)
grs2 <- grs1+1
ngr1 <- outerfold*grs2 - nsamnonev
foldsnonev <- lapply(1:outerfold,function(xg) {
  if(xg <= ngr1) els <- randnonev[(1+(xg-1)*grs1):(xg*grs1)] else els <- randnonev[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
  return(els)
}
)
folds <- lapply(1:outerfold,function(i) c(foldsev[[i]],foldsnonev[[i]]))
}


if(!(is.null(grr$arguments$dataunpen)) & recalibrate==TRUE) {
  recalibrate <- FALSE
  print("Recalibration currently only feasible for designs without unpenalized covariates. Recalibrate set to FALSE")
}

if(model=="survival" & recalibrate==TRUE) {
  recalibrate <- FALSE
  print("Recalibration currently only feasible for linear and binary response. Recalibrate set to FALSE")
}

ntest <- length(folds[[1]])
if(ntest < 25 & recalibrate==TRUE) {
  recalibrate <- FALSE
  print("Test sample size too small for recalibration. Need at least 25 test samples. Recalibrate set to FALSE")
}

whichfold <- rep(NA,nsam)
linpreds <- c()
preds <- c()
responsesam <- c()
predsbres <- c()
for(k in 1:length(folds)){
#k<-1
print(paste("Fold nr:",k))
samout <- folds[[k]]
whichfold[samout] <- k
nout <- length(samout)
highdimdatamin <- highdimdata[,-samout]
responsemin <- response[-samout]
arg <- grr$arguments
dataunpenmin <- arg$dataunpen[-samout,,drop=FALSE] 
offsarg  <- arg$offset
if(length(offsarg)>1) offsmin <- offsarg[-samout] else offsmin <- offsarg
#NEW
grmin <- grridge(highdimdatamin,responsemin,partitions=arg$partitions,unpenal=arg$unpenal, offset=offsmin, 
method=arg$method,niter=arg$niter, monotone=arg$monotone, optl=arg$optl, innfold=arg$innfold,
fixedfoldsinn=arg$fixedfoldsinn, selectionForward=arg$selectionForward, maxsel=arg$maxsel,
stepsel=arg$stepsel,cvlmarg=arg$cvlmarg, dataunpen=dataunpenmin,savepredobj=arg$savepredobj,ord=arg$ord,
comparelasso = arg$comparelasso, optllasso=arg$optllasso, selectionEN = arg$selectionEN, 
compareunpenal = arg$compareunpenal,modus=arg$modus)
penobj <- grmin$predobj
dataunpen <- arg$dataunpen
Xsam <- t(highdimdata[,samout,drop=FALSE])
dataunpensam <- dataunpen[samout,,drop=FALSE]
responsesam <- c(responsesam,samout)
responseout <- response[samout]
if(is.null(dataunpensam)) dataunpensam <- data.frame(fake=rep(NA,length(samout)))
predellall <- c()
predellall2 <- c()
linpredall <- c()
if(is.null(penobj)) {
cat("No prediction objection available. Run grridge using either savepredobj=\"last\" or savepredobj=\"all\"\n")
return(NULL)
}

unpenal<-arg$unpenal
#offsarg <- log(1/9);unpenal= ~0;nout=1
if(!is.null(offsarg)){
noffs <- length(offsarg)
offsargs <-"c("
if(nout>1){
    if(noffs==1) {for(i in 1:(nout-1)) offsargs <- paste(offsargs,offsarg,",",sep="")} else {
    for(i in 1:(nout-1)) offsargs <- paste(offsargs,offsarg[i],",",sep="")}
}
if(noffs==1) offsargs <- paste(offsargs,offsarg,")",sep="")  else offsargs <- paste(offsargs,offsarg[nout],")",sep="")
if((unpenal != ~0) & (unpenal != ~1)){
unpenal <- formula(paste(deparse(unpenal),"+ offset(",offsargs,")",sep=""))
} else {
unpenal <- formula(paste("~","offset(",offsargs,")",sep=""))
}
}

nmp <- names(penobj)

if(model=="logistic" | model =="linear"){
    #NEW
    if(arg$selectionForward) npreds <- length(penobj) -1 else npreds <- length(penobj)
    if(arg$comparelasso) npreds <- npreds-1
    if(arg$compareunpenal) npreds <- npreds-1
    if(arg$selectionEN) npreds <- npreds-1
    if(npreds>0){
        for(ell in 1:npreds){
        #ell<-1
        predobj <- penobj[[ell]]
        lmvecall <-  grmin$lambdamultvec[,ell]
        Xsamw <- t(t(Xsam)/sqrt(lmvecall))
        predell <- predict(predobj,Xsamw,unpenalized=unpenal,data=dataunpensam)[1:nout]
        predellall <- cbind(predellall,predell)
        
        if(recalibrate){
          datlp <- Xsamw %*% predobj@penalized
          if(model=="logistic") refitmod <- glm(responseout ~  1 + datlp, family="binomial") else
            refitmod <- glm(responseout ~  1 + datlp, family="gaussian")
          print(paste("Recalibration intercept and slope for (gr)-ridge model",ell,":",round(refitmod$coefficients[1],3),";",
                      round(refitmod$coefficients[2],3)))
          predell2 <- predict(refitmod,type="response")
          predellall2 <- cbind(predellall2, predell2)
        }
        }
    }
  if(arg$selectionForward){
      predobj <- penobj[[npreds+1]]
      whsel <- grmin$whichsel
      nc <- ncol(grmin$lambdamultvec)
      lmvecall <-  grmin$lambdamultvec[,nc]
      Xsamw <- t(t(Xsam)/sqrt(lmvecall))
      Xsamw <- Xsamw[,whsel,drop=FALSE]
      predell <- predict(predobj,Xsamw,unpenalized=unpenal,data=dataunpensam)[1:nout]
      predellall <- cbind(predellall,predell)
      
      if(recalibrate){
        datlp <- Xsamw %*% predobj@penalized
        if(model=="logistic") refitmod <- glm(responseout ~  1 + datlp, family="binomial") else
          refitmod <- glm(responseout ~  1 + datlp, family="gaussian")
        print(paste("Recalibration intercept and slope for model selected variables:",
                    round(refitmod$coefficients[1],3),";",
                    round(refitmod$coefficients[2],3)))
        predell2 <- predict(refitmod,type="response")
        predellall2 <- cbind(predellall2, predell2)
      }
  }
  if(arg$comparelasso){ 
      if(arg$selectionForward) take <- npreds+2 else take <- npreds + 1
      predobj <- penobj[[take]]
      predell <- predict(predobj,Xsam,unpenalized=unpenal,data=dataunpensam)[1:nout]
      predellall <- cbind(predellall,predell)
      
      if(recalibrate){
        datlp <- Xsam %*% predobj@penalized
        if(model=="logistic") refitmod <- glm(responseout ~  1 + datlp, family="binomial") else
          refitmod <- glm(responseout ~  1 + datlp, family="gaussian")
        print(paste("Recalibration intercept and slope for lasso model:",round(refitmod$coefficients[1],3),";",
                    round(refitmod$coefficients[2],3)))
        predell2 <- predict(refitmod,type="response")
        predellall2 <- cbind(predellall2, predell2)
      }
  }
    
    if(arg$selectionEN){
      nadd <- arg$selectionForward + arg$comparelasso
      take <- npreds+nadd + 1
      predobj <- penobj[[take]]
      nc <- ncol(grmin$lambdamultvec)
      lmvecall <-  grmin$lambdamultvec[,nc]
      Xsamw <- t(t(Xsam)/sqrt(lmvecall))
      whEN <- grmin$resEN$whichEN
      Xsamw <- Xsamw[,whEN,drop=FALSE]
      predell <- predict(predobj,Xsamw,unpenalized=unpenal,data=dataunpensam)[1:nout]
      predellall <- cbind(predellall,predell)
      
      if(recalibrate){
        datlp <- Xsamw %*% predobj@penalized
        if(model=="logistic") refitmod <- glm(responseout ~  1 + datlp, family="binomial") else
          refitmod <- glm(responseout ~  1 + datlp, family="gaussian")
        print(paste("Recalibration intercept and slope for EN model",ell,":",round(refitmod$coefficients[1],3),";",
                    round(refitmod$coefficients[2],3)))
        predell2 <- predict(refitmod,type="response")
        predellall2 <- cbind(predellall2, predell2)
      }
    }
  
  if(arg$compareunpenal){
      nadd <- arg$selectionForward + arg$comparelasso + arg$selectionEN
      take <- npreds+nadd + 1
      predobj <- penobj[[take]]
      predell <- predict(predobj,data=dataunpensam,type="response")[1:nout]
      predellall <- cbind(predellall,predell)
  }
  
      print(paste("Sample(s) left out:",samout))
      print("True:")
      print(response[samout])
      print("Prediction(s):")
      colnames(predellall)<-nmp
      if(recalibrate) {
        colnames(predellall2) <- paste(nmp,"_recalib",sep="") 
        predellall <- cbind(predellall,predellall2)
      }
      print(data.frame(response=response[samout],predellall))
  preds <- rbind(preds,predellall)
  
} else {  #model="survival"
  
  mmunpen <- as.matrix(model.matrix(unpenal,dataunpensam)[,-1])
  
  if(arg$selectionForward) npreds <- length(penobj) -1 else npreds <- length(penobj)
  if(arg$comparelasso) npreds <- npreds-1
  if(arg$compareunpenal) npreds <- npreds-1
  if(arg$selectionEN) npreds <- npreds-1
  if(npreds>0){
    for(ell in 1:npreds){
      predobj <- penobj[[ell]]
      lmvecall <-  grmin$lambdamultvec[,ell]
      Xsamw <- t(t(Xsam)/sqrt(lmvecall))
      predell0 <- predict(predobj,Xsamw,unpenalized=unpenal,data=dataunpensam)
      predell <- sapply(allobstimes,function(tim) return(survival(predell0,time=tim)))
      predellall <- cbind(predellall,predell)
      coeffpen <- matrix(predobj@penalized,ncol=1)
      coeffunpen <- matrix(predobj@unpenalized,ncol=1)
      if(length(coeffunpen) == 0) linpredall <- cbind(linpredall,Xsamw %*% coeffpen) else 
        linpredall <- cbind(linpredall,Xsamw %*% coeffpen + mmunpen %*% coeffunpen) 
    }
  }
  if(arg$selectionForward){
    predobj <- penobj[[npreds+1]]
    whsel <- grmin$whichsel
    nc <- ncol(grmin$lambdamultvec)
    lmvecall <-  grmin$lambdamultvec[,nc]
    Xsamw <- t(t(Xsam)/sqrt(lmvecall))
    Xsamw <- Xsamw[,whsel,drop=FALSE]
    predell0 <- predict(predobj,Xsamw,unpenalized=unpenal,data=dataunpensam)
    predell <- sapply(allobstimes,function(tim) return(survival(predell0,time=tim)))  
    predellall <- cbind(predellall,predell) 
    coeffpen <- matrix(predobj@penalized,ncol=1)
    coeffunpen <- matrix(predobj@unpenalized,ncol=1)
    if(length(coeffunpen) == 0) linpredall <- cbind(linpredall,Xsamw %*% coeffpen) else 
      linpredall <- cbind(linpredall,Xsamw %*% coeffpen +  mmunpen %*% coeffunpen)     
  }
  
  if(arg$comparelasso){
    if(arg$selectionForward) take <- npreds+2 else take <- npreds + 1
    predobj <- penobj[[take]]
    predell0 <- predict(predobj,Xsam,unpenalized=unpenal,data=dataunpensam)
    predell <- sapply(allobstimes,function(tim) return(survival(predell0,time=tim)))
    predellall <- cbind(predellall,predell)
    coeffpen <- matrix(predobj@penalized,ncol=1)
    coeffunpen <- matrix(predobj@unpenalized,ncol=1)
    if(length(coeffunpen) == 0) linpredall <- cbind(linpredall,Xsam %*% coeffpen) else 
      linpredall <- cbind(linpredall,Xsam %*% coeffpen +  mmunpen %*% coeffunpen) 
  }
  
  if(arg$selectionEN){ #EN
    nadd <- arg$selectionForward + arg$comparelasso
    take <- npreds+nadd + 1
    predobj <- penobj[[take]]
    nc <- ncol(grmin$lambdamultvec)
    lmvecall <-  grmin$lambdamultvec[,nc]
    Xsamw <- t(t(Xsam)/sqrt(lmvecall))
    whEN <- grmin$resEN$whichEN
    Xsamw <- Xsamw[,whEN,drop=FALSE]
    predell0 <- predict(predobj,Xsamw,unpenalized=unpenal,data=dataunpensam)
    predell <- sapply(allobstimes,function(tim) return(survival(predell0,time=tim)))  
    predellall <- cbind(predellall,predell) 
    coeffpen <- matrix(predobj@penalized,ncol=1)
    coeffunpen <- matrix(predobj@unpenalized,ncol=1)
    if(length(coeffunpen) == 0) linpredall <- cbind(linpredall,Xsamw %*% coeffpen) else 
      linpredall <- cbind(linpredall,Xsamw %*% coeffpen +  mmunpen %*% coeffunpen)    
  }
  
  if(arg$compareunpenal){
    nadd <- arg$selectionForward + arg$comparelasso
    take <- npreds+nadd + 1
    predobj <- penobj[[take]]
    bogus <- rep(0,nout)
    predell0 <- predict(predobj,penalized=~bogus,unpenalized=unpenal,data=dataunpensam)
    predell <- sapply(allobstimes,function(tim) return(survival(predell0,time=tim)))
    predellall <- cbind(predellall,predell)
    coeffunpen <- matrix(predobj@unpenalized,ncol=1) 
    linpredall <- cbind(linpredall,mmunpen %*% coeffunpen) 
  }

  print(paste("Sample(s) left out:",samout))
  print("True:")
  print(response[samout,])
  print("True and Prediction(s):")
  names(predellall)<-nmp
  print(cbind(response[samout,1],response[samout,2],linpredall))
  linpreds <- rbind(linpreds,linpredall)
  predsbres <- rbind(predsbres,predellall)
}
}

if(model=="logistic" | model =="linear"){
sams <- unlist(folds)
od <- order(sams)
mat <- data.frame(response[responsesam],preds)[od,]
mat <- data.frame(mat,whichfold)
colnames(mat) <- c("TrueResponse",colnames(preds),"whichfold")
return(mat)
} else { 
   
  sams <- unlist(folds)
  od <- order(sams)
  mat <- data.frame(response[responsesam,1],response[responsesam,2],linpreds)[od,]
  mat <- data.frame(mat,whichfold)
  
  predsbresord <- predsbres[od,]
  predsbresordlist <- list()
  npred <- ncol(mat)-2
  nc <- ncol(predsbresord)
  ncmat <- nc/npred
  
  for(j in 1:npred)  predsbresordlist <- c(predsbresordlist,list(predsbresord[,((j-1)*ncmat+1):(j*ncmat)])) 
  names(predsbresordlist) <- nmp
  
   colnames(mat) <- c("TrueSurv","Status",nmp,"whichfold")
  return(list(linpredmat = mat,breslowmats=predsbresordlist,obstimes=allobstimes)) 
}
}
