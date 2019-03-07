predict.grridge <- function (object, datanew, printpred = FALSE,
                             dataunpennew=NULL,responsetest=NULL,
                             recalibrate=FALSE,...) 
{
#object=grMaarten1part;datanew=mirnormcen_resp;dataunpennew=datfr;responsetest=NULL;recalibrate=FALSE
  grr <- object
  model <- grr$model
  
  if (!(is.null(grr$arguments$dataunpen)) & recalibrate == 
        TRUE) {
    recalibrate <- FALSE
    print("Recalibration currently only feasible for designs without unpenalized covariates. Recalibrate set to FALSE")
  }
  if (grr$model == "survival" & recalibrate == TRUE) {
    recalibrate <- FALSE
    print("Recalibration currently only feasible for linear and binary response. Recalibrate set to FALSE")
  }
  ntest <- ncol(datanew)
  if (ntest < 25 & recalibrate == TRUE) {
    recalibrate <- FALSE
    print("Test sample size too small for recalibration. Need at least 25 test samples. Recalibrate set to FALSE")
  }
  penobj <- grr$predobj
  arg <- grr$arguments
  
  if(arg$standardizeX) {
    cat("Model has been fit on standardized features. Make sure datanew contains standardized features. 
        Use either the training data or datanew (if sufficiently large, say N>20) to standardize. See 
        ?predict.grridge.")
  }
  
  offsarg <- arg$offset
  if (class(datanew) != "data.frame") 
    datanew <- data.frame(datanew)
  npred <- ncol(datanew)
  Xsam <- t(datanew)
  
  optllasso <- grr$arguments$optllasso
  if (length(offsarg) > 1) offsout <- offsarg[samout] else offsout <- rep(0,npred)
  
  if (is.null(dataunpennew)) {
    dataunpensam <- data.frame(fake = rep(NA, npred))} else{
      if (class(dataunpennew) != "data.frame") {dataunpensam <- data.frame(dataunpennew)}else {
        dataunpensam <- dataunpennew}
    }
  predellall <- c()
  predellall2 <- c()
  if (is.null(penobj)) {
    cat("No prediction objection available. Run grridge using either savepredobj=\"last\" or savepredobj=\"all\"\n")
    return(NULL)
  }
    
  unpenal <- arg$unpenal
  if (!is.null(offsarg)) {
    noffs <- length(offsarg)
    offsargs <- "c("
    if (npred > 1) {
      if (noffs == 1) {
        for (i in 1:(npred - 1)) offsargs <- paste(offsargs, 
                                                   offsarg, ",", sep = "")
      }
      else {
        for (i in 1:(npred - 1)) offsargs <- paste(offsargs, 
                                                   offsarg[i], ",", sep = "")
      }
    }
    if (noffs == 1) 
      offsargs <- paste(offsargs, offsarg, ")", sep = "")
    else offsargs <- paste(offsargs, offsarg[npred], ")", 
                           sep = "")
    if ((unpenal != ~0) & (unpenal != ~1)) {
      unpenal <- formula(paste(deparse(unpenal), "+ offset(", 
                               offsargs, ")", sep = ""))
    }
    else {
      unpenal <- formula(paste("~", "offset(", offsargs, 
                               ")", sep = ""))
    }
  }
  
  #preparing for glmnet (lasso)
  if((is.null(dataunpennew)) | (unpenal == ~0) | (unpenal == ~1)){
    mmout <- NULL
  } else  {
    mmout <- model.matrix(unpenal,dataunpennew)
    if(prod(mmout[,1]==rep(1,npred))==1) {
      mmout <- mmout[,-1,drop=FALSE]
    }
  }
  
  nmp <- names(penobj)
  
  if (arg$selectionEN) {npreds <- length(penobj) - length(arg$maxsel)} else {npreds <- length(penobj)}
   if (arg$comparelasso)
    npreds <- npreds - 1
  if (arg$compareunpenal) 
    npreds <- npreds - 1
  if (npreds > 0) {
    for (ell in 1:npreds) {
      predobj <- penobj[[ell]]
      lmvecall <- grr$lambdamultvec[, ell]
      Xsamw <- t(t(Xsam)/sqrt(lmvecall))
      #predell <- predict(predobj, Xsamw, data = dataunpensam)[1:npred]
      #UPDATE 7-3-2019
      if(model != "survival")  predell <- predict(predobj, Xsamw, unpenalized = unpenal,  data = dataunpensam)[1:npred]
      else { #survival (PENALIZED DOES NOT SUPPLY linear predictor...
        mmunpen <- as.matrix(model.matrix(unpenal, dataunpensam))
        wh <- match(names(predobj@unpenalized),colnames(mmunpen))
        coeffpen <- matrix(predobj@penalized, ncol = 1)
        coeffunpen <- matrix(predobj@unpenalized, ncol = 1)
        if (length(coeffunpen) == 0) predell <- exp(Xsamw %*% coeffpen) else predell <- exp(Xsamw %*%  coeffpen + 
                                                                                              mmunpen[,wh,drop=FALSE] %*% coeffunpen)
      }
      predellall <- cbind(predellall, predell)
      if (recalibrate) {
        datlp <- Xsamw %*% predobj@penalized
        refitmod <- glm(responsetest ~ 1 + datlp, family = "binomial")
        predell2 <- predict(refitmod, type = "response")
        predellall2 <- cbind(predellall2, predell2)
      }
    }
  }
  
  
  if (arg$selectionEN) {
    nadd =  arg$comparelasso
    toadd <- npreds + nadd
    for (ell in 1:length(arg$maxsel)) {
    predobj <- penobj[[toadd+ell]]
    #whsel <- grr$resEN$whichEN
    whsel <- grr$resEN[[ell]]$whichEN
    nc <- ncol(grr$lambdamultvec)
    lmvecall <- grr$lambdamultvec[, nc]
    Xsamw <- t(t(Xsam)/sqrt(lmvecall))
    Xsamw <- Xsamw[, whsel, drop = FALSE]
    #predell <- predict(predobj, Xsamw, data = dataunpensam)[1:npred]
    #UPDATE 7-3-2019
    if(model != "survival")  predell <- predict(predobj, Xsamw, unpenalized = unpenal,  data = dataunpensam)[1:npred]
    else { #survival (PENALIZED DOES NOT SUPPLY linear predictor...
      mmunpen <- as.matrix(model.matrix(unpenal, dataunpensam))
      wh <- match(names(predobj@unpenalized),colnames(mmunpen))
      coeffpen <- matrix(predobj@penalized, ncol = 1)
      coeffunpen <- matrix(predobj@unpenalized, ncol = 1)
      if (length(coeffunpen) == 0) predell <- exp(Xsamw %*% coeffpen) else predell <- exp(Xsamw %*%  coeffpen +                                                                                   mmunpen[,wh,drop=FALSE] %*% coeffunpen)
    }
    predellall <- cbind(predellall, predell)
    if (recalibrate) {
      datlp <- Xsamw %*% predobj@penalized
      refitmod <- glm(responsetest ~ 1 + datlp, family = "binomial")
      predell2 <- predict(refitmod, type = "response")
      predellall2 <- cbind(predellall2, predell2)
    }
  } #end for ell
  }
  
  if (arg$comparelasso) {
    
    take <- npreds+1
    
    print(paste("takelasso",take))
    predobj <- penobj[[take]]
    
    #UPDATE 7-3-2019
    predell <- try(as.numeric(predict(predobj,cbind(Xsam,mmout),s=c(optllasso),offset=offsout, type="response")), silent=T)
    if(class(predell) == "try-error") predell <- as.numeric(predict(predobj,cbind(Xsam,mmout),s=c(optllasso),newoffset=offsout,type="response"))
    
    
    predellall <- cbind(predellall, predell)
    if (recalibrate) {
      datlp <- Xsam %*% predobj@penalized
      if (model == "logistic") 
        refitmod <- glm(responsetest ~ 1 + datlp, family = "binomial")
      else refitmod <- glm(responsetest ~ 1 + datlp, family = "gaussian")
      predell2 <- predict(refitmod, type = "response")
      predellall2 <- cbind(predellall2, predell2)
    }
  }
  if (arg$compareunpenal) {
    nadd <- arg$selectionEN + arg$comparelasso
    take <- npreds + nadd + 1
    predobj <- penobj[[take]]
    predell <- predict(predobj, data = dataunpensam, type = "response")[1:npred]
    predellall <- cbind(predellall, predell)
  }
  colnames(predellall) <- nmp
  if (recalibrate) {
    colnames(predellall2) <- paste(nmp, "_recalib", sep = "")
    predellall <- cbind(predellall, predellall2)
  }
  
  return(predellall)
}
