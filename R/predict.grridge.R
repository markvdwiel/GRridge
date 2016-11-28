predict.grridge <- function (object, datanew, printpred = FALSE,
                             dataunpennew=NULL,responsetest=NULL,
                             recalibrate=FALSE,...) 
{
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
  offsarg <- arg$offset
  if (class(datanew) != "data.frame") 
    datanew <- data.frame(datanew)
  npred <- ncol(datanew)
  Xsam <- t(datanew)
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
  nmp <- names(penobj)
  if (arg$selectionForward || arg$compareEN || arg$selectionEN) {npreds <- length(penobj) - 1}else {npreds <- length(penobj)}
  if((arg$selectionForward)&&(arg$compareEN) || (arg$selectionForward)&&(arg$selectionEN)){npreds <- length(penobj) - 2}
  if (arg$comparelasso)
    npreds <- npreds - 1
  if (arg$compareunpenal) 
    npreds <- npreds - 1
  if (npreds > 0) {
    for (ell in 1:npreds) {
      predobj <- penobj[[ell]]
      lmvecall <- grr$lambdamultvec[, ell]
      Xsamw <- t(t(Xsam)/sqrt(lmvecall))
      predell <- predict(predobj, Xsamw, data = dataunpensam)[1:npred]
      predellall <- cbind(predellall, predell)
      if (recalibrate) {
        datlp <- Xsamw %*% predobj@penalized
        refitmod <- glm(responsetest ~ 1 + datlp, family = "binomial")
        predell2 <- predict(refitmod, type = "response")
        predellall2 <- cbind(predellall2, predell2)
      }
    }
  }
  
  if (arg$selectionForward) {
    predobj <- penobj[[npreds+1]]
    whsel <- grr$whichsel
    nc <- ncol(grr$lambdamultvec)
    lmvecall <- grr$lambdamultvec[, nc]
    Xsamw <- t(t(Xsam)/sqrt(lmvecall))
    Xsamw <- Xsamw[, whsel, drop = FALSE]
    predell <- predict(predobj, Xsamw, data = dataunpensam)[1:npred]
    predellall <- cbind(predellall, predell)
    if (recalibrate) {
      datlp <- Xsamw %*% predobj@penalized
      refitmod <- glm(responsetest ~ 1 + datlp, family = "binomial")
      predell2 <- predict(refitmod, type = "response")
      predellall2 <- cbind(predellall2, predell2)
    }
  }
  
  if (arg$compareEN || arg$selectionEN) {
    add = arg$selectionForward + arg$comparelasso
    take = npreds + add + 1
    predobj <- penobj[[take]]
    whsel <- grr$resEN$whichEN
    nc <- ncol(grr$lambdamultvec)
    lmvecall <- grr$lambdamultvec[, nc]
    Xsamw <- t(t(Xsam)/sqrt(lmvecall))
    Xsamw <- Xsamw[, whsel, drop = FALSE]
    predell <- predict(predobj, Xsamw, data = dataunpensam)[1:npred]
    predellall <- cbind(predellall, predell)
    if (recalibrate) {
      datlp <- Xsamw %*% predobj@penalized
      refitmod <- glm(responsetest ~ 1 + datlp, family = "binomial")
      predell2 <- predict(refitmod, type = "response")
      predellall2 <- cbind(predellall2, predell2)
    }
  }
  
  if (arg$comparelasso) {
    
    if (arg$selectionForward)
      take <- npreds + 2
    else take <- npreds
    
    print(paste("takelasso",take))
    predobj <- penobj[[take]]
    predell <- predict(predobj, Xsam, unpenalized = unpenal, 
                       data = dataunpensam)[1:npred]
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
    nadd <- arg$selectionForward + arg$compareEN + arg$comparelasso
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
