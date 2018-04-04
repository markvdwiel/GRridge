grridgeCV <- function (grr, highdimdata, response, outerfold = length(response), 
          fixedfolds = TRUE, recalibrate = FALSE) {
  model <- grr$model
  if(model=="linear") return(.grridgeCVlin(grr=grr, highdimdata=highdimdata, response=response,  outerfold = outerfold, fixedfolds = fixedfolds, recalibrate = recalibrate))
  
  if(grr$arg$standardizeX=TRUE){
    print("Covariates are standardized")
    sds <- apply(highdimdata,1,sd)
    sds2 <- sapply(sds,function(x) max(x,10^{-5}))
    highdimdata <- (highdimdata-apply(highdimdata,1,mean))/sds2 
  }
  
  if (model == "survival") 
    allobstimes <- sort(response[, 1])
  balance <- FALSE
  nsam <- ncol(highdimdata)
  if (is.null(colnames(highdimdata))) {
    print("No sample names available. Creating names S1, ..., Sn.")
    cnames <- sapply(1:nsam, function(i) paste("S", i, sep = ""))
    colnames(highdimdata) <- cnames
  }
  if (outerfold < 10 & model != "linear") {
    print("Stratified splits for events applied")
    balance <- TRUE
  }
  if (fixedfolds) 
    set.seed(3534)
  else set.seed(NULL)
  if (!balance) {
    rand <- sample(1:nsam)
    grs1 <- floor(nsam/outerfold)
    grs2 <- grs1 + 1
    ngr1 <- outerfold * grs2 - nsam
    folds <- lapply(1:outerfold, function(xg) {
      if (xg <= ngr1) 
        els <- rand[(1 + (xg - 1) * grs1):(xg * grs1)]
      else els <- rand[(ngr1 * grs1 + 1 + (xg - ngr1 - 
                                             1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
      return(els)
    })
  }
  else {
    if (model == "logistic") 
      if (class(response) == "factor") 
        nev <- which((as.numeric(response) - 1) == 1)
      else nev <- which(response == 1)
      if (model == "survival") 
        nev <- which(response[, 1] == 1)
      nsamev <- length(nev)
      randev <- sample(nev)
      grs1 <- floor(nsamev/outerfold)
      grs2 <- grs1 + 1
      ngr1 <- outerfold * grs2 - nsamev
      foldsev <- lapply(1:outerfold, function(xg) {
        if (xg <= ngr1) 
          els <- randev[(1 + (xg - 1) * grs1):(xg * grs1)]
        else els <- randev[(ngr1 * grs1 + 1 + (xg - ngr1 - 
                                                 1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
        return(els)
      })
      nonev <- setdiff(1:nsam, nev)
      nsamnonev <- length(nonev)
      randnonev <- sample(nonev)
      grs1 <- floor(nsamnonev/outerfold)
      grs2 <- grs1 + 1
      ngr1 <- outerfold * grs2 - nsamnonev
      foldsnonev <- lapply(1:outerfold, function(xg) {
        if (xg <= ngr1) 
          els <- randnonev[(1 + (xg - 1) * grs1):(xg * 
                                                    grs1)]
        else els <- randnonev[(ngr1 * grs1 + 1 + (xg - ngr1 - 
                                                    1) * grs2):(ngr1 * grs1 + (xg - ngr1) * grs2)]
        return(els)
      })
      folds <- lapply(1:outerfold, function(i) c(foldsev[[i]], 
                                                 foldsnonev[[i]]))
  }
  if (!(is.null(grr$arguments$dataunpen)) & recalibrate == 
      TRUE) {
    recalibrate <- FALSE
    print("Recalibration currently only feasible for designs without unpenalized covariates. Recalibrate set to FALSE")
  }
  if (model == "survival" & recalibrate == TRUE) {
    recalibrate <- FALSE
    print("Recalibration currently only feasible for linear and binary response. Recalibrate set to FALSE")
  }
  ntest <- length(folds[[1]])
  if (ntest < 25 & recalibrate == TRUE) {
    recalibrate <- FALSE
    print("Test sample size too small for recalibration. Need at least 25 test samples. Recalibrate set to FALSE")
  }
  whichfold <- rep(NA, nsam)
  linpreds <- c()
  preds <- c()
  responsesam <- c()
  predsbres <- c()
  for (k in 1:length(folds)) {
    #k<-1
    print(paste("Fold nr:", k))
    samout <- folds[[k]]
    whichfold[samout] <- k
    nout <- length(samout)
    highdimdatamin <- highdimdata[, -samout]
    responsemin <- response[-samout]
    arg <- grr$arguments
    dataunpenmin <- arg$dataunpen[-samout, , drop = FALSE]
    offsarg <- arg$offset
    unpenal = arg$unpenal
    dataunpenout <- arg$dataunpen[samout, , drop = FALSE]
    if((is.null(dataunpenout)) | (unpenal == ~0) | (unpenal == ~1)) {
      mmout <- NULL
    } else  {
      mmout <- model.matrix(unpenal,dataunpenout)
      if(prod(mmout[,1]==rep(1,nout))==1) {
        mmout <- mmout[,-1,drop=FALSE]
      }
    }
    
    if (length(offsarg) > 1) offsmin <- offsarg[-samout] else offsmin <- offsarg
    #switch of standardization, as this is done before fitting
    grmin <- grridge(highdimdatamin, responsemin, partitions = arg$partitions, 
                     unpenal = arg$unpenal, offset = offsmin, method = arg$method, 
                     niter = arg$niter, monotone = arg$monotone, optl = arg$optl, 
                     innfold = arg$innfold, fixedfoldsinn = arg$fixedfoldsinn, 
                     maxsel = arg$maxsel,  cvlmarg = arg$cvlmarg, dataunpen = dataunpenmin, 
                     savepredobj = arg$savepredobj, ord = arg$ord, comparelasso = arg$comparelasso, 
                     optllasso = arg$optllasso, selectionEN = arg$selectionEN, 
                     compareunpenal = arg$compareunpenal, modus = arg$modus,EBlambda=arg$EBlambda, standardizeX = FALSE)
    penobj <- grmin$predobj
    dataunpen <- arg$dataunpen
    Xsam <- t(highdimdata[, samout, drop = FALSE])
    optllasso <- grmin$arguments$optllasso
    if (length(offsarg) > 1) offsout <- offsarg[samout] else offsout <- rep(0,nout)
    
    dataunpensam <- dataunpen[samout, , drop = FALSE]
    responsesam <- c(responsesam, samout)
    responseout <- response[samout]
    if (is.null(dataunpensam)) dataunpensam <- data.frame(fake = rep(NA, length(samout)))
    predellall <- c()
    predellall2 <- c()
    linpredall <- c()
    if (is.null(penobj)) {
      cat("No prediction objection available. Run grridge using either savepredobj=\"last\" or savepredobj=\"all\"\n")
      return(NULL)
    }
    unpenal <- arg$unpenal
    if (!is.null(offsarg)) {
      noffs <- length(offsarg)
      offsargs <- "c("
      if (nout > 1) {
        if (noffs == 1) {
          for (i in 1:(nout - 1)) offsargs <- paste(offsargs, 
                                                    offsarg, ",", sep = "")
        }
        else {
          for (i in 1:(nout - 1)) offsargs <- paste(offsargs, 
                                                    offsarg[i], ",", sep = "")
        }
      }
      if (noffs == 1) 
        offsargs <- paste(offsargs, offsarg, ")", sep = "")
      else offsargs <- paste(offsargs, offsarg[nout], ")", 
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
    
    npreds <- length(penobj)
    if (arg$comparelasso) 
      npreds <- npreds - 1
    if (arg$compareunpenal) 
      npreds <- npreds - 1
    if (arg$selectionEN) 
      npreds <- npreds - 1
    if (npreds > 0) {
      for (ell in 1:npreds) {
        predobj <- penobj[[ell]]
        lmvecall <- grmin$lambdamultvec[, ell]
        Xsamw <- t(t(Xsam)/sqrt(lmvecall))
        
        if(model != "survival")  predell <- predict(predobj, Xsamw, unpenalized = unpenal,  data = dataunpensam)[1:nout]
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
          if (model == "logistic") 
            refitmod <- glm(responseout ~ 1 + datlp, 
                            family = "binomial")
          else refitmod <- glm(responseout ~ 1 + datlp, 
                               family = "gaussian")
          print(paste("Recalibration intercept and slope for (gr)-ridge model", 
                      ell, ":", round(refitmod$coefficients[1], 
                                      3), ";", round(refitmod$coefficients[2], 
                                                     3)))
          predell2 <- predict(refitmod, type = "response")
          predellall2 <- cbind(predellall2, predell2)
        }
      }
    }
    
    if (arg$comparelasso) {
      take <- npreds + 1
      predobj <- penobj[[take]]
      
      predell <- as.numeric(predict(predobj,cbind(Xsam,mmout),s=c(optllasso),offset=offsout,type="response"))
      
      predellall <- cbind(predellall, predell)
      if (recalibrate) {
        datlp <- Xsam %*% predobj$beta
        if (model == "logistic") 
          refitmod <- glm(responseout ~ 1 + datlp, 
                          family = "binomial")
        else refitmod <- glm(responseout ~ 1 + datlp, 
                             family = "gaussian")
        print(paste("Recalibration intercept and slope for lasso model:", 
                    round(refitmod$coefficients[1], 3), ";", 
                    round(refitmod$coefficients[2], 3)))
        predell2 <- predict(refitmod, type = "response")
        predellall2 <- cbind(predellall2, predell2)
      }
    }
    if (arg$selectionEN) {
      nadd <- arg$comparelasso
      take <- npreds + nadd + 1
      predobj <- penobj[[take]]
      nc <- ncol(grmin$lambdamultvec)
      lmvecall <- grmin$lambdamultvec[, nc]
      Xsamw <- t(t(Xsam)/sqrt(lmvecall))
      whEN <- grmin$resEN$whichEN
      Xsamw <- Xsamw[, whEN, drop = FALSE]
      
      if(model != "survival")  predell <- predict(predobj, Xsamw, unpenalized = unpenal,  data = dataunpensam)[1:nout]
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
        if (model == "logistic") 
          refitmod <- glm(responseout ~ 1 + datlp, 
                          family = "binomial")
        else refitmod <- glm(responseout ~ 1 + datlp, 
                             family = "gaussian")
        print(paste("Recalibration intercept and slope for EN model", 
                    ell, ":", round(refitmod$coefficients[1], 
                                    3), ";", round(refitmod$coefficients[2], 
                                                   3)))
        predell2 <- predict(refitmod, type = "response")
        predellall2 <- cbind(predellall2, predell2)
      }
    }
    if (arg$compareunpenal) {
      nadd <- arg$comparelasso + arg$selectionEN
      take <- npreds + nadd + 1
      predobj <- penobj[[take]]
      
      if(model != "survival")  predell <- predict(predobj, newdata = dataunpensam, #bug repaired 19/12/2016: data -> newdata
                                                  type = "response")[1:nout]
      else { #survival (PENALIZED DOES NOT SUPPLY linear predictor...
        mmunpen <- as.matrix(model.matrix(unpenal, dataunpensam))
        wh <- match(names(predobj@unpenalized),colnames(mmunpen))
        coeffunpen <- matrix(predobj@unpenalized, ncol = 1)
        predell <- exp(mmunpen[,wh,drop=FALSE] %*% coeffunpen)  
      }
      predellall <- cbind(predellall, predell)
    }
    print(paste("Sample(s) left out:", samout))
    print("True:")
    print(response[samout])
    print("Prediction(s):")
    colnames(predellall) <- nmp
    if (recalibrate) {
      colnames(predellall2) <- paste(nmp, "_recalib", 
                                     sep = "")
      predellall <- cbind(predellall, predellall2)
    }
    print(data.frame(response = response[samout], predellall))
    preds <- rbind(preds, predellall)
  } #end outer CV for loop
  sams <- unlist(folds)
  od <- order(sams)
  mat <- data.frame(response[responsesam], preds)[od, ]
  mat <- data.frame(mat, whichfold)
  colnames(mat) <- c("TrueResponse", colnames(preds), "whichfold")
  return(mat)
} #end non-linear
