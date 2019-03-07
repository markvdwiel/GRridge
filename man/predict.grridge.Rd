\name{predict.grridge}
\alias{predict.grridge}
\alias{predict}
\title{
Predictions for new samples
}
\description{
Returns predictions for new samples from a \code{grridge} object
}

\usage{
\method{predict}{grridge}(object, datanew, printpred = FALSE, dataunpennew=NULL,
        responsetest=NULL, recalibrate=FALSE, \dots)
}

\arguments{
 \item{object}{
A model object resulted from the grridge function.
}

  \item{datanew}{
Vector or data frame. Contains the new data. For a data frame: columns are samples, rows are variables (features).
}

\item{printpred}{
Boolean. Should the predictions be printed on the screen?
}

\item{dataunpennew}{
Vector or data frame. Optional new data for unpenalized variables. NOTE: columns are covariates, rows are samples.
}

\item{responsetest}{
Factor, numeric, binary or survival. Response values of test samples. The number of response values should equal 
\code{ncol(datanew)}. Only relevant if \code{recalibrate=TRUE}.
}

\item{recalibrate}{
Boolean. Should the prediction model be recalibrated on the test samples? Only implemented for 
logistic and linear regression with only penalized covariates.
}

\item{...}{
There is no further argument is used.
}

}
\details{
This function returns predictions of the response using the \code{\link{grridge}} output. It should be applied to samples NOT used for fitting of the models. About \code{recalibrate}: we noticed that 
 recalibration of the linear predictor using a simple regression (with intercept and slope) can improve predictive performance, e.g. in terms of brier score or mean square error (not in terms of AUC, which is rank-based). We recommend to use it for large enough test sets (say >= 25 samples), in particular when one suspects that the test set could have somewhat different properties than the training set. 
For survival, the proportional hazards are returned, i.e. exp(linear predictor). If survival time predictions are desired, these may be obtained by applying \code{predict} to fit objects, which are stored as output of \code{\link{grridge}}.
}
\value{
A matrix containing the predictions from all models available in \code{grridge}. 
}

\references{
Mark van de Wiel, Tonje Lien, Wina Verlaat, Wessel van Wieringen, Saskia Wilting. (2016). 
Better prediction by use of co-data: adaptive group-regularized ridge regression.
Statistics in Medicine, 35(3), 368-81.
}

\author{
Mark A. van de Wiel
}

\seealso{
Cross-validated predictions: \code{\link{grridgeCV}}. 
Examples: \code{\link{grridge}}.
}


\examples{
#data(dataFarkas)

#firstPartition <- CreatePartition(CpGannFarkas)

#sdsF <- apply(datcenFarkas,1,sd)
#secondPartition <- CreatePartition(sdsF,decreasing=FALSE, uniform=TRUE, grsize=5000)

## Concatenate two partitions
#partitionsFarkas <- list(cpg=firstPartition, sds=secondPartition)

## A list of monotone functions from the corresponding partition
#monotoneFarkas <- c(FALSE,TRUE)

#Hold-out first two samples
#testset <- datcenFarkas[,1:2]
#resptest <- respFarkas[1:2]

#trainingset <- datcenFarkas[,-(1:2)]
#resptraining <- respFarkas[-(1:2)]

#grFarkas <- grridge(trainingset,resptraining,optl=5.680087,
#                    partitionsFarkas,monotone=monotoneFarkas)

## Prediction of the grridge model to the test samples

#Standardize variables in testset, becasuse by default standardizeX = TRUE in grridge function.
#Here, test set is small so we use the summaries of the training set. If the test set
#is large one may opt to standardize wrt to the test set itself. 

#sds <- apply(trainingset,1,sd)
#sds2 <- sapply(sds,function(x) max(x,10^{-5}))
#teststd <- (testset-apply(trainingset,1,mean))/sds2


#yhat <- predict.grridge(grFarkas,teststd)
}

