\name{grridgeCV}
\alias{grridgeCV}
\title{
Returns the cross-validated predictions
}
\description{
Returns the cross-validated predictions for a \code{grridge} logistic, linear or Cox regression. 
}
\usage{
grridgeCV(grr, highdimdata, response, outerfold = length(response),
           fixedfolds = TRUE, recalibrate=FALSE)}
\arguments{
\item{grr}{
List. Output of GRridge function.
}
  \item{highdimdata}{
Matrix or numerical data frame. Contains the primary data of the study. Columns are samples, rows are variables (features).
}
  \item{response}{
Factor, numeric, binary or survival. Response values. The number of response values should equal \code{ncol(highdimdata)}.
}
  \item{outerfold}{
Integer. Fold used for cross-validation loop.
}
  \item{fixedfolds}{
Boolean. Use fixed folds for cross-validation?
}

\item{recalibrate}{
Boolean. Should the prediction model be recalibrated on the test samples? Only implemented for 
logistic and linear regression with only penalized covariates.
}
}
\details{
This convenience function returns cross-validated predictions from \code{\link{grridge}}, including those from 
ordinary (logistic) ridge regression. It can be used to compute ROC-curves. About \code{recalibrate}: this option allows to compare recalibrated models, but only if the test sample size is large enough. See 
\code{\link[GRridge]{predict.grridge}}
}
\value{
For linear and logistic regression: A matrix containing the predictions. The first column contains the sample indices, 
the second the prediction by ordinary ridge, the third the predictions by group-regularized ridge, the 
fourth (optional) the predictions by group-regularized ridge plus selection. Finally, it may contain
predictions by lasso and/or a regression model with unpenalized covariates only. 
For Cox regression: a list with three components. First, a matrix with the linear predictions, in the same format as above.
Second, a list of survival probability matrices (one for each prediction model), the rows of which indicate the individuals and the columns indicate the observed times. Third, a vector with the observed times. 
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
For logistic regression: ROC-curves: \code{\link{roc}}. AUC: \code{\link{auc}}.
GRridge: \code{link{grridge}}.
}

\examples{

# load data objects
data(dataFarkas)

# In this example, we provide one partition only
# see "CreatePartition" for examples in creating multiple partitions
firstPartition <- CreatePartition(CpGannFarkas)

# the optimum lambda2 is provided in this example
# worth to try:
# grFarkas <- grridge(datcenFarkas,respFarkas,firstPartition,monotone=FALSE)
# grFarkas$optl
grFarkas <- grridge(datcenFarkas,respFarkas, optl=5.680087,
                      firstPartition,monotone=FALSE)

###
# grFarkasCV <- grridgeCV(grFarkas,datcenFarkas,respFarkas,outerfold=10)

}
