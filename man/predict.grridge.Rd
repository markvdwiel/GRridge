\name{predict.grridge}
\alias{predict.grridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions for new samples
}
\description{
Returns predictions for new samples for a \code{\link[GRridge]{grridge}} object
}
\usage{
predict.grridge(grr, datanew, printpred = FALSE, dataunpennew=NULL, responsetest=NULL, recalibrate=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{grr}{
List. Output of GRridge function.
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
}
\details{
This function returns predictions of the response using the \code{\link[GRridge]{grridge}} output. It should be applied to samples NOT used for fitting of the models. About \code{recalibrate}: we noticed that 
 recalibration of the linear predictor using a simple regression (with intercept and slope) can improve predictive performance, e.g. in terms of brier score or mean square error (not in terms of AUC, which is rank-based). We recommend to use it for large enough test sets (say >= 25 samples), in particular when one suspects that the test set could have somewhat different properties than the training set. 
}
\value{
A matrix containing the predictions from all models available in \code{grr}. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mark A. van de Wiel
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
Cross-validated predictions: \code{\link[GRridge]{grridge.cv}}. Examples: \code{\link[GRridge]{grridge}}.
}

