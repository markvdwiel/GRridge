\name{auc}
\alias{auc}
\title{
Area under the ROC curve
}
\description{
Computes Area-under-the-ROC-curve}
\usage{
auc(rocout)
}
\arguments{
  \item{rocout}{
Matrix with two rows
}
}
\details{
Operates on the output of \code{\link{roc}}. The rows of the input matrix represent the False Positive Rates (FPR) and 
corresponding True Positive Rates (TPR) at fixed thresholds. 
}
\value{
Numeric. Interpolated Area-under-the-ROC-curve.
}

\author{
Mark A. van de Wiel
}

\seealso{
ROC-curves: \code{\link{roc}}. creating multiple partitions: \code{\link{CreatePartition}}. 
Examples: \code{\link{grridge}}.
}


\examples{
# Load data objects
data(dataFarkas)

## In this example, we provide one partition only
## see "CreatePartition" for examples in creating multiple partitions
firstPartition <- CreatePartition(CpGannFarkas)

# grFarkas <- grridge(datcenFarkas,respFarkas, firstPartition,monotone=FALSE)

## Prediction of the grridge model to the training samples
#cutoffs <- rev(seq(0,1,by=0.1))
#fakenew <- datcenFarkas
#yhat <- predict.grridge(grFarkas,fakenew)

#rocridgeF <- roc(probs=as.numeric(yhat[,2]),true=respFarkas[1:30],cutoffs=cutoffs)
#auc(rocridgeF)
}

