\name{roc}
\alias{roc}
\title{
Produces ROC curve for probabilistic classifiers (e.g. logistic regression)
}
\description{
Computes an ROC-curve for probabilistic classifiers.}

\usage{
roc(probs, true, cutoffs)
}
\arguments{
  \item{probs}{
Numeric vector, with values between 0 and 1
}
  \item{true}{
Binary vector.
}
  \item{cutoffs}{
Numeric vector, with increasing values between 0 and 1.
}
}
\details{
The vector probs contains predicted probabilities for the response to equal 1, as produced by 
a probabilistic classifier like logistic regression. The cutoffs can simply represent a grid of values between 0 and 1. 
}

\value{
A matrix with two rows which contain corresponding False Positive and True Positive Rates for all cutoffs.
}

\author{
Mark A. van de Wiel
}

\seealso{
For area-under-the ROC-curve: \code{\link{auc}}. Examples: \code{\link{grridge}}.
}


\examples{
# Load data objects
data(dataFarkas)

firstPartition <- CreatePartition(CpGannFarkas)

sdsF <- apply(datcenFarkas,1,sd)
secondPartition <- CreatePartition(sdsF,decreasing=FALSE, uniform=TRUE, grsize=5000)

# Concatenate two partitions
partitionsFarkas <- list(cpg=firstPartition, sds=secondPartition)

# A list of monotone functions from the corresponding partition
monotoneFarkas <- c(FALSE,TRUE)

#grFarkas <- grridge(datcenFarkas,respFarkas,optl=5.680087,partitionsFarkas,monotone=monotoneFarkas)
#grFarkascv <- grridgeCV(grFarkas,datcenFarkas,respFarkas,outerfold=10)

#cutoffs <- rev(seq(0,1,by=0.01))
#rocgrridgeF <- roc(probs=grFarkascv[,3],true=grFarkascv[,1],cutoffs=cutoffs)
#rocridgeF <- roc(probs=grFarkascv[,2],true=grFarkascv[,1],cutoffs=cutoffs)
#plot(rocridgeF[1,],rocridgeF[2,],type="l",lty=1,ann=FALSE,col="grey")
#points(rocgrridgeF[1,],rocgrridgeF[2,],type="l",lty=1,col="black")
#legend(0.75,0.1, legend=c("ridge","GRridge"),
#       lty=c(1,1), lwd=c(1,1),col=c("grey","black"))

}

