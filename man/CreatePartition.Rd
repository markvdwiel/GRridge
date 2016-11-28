\name{CreatePartition}
\alias{CreatePartition}
\title{
Creates a partition (groups) of variables
}
\description{
Creates a partition (groups) of variables from nominal (factor) or numeric input
}
\usage{
CreatePartition(vec,varnamesdata=NULL,subset=NULL,grsize=NULL,
                decreasing=TRUE,uniform=FALSE,ngroup=10,mingr=25)
}
\arguments{
  \item{vec}{
Factor, numeric vector or character vector. 
}

\item{subset}{
Character vector. Names of variables (features) that correspond to the values in \code{vec}. 
Allows to make a partition on a subset of all variables. Requires \code{varnamesdata}. 
}

\item{varnamesdata}{
Character vector. Names of the variables (features). Only relevant when \code{vec} is a character vector OR when \code{subset} is specified.
}

  \item{grsize}{
Numeric. Size of the groups. Only relevant when \code{vec} is a numeric vector and \code{uniform=TRUE}.
}
  \item{decreasing}{
Boolean. If \code{TRUE} then \code{vec} is sorted in decreasing order.
}

\item{uniform}{
Boolean. If \code{TRUE} the group sizes are as equal as possible.
}

\item{ngroup}{
Numeric. Number of the groups to create. Only relevant when \code{vec} is a numeric vector.
}

\item{mingr}{
Numeric. Minimum group size. Only relevant when \code{vec} is a numeric vector and \code{uniform=FALSE}.
}
}

\details{
A convenience function to create partitions of variables from external information that is stored in \code{vec}. If \code{vec} is 
a factor then the levels of the factor define the groups. If \code{vec} is a character vector,
then \code{varnamesdata} need to be specified (\code{vec} is supposed to be a subset of \code{varnamesdata}, e.g. a published gene list). 
In this case a partition of two groups is created: one with those variables of \code{varnamesdata} that also appear in \code{vec} and 
one which do not appear in \code{vec}. If \code{vec} is a numeric vector, then groups contain the variables corresponding to 
\code{grsize} consecutive values of the values in \code{vec}. Alternatively, the group size 
is determined automatically from \code{ngroup}. If \code{uniform=FALSE}, a group with rank $r$ is  
of approximate size \code{mingr*(r^f)}, where \code{f>1} is determined such that the total number of groups equals \code{ngroup}.
Such unequal group sizes enable the use of fewer groups (and hence faster computations) while still maintaining a 
good `resolution' for the extreme values in \code{vec}. About \code{decreasing}: if smaller values of components of \code{vec}
mean `less relevant' (e.g. test statistics, absolute regression coefficients) use \code{decreasing=TRUE}, else use \code{decreasing=FALSE}, e.g. for p-values. If \code{subset} is defined, then \code{varnamesdata} should be specified as well. The parition will then only be applied to variables in \code{subset} and in 
\code{varnamesdata}.
}
\value{
A list the components of which contain the indices of the variables belonging to each of the groups. 
}

\author{
Mark A. van de Wiel
}

\seealso{
For gene sets (overlapping groups): \code{\link{matchGeneSets}}. 
Further example in real life dataset: \code{\link{grridge}}.

}

\examples{
#SOME EXAMPLES ON SMALL NR OF VARIABLES

#EXAMPLE 1: partition based on known gene signature
genset <- sapply(1:100,function(x) paste("Gene",x))
signature <- sapply(seq(1,100,by=2),function(x) paste("Gene",x))
SignaturePartition <- CreatePartition(signature,varnamesdata=genset)

#EXAMPLE 2: partition based on factor variable
Genetype <- factor(sapply(rep(1:4,25),function(x) paste("Type",x)))
TypePartition <- CreatePartition(Genetype)

#EXAMPLE 3: partition based on continuous variable, e.g. p-value
pvals <- rbeta(100,1,4)

#Creating a partition of 10 equally-sized groups, corresponding to increasing p-values.
PvPartition <- CreatePartition(pvals, decreasing=FALSE,uniform=TRUE,ngroup=10)

#Alternatively, create a partition of 5 unequally-sized groups,
#with minimal size at least 10. Group size
#increases with less relevant p-values.
# Recommended when nr of variables is large.
PvPartition2 <- CreatePartition(pvals, decreasing=FALSE,uniform=FALSE,ngroup=5,mingr=10)

#EXAMPLE 4: partition based on subset of variables,
#e.g. p-values only available for 50 genes. 
genset <- sapply(1:100,function(x) paste("Gene",x))

subsetgenes <- sort(sapply(sample(1:100,50),function(x) paste("Gene",x)))

pvals50 <- rbeta(50,1,6)

#Returns the partition for the subset based on the indices of 
#the variables in entire genset. Variables not
#present in subset will receive group-penalty = 1 for this partition. 

PvPartitionSubset <- CreatePartition(pvals50, varnamesdata = genset,subset = subsetgenes,
                                     decreasing=FALSE,uniform=TRUE, ngroup=5)

#EXAMPLE 5: COMBINING PARTITIONS

#Combines partitions into one list with named components. 
#This can be use as input for the grridge() #function.
#NOTE: if one aims to use one partition only, then this can be directly used in grridge(). 

MyPart <- list(signature=SignaturePartition, type = TypePartition,
               pval = PvPartition, pvalsubset=PvPartitionSubset)
}


