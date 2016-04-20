\name{CreatePartition}
\alias{CreatePartition}
\title{
Creates a partition (groups) of variables
}
\description{
Creates a partition (groups) of variables from nominal (factor) or numeric input
}
\usage{
CreatePartition(vec,varNames=NULL,grsize=NULL,decreasing=TRUE,uniform=FALSE,ngroup=100,mingr=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
Factor, numeric vector or character vector. 
}

\item{varNames}{
Character vector. Names of the variables (features). Only relevant when \code{vec} is a character vector.
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
then \code{varNames} need to be specified (\code{vec} is supposed to be a subset of \code{varNames}, e.g. a published gene list). 
In this case a partition of two groups is created: one with those variables of \code{varNames} that also appear in \code{vec} and 
one which do not appear in \code{vec}. If \code{vec} is a numeric vector groups, contain the variables corresponding to 
\code{grsize} consecutive values of the values in \code{vec}. Alternatively, the group size 
is determined automatically from \code{ngroup}. If \code{uniform=FALSE}, a group with rank $r$ is  
of approximate size \code{mingr*(r^f)}, where \code{f>1} is determined such that the total number of groups equals \code{ngroup}.
Such unequal group sizes enable the use of fewer groups (and hence faster computations) while still maintaining a 
good `resolution' for the extreme values in \code{vec}. About \code{decreasing}: if smaller values of components of \code{vec}
mean `less relevant' (e.g. test statistics, absolute regression coefficients) use \code{decreasing=TRUE}, else use 
\code{decreasing=FALSE}, e.g. for p-values.
}
\value{
A list the components of which contain the indices of the variables belonging to each of the groups. 
}

\author{
Mark A. van de Wiel
}

\seealso{
Examples: \code{\link[GRridge]{grridge}}.
}

