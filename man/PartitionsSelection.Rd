\name{PartitionsSelection}
\alias{PartitionsSelection}
\title{
 Co-data selection in a Group-regularized ridge regression model
}

\description{
This function implements a procedure to optimize the use of co-data in a GRridge model. Although there is no harm to include as much as co-data in a GRridge model, ordering and selecting co-data can optimized the performance of a GRridge model. This procedure is similar with forward feature selection in a classical regression model
}

\usage{
PartitionsSelection(highdimdata, response, partitions, 
                    monotoneFunctions, optl=NULL, innfold=NULL) 
}

\arguments{
  \item{highdimdata}{
Matrix or numerical data frame. Contains the primary data of the study. Columns are samples, rows are features.
}
  \item{response}{
Factor, numeric, binary or survival. Response values. The number of response values should equal \code{ncol(highdimdata)}.
}
  \item{partitions}{
List of lists. Each list component contains a partition of the variables, which is again a list.
}
  \item{monotoneFunctions}{
  Vector. Monotone functions from each partition. This argument is necesarily specified. If the jth component of monotone equals TRUE, then the group-penalties are forced to be monotone. If \code{monotone=NULL} monotony is not imposed for any partition.
}
  \item{optl}{
  Global penalty parameter lambda
  }
  \item{innfold}{
  Integer. The fold for cross-validating the global regularization parameter lambda and for computing cross-validated likelihoods. Defaults  LOOCV.
}
}

\value{
A list containing (i) the indeces of the selected and ordered partitions and (ii) the optimum lambda penalty from the ridge regression. 
}

\author{
Putri W. Novianti
}

\seealso{
Creating partitions: \code{\link{CreatePartition}}. 
Group-regularized ridge regression: \code{\link{grridge}}. 
}

\examples{

# # Load data objects
# data(dataWurdinger)
# 
# # Transform the data set to the square root scale
# dataSqrtWurdinger <- sqrt(datWurdinger_BC)
# 
# #Standardize the transformed data
# datStdWurdinger <- t(apply(dataSqrtWurdinger,1,function(x){(x-mean(x))/sd(x)}))
# 
# # A list of gene names in the primary RNAseq data
# genesWurdinger <- annotationWurdinger$geneSymbol
# 
# # co-data 1: a partition based on immunologic signature pathway
# # The initial gene sets (groups) are merged into five new groups, using the "mergeGroups" function
# immunPathway <- coDataWurdinger$immunologicPathway
# parImmun <- immunPathway$newClust
# 
# # co-data 2: a partition based on chromosomal location
# ChromosomeWur0 <- as.vector(annotationWurdinger$chromosome)
# ChromosomeWur <- ChromosomeWur0
# idC <- which(ChromosomeWur0=="MT" | ChromosomeWur0=="notBiomart" | ChromosomeWur0=="Un")
# ChromosomeWur[idC] <- "notMapped"
# parChromosome <- CreatePartition(as.factor(ChromosomeWur))
# 
# # co-data 3: a partition based on a list of platelets expressed genes
# plateletsExprGenes <- coDataWurdinger$plateletgenes
# # Group genes in the primary data based on the list
# # The genes are grouped into
# # either "NormalGenes" or "Non-overlapGenes"
# is <- intersect(plateletsExprGenes,genesWurdinger)
# im <- match(is, genesWurdinger)
# plateletsGenes <- replicate(length(genesWurdinger),"Non-overlapGenes")
# plateletsGenes[im] <- "NormalGenes"
# plateletsGenes <- as.factor(plateletsGenes)
# parPlateletGenes <- CreatePartition(plateletsGenes)
# 
# # co-data 4: a partition based on chromosomal location. 
# # A list of chromosomal location based on {\tt biomaRt} data bases.
# ChromosomeWur0 <- as.vector(annotationWurdinger$chromosome)
# ChromosomeWur <- ChromosomeWur0
# idC <- which(ChromosomeWur0=="MT" | ChromosomeWur0=="notBiomart" |
#                ChromosomeWur0=="Un")
# ChromosomeWur[idC] <- "notMapped"
# table(ChromosomeWur)
# parChromosome <- CreatePartition(as.factor(ChromosomeWur))
# 
# partitionsWurdinger <- list(immunPathway=parImmun,
#                             transcriptionFactor=parTranscriptFactor,
#                             plateletsGenes=parPlateletGenes,
#                             chromosome=parChromosome)
# 
# #A list of monotone functions from the corresponding partitions
# monotoneWurdinger <- c(FALSE,FALSE,FALSE,FALSE)
# 
# # Start ordering and selecting partitions
# optPartitions <- PartitionsSelection(datStdWurdinger, respWurdinger, 
#                                      partitions=partitionsWurdinger,
#                                      monotoneFunctions=monotoneWurdinger)
}


