\name{grridge}
\alias{grridge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Group-regularized (logistic) ridge regression 
}
\description{
This function implements adaptive group-regularized (logistic) ridge regression by use of co-data. It uses co-data to improve 
predictions of binary and continuous response from high-dimension (e.g. genomics) data. Here, co-data is 
auxiliary information on variables (e.g. genes), such as annotation or p-values from other studies.
}
\usage{
grridge(highdimdata, response, partitions, unpenal = ~1, offset=NULL, method="stable", niter=10, monotone=NULL, optl=NULL, innfold=NULL, 
fixedfoldsinn=TRUE, selection=FALSE,maxsel=100,stepsel=1,cvlmarg=1,savepredobj="last", dataunpen=NULL, ord = 1:length(partitions),
comparelasso=FALSE,optllasso=NULL,compareunpenal=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{highdimdata}{
Matrix or numerical data frame. Contains the primary data of the study. Columns are samples, rows are variables (features).
}
  \item{response}{
Factor, numeric or binary. Response values. The \code{length(response)} should equal \code{ncol(highdimdata)}.
}
  \item{partitions}{
List of lists. Each list component contains a partition of the variables, which is again a list. See details.
}
  \item{unpenal}{
Formula. Includes unpenalized variables. Set to \code{unpenal = ~0} if an intercept is not desired.
}

\item{offset}{
Numeric (vector). Optional offset, either one constant or sample-specific, in which case \code{length(offset)=ncol(highdimdata)}
}

\item{method}{
Character. Equal to \code{"stable"}: the iterative method, \code{"exact"}: the non-iterative, systems-based method, \code{"adaptridge"}:
adaptive ridge (not recommended). 
}

  \item{niter}{
Integer. Maximum number of re-penalization iterations.
}
  \item{monotone}{
Vector of booleans. If the jth component of \code{monotone} equals \code{TRUE}, then 
the group-penalties are forced to be monotone. If \code{monotone=NULL} monotony is not imposed for any partition. 
}
  \item{optl}{
Numeric. Value of the global regularization parameter (lambda). If specified, it skips optimization by cross-validation.
}
  \item{innfold}{
Integer. The fold for cross-validating the global regularization parameter lambda and for computing cross-validated likelihoods. 
Defaults too LOOCV.
}

  \item{fixedfoldsinn}{
Boolean. Use fixed folds for inner cross-validation? 
}

  \item{selection}{
Boolean. If \code{selection=TRUE} post-hoc variable selection is performed.
}
  \item{maxsel}{
Integer. The maximum number of selected variables.
}
  \item{stepsel}{
Integer. Step-size for variable selection.
}
  \item{cvlmarg}{
Numeric. Maximum margin (in percentage) that the cross-validated likelihood of the model with selected variables may deviate from the optimum one.
}
  \item{savepredobj}{
Character. If \code{savepredobj="last"}, only the last penalized prediction object is saved;
if \code{savepredobj="all"} all are saved; if \code{savepredobj="none"}, none are saved. 
}
  \item{dataunpen}{
Data frame. Optional data for unpenalized variables.
}
  \item{ord}{
Integer vector. The order in which the partitions in \code{partitions} are used. 
}

  \item{comparelasso}{
Boolean. If \code{comparelasso=TRUE} the results of lasso regression are included.
}

 \item{optl}{
Numeric. Value of the global regularization parameter (lambda) in the lasso. If specified, 
optimization by cross-validation is skipped.
}

  \item{comparelasso}{
Boolean. If \code{compareunpenal=TRUE} the results of regression with unpenalized covariates only are included.
Only relevant when \code{dataunpenal} is specified. 
}
}
\details{
About \code{partitions}: this is a list of partitions.  
Each partition is a (named) list that contains the indices (row numbers) of the variables in the concerning group. Such a partition is usually created by 
\code{\link[GRridge]{CreatePartition}}. All indices should be represented in the partition. 
NOTE: the length of \code{partitions} should equal the number of partitions. If you use only one partition, then use for example
\code{mypartitions <- list(mypartition=onepartition)}. 

About \code{savepredobj}: use \code{savepredobj="all"} if you want to compare performances of the various predictors (e.g. ordinary ridge, 
group-regularized ridge, group-regularized ridge + selection) using \code{\link[GRridge]{grridge.cv}}.
About \code{monotone}: We recommend to set the jth component of \code{monotone} to \code{TRUE} when the jth partition 
is based on external p-values, test statistics or regression coeeficients. This increases stability of the predictions. 
About \code{cvlmarg}: We recommended to use values between 0 and 2. A larger value will generally result in fewer selected variables. 
About \code{innfold}: for large data sets considerable computing time may be saved when setting \code{innfold=10} instead of default 
leave-one-out-cross-validation (LOOCV). About \code{method}: \code{"stable"} is recommended. If the number of variables is
not very large, say \code{<2000}, the faster non-iterative \code{"exact"} method can be used as an alternative. 
}
\value{
A list object containing:
 \item{true}{True values of the response}
  \item{cvls}{Cross-validated likelihoods from the iterations}
  \item{lambdamults}{List of lists object containing the penalty multipliers per group per partition}
  \item{optl}{Global penalty parameter lambda}
  \item{lambdamultvec}{Vector with penalty multipliers per variable}
  \item{predobj}{List of prediction objects}
  \item{betas}{Estimated regression coefficients}
  \item{whichsel}{Indices of selected variables}
  \item{cvlssel}{Trace of cross-validated likelihoods for the number of selected variables}
  \item{reslasso}{Results of the lasso. \code{NULL} when \code{comparelasso=FALSE}}
  \item{arguments}{Arguments used to call the function}
}

\references{
Mark A. van de Wiel, Tonje G. Lien, Wina Verlaat, Wessel N. van Wieringen, Saskia M. Wilting (2015).  
Better prediction by use of co-data: Adaptive group-regularized ridge regression. 
Statistics in Medicine. Preliminary version: http://arxiv.org/abs/1411.3496.
}

\author{
Mark A. van de Wiel
}

\seealso{
Creating partitions: \code{\link[GRridge]{CreatePartition}};
Cross-validation for assessing predictive performance: \code{\link[GRridge]{grridge.cv}}.

}
\examples{
\dontrun{
#NOTE: EXAMPLE DEVIATES SOMEWHAT FROM THE EXAMPLE IN THE MANUSCRIPT IN ORDER TO SHOW SOME
#OTHER FUNCTIONALITIES.

# 1ST EXAMPLE: Farkas DATA, USING ANNOTATION: DISTANCE TO CpG

#load data objects: datcenFarkas: methylation data for cervix samples (arcsine-transformed beta values), 
#respFarkas: binary response (Normal and Precursor) and CpGannFarkas: annotation of probes according to location 
#(CpG-Island, North-Shelf, South-Shelf, North-Shore, South-Shore, Distant) 
data(dataFarkas)

#Create list of partition(s), here only one
partitionFarkas <- list(cpg=CreatePartition(CpGannFarkas))

#Group-regularized ridge applied to data datcenFarkas, response respFarkas and partition partitionFarkas. 
#Saves the prediction objects from ordinary and group-regularized ridge.
#Includes unpenalized intercept by default.
grFarkas <- grridge(datcenFarkas,respFarkas,partitionFarkas,savepredobj="all") 

#Predictions from all models for NEW samples. Here, illustrated on 2 samples of the training data. 
fakenew <- datcenFarkas[,3:4]
predict.grridge(grFarkas,fakenew)

#10-fold cross-validation to assess predictive performances of the predictors saved in the grFarkas object.
#grridge.cv invokes the grridge function (on the smaller training sets) using the same arguments as used 
#to create grFarkas. MAY TAKE CONSIDERABLE TIME!
grFarkascv <- grridge.cv(grFarkas, datcenFarkas,respFarkas,outerfold=10)

#Computes the ROC curve and AUC for probabilistic classifiers
cutoffs <- rev(seq(0,1,by=0.01))
rocridgeF <- roc(probs=grFarkascv[,2],true=grFarkascv[,1],cutoffs=cutoffs)
auc(rocridgeF)
rocgrridgeF <- roc(probs=grFarkascv[,3],true=grFarkascv[,1],cutoffs=cutoffs)
auc(rocgrridgeF)

plot(rocridgeF[1,],rocridgeF[2,],type="l",lty=1,ann=F,col="black") 
points(rocgrridgeF[1,],rocgrridgeF[2,],type="l",lty=1,col="grey")


# 2ND EXAMPLE: Verlaat DATA, USING P-VALUES AND SIGN OF EFFECT FROM FARKAS DATA

#Load objects: datcenVerlaat (data), respVerlaat (response), pvalFarkas (p-values from external study), 
#diffmeanFarkas (external effect sizes, cases minus controls), CpGann: annotation as factor
data(dataVerlaat)

#Creates a partition of variables with 100 groups, sorted in increasing order of values of vec. 
#Here, minimal group size equals 10, group sizes gradually increase. Use decreasing = TRUE, if a large 
#(summary) value from the external data indicates more relevance (eg test statistic).
firstPartition <- CreatePartition(vec=pvalFarkas,decreasing=FALSE,mingr=10,ngroup=100)

#Create partition of variables based on sign of the effect size in external data. Positive means hypermethylated in cases.
whpos <- which(diffmeanFarkas>0)
secondPartition <- list(Pos=whpos,Neg=(1:length(diffmeanFarkas))[-whpos])

#Create annotation-based classification
thirdPartition <- CreatePartition(CpGann)

#create a fake prior probe list of size 500 to illustrate CreatePartition() on a character vector.
set.seed(252345);rn <- rownames(datcenVerlaat); fakeprobes <- sample(rn,500); 
fourthPartition <- CreatePartition(fakeprobes,varNames=rn)


#Concatenate four partitions 
partitionsVerlaat <- list(pFarkas=firstPartition,posneg=secondPartition,cpg=thirdPartition, fake=fourthPartition)

#Group-regularized ridge; unpenal=~0 implies no intercept. In this study, cases are deliberately oversampled which may bias the intercept.
#monotone =c(TRUE,FALSE,FALSE,FALSE) indicates that monotonely decreasing group-penalties are desired for the pvalFarkas-based partition (and not for 
#other partitions). Selection = TRUE indicates that CV-likelihood-based a posteriori variable selection is applied. A comparison with the lasso 
#is included.
grVerlaat <- grridge(datcenVerlaat,respVerlaat,unpenal=~0,partitionsVerlaat,monotone = c(TRUE,FALSE,FALSE,FALSE),selection=TRUE,savepredobj="all",
comparelasso=TRUE)


#Displays number of selected markers
length(grVerlaat$whichsel)

#Leave-one-out-cross-validation (default) for assessing predictive performance. 
#WILL TAKE SOME TIME! USE E.G. 5-FOLD CV (outerfold = 5) TO SPEED UP
grVerlaatcv <- grridge.cv(grVerlaat, datcenVerlaat,respVerlaat) 

#Displays true response [,1]; ridge prediction [,2]; grridge prediction[,3]; 
#grridge prediction after variable selection [,4]; lasso [,5]
grVerlaatcv

#Computes ROC curves and AUCs for the probabilistic classifiers 
cutoffs <- rev(seq(0,1,by=0.01))
rocridgeV <- roc(probs=grVerlaatcv[,2],true=grVerlaatcv[,1],cutoffs)
auc(rocridgeV)
rocgrridgeV <- roc(probs=grVerlaatcv[,3],true=grVerlaatcv[,1],cutoffs)
auc(rocgrridgeV)
rocgrridgeSelV <- roc(probs=grVerlaatcv[,4],true=grVerlaatcv[,1],cutoffs) 
auc(rocgrridgeSelV)
roclasso <- roc(probs=grVerlaatcv[,5],true=grVerlaatcv[,1],cutoffs) 
auc(roclasso)


plot(rocridgeV[1,],rocridgeV[2,],type="l",lty=1,ann=F,col="black") 
points(rocgrridgeV[1,],rocgrridgeV[2,],type="l",lty=1,col="grey")
points(rocgrridgeSelV[1,],rocgrridgeSelV[2,],type="l",lty=2,col="grey")
points(roclasso[1,],roclasso[2,],type="l",lty=2,col="red")
}
}

