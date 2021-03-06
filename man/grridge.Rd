\name{grridge}
\alias{grridge}
\title{
Group-regularized (logistic) ridge regression 
}
\description{
This function implements adaptive group-regularized (logistic) ridge regression by use of co-data. It uses co-data to improve predictions of binary and continuous response from high-dimension (e.g. genomics) data. Here, co-data is auxiliary information on variables (e.g. genes), such as annotation or p-values from other studies.
}
\usage{
grridge(highdimdata, response, partitions, unpenal = ~1, 
        offset=NULL, method="exactstable",
        niter=10, monotone=NULL, optl=NULL, innfold=NULL, 
        fixedfoldsinn=TRUE, maxsel=c(25,100),selectionEN=FALSE,cvlmarg=1,
        savepredobj="all", dataunpen=NULL, ord = 1:length(partitions),
        comparelasso=FALSE,optllasso=NULL,cvllasso=TRUE,
        compareunpenal=FALSE,trace=FALSE,modus=1,
        EBlambda=FALSE,standardizeX = TRUE)       
}
\arguments{
  \item{highdimdata}{
Matrix or numerical data frame. Contains the primary data of the study. Columns are samples, rows are variables (features).
}
  \item{response}{
Factor, numeric, binary or survival. Response values. The number of response values should equal \code{ncol(highdimdata)}.
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
Character. Equal to \code{"exactstable"}: the stable iterative, systems-based method, \code{"stable"}: the iterative non-systems-based method, \code{"exact"}: the non-iterative, systems-based method, \code{"adaptridge"}:
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

  \item{selectionEN}{
Boolean. If \code{selectionEN=TRUE} then post-hoc variable selection by weighted elastic net is performed.
}
  \item{maxsel}{
Vector of integers. The maximum number of selected variables. Can be multiple to allow comparing models of various sizes.
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

 \item{optllasso}{
Numeric. Value of the global regularization parameter (lambda) in the lasso. If specified, 
optimization by cross-validation is skipped.
}

 \item{cvllasso}{
Boolean. If \code{cvllasso=TRUE} it returns the cross-validated likelihood for lasso when \code{comparelasso=TRUE}.
}

  \item{compareunpenal}{
Boolean. If \code{compareunpenal=TRUE} the results of regression with unpenalized covariates only are included.
Only relevant when \code{dataunpenal} is specified. 
}

\item{trace}{
Boolean. If \code{trace=TRUE} the results of the cross-validation for parameter (lambda) tuning are shown.
}

\item{modus}{
Integer. Please use \code{modus=1}. Only use \code{modus=2} when backward compatibility with versions <= 1.6 is desired.
}

\item{EBlambda}{
Boolean. If \code{EBlambda=TRUE} global lambda is estimated by empirical Bayes (currently only available for linear model).
}
\item{standardizeX}{
Boolean. If \code{standardizeX=TRUE} variables in X are standardized prior to the analysis.
}

}
\details{
About \code{partitions}: this is a list of partitions or one partition represented as a simple list.  
Each partition is a (named) list that contains the indices (row numbers) of the variables in the concerning group. Such a partition is usually created by 
\code{\link{CreatePartition}}. 
About \code{savepredobj}: use \code{savepredobj="all"} if you want to compare performances of the various predictors (e.g. ordinary ridge, 
group-regularized ridge, group-regularized ridge + selection) using \code{\link{grridgeCV}}.
About \code{monotone}: We recommend to set the jth component of \code{monotone} to \code{TRUE} when the jth partition 
is based on external p-values, test statistics or regression coeeficients. This increases stability of the predictions. If \code{selectionEN=TRUE}, EN selection will, for all elements m of \code{maxsel}, select exactly m or fewer variables. Note that EN is only used for selection; 
the final predictive model is a group-ridge model fitted only on the selected variables using the penalties estimated by GRridge. Using multiple values for 
\code{maxsel} allows comparing models of various sizes, also in terms of cross-validated performance when using \code{\link{grridgeCV}}.
About \code{cvlmarg}: We recommended to use values between 0 and 2. A larger value will generally result in fewer selected variables by forward selection. 
About \code{innfold}: for large data sets considerable computing time may be saved when setting \code{innfold=10} instead of default leave-one-out-cross-validation (LOOCV). About \code{method}: \code{"exactstable"} is recommended. If the number of variables is not very large, say \code{<2000}, the faster non-iterative \code{"exact"} method can be used as an alternative. \code{grridge} uses the \code{penalized} package to fit logistic and survival ridge models; \code{glmnet} is used for linear response and for fitting lasso when \code{comparelasso=TRUE}. 
}
\value{
A list object containing:
 \item{true}{True values of the response}
  \item{cvfit}{Measure of fit. Cross-validated likelihoods from the iterations for linear and survival model; minus CV error for linear model}
  \item{lambdamults}{List of lists object containing the penalty multipliers per group per partition}
  \item{optl}{Global penalty parameter lambda}
  \item{lambdamultvec}{Vector with penalty multipliers per variable}
  \item{predobj}{List of prediction objects}
  \item{betas}{Estimated regression coefficients}
  \item{reslasso}{Results of the lasso. \code{NULL} when \code{comparelasso=FALSE}}
  \item{resEN}{Results of the Elastic Net selection for all elements of \code{maxsel}. \code{list()} when \code{selectionEN=FALSE}}
  \item{model}{Model used for fitting: logistic, linear or survival}
  \item{arguments}{Arguments used to call the function}
  \item{allpreds}{Predictions on the same data}
}

\references{
Mark van de Wiel, Tonje Lien, Wina Verlaat, Wessel van Wieringen, Saskia Wilting. (2016). 
Better prediction by use of co-data: adaptive group-regularized ridge regression.
Statistics in Medicine, 35(3), 368-81.

Novianti PW, Snoek B, Wilting SM, van de Wiel MA (2017). 
Better diagnostic signatures from RNAseq data through use of auxiliary co-data. 
Bioinformatics, 33, 1572-1574.
}

\author{
Mark A. van de Wiel
}

\seealso{
Creating partitions: \code{\link{CreatePartition}};
Cross-validation for assessing predictive performance: \code{\link{grridgeCV}}.

}
\examples{

## NOTE: 
## 1. EXAMPLE DEVIATES SOMEWHAT FROM THE EXAMPLE IN THE MANUSCRIPT IN ORDER TO SHOW SOME
##    OTHER FUNCTIONALITIES.
## 2. HERE WE SHOW A SIMPLE EXAMPLE FROM THE FARKAS DATA SET 
## MORE EXTENSIVE EXAMPLES OF FUNCTIONALITIES IN THE GRRIGDE PACKAGE ARE PROVIDED IN 
## VIGNETTE DOCUMENTATION FILE


## 1ST EXAMPLE: Farkas DATA, USING ANNOTATION: DISTANCE TO CpG

##load data objects:
##datcenFarkas: methylation data for cervix samples (arcsine-transformed beta values)
##respFarkas: binary response (Normal and Precursor)
##CpGannFarkas: annotation of probes according to location
##(CpG-Island, North-Shelf, South-Shelf, North-Shore, South-Shore, Distant) 
data(dataFarkas)

##Create list of partition(s), here only one partition included
partitionFarkas <- list(cpg=CreatePartition(CpGannFarkas))

##Group-regularized ridge applied to data datcenFarkas, 
##response respFarkas and partition partitionFarkas. 
##Saves the prediction objects from ordinary and group-regularized ridge.
##Includes unpenalized intercept by default.

#grFarkas <- grridge(datcenFarkas,respFarkas, optl=5.680087,
#                      partitionFarkas,monotone=FALSE)

## 2ND EXAMPLE: Verlaat DATA, USING P-VALUES AND SIGN OF EFFECT FROM FARKAS DATA
## see vignette documentation file!
}

