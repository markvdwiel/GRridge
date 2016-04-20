grridge <-
function(highdimdata, response, partitions, unpenal = ~1, offset=NULL, method="stable", niter=10, monotone=NULL, optl=NULL, innfold=NULL, 
fixedfoldsinn=TRUE, selection=FALSE,maxsel=100,stepsel=1,cvlmarg=1,savepredobj="last", dataunpen=NULL, ord = 1:length(partitions),
comparelasso=FALSE,optllasso=NULL,compareunpenal=FALSE){#
#method: "stable", "exact","adaptridge"
#highdimdata<-betasfcen;response<- fl;unpenal=~0;offset=NULL;partitions<-partitionsVerlaat;monotone = c(TRUE,FALSE);
#selection=TRUE;savepredobj="all"; optl=3;innfold=NULL; selection=TRUE;maxsel=100;stepsel=1;cvlmarg=1
#dataunpen=NULL;niter<-10;ord <- 1:2
#partitions <- list(partitionsVerlaat[[2]],partitionsVerlaat[[1]]);monotone=c(FALSE,TRUE)
#highdimdata <- t(X);response <- y;calssif <- partitions; unpenal = ~0;printl=FALSE; niter=10; monotone=NULL; optl=NULL; innfold=NULL;selection=FALSE;maxsel=100;stepsel=1;savepredobj="last"; dataunpen=NULL
#highdimdata <- gesel;response <- nodestat;partitions <- clroepman; printl=T;printp=T;unpenal = ~0; niter=10; monotone=NULL; optl=NULL; innfold=10;selection=FALSE;maxsel=100;stepsel=1;savepredobj="last"; dataunpen=NULL
#response: factor with two levels (logistic ridge regression) or numeric (linear ridge regression, except when only 0-1 values appear). 
#savepredobj = "last": save only last pred obcject; "all": save all; "none": save none
#highdimdata <- betasfcen;response<-fl;monotone <- c(TRUE);niter <- 10;partitions <- partitionwina;innfold=NULL;
#unpenal=~0+hpv16;selection=TRUE;maxsel=100;stepsel=1;optl=37.32;dataunpen=hpv16dat
#addintercept=TRUE;addunpenal=TRUE
#... futher arguments for penalized, optL2 and cvl, such as unpenalized
#model is determined form the response
#cvlmarg: Margin for CVL in percentages when selecting
#printl: print update for lambda multiplier
#niter: setting it to zero renders results for ordinary ridge regression
#partitions<-partitionsVerlaat;monotone = c(TRUE,FALSE,FALSE);ord <- c(1,2,3);highdimdata<-datcenVerlaat
#highdimdata<-datcenFarkas;response<-respFarkas;partitions<- partitionFarkas;offset=NULL;monotone = NULL;
#selection=TRUE;savepredobj="all";optl=5.68;innfold=NULL; selection=TRUE;maxsel=100;stepsel=1;cvlmarg=1
#dataunpen=NULL;niter<-10;ord <- 1;unpenal = ~1;fixedfoldsinn=TRUE;method<-"adapt ridge";controlbound<-5
#highdimdata<-mirnormcen;response<-resp;partitions<- partkeep;offset=NULL;monotone = c(TRUE,TRUE);
#selection=FALSE;savepredobj="all";optl=192;innfold=NULL; maxsel=100;stepsel=1;cvlmarg=1
#dataunpen=NULL;niter<-10;ord <- 1:2;unpenal = ~1;fixedfoldsinn=TRUE;comparelasso=TRUE;optllasso=NULL
#highdimdata<-mirnormcen;response<-resp;partitions<- partkeep;offset=NULL;monotone = c(TRUE);
#selection=FALSE;savepredobj="all";optl=664;innfold=10; maxsel=100;stepsel=1;cvlmarg=1;compareunpenal=T;
#dataunpen=data.frame(therapy=therapy);niter<-10;ord <- 1;unpenal = ~1+therapy;fixedfoldsinn=TRUE;comparelasso=TRUE;optllasso=7.06;

#see grridgeExample for creation of data

#sizes = c(5000,2000,1300,30,50,300)
#selectuneq <- sapply(1:6,function(i){part <- partitionFarkas[[1]][[i]]; return(sample(part,size=sizes[i]))})
#wh <- unlist(selectuneq)
#datcenF <- datcenF
#CpGF<- CpGannFarkas
#partitionF <- list(cpg=CreatePartition(CpGF))
#
#partitionFarkas <- list(cpg=CreatePartition(CpGannFarkas))

#highdimdata<-datcenFarkas;response<-respFarkas;partitions<- partitionFarkas;offset=NULL;monotone = NULL;
#selection=TRUE;savepredobj="all";innfold=10; selection=FALSE;maxsel=100;stepsel=1;cvlmarg=1
#dataunpen=NULL;niter<-10;ord <- 1;unpenal = ~0;fixedfoldsinn=TRUE;compareunpenal=F;comparelasso=F;optllasso=NULL;
#controlbound<-10;optl=5.68
#

if(method=="adaptridge" | method=="exact") niter <- 1
nclass <- length(partitions)
partitions <- partitions[ord]
monotone <- monotone[ord]
nr <- nrow(highdimdata)

#check whether each partition is correct
for(ncl in 1:nclass){
indexset <- unlist(partitions[[ncl]])
if(length(indexset) != nr){
print(paste("Error: partition",ncl,"does not contain all row indices of the data"))
return(NULL)
}
}

#NEW
arguments <- list(partitions=partitions,unpenal=unpenal, offset=offset, method=method, niter=niter, monotone=monotone, optl=optl, innfold=innfold, fixedfoldsinn=fixedfoldsinn, 
selection=selection, maxsel=maxsel,stepsel=stepsel, cvlmarg=cvlmarg, dataunpen=dataunpen,savepredobj=savepredobj, ord=ord, 
comparelasso=comparelasso, optllasso=optllasso, compareunpenal=compareunpenal)


if(nr > 10000 & is.null(innfold)) print("NOTE: consider setting innfold=10 to save computing time")

nmp0 <- names(partitions)
if(is.null(nmp0)) nmp0 <- sapply(1:length(partitions),function(i) paste("Grouping",i))  
nmp0 <- sapply(1:length(partitions),function(i) {
if(nmp0[i]=="") return(paste("Grouping",i)) else return(nmp0[i])})


nmp <- c("NoGroups","GroupRegul")
nmpweight <- nmp #to be used later
if(selection) nmp <- c(nmp,"ForwSel") 
if(comparelasso) nmp <- c(nmp,"lasso") #NEW
if(compareunpenal) nmp <- c(nmp,"modelunpen") #NEW


if(class(response) =="factor") {
nlevel <- length(levels(response))
    if(nlevel != 2){
    print("Response is not binary, so not suitable for two-class classification.")
    return(NULL)
    } else {
    model = "logistic"
    print("Binary response, executing logistic ridge regression")
    lev <- levels(response)
    print(paste("Predicting probability on factor level",lev[2]))
    }} else {
    if(class(response) == "numeric"){
    valresp <- sort(unique(response)) 
    if(length(valresp)==2 & valresp[1]==0 & valresp[2]==1) {
        model = "logistic"
        print("Binary response, executing logistic ridge regression")
        } else {
        model = "linear"
        print("Numeric continuous response, executing linear ridge regression")
        }
    } else {
    print("Non-valid response. Should be binary or numeric.")
    return(NULL)
    }
}

if((unpenal != ~0) & (unpenal != ~1)) {
if(is.null(dataunpen)) {print("If unpenal contains variables, data of the unpenalized variables should be specified in the data slot!")
return(NULL)
} 
}


nsam <- ncol(highdimdata)
#offset <- log(1/9);unpenal= ~0;nsam=10
if(!is.null(offset)){
noffs <- length(offset)
offsets <-"c("
if(noffs==1) {for(i in 1:(nsam-1)) offsets <- paste(offsets,offset,",",sep="")} else {
for(i in 1:(nsam-1)) offsets <- paste(offsets,offset[i],",",sep="")}
if(noffs==1) offsets <- paste(offsets,offset,")",sep="")  else offsets <- paste(offsets,offset[nsam],")",sep="")
if((unpenal != ~0) & (unpenal != ~1)){
unpenal <- formula(paste(deparse(unpenal),"+ offset(",offsets,")",sep=""))
} else {
unpenal <- formula(paste("~","offset(",offsets,")",sep=""))
}
}


if(is.null(dataunpen)) datapred <- data.frame(fake=rep(NA,ncol(highdimdata))) else datapred <- dataunpen
nopen <- unpenal

if(is.null(innfold)) foldinit <- nsam else foldinit <- innfold
pmt0<- proc.time()
optl0 <- optl #optl0 for later use
if(is.null(optl)){
print("Finding lambda for initial ridge regression")
if(fixedfoldsinn) set.seed(346477)
opt <- optL2(response, penalized = t(highdimdata),fold=foldinit,trace=T,unpenalized=nopen,data=datapred)
#opt <- optL2(response, penalized = t(highdimdata),fold=foldinit,trace=T,unpenalized=nopen,data=datapred,...)
time1 <- proc.time()-pmt0
print(opt$cv)
print(paste("Computation time for cross-validating main penalty parameter:",time1[3]))
optl <- opt$lambda
arguments$optl <- optl
}




pmt <- proc.time()
nsam <- ncol(highdimdata)
nfeat <- nrow(highdimdata)
XM0 <- t(highdimdata)
response0 <- response

if(is.null(monotone)) monotone <- rep(FALSE,nclass)
if(!selection) cvlnstot <- rep(0,(nclass+1)) else cvlnstot <- rep(0,(nclass+2))
allpreds <- c()
whsam <- 1:nsam


#whsam <- c(11,12,31,35)

responsemin <- response0
pen0 <- penalized(responsemin, penalized = XM0, lambda2 = optl, unpenalized=nopen,data=cbind(XM0,datapred))
nmunpen <- names(pen0@unpenalized)
if(is.element("(Intercept)",nmunpen)) addintercept <- TRUE else addintercept <- FALSE #needed later for matrix operations
#pen0 <- penalized(responsemin, penalized = XM, lambda2 = optl, unpenalized=nopen,data=datapred,...)
if(is.null(innfold)) {nf <- nrow(XM0)} else {
    if(!is.null(optl0)) {
                nf <- innfold
                if(fixedfoldsinn) set.seed(346477)
                } else {nf <- opt$fold}
    }
opt2 <- cvl(responsemin, penalized = XM0,fold=nf, lambda2 = optl,unpenalized=nopen,data=datapred)
nf <- opt2$fold  #folds are fixed for all cvs hereafter
#opt2 <- cvl(responsemin, penalized = XM,fold=nf, lambda2 = optl, unpenalized=nopen,data=datapred,...)
cvln0 <- opt2$cvl





cvlnprev <- cvln0
penprev <- pen0
pen <- pen0
print(cvln0)
XMw0 <- XM0
XMw0prev <- XM0
converged <- FALSE
conv <- rep(FALSE,nclass)
controlbound <-10 #maximal relative deviation of the estimates

almvecall <- rep(1,nfeat) #lambda multipliers for unweighted ridge 
lambdas <- lapply(partitions, function(cla) {
ngroup <- length(cla)
return(rep(1,ngroup))
})
lmvec <- lmvecprev <- array(1,nfeat)
i <- 1

while(!converged & i <= niter){
cl <- 1
if(method=="adaptridge") cl <- nclass
while(cl <= nclass){
convcl <- conv[cl]
if(!convcl){
whgr <- partitions[[cl]]
ngroup1 <- length(whgr)
names(lambdas[[cl]]) <- names(whgr)


coeff <- penprev@penalized
preds <- predict(penprev,XMw0,data=datapred)[1:nsam]  #needed because linear model returns matrix 
   

coeffsq <- coeff^2
if(model=="logistic") {
    Wi <- sqrt(preds*(1-preds)) 
    constlam <- 2
    } else { #model = "linear"
    Wi <- rep(1,length(preds))
    constlam <- 1
    }

if(!is.null(dataunpen)) {
    mm <- model.matrix(nopen,dataunpen) #model.matrix automatically adds the intercept
    XMW <- t(t(cbind(XMw0,10^5*mm)) %*% diag(Wi)) #10^8: large number that 'undoes' the generic penalty 
    } else {
    if(addintercept) XMW <- t(t(cbind(XMw0,rep(10^5,nsam))) %*% diag(Wi)) else XMW <- t(t(XMw0) %*% diag(Wi))
    }
    
SVD <- svd(XMW)

leftmat <- SVD$v %*% diag(1/((SVD$d)^2+constlam*optl)) %*% diag(SVD$d) %*% t(SVD$u)

if(model=="linear"){
vars3 <- VarRes*rowSums(leftmat^2)
} else {
vars3 <- rowSums(leftmat^2)
}

#XMW <- XMW  #checked; coincides with svd solution above
#p1 <- t(XMW) %*% XMW + constlam*optl*diag(x=1,nrow=P)
#inv1 <- solve(a=p1,b=diag(x=1,nrow=P))
#vars1 <- inv1 %*% t(XMW) %*% XMW %*% inv1
#vs <- VarRes*diag(vars1)

if(model=="linear"){
mycoeff2svd <- (leftmat %*% response)^2 #for linear ridge the formula is exact
} else {
respnum <- as.numeric(response)-1
z <- matrix(log(preds/(1-preds))+(respnum - preds)/(preds*(1-preds)),ncol=1) 
mycoeff2svd <- (leftmat %*% z)^2 
}
#
#coefmatrix3 <- (vars2left %*% XMW)^2
#coefmatrix2 <- (inv1 %*% t(XMW) %*% XMW)^2

#mycoeff2 <- (inv1 %*% t(XMW) %*% z)^2 #checked: coincides with svd solution

#cii <- diag(coefmatrix2)
cii2 <- (rowSums(leftmat * t(XMW)))^2  #checked


leftmat <- leftmat/sqrt(vars3)

#lefts <- function(group){
###group <- 1
#ind <- whgr[[group]]
#coefftau2 <- sum(mycoeff2svd[ind]-vars3[ind])
#return(coefftau2)
#}
#leftside <- sapply(1:length(whgr),lefts)

#random:
#whgr <- lapply(1:ngroup1,function(i){ni <- length(whgr[[i]]);return(sample(1:nr,ni))})

lefts2 <- function(group){
##group <- 1
ind <- whgr[[group]]
coefftau2 <- sum(mycoeff2svd[ind]/vars3[ind])-length(ind)
return(coefftau2)
}
leftside <- sapply(1:length(whgr),lefts2)

#lefts2 <- function(group){
###group <- 1
#ind <- whgr[[group]]
#coefftau2 <- sum(coeffsq[ind]/vars3[ind])-length(ind)
#return(coefftau2)
#}
#leftside <- sapply(1:length(whgr),lefts2)
#

#leftmat of dimension p*n
rightmat = XMW  #n*p

#store the rightmat products (dim n*n) to avoid recomputing
rightmats <- lapply(1:length(whgr),function(j){
rightj2 <- rightmat[,whgr[[j]]]
rcp <- rightj2 %*% t(rightj2)
return(rcp)
})


coefmatfast <- t(apply(matrix(1:length(whgr),nrow=length(whgr)),1,  
function(i){
#i<-1;j<-1
lefti2  <- leftmat[whgr[[i]],]
#lefti2 <- matrix(c(1,2,3,4,3,8),nrow=3)
lcp <- t(lefti2) %*% lefti2
#lcp <- leftcrosspr[upper.tri(leftcrosspr)]
ckls <- sapply(1:length(whgr),function(j){
#rightj2 <- rightmat[,whgr[[j]]]
##rightj2 <- matrix(c(1,0,2,4,4,1),nrow=2)
#rcp <- rightj2 %*% t(rightj2)
rcp <- rightmats[[j]]
return(sum(lcp*rcp))
})
return(ckls)
}))


#sum((lefti2 %*% rightj2)^2)

#coefmatFast <- function(L,R){
#lcp <- t(L) %*% L
#rcp <- R %*% t(R)
#return(sum(lcp*rcp))
#}


#Exact
if(method=="exact"){
soltau = solve(sum(coefmatfast),sum(leftside))
sol = solve(coefmatfast,leftside)
low <- soltau/controlbound;up = soltau*controlbound
parinint <- sapply(sol,function(x) min(max(low,x),up))
minopt <- optim(par = parinint,fn = function(pars=c(parinint)) sum(leftside-coefmatfast %*% pars)^2,method="L-BFGS-B",
lower=rep(low,ngroup1),upper=rep(up,ngroup1))
tausqest0 <- minopt$par
tausqest0
}

#controlbound <- 10
#soltau = solve(sum(coefmatfast),sum(leftside))
#sol = solve(coefmatfast,leftside)
#low <- soltau/controlbound;up = soltau*controlbound
#parinint <- sapply(sol,function(x) min(max(low,x),up))
#ir <- isoreg(1/parinint)
#parinint2 <- 1/(ir$yf)
#startup <- parinint2 - c(parinint2[-1],0)
#coefmatfastcs <- t(apply(coefmatfast,1,function(row) rev(cumsum(rev(row)))))
#coefmatfastcs[which(lower.tri(coefmatfastcs))]<-0
#minopt <- optim(par = parinint,fn = function(pars=c(parinint)) sum(leftside-coefmatfastcs %*% pars)^2,method="L-BFGS-B",
#lower=c(rep(0,ngroup1-1),low))
#tausqest0 <- rev(cumsum(rev(minopt$par)))
#tausqest0



#hybrid
#i<-2
if(method=="stable"){
soltau = solve(sum(coefmatfast),sum(leftside))
solhyb <- sapply(1:ngroup1,function(i){
leftsidei <- leftside[i] 
rightsidei <- c(coefmatfast[i,i], sum(coefmatfast[i,-i]*soltau))
soli <- (leftsidei-rightsidei[2])/rightsidei[1]
return(max(min(soltau*controlbound,soli),soltau/controlbound))
})
}

#if(method=="stable2"){
#soltau = solve(sum(coefmatfast),sum(leftside))
#solhyb <- sapply(1:ngroup1,function(i){
##i<-4
#leftsidei <- leftside[i]
#leftsidei2 <- sum(leftside[-i]) 
#rightsidei <- c(coefmatfast[i,i], sum(coefmatfast[i,-i]))
#rightsidei2 <- c(sum(coefmatfast[-i,i]),sum((coefmatfast[,-i])[-i,]))
#leftsi <- c(leftsidei,leftsidei2)
#rightsi <- rbind(rightsidei,rightsidei2)
#solve(rightsi,leftsi)
#return(max(min(soltau*controlbound,soli),soltau/controlbound))
#})
#}

#rightsides <- t(sapply(1:ngroup1,function(i){
#leftsidei <- leftside[i] 
#alli <- c(leftsidei,coefmatfast[i,i], sum(coefmatfast[i,-i]*soltau))
#return(alli)
#}))
#
#

#solhybrec <- c(solhyb[ngroup1])
#solhybprev <- solhyb[ngroup1]
#for(i in (ngroup1-1):1){
#shb <- max(min(solhybprev*3,solhyb[i]),solhybprev/3)
#solhybrec <- c(shb,solhybrec)
#solhybprev <- shb
#}



#groupeq2 <- function(group){
##group <- 1
#ind <- whgr[[group]]
#coefftau2 <- (sum(coeffsq[ind]/vars3[ind])-length(ind))/sum(cii2[ind]/(vars3[ind]))
#return(coefftau2)
#}
#
##groupeq <- function(group){
###group <- 1
##ind <- whgr[[group]]
##coefftau2 <- mean(coeffsq[ind])
##return(coefftau2)
##}
##
##tausqest0 <- sapply(1:length(whgr),groupeq)
#tausqest1 <- sapply(1:length(whgr),groupeq2)
#tausqest1




if(method=="exact") {
tausqest <- tausqest0
print("Exact")
} 
if(method=="stable") {
print("stable")
#tausqest <- solhybrec
tausqest <- solhyb
} 
if(method=="adaptridge") {
print("adaptive ridge")
#tausqest <- solhybrec
} 

if(method=="stable" | method=="exact"){

lambdanoncal <- 1/tausqest

if(monotone[cl]){
weigh = unlist(lapply(whgr,length))
lambdamultnoncal <- pava(lambdanoncal,w=weigh)
} else lambdamultnoncal <- lambdanoncal

#lambdamultnoncal<-1/rev(cumsum(rev(log(1/lambdanoncal))))/(ngroup1:1)

#logtausqest<-log(1/lambdamultnoncal)
#con3 <- sum(sapply(1:length(whgr),function(gr){return(length(whgr[[gr]])*logtausqest[gr])}))/nfeat
#loglambdamult<-  - logtausqest + con3
#lambdamult <- exp(loglambdamult)

tausqest<-1/lambdamultnoncal
con3 <- sum(sapply(1:length(whgr),function(gr){return(length(whgr[[gr]])*tausqest[gr])}))
tausqestcal<-  nfeat/con3*tausqest
lambdamult <- 1/tausqestcal
print(lambdamult)



for(k in 1:length(whgr)){
wh <- whgr[[k]]
XMw0[,wh] <- XMw0[,wh]/sqrt(lambdamult[k])
}
} else {
tausqest <- coeffsq
con3 <- sum(tausqest)
tausqestcal<-  nfeat/con3*tausqest
lambdamult <- 1/tausqestcal
XMw0 <- t(t(XMw0)/sqrt(lambdamult))
}

opt2w <- cvl(responsemin, penalized = XMw0,fold=nf,lambda2=optl,unpenalized=nopen,data=datapred)
#opt2w <- cvl(responsemin, penalized = XMw,fold=nf,lambda2=optl, unpenalized=nopen,data=datapred,...)
#pred <- opt2w$predictions
cvln1 <- opt2w$cvl
print(cvln1)


if((cvln1 - cvlnprev)/abs(cvlnprev) > 0.1/100 | ((cvln1 - cvlnprev)/abs(cvlnprev) >= -0.002 & i==1)) { 
#if((cvln1 - cvlnprev)/abs(cvlnprev) >= cvlmarg) { 
pen <- penalized(responsemin, penalized = XMw0, lambda2 = optl,unpenalized=nopen,data=datapred)
if(method == "exact" | method =="stable"){
for(group in 1:ngroup1){
ind <- whgr[[group]]
lmvec[ind] <- lmvec[ind]*lambdamult[group]
}
lambdas[[cl]] <- lambdas[[cl]]*lambdamult
}
cvlnprev <- cvln1
penprev <- pen
XMw0prev <- XMw0
print(paste("Partition",nmp0[cl],"improved results"))
} else {
if(method=="stable") print(paste("Partition",nmp0[cl],"CONVERGED after",i,"iterations")) 
  else print(paste("Partition",nmp0[cl],"did not improve results"))
  
conv[cl] <- TRUE
cvln1 <- cvlnprev
XMw0 <- XMw0prev
pen <- penprev
}
}
cl <- cl+1
}
if(sum(conv)==nclass){
    converged <- TRUE
if(method=="stable") print(paste("All partitions CONVERGED after",i,"iterations"))
}
i <- i+1
}



#if(reCV){
#print("Finding NEW lambda")
#if(fixedfoldsinn) set.seed(346477)
#optnew <- optL2(response, penalized = XMw0,fold=foldinit,trace=T,unpenalized=nopen,data=datapred)
##opt <- optL2(response, penalized = t(highdimdata),fold=foldinit,trace=T,unpenalized=nopen,data=datapred,...)
#time1 <- proc.time()-pmt0
#print(optnew$cv)
#print(paste("Computation time for re-cross-validating main penalty parameter:",time1[3]))
#optlnew <- optnew$lambda
#pen <- penalized(responsemin, penalized = XMw0, lambda2 = optlnew, unpenalized=nopen,data=datapred)
#}

if(niter==0) {pen <- pen0;XMw0<-XM0;cvln1<-cvln0;soltau <- NULL}
pred0 <- predict(pen0,XM0,unpenalized=nopen,data=datapred)[1:nsam]  
predw <- predict(pen,XMw0,unpenalized=nopen,data=datapred)[1:nsam]
predshere <- cbind(pred0,predw)
cvlnssam <- c(cvln0,cvln1)
lmvecall <- lmvec
almvecall <- cbind(almvecall,lmvecall)
predobj <- c(pen0,pen)


allpreds <- predshere
whichsel <- NULL
betassel <- NULL

npr <- length(predobj)
pred2 <-predobj[[npr]]
lambs <- lmvecall
oldbeta <- pred2@penalized
newbeta <- oldbeta/sqrt(lambs)
time2 <- proc.time()-pmt
print(paste("Computation time for adaptive weigthing:",time2[3]))




if(selection){
pmt <- proc.time()
ord <- order(abs(newbeta),decreasing=T)
sequ <- seq(0,maxsel,by=stepsel)
cvlsels <- c()
for(nsel in sequ){
whsel <- ord[1:nsel]
datwsel <- XMw0[,whsel]
pensel <- penalized(responsemin,datwsel,lambda2=optl, unpenalized=nopen,data=datapred)
#pensel <- penalized(responsemin,datwsel,lambda2=optl, unpenalized=nopen,data=datapred,...)
optsel <- cvl(responsemin,datwsel,fold=nf,lambda2=optl,unpenalized=nopen,data=datapred)
#optsel <- cvl(responsemin,datwsel,fold=nf,lambda2=optl, unpenalized=nopen,data=datapred,...)
cvlsel <- optsel$cvl
cvlsels <- c(cvlsels,cvlsel)
}
whbest <- which.max(cvlsels)
cvmax <- cvlsels[whbest]
nsel2 <- sequ[(which((cvlsels - cvmax) >= cvlmarg/100*(cvmax)))[1]]
whsel2 <- ord[1:nsel2]
whichsel <- whsel2
betassel <- newbeta[whsel2]
datwsel2 <- XMw0[,whsel2]
pensel2 <- penalized(responsemin,datwsel2,lambda2=optl,unpenalized=nopen,data=datapred)
optsel2 <- cvl(responsemin,datwsel2,fold=nf,lambda2=optl,unpenalized=nopen,data=datapred)
#pensel2 <- penalized(responsemin,datwsel2,lambda2=optl, unpenalized=nopen,data=datapred,...)
#optsel2 <- cvl(responsemin,datwsel2,fold=nf,lambda2=optl, unpenalized=nopen,data=datapred,...)
cvlsel2 <- optsel2$cvl
#predsel <- predict(pensel2,penalized=XMw0[,whsel2,drop=FALSE], unpenalized=nopen)
#predsel <- predict(pensel2,penalized=XMw0[,whsel2,drop=FALSE], unpenalized=nopen)
predsel <- predict(pensel2,penalized=XMw0[,whsel2,drop=FALSE],unpenalized=nopen,data=datapred)[1:nsam]

allpreds <- cbind(allpreds,predsel)
predobj <- c(predobj,pensel2)
cvlnssam <- c(cvlnssam,cvlsel2)
print(paste("Number of selected markers by forward selection:",nsel2))
time3 <- proc.time()-pmt
print(paste("Computation time for feature selection by forward selection:",time3[3]))
}
if(selection) cvlssel <- data.frame(nsel=sequ,cvl=cvlsels) else cvlssel <- NULL

cvlnstot <- cvlnssam
mat <- cbind(response0,allpreds) 
#if(printp){
#print("True:")
#print(response0)
#print(round(allpreds,3))
#}

#NEW
reslasso <- NULL
if(comparelasso){
print("Starting lasso")
    if(is.null(optllasso)){
    print("Finding lambda for lasso regression")
    opt <- optL1(response, penalized = t(highdimdata),fold=nf,trace=T,unpenalized=nopen,data=datapred)
    print(opt$cv)
    optllasso <- opt$lambda
    arguments$optllasso <- optllasso
    cvllasso <- opt$cv
    } else {
    cvliklasso <- cvl(response,penalized = t(highdimdata),lambda1 = optllasso,fold=nf,unpenalized=nopen,data=datapred)
#optsel <- cvl(responsemin,datwsel,fold=nf,lambda2=optl, unpenalized=nopen,data=datapred,...)
    cvllasso <- cvliklasso$cvl
    }
cvlnstot <- c(cvlnstot,cvllasso)
penlasso <- penalized(response, penalized = t(highdimdata), lambda1 = optllasso, unpenalized=nopen,data=cbind(XM0,datapred))
whichlasso <- which(penlasso@penalized != 0)
betaslasso <- penlasso@penalized[whichlasso]
predobj <- c(predobj,penlasso)
reslasso <- list(cvllasso=cvllasso,whichlasso=whichlasso,betaslasso=betaslasso)
}

if(compareunpenal){
if(model == "logistic") famglm <- "binomial" else famglm <- "gaussian"
print("Starting unpenalized glm")
form <- formula(paste("response","~",as.character(unpenal)[2]))
modelglm <- glm(form,family=famglm,data=dataunpen)
predobj <- c(predobj,list(modelglm))
}

printlam <- function(lambs) {if(length(lambs)<=10) return(lambs) else return(summary(lambs))}

suml <- lapply(lambdas,printlam)
print("Final lambda multipliers (summary):")
print(suml)
print(paste("CVLs",cvlnstot))
timetot <- proc.time()-pmt0
print(paste("Total computation time:",timetot[3]))
names(predobj) <- nmp
colnames(almvecall) <- nmpweight
if(savepredobj=="last") {
predobj <- predobj[length(predobj)]
almvecall <- matrix(lmvecall,ncol=1)
}
if(savepredobj=="none") predobj <- NULL
#colnames(allpreds) <- nmp
#NEW
return(list(true=response,cvls = cvlnstot,lambdamults = lambdas, optl=optl, lambdamultvec = almvecall, predobj=predobj,betas=newbeta, 
whichsel = whichsel,cvlssel = cvlssel,reslasso=reslasso,arguments=arguments))
}
