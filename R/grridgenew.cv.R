grridge.cv <-
function(grr,highdimdata,response,outerfold=length(response),fixedfolds=TRUE){
#highdimdata<-datcenFarkas;response<-respFarkas;grr <- grFarkas;outerfold=10;fixedfolds<-T
#highdimdata<-datcenVerlaat;response<-respVerlaat;grr <- grV;outerfold=5;fixedfolds<-T
#highdimdata<-betasfcen;response<-fl;grr<-grwinacomploffset;outerfold=length(fl);fixedfolds<-T
#highdimdata<-ncftrcen;response<-resp;grr<-grBarbara;outerfold=10;fixedfolds<-T
#highdimdata<-mirnormcen;response<-resp;grr<-grMaarten;outerfold=10;fixedfolds<-T



nsam <- ncol(highdimdata)
if(fixedfolds) set.seed(3534)
rand<-sample(1:nsam)
grs1 <- floor(nsam/outerfold)
grs2 <- grs1+1
ngr1 <- outerfold*grs2 - nsam
folds <- lapply(1:outerfold,function(xg) {
if(xg <= ngr1) els <- rand[(1+(xg-1)*grs1):(xg*grs1)] else els <- rand[(ngr1*grs1 + 1+(xg-ngr1-1)*grs2):(ngr1*grs1 + (xg-ngr1)*grs2)]
return(els)
}
)

preds <- c()
responsesam <- c()
for(k in 1:length(folds)){
#k<-1
print(paste("Fold nr:",k))
samout <- folds[[k]]
nout <- length(samout)
highdimdatamin <- highdimdata[,-samout]
responsemin <- response[-samout]
arg <- grr$arguments
dataunpenmin <- arg$dataunpen[-samout,,drop=FALSE] 
offsarg  <- arg$offset
if(length(offsarg)>1) offsmin <- offsarg[-samout] else offsmin <- offsarg
#NEW
grmin <- grridge(highdimdatamin,responsemin,partitions=arg$partitions,unpenal=arg$unpenal, offset=offsmin, 
method=arg$method,niter=arg$niter, monotone=arg$monotone, optl=arg$optl, innfold=arg$innfold,
fixedfoldsinn=arg$fixedfoldsinn, selection=arg$selection,maxsel=arg$maxsel,stepsel=arg$stepsel,cvlmarg=arg$cvlmarg, 
dataunpen=dataunpenmin,savepredobj=arg$savepredobj,ord=arg$ord,
comparelasso = arg$comparelasso, optllasso=arg$optllasso, compareunpenal = arg$compareunpenal)
penobj <- grmin$predobj
dataunpen <- arg$dataunpen
Xsam <- t(highdimdata[,samout,drop=FALSE])
dataunpensam <- dataunpen[samout,,drop=FALSE]
responsesam <- c(responsesam,samout)
if(is.null(dataunpensam)) dataunpensam <- data.frame(fake=rep(NA,length(samout)))
predellall <- c()
if(is.null(penobj)) {
cat("No prediction objection available. Run grridge using either savepredobj=\"last\" or savepredobj=\"all\"\n")
return(NULL)
}

unpenal<-arg$unpenal
#offsarg <- log(1/9);unpenal= ~0;nout=1
if(!is.null(offsarg)){
noffs <- length(offsarg)
offsargs <-"c("
if(nout>1){
    if(noffs==1) {for(i in 1:(nout-1)) offsargs <- paste(offsargs,offsarg,",",sep="")} else {
    for(i in 1:(nout-1)) offsargs <- paste(offsargs,offsarg[i],",",sep="")}
}
if(noffs==1) offsargs <- paste(offsargs,offsarg,")",sep="")  else offsargs <- paste(offsargs,offsarg[nout],")",sep="")
if((unpenal != ~0) & (unpenal != ~1)){
unpenal <- formula(paste(deparse(unpenal),"+ offset(",offsargs,")",sep=""))
} else {
unpenal <- formula(paste("~","offset(",offsargs,")",sep=""))
}
}

nmp <- names(penobj)
    #NEW
    if(arg$selection) npreds <- length(penobj) -1 else npreds <- length(penobj)
    if(arg$comparelasso) npreds <- npreds-1
    if(arg$compareunpenal) npreds <- npreds-1
    if(npreds>0){
        for(ell in 1:npreds){
        #ell<-1
        predobj <- penobj[[ell]]
        lmvecall <-  grmin$lambdamultvec[,ell]
        Xsamw <- t(t(Xsam)/sqrt(lmvecall))
        predell <- predict(predobj,Xsamw,unpenalized=unpenal,data=dataunpensam)[1:nout]
        predellall <- cbind(predellall,predell)
        }
    }
if(arg$selection){
    predobj <- penobj[[npreds+1]]
    whsel <- grmin$whichsel
    nc <- ncol(grmin$lambdamultvec)
    lmvecall <-  grmin$lambdamultvec[,nc]
    Xsamw <- t(t(Xsam)/sqrt(lmvecall))
    Xsamw <- Xsamw[,whsel,drop=FALSE]
    predell <- predict(predobj,Xsamw,unpenalized=unpenal,data=dataunpensam)[1:nout]
    predellall <- cbind(predellall,predell)
}
#NEW
if(arg$comparelasso){
    if(arg$selection) take <- npreds+2 else take <- npreds + 1
    predobj <- penobj[[take]]
    predell <- predict(predobj,Xsam,unpenalized=unpenal,data=dataunpensam)[1:nout]
    predellall <- cbind(predellall,predell)
}

if(arg$compareunpenal){
    nadd <- arg$selection + arg$comparelasso
    take <- npreds+nadd + 1
    predobj <- penobj[[take]]
    predell <- predict(predobj,data=dataunpensam,type="response")[1:nout]
    predellall <- cbind(predellall,predell)
}

    print(paste("Sample(s) left out:",samout))
    print("True:")
    print(response[samout])
    print("Prediction(s):")
    colnames(predellall)<-nmp
    print(data.frame(response=response[samout],predellall))
preds <- rbind(preds,predellall)
}
sams <- unlist(folds)
od <- order(sams)
mat <- data.frame(response[responsesam],preds)[od,]
colnames(mat) <- c("TrueResponse",nmp)
return(mat)
}
