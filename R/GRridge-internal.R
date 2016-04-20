.preddich <-
function(probs,cutoff){
        sapply(probs,function(x) if(x<=cutoff) 0 else 1)
    }
    
.rocaux <-
function(probs,true,cutoff){  
    predlab <- .preddich(probs,cutoff)
    if(class(true)=="factor") true <- as.numeric(true)-1
    NN <- length(true[true==0])
    NP <- length(true[true==1])
    predtrue <- cbind(predlab,true)
    FP <- length(intersect(which(true==0),which(predlab==1)))
    FPR <- FP/NN
    TP <- length(intersect(which(true==1),which(predlab==1)))
    TPR <- TP/NP
    Acc <- (length(intersect(which(true==0),which(predlab==0)))+length(intersect(which(true==1),which(predlab==1))))/length(true)
    roc <- c(FPR=FPR,TPR=TPR,Accuracy=Acc)
    return(roc)
}

