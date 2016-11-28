auc <-
function(rocout){
    TPRthr <- rocout[2,]
    FPRthr <- rocout[1,]
    FPRdif <- FPRthr - c(0, FPRthr[-length(FPRthr)])
    TPRinterpolated <- (TPRthr + c(0, TPRthr[-length(TPRthr)]))/2
    AUC <- TPRinterpolated %*% FPRdif
    return(AUC)
}
