roc <-
function(probs,true,cutoffs){  
    rocs <- sapply(cutoffs, .rocaux, probs=probs, true=true)
    return(rocs)
}
