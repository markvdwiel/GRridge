roc <-
function(probs,true,cutoffs){  
    if(cutoffs[2] >= cutoffs[1]){
      print("ERROR: cut-offs needs to be a DECREASING vector")
      return(NULL)
    } else {
    rocs <- sapply(cutoffs, .rocaux, probs=probs, true=true)
    return(rocs)
    }
}
