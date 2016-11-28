matchGeneSets <- function(GeneIds,GeneSets,minlen=25,remain=TRUE){
  sets <- lapply(GeneSets, function(gs) {
    gs2 <- gs@geneIds
    ma <- match(gs2, GeneIds)
    ma <- ma[!is.na(ma)]
    el <- length(ma)
    if (el < minlen) 
      return(NULL)
    else return(ma)
    return(ma)
  })
  names(sets) <- names(GeneSets)
  whnull <- which(sapply(sets, is.null))
  if(length(whnull)!=0) {sets2 <- sets[-whnull]} else {sets2 <- sets}
  if (remain) {
    ngene <- length(GeneIds)
    remainder <- setdiff(1:ngene, unique(unlist(sets)))
    if (length(remainder) >= minlen) {
      sets2 <- c(sets2, list(Remainder = remainder))
    }
  }
  return(sets2)
}


