matchGeneSets <- function(GeneIds,GeneSets,minlen=25,remain=TRUE){
  #seq <- read.table(file="C:\\ExternData\\NeuroBlastoma\\GSE49711_SEQC_NB_MAV_G_log2.20121127.txt",header=T,sep="\t");
  #GeneSets <- oncocursym; remain=FALSE;minlen=25; GeneIds <- seq[,1]
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
  
  #the modified part
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


