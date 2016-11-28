mergeGroups = function(highdimdata=highdimdata, initGroups=initGroups,maxGroups=maxGroups,
                       methodDistance="manhattan", methodClust="complete"){
   
  if (is.null(maxGroups) | is.null(highdimdata) | is.null(initGroups)) {
    print("Please check these following arguments: highdimdata, initiGroups and maxGroups.")
    return(NULL)
  }
  
  if (maxGroups<2) {
    print("The number of new groups is too few. Please consider higher number of groups.")
    return(NULL)
  }
  
  if (length(initGroups)<=maxGroups) {
    print("The number of new groups cannot be greater than or equal to the number of initial groups.")
    return(NULL)
  }  
  
  clustNames = names(initGroups)
  nClust = length(clustNames)
  nSamp = ncol(highdimdata)
  nVar = nrow(highdimdata)
    
  eigenClust = matrix(NA,nSamp,(nClust-1))
  for(i in 1:(nClust-1)){
    idClust = initGroups[[i]]                          
    datClust = t(highdimdata[idClust,])               
    covDat = cov(datClust)                             
    eigClust = eigen(covDat)                           
    eigVec1 = eigen(covDat)$vectors[1,]                
    eigenClust[,i] = t(t(eigVec1) %*% t(datClust))     
  }
  colnames(eigenClust) = clustNames[-nClust]
  
  dist_eigenClust = dist(t(eigenClust), method = methodDistance)
  clust = hclust(dist_eigenClust, method=methodClust) 
  newgroups = cutree(clust, k=maxGroups-1)             
  newgroups = as.numeric(c(newgroups,maxGroups))       
  
  GroupNew = newClustmembers = list()
  clustmembers = unique(newgroups)
  for(i in clustmembers){
    idClust.i = which(newgroups==i) 
    newClustmembers[[i]] = clustNames[idClust.i]
    GroupNew[[i]] = unique(as.numeric(unlist(initGroups[idClust.i])))  
  }
  names(GroupNew) = names(newClustmembers) = clustmembers
  return(list(newGroups=GroupNew,newGroupMembers=newClustmembers))
}
