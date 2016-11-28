PartitionsSelection = function(highdimdata, response, partitions, 
                               monotoneFunctions, optl=NULL, innfold=NULL){ 


  if(length(partitions)==1){
    print("There is only one partition available. There is no need to do partition selection.")
    return(NULL)
  }
  
  if(length(partitions)!=length(monotoneFunctions)){
    print("ERROR: length 'monotone' unequal to length 'partitions' ")
    return(NULL)
  }

  PartitionsNames = names(partitions)
  Partitions.i = partitions
  nPar.i = length(Partitions.i)
  monotone.i = monotoneFunctions
  parIn = parIn2 = im = ordPar = c()
  lambdamults = vector("list",nPar.i)
  names(lambdamults) = PartitionsNames
  cvlmarg = 1
  cvlmodel = c()

  while(nPar.i != 0){
    cvl.i = c()
    for(i in 1:nPar.i){
      idTemp = c(im,which(PartitionsNames==names(Partitions.i)[i]))
      if(is.null(optl)==TRUE){
      grMod.i =  grridge(highdimdata,response,partitions=partitions[idTemp],
                          monotone=monotoneFunctions[idTemp],innfold=innfold) 
      optl = grMod.i$optl
      }else{
      grMod.i =  grridge(highdimdata,response,partitions=partitions[idTemp],
                          monotone=monotoneFunctions[idTemp],optl= optl, innfold=innfold) 
      }
      cvl.i[i] =  grMod.i$cvls[[length(grMod.i$cvls)]]
    }
    names(cvl.i) = names(Partitions.i)
    idMax = which(cvl.i==max(cvl.i))  
    if(nPar.i == length(monotoneFunctions)){cvlmodel=-300}else{cvlmodel=cvlmodel}
    if((cvlmodel[length(cvlmodel)] - cvl.i[idMax]) >= (cvlmarg/100 * cvl.i[idMax])){break}  
    idNext = which(PartitionsNames==names(cvl.i)[idMax])
    ordPar = c(ordPar,idNext)
    parIn = c(parIn,PartitionsNames[idNext])  
    parIn2 = c(parIn2, paste(parIn2,PartitionsNames[idNext],sep=";"))
    is =  intersect(PartitionsNames,parIn)
    im = match(parIn,PartitionsNames)
    
    grMod = grridge(highdimdata,response,partitions=partitions[im],monotone=monotoneFunctions[im],
                    optl= optl,innfold=innfold) 
    cvlmodel = c(cvlmodel,grMod$cvls[[length(grMod$cvls)]])
    
    Partitions.i = Partitions.i[-idMax]
    monotone.i = monotone.i[-idMax]
    nPar.i=length(Partitions.i)  
  }
  
  resMat = data.frame(Partitions_in_theModel=parIn2, cvl=cvlmodel[-1], 
                      gainCvl= if(length(cvlmodel)==2){cvlmodel[-1]}else{diff(cvlmodel[-1])})
  print(resMat)
  
  return(list(ordPar=ordPar,optl=optl))
}

