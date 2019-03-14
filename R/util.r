ndc <- function(sqns, y, stop){
  cTable = arg_check(sqns, 'KZsqns', y)
  temp1 = table(cTable[,2])
  temp2 = names(temp1[temp1==1])
  bin = numeric()
  if(length(temp2)!=0){
    bin = sapply(temp2,function(x){which(cTable[,2]==x)})
  }
  if(!stop){
    bin = c(bin, which(cTable[,2]=='*'))
  }
  if(length(bin)!=0){
    x = table(sqns)[-bin]
    sqns = unlist(sapply(1:length(x), function(i){rep(names(x[i]), x[i])}))
    attr(sqns, 'class') <- 'KZsqns'
  }
  return(sqns)
}

scu <- function(tx, cTable, stop = FALSE){

  if(method[1] == "sharp"){
    tx[tx==0] <- 0.5
  }
  w = tx[cTable[,1]]
  names(w) <- cTable[,1]
  absentAA = which(is.na(w))
  if(length(absentAA)!=0){
    w[absentAA] = 0
  }

  temp = table(cTable[,2])
  degv = numeric(0)
  for(j in names(temp)){
    dx = which(cTable[,2]==j)
    if(sum(w[dx])==0){
      degv[dx] = NA
      cat('An animo acid is missing, and NAs will be produced for this animo acid\n')
    } else {
      degv[dx] = w[dx]/max(w[dx])
    }
  }
  names(degv) = cTable[,1]
  if(!stop) degv[which(cTable[,2]=='*')] <- NA
  if(method[1]=="bulmer"){
    degv[degv==0] <- 0.01
  }
  return(degv)
}

rscu_ <- function(x, y, stop, supplementory.info){
  cTable = arg_check(x, 'KZsqns', y)

  w = table(x)[cTable[,1]]
  names(w) <- cTable[,1]
  absentAA = which(is.na(w))
  if(length(absentAA)!=0){
    w[absentAA] = 0
  }

  temp = table(cTable[,2])
  degv = numeric(0)
  for(j in names(temp)){
    dx = which(cTable[,2]==j)
    if(sum(w[dx])==0){
      degv[dx] = NA
      cat('An animo acid is missing, and NAs will be produced for this animo acid\n')
    } else {
      degv[dx] = w[dx]/sum(w[dx])*temp[j]
    }
  }
  names(degv) = cTable[,1]
  if(!stop) degv[which(cTable[,2]=='*')] <- NA
  if (supplementory.info){
    return(list(data=degv, top5 = names(degv[order(degv, decreasing = TRUE)[1:5]]), bottom5 = names(degv[order(degv)[1:5]])))
  } else {
    return(degv)
  }
}

nc_fuglsang4m <- function(tx, cTable){
  w = tx[cTable[,1]]
  names(w) <- cTable[,1]
  absentAA = which(is.na(w))
  if(length(absentAA)!=0){
    w[absentAA] = 0
  }

  temp = table(cTable[,2])
  degv = numeric(0)
  for(j in names(temp[-1])){
    xd = which(cTable[,2]==j)
    if(sum(w[xd])==0){
      degv[j] = NA
    } else {
      degv[j] = (sum(w[xd])*sum((w[xd]/sum(w[xd]))^2)-1)/(sum(w[xd])-1)
    }
  }
  du = mean((1/degv)/(temp[-1]), na.rm=T)
  for(j in which(is.na(degv))){
    degv[j] = du*(temp[names(degv[j])] - 1) + 1
  }

  und = names(table((temp[-1])))
  slots = rep(0, length(und))
  dx = 1
  for(j in und){
    if(j==1){
      slots[dx] = 1
    } else {
      slots[dx] = mean(degv[which(names(degv) %in% names(which(temp==j)))])
    }
    dx=dx+1
  }
  return(sum(table(temp[-1])/slots))
}

nc_wright <- function(tx, cTable){
  w = tx[cTable[,1]]
  names(w) <- cTable[,1]
  absentAA = which(is.na(w))
  if(length(absentAA)!=0){
    w[absentAA] = 0
  }

  temp = table(cTable[,2])
  degv = numeric(0)
  for(j in names(temp[-1])){
    xd = which(cTable[,2]==j)
    if(sum(w[xd])==0){
      degv[j] = NA
    } else {
      degv[j] = (sum(w[xd])*sum((w[xd]/sum(w[xd]))^2)-1)/(sum(w[xd])-1)
    }
  }
  und = names(table((temp[-1])))
  slots = rep(0, length(und))
  dx = 1
  for(j in und){
    if(j==1){
      slots[dx] = 1
    } else {
      slots[dx] = mean(degv[which(names(degv) %in% names(which(temp==j)))] ,na.rm=TRUE)
    }
    dx=dx+1
  }

  if(NaN %in% slots) cat("At least one degenerate class of AA is missing, NaN will be produced. Try Fuglsang's formula to get an estimate.\n")
  return(sum(table(temp[-1])/slots))
}

dna2aa_ <- function(dnaSequence, aaNoLetters, codonTable){
  cTable = arg_check(dnaSequence, 'KZsqns', codonTable)
  if(aaNoLetters==1){
    cTable = setNames(cTable[,2],cTable[,1])
  } else {
    cTable = setNames(cTable[,3],cTable[,1])
  }
  return(unname(cTable[dnaSequence]))
}

scuo_ = function(x, y){
  cTable = arg_check(x, 'KZsqns', y)

  temp = table(x)[cTable[,1]]
  exclude = which(cTable[,2]=='*')
  degv = table(cTable[,2])[-1]
  degv = degv[degv!=1]
  haa = numeric(length(degv))
  hmax = haa
  for(i in 1:length(degv)){
    tx = temp[cTable[,2] == names(degv[i])]
    haa[i] = - crossprod(tx/(sum(tx)), log(tx/(sum(tx))))
    hmax[i] = -log(1/length(tx))
  }
  haa[is.na(haa)] = 0
  scuo_aa = (hmax-haa)/hmax
  aa = table(dna2aa(x, 1, y))[names(degv)]
  aa[is.na(aa)] = 0
  return(crossprod(scuo_aa, aa)/sum(aa))
}

