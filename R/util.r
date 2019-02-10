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
      warning('An animo acid is missing, and NAs will be produced for this animo acid')
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

bases = strsplit("ACGTUWSMKRYBDHVNZ", '')[[1]]
