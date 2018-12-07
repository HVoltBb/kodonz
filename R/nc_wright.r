#' Effective number of codons based on Wright's formula
#'

nc_wright <- function(x, cTable){
  w = table(x)[cTable[,1]]
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

  if(NaN %in% slots) warning("At least one degenerate class of AA is missing, NaN will be produced. Try Fuglsang's formula to get an estimate.")
  return(sum(table(temp[-1])/slots))
}
