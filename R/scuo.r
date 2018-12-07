#' Synonymous codon usage order
#' 
#' Calculates the synonymous codon usage order of a codon sequence based on a codon table.
#' 
#' @param x a char sequence of codons
#' @param y a number to assign the type of codon table to be used for the translation. 0 = standard codon table, 2 = vertibrate mitchodrial, . If no number is specified or an option other than those provided above is specified, there will be a warning and the standard codon table will be used
#' @examples 
#' 
#' @export

scuo = function(x, y=0){
  cTable = arg_check(x, 'character', y)
  
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
  scuo_aa = (hmax-haa)/hmax
  aa = table(dna2aa(x, 1, y))[names(degv)]
  return(crossprod(scuo_aa, aa)/sum(aa))
}