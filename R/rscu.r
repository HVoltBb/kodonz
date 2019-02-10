#' Relative synonymous codon usage
#'
#' Calculates the relative synomynous codon usage (RSCU) for each codon.
#'
#' @details The RSCU is calculated for each animo acid coding codon present in the provided codon sequence. When any animo acid is missing from the sequence, each coding codon for that animo acid will have NA returned as their RSCU. The RSCU of stop codons also have NAs in their position as the default setting. To include RSCUs for the stop codons in the returned array, set \emph{stop} as TRUE.
#' @param x a char array of codons
#' @param y an optional number denoting the codon table to be used. 0 = standard codon table (default), 2 = vertibrate mitchodrial (See dna2aa for additional options). When an option other than those mentioned above is provided, the standard codon table will be used
#' @param stop a optional boolean denoting whether the RSCU values for the stop codons should be returned. FALSE: NAs will be returned (default); TRUE: the RSCUs will be returned
#' @return a numeric array of RSCUs for each animo acid
#' @examples
#'
#' @export

rscu <- function(x, y=0, stop = FALSE){
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
      warning('An animo acid is missing, and NAs will be produced for this animo acid')
    } else {
      degv[dx] = w[dx]/sum(w[dx])*temp[j]
    }
  }
  names(degv) = cTable[,1]
  if(!stop) degv[which(cTable[,2]=='*')] <- NA
  return(degv)
}
