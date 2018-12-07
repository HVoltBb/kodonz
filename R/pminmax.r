#' %MinMax of a coding DNA sequence
#'
#' Calculates the %MinMax (Clarke and Clark, 2008) of a DNA sequence coding a protein of interest.
#'
#' @param x a char array of the codons for which %MinMax will be calculated
#' @param z the width of the sliding window.
#' @param cut a named numeric array of the codon usage table of the host species. Each entry is the number of codons of genes from the reference species, and the name is the case insensitive three letter base char string of the codon.
#' @param y a number denoting the codon table to be used. When \emph{cut} is specified, \emph{y} also needs to be specified. 0 = standard codon table (default), 2 = vertibrate mitchodrial, 3 = ?, 4 = ?. When an option other than those mentioned above is provided, the standard codon table will be used
#' @param spp a character string denoting the host species name on which the codon usage table is based. Four reference tables are available: 'ecoli', E. coli; 'insect', insect; 'yeast', yeast; 'celegans', C. elegans. When specified, this option overides \emph{cut}. These build-in codon usage tables are taken from https://www.genscript.com/tools/codon-frequency-table.
#' @examples
#'
#' @export

pminmax <- function(x, z, cut, y=0, spp=''){
  switch(as.character(spp),
         'ecoli' = {
           cft = cft_ecoli[,5]
           names(cft) = cft_ecoli[,1]
           y=0
         },
         'insect' = {
           cft = cft_insect[,5]
           names(cft) = cft_insect[,1]
           y=0
         },
         'yeast' = {
           cft = cft_yeast[,5]
           names(cft) = cft_yeast[,1]
           y=0
         },
         'celegans' = {
           cft = cft_celegans[,5]
           names(cft) = cft_celegans[,1]
           y=0
         },
         {
           cft = cut
         })
  cTable = arg_check(x, 'character', y)

  temp = table(cTable[,2])[-1]
  xavg = rep(0, 64)
  xmax = rep(0, 64)
  xmin = rep(0, 64)
  aa = names(temp)
  names(xavg) = cTable[,1]
  names(xmax) = cTable[,1]
  names(xmin) = cTable[,1]
  for(i in 1:20){
    codonlist = cft[cTable[cTable[,2]==aa[i],1]]
    xavg[cTable[,2]==aa[i]] = mean(codonlist)
    xmax[cTable[,2]==aa[i]] = max(codonlist)
    xmin[cTable[,2]==aa[i]] = min(codonlist)
  }

  xn = length(x)
  if(xn < z) stop('The width of the sliding window is too large with respect to the DNA sequence provided.')
  ans = rep(0, xn-z+1)
  for(i in 1:(xn-z+1)){
    span = (1:z)+i-1
    num = sum(cft[x[span]]-xavg[x[span]])
    if(num>0){
      den = sum(xmax[x[span]]-xavg[x[span]])
    } else {
      den = sum(xavg[x[span]]-xmin[x[span]])
    }
    ans[i] = num/den
  }
  return(ans)
}
