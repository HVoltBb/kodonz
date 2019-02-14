#' \%MinMax of a coding DNA sequence
#'
#' Calculates the \%MinMax (Clarke and Clark, 2008) of a DNA sequence coding a protein of interest.
#'
#' @param x a KZsqns object of the codons for which \%MinMax will be calculated; when a list of such objects are provided, \%MinMax will be calculated for the first object only
#' @param z the width of the sliding window. Default 10
#' @param cut a named numeric array of the codon usage table of the host species. Each entry is the number of codons of genes from the reference species, and the name is the case insensitive three letter base char string of the codon. Optional if \emph{spp} is specified
#' @param y a number denoting the codon table to be used. When \emph{cut} is specified, {y} also needs to be specified. 0 = standard codon table (default), 2 = vertibrate mitchodrial (See dna2aa for additional options). When an option other than those mentioned above is provided, the standard codon table will be used
#' @param spp a character string denoting the host species name on which the codon usage table is based. Sixteen build-in reference tables are available: 'ecoli', E. coli (standard code), default option; 'insect', insect (standard code); 'yeast', yeast (standard code); 'celegans', C. elegans (standard code); 'drosophila', Drosophila melanogaster (standard code); 'human', human (standard code); 'mouse', mouse (standard code); 'rat', rat (standard code); 'pig', pig (standard code); 'pichia', Pichia pastoris (standard code); 'arabidopsis', Arabidopsis thaliana (standard code); 'streptomyces', Streptomyces (standard code); 'zea', Zea mays (standard code); 'nicotiana', Nicotiana tabacum (standard code); 'saccharomyces', Saccharomyces cerevisiae (standard code); 'cricetulus', Cricetulus griseus (standard code). When specified, this option overides {cut}.
#' @return a numerical array of \%MinMax values
#' @details Frequently used codon frequency table in 16 expression host organisms were built-in based on Genscript [1]. In addition, when your host organism is not one of those built-in, you can supply a custom array of codon usage table and the corresponding codon table can be supplied using {cut} and {y} to calculate \%MinMax statistic for additional host organisms.
#' @references [1] https://www.genscript.com/tools/codon-frequency-table [retreaved 2/4/2019]
#' @examples
#'
#' @export

pminmax <- function(x, z=10, cut='n', y=0, spp='ecoli'){
  if(is.list(x)){
    x = x[[1]]
    cat("Expecting a KZsqns object and a list is given\nOnly the first element of the list will be analyzed\n")
  }
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
         'drosophila' = {
           cft = cft_drosophila[,5]
           names(cft) = cft_drosophila[,1]
           y=0
         },
         'human' = {
           cft = cft_human[,5]
           names(cft) = cft_human[,1]
           y=0
         },
         'mouse' = {
           cft = cft_mouse[,5]
           names(cft) = cft_mouse[,1]
           y=0
         },
         'rat' = {
           cft = cft_rat[,5]
           names(cft) = cft_rat[,1]
           y=0
         },
         'pig' = {
           cft = cft_pig[,5]
           names(cft) = cft_pig[,1]
           y=0
         },
         'pichia' = {
           cft = cft_pichia[,5]
           names(cft) = cft_pichia[,1]
           y=0
         },
         'arabidopsis' = {
           cft = cft_arabidopsis[,5]
           names(cft) = cft_arabidopsis[,1]
           y=0
         },
         'streptomyces' = {
           cft = cft_streptomyces[,5]
           names(cft) = cft_streptomyces[,1]
           y=0
         },
         'zea' = {
           cft = cft_zea[,5]
           names(cft) = cft_zea[,1]
           y=0
         },
         'nicotiana' = {
           cft = cft_nicotiana[,5]
           names(cft) = cft_nicotiana[,1]
           y=0
         },
         'saccharomyces' = {
           cft = cft_saccharomyces[,5]
           names(cft) = cft_saccharomyces[,1]
           y=0
         },
         'cricetulus' = {
           cft = cft_cricetulus[,5]
           names(cft) = cft_cricetulus[,1]
           y=0
         },
         {
           cft = cut
         })
  cTable = arg_check(x, 'KZsqns', y)
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
