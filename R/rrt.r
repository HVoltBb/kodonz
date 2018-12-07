#' Random reverse translation
#' 
#' Randomly reverse translate an animo acid sequence into a DNA sequence according to a codon usage table and a codon table
#' 
#' @param x a char array of animo acids (e.g., either 'G' or 'Gly' is acceptable but do not mix one letter symbols with three letter symbols in the same array) or codons (e.g., 'GGG'). The first element of the array will be used to determine the type of the input
#' @param cut a named numerical array of codon usage table. Case-insensitive three letter codon names are used. Numbers could be the actuall codon counts or scaled ratios
#' @param y a number denoting the codon table to be used. When \emph{cut} is specified, \emph{y} also needs to be specified. 0 = standard codon table (default), 2 = vertibrate mitchodrial, 3 = ?, 4 = ?. When an option other than those mentioned above is provided, the standard codon table will be used
#' @param spp a char string of the host species name for which the codon usage table will be used. Four build-in reference tables are available: 'ecoli', E. coli; 'insect', insect; 'yeast', yeast; 'celegans', C. elegans. When specified, this option overides \emph{cut}. These build-in codon usage tables are taken from https://www.genscript.com/tools/codon-frequency-table.
#' @return a char array of random reverse translated codons
#' @examples 
#' 
#' @export

rrt <- function(x, cut, y=0, spp){
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
  
  if(x[1]%in%cTable[,1]){
    temp = dna2aa(x,1,y)
  } else if (x[1]%in%cTable[,2]){
    temp = x
  } else {
    tabx = unique(cTable[,2:3])
    tab = tabx[,1]
    names(tab) = tabx[,2]
    temp = unname(tab[x])
  }
  
  nx = length(x)
  ans = character(nx)
  types = table(cTable[,2])
  typename = names(types)
  for(i in 1:length(types)){
    indx = which(temp==typename[i])
    codons = cTable[cTable[,2]==typename[i],1]
    temp[indx] <- apply(rmultinom(length(indx), 1, cft[codons]), 2, function(x){codons[as.logical(x)]})
  }
  return(temp)
}