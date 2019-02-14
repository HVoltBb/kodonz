#' Random reverse translation
#'
#' Randomly reverse translate an animo acid sequence into a DNA sequence according to a codon usage table and a codon table
#'
#' @param x a char array of animo acids (e.g., either 'G' or 'Gly' is acceptable but do not mix one letter symbols with three letter symbols in the same array) or codons (e.g., 'GGG'). The first element of the array will be used to determine the type of the input
#' @param cut a named numerical array of codon usage table. Case-insensitive three letter codon names are used. Numbers could be the actuall codon counts or scaled ratios
#' @param y a number denoting the codon table to be used. When \emph{cut} is specified, \emph{y} also needs to be specified. 0 = standard codon table (default), 2 = vertibrate mitchodrial and see dna2aa for more options. When an option other than those mentioned above is provided, the standard codon table will be used
#' @param spp a char string of the host species name for which the codon usage table will be used: 'ecoli', E. coli; 'insect', insect; 'yeast', yeast; 'celegans', C. elegans, and see pminmax for more options. When specified, this option overides \emph{cut}. These build-in codon usage tables are taken from https://www.genscript.com/tools/codon-frequency-table.
#' @return an KZsqns object of random reverse translated codons
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
  attr(temp, 'class') <- "KZsqns"
  return(temp)
}
