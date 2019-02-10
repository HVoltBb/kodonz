#' Effective number of codons
#'
#' This function calculate the effective number of codons (ENC) using Wright's formula.
#'
#' @param x a KZsqns object
#' @param y a number to assign the type of codon table to be used for the translation. 0 = standard codon table, 2 = vertibrate mitchodrial (See dna2aa for additional options). If no number is specified or an option other than those provided above is specified, there will be a warning and the standard codon table will be used
#' @param formula a char array to denote the formula for calculating ENC. Options include: 'w': Wright's formula, 'f4-': Fuglsang's Ncf4- formula.
#' @return an numeric array of the ENCs
#' @references
#' Fuglsang, A., 2006 Estimating the "Effective Number of Codons": The Wright way of determining codon homozygosity leads to superior estimates. Genetics 172:1301-1307
#' Wright, F., 1990 The 'effective number of codons' used in a gene. Gene 87: 23-29.
#' @examples
#' dnalist = vector('list', 2)
#' dnalist[[1]] = CodonTable2[sample(1:64, 1e4, TRUE),1]
#' attr(dnalist[[1]], 'class') <- 'KZsqns'
#' dnalist[[2]] = CodonTable2[sample(1:64, 1e4, TRUE),1]
#' attr(dnalist[[2]], 'class') <- 'KZsqns'
#' enc(x = dnalist, y = 2) # Two numbers should be around 60 as there are 60 animo acid coding codons and 4 stop codons for codon table 2
#'
#' dnalist = vector('list', 2)
#' dnalist[[1]] = CodonTable2[sample(1:64, 1e4, TRUE),1]
#' attr(dnalist[[1]], 'class') <- 'KZsqns'
#' dnalist[[2]] = CodonTable2[sample(1:64, 1e4, TRUE),1]
#' attr(dnalist[[2]], 'class') <- 'KZsqns'
#' enc(x = dnalist, y = 4) # Two numbers should be around 62 as there are 60 animo acid coding codons and 2 stop codons for codon table 4
#' @export

enc <- function(x, y=0, formula = 'w'){
  cTable = arg_check(x[[1]], "KZsqns", y)
  calc = nc_wright
  switch(formula,
         'f4-' = {calc = nc_fuglsang4m},
         'w'   = {},
         {
           warning("Not a valid formula, and Wright's formula will be used.")
         })
  return(sapply(x, calc, cTable))
}
