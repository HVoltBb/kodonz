#' Effective number of codons
#'
#' This function calculate the effective number of codons (ENC) using Wright's formula.
#'
#' @param x a list of KZsqns objects
#' @param spp.list an array indicating the the species identity of each KZsqns object in the list \emph{x}. The default option is to calculate an ENC for each KZsqns object. When \emph{spp.list} is specified, an ENC is calculated for each unique element of \emph{spp.list}
#' @param y a number to assign the type of codon table to be used for the translation. 0 = standard codon table, 2 = vertibrate mitchodrial (See dna2aa for additional options). If no number is specified or an option other than those provided above is specified, there will be a warning and the standard codon table will be used
#' @param formula a char array to denote the formula for calculating ENC. Options include: 'w': Wright's formula, 'f4-': Fuglsang's Ncf4- formula.
#' @return an numeric array of the ENCs
#' @section Warning:
#' According to Wright's formula, the ENC is the sum of the heterogeneity of all the degeneracy classes. When all the codons belonging to a degeneracy class is missing, the ENC can not be calculated. A warning will result when it happens. The most likely situation leading to such a condition is that the sequence of codons is too short. When the sequence is short, the calculated ENC is likely to deviate much from the underlying value. It is recommended to increase the sample size to solve this issue when Wright's formula is prefered. Alternatively, Fuglsang's Ncf4- formula could be used instead.
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

enc <- function(x, y=0, spp.list='n', formula = 'w'){
  if(!is.list(x)) x = list(x)
  cTable = arg_check(x[[1]], "KZsqns", y)
  calc = nc_wright
  switch(formula,
         'f4-' = {calc = nc_fuglsang4m},
         'w'   = {},
         {
           warning("Not a valid formula, and Wright's formula will be used.")
         })
  if(length(spp.list)==1 && spp.list=='n'){
    cat("No species list is provided\nCalculating ENCs for each KZsqns object.\n")
    return(sapply(x, function(i){calc(table(i), cTable)}))
  } else {
    if (length(spp.list)!=length(x)) stop("The lengths of the species list and the list of KZsqns do not match")
    spp = unique(spp.list)
    fixx = function(i){
      jobs = which(spp.list==spp[i])
      ans = rep(0, 64)
      for(j in jobs){
        ans = rbind(ans, table(x[[j]])[cTable[,1]])
      }
      ans = colSums(ans, na.rm = TRUE)
      names(ans) <- cTable[,1]
      return(ans)
    }
    ans = sapply(1:length(spp), function(i){calc(fixx(i), cTable)})
    names(ans) <- spp
    return(ans)
  }
}
