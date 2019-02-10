#' Custom CAI reference table
#'
#' @param x a list of KZsqns object, from which the reference table will be calculated
#' @param y y an optional number denoting the codon table to be used. 0 = standard codon table (default), 2 = vertibrate mitchodrial (See dna2aa for additional options). When an option other than those mentioned above is provided, the standard codon table will be used
#' @param method a charactor string denoting which method to use for codons absent from the reference set: "bulmer", set the adaptiveness of the codon to 0.01; "sharp", in calculating the adaptiveness of an absent codon, use a frequency of 0.5 instead
#' @details If any codon is missing from the reference set, the adaptiveness value of that codon need to be adjusted in order to obtain valid CAI for sequences that have that codon in their sequence. The sharp's method is the default option [1], but the user may choose to use bulmer's method [2].
#'
#' @description A helper method to calculate the codon adaptiveness index reference table based on a list of highly expressed genes
#' @references [1] Sharp, P. M., & Li, W.-H. (1987). The codon adaptation index-a measure of directional synonymous codon usage bias, and its potential applications. Nucleic acids research, 15(3), 1281-1295.
#' [2] Bulmer, M. (1988). Are codon usage patterns in unicellular organisms determined by selection-mutation balance? Journal of Evolutionary Biology, 1(1), 15-26.
#' @export

cai_crt <- function(x, y, method=c("sharp", "bulmer")){
  cTable = arg_check(x[[1]], 'KZsqns',y)
  ttx = lapply(x, table)

  tx = rep(0, 64)
  names(tx) <- cTable[,1]

  for (i in 1:length(ttx)){
    tx = colSums(rbind(tx+ttx[[i]][cTable[,1]]), na.rm = T)
  }
  ## cat(which(tx==0))

  if(method[1] == "sharp"){
    tx[tx==0] <- 0.5
  }

  temp = table(cTable[,2])
  degv = numeric(0)
  for(j in names(temp)){
    dx = which(cTable[,2]==j)
    if(sum(tx[dx])==0){
      degv[dx] = NA
      warning('An animo acid is missing, and NAs will be produced for this animo acid')
    } else {
      degv[dx] = tx[dx]/max(tx[dx])
    }
  }
  names(degv) = cTable[,1]
  degv[which(cTable[,2]=='*')] <- NA
  if(method[1]=="bulmer"){
    degv[degv==0] <- 0.01
  }
  return(degv)
}
