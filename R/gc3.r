#' Calculate GC3 value of a codon sequence
#'
#' Outputs an array of numbers, representing the proportion of codons with G/C bases in the third position
#'
#' @param x a list of KZsqns objects.
#' @return a number array of GC3s.
#'
#' @examples
#' data('CodonTable0')
#' x = vector('list', 5) # Creating an empty list of length 5
#' for(i in 1:5){
#'   x[[i]] = CodonTable0[sample(1:64, 10*i, TRUE),1]
#'   attr(x[[i]], 'class') = 'KZsqns'
#' }
#'
#' cat(gc3(x))
#' @export

gc3 <- function(x){
  ans = matrix(0, length(x), 1)
  colnames(ans) <- 'GC3'
  rownames(ans) <- paste0('s_', 1:length(x))
  temp = n3_freq(x)
  ans[,1] = rowSums(temp[,c('G3', 'C3'), drop=F])
  return(ans)
}
