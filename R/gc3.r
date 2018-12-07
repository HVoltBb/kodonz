#' Calculate GC3 value of a codon sequence
#'
#' Outputs an array of numbers, representing the proportion of codons with G/C bases in the third position
#'
#' @param x a list of objects of class KZsqns.
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
  if(class(x)!='list') stop("A list of KZsqns is required. Use load.fasta() for importing DNA sequences.")

  ans = numeric(length(x))
  for(i in 1:length(x)){
    if(class(x[[i]]!='KZsqns')) warning("KZsqns objects are expected")
    gc = 0
    for(j in x[[i]]){
      if (substr(j, 3,3)=="G" || substr(j,3,3)=='C') gc = gc + 1
    }
    ans[i] = gc/length(x[[i]])
  }
  return(ans)
}
