#' Nucleotide at the third position of the codon
#'
#' Calculates the nucleotide frequency at the third position of the codon
#' @param x a list of KZsqns objects.
#' @return a matrix of nucleotide composition at the third position
#' @export

n3_freq <- function(x){
  if(!is.list(x)){
    cat("Just one perhap very long sequence?\n")
    x = list(x)
  }

  ans = matrix(0, length(x), 4)
  for(i in 1:length(x)){
    if(class(x[[i]])!='KZsqns') warning("KZsqns objects are expected")
    freq_table = table(c('A', 'C', 'G', 'T', strsplit(paste0(x[[i]], collapse = ''),'')[[1]][c(F, F, T)]))
    ans[i,] = (freq_table-rep(1,4))/length(x[[i]])
  }

  colnames(ans) <- c('A3', 'C3', 'G3', 'T3')
  rownames(ans) <- paste0('s_', 1:length(x), sep='')
  return(ans)
}
