#' Overall nucleotides frequency
#'
#' Calculate the overall nucleotides frequency composition
#'
#' @param x a list of KZsqns objects.
#' @return a n-by-4 numeric matrix, where n is the length of list x
#'
#' @examples
#' data(CodonTable0)
#' x = vector('list', 5) # Creating an empty list of length 5
#'
#' for(i in 1:5){
#' x[[i]] = CodonTable0[sample(1:64, 10*i, TRUE),1]
#' attr(x[[i]], 'class') = 'KZsqns'
#' }
#'
#' cat(base_freq(x))
#' @export

base_freq<- function(x){
  if(!is.list(x)){
    cat("Just one perhap very long sequence?\n")
    x = list(x)
  }

  ans = matrix(0, length(x), 4)
  for(i in 1:length(x)){
    if(class(x[[i]])!='KZsqns') warning("KZsqns objects are expected")
    freq_table = table(strsplit(paste0(paste0(x[[i]], collapse = ''),'ATGC', sep=''),''))
    ans[i,] = (freq_table-rep(1,4))/length(x[[i]])/3
  }

  colnames(ans) <- c('A', 'C', 'G', 'T')
  rownames(ans) <- paste0('s_', 1:length(x), sep='')
  return(ans)
}
