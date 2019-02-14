#' Overall AT and GC content
#'
#' Calculate the overall AT and GC content
#' @param x a list of KZsqns objects
#' @return a numeric matrix of AU and GC content
#' @export

augc_freq <- function(x){
  if(!is.list(x)){
    cat("Just one sequence?\n")
    x = list(x)
  }
  temp = base_freq(x)
  ans = cbind(temp[,'A',drop=F]+temp[,'T',drop=F], temp[, 'G', drop=F]+temp[, 'C', drop=F])
  colnames(ans) <- c('AT', 'GC')
  return(ans)
}
