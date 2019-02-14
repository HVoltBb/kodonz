#' Composition analysis
#'
#' Calculates nucleotide composition of a list of DNA sequenses
#'
#' @export

comp_ly <- function(x){
  if(!is.list(x)){
    cat("Just one perhap very long sequence?\n")
    x = list(x)
  }
  return(cbind(base_freq(x), n3_freq(x), augc_freq(x), gcx(x,1), gcx(x,2), gcx(x,3)))
}
