#' Synonymous codon usage order
#'
#' Calculates the synonymous codon usage order of a codon sequence based on a codon table.
#'
#' @param x a KZsqns object of codon sequence or a list of such objects
#' @param y a number to assign the type of codon table to be used for the translation. 0 = standard codon table, 2 = vertibrate mitchodrial, . If no number is specified or an option other than those provided above is specified, there will be a warning and the standard codon table will be used
#' @return a numeric array of synonymous codon usage order, or a list of such arrays if \emph{x} is a list
#' @examples
#'
#' @export

scuo = function(x, y=0){
  if(is.list(x)){
    return(lapply(x, scuo_, y))
  } else {
    return(scuo_(x, y))
  }
}
