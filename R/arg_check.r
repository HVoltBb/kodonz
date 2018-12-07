#' Arguments checks
#'
#' Performs the usual checks of arguments consistency of codon sequences and codon table to be used
#'
#' @param x a char array of codons
#' @param mode a char array denoting the expected mode of the argument x
#' @param y a number denoting the codon table to be used
#' @return a numeric matrix of appropriate codon table
#'

arg_check <- function(x, mode, y){
  if(class(x)!=mode) stop(paste0("A ", mode, " object is expected."))

  cTable = CodonTable0
  switch(as.character(y),
         '0' = {},
         '2' = {cTable = CodonTable2},
         '3' = {cTable = CodonTable3},
         '4' = {cTable = CodonTable4},
         {
           warning("Not a valid type of codon table, and the standard codon table will be used.")
         }
  )
  return(cTable)
}
