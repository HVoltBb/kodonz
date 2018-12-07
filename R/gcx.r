#' GCx of a gene
#'
#' Calculates the GC1,2,3 and the overall GC content of a gene
#'
#' @param x a char array of codons
#' @param pos the position of the codon for which GCx is calculated. 1: GC1; 2:GC2; 3: GC3; for any other values, the overall GC content will be calculated with a warning
#' @return a number of the GCx of the gene
#' @examples
#'
#' @export

gcx <- function(x, pos=1){
  if(class(x)!='character') stop("A character array of codons is expected.")

  indx = c(T,T,T)
  switch(as.character(pos),
         '1' = {
           indx[2:3] = F
         },
         '2' = {
           indx[c(1,3)] = F
         },
         '3' = {
           indx[1:2] = F
         },
         {
           warning("Not a valid pos option. Overall G/C content will be calculated.")
         })
  w = strsplit(paste(x,collapse = ''), split = '')[[1]][indx]
  nogc = table(w)[c("G", "C")]
  nogc[is.na(nogc)] <- 0
  return(unname(nogc[1]+nogc[2])/length(w))
}
