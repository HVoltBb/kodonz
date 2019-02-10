#' Remove non-degenerate and stop codons
#'
#' @description Removes non-redundant and stop codons based on the specified codon table. For some compositional analysis, it is advisable to remove non-redundant amino acid coding codons and stop codons.
#' @param x a list of KZsqns objects
#' @param y a vector of numbers denoting the codon table to be used for each KZsqns object. 0 = standard codon table (default), 2 = vertibrate mitchodrial (See dna2aa for additional options). When an option other than those mentioned above is provided, the standard codon table will be used. When the length of this vector is smaller than the number of KZsqns objects supplied, standard codon tables will be used for those remaining KZsqns objects
#' @param keep.stop an option indication whether stop codons should also be kept. TRUE: keep stop codons if they are degenerate (default); FALSE: remove stop codons
#' @export

rm_ndc <- function(x, y=0, keep.stop=TRUE){
  if(class(x)!='list') stop("A list of KZsqns is required. Use load.fasta() to import DNA sequences.")
  return(lapply(x, ndc, y, keep.stop))
}
