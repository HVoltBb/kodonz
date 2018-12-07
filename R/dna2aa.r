#' Translating a codon sequence to into animo-acids
#'
#' This function translate an array of codons into the corresponding amino acid sequence using the codon table of choice.
#'
#' @param dnaSequence a character array of codons
#' @param aaNoLetters a number to assign the format of animo acid. 1 = one-letter abbreviation, 3 = three-letter abbreviation. If it defaults into the 3-letter abbreviation when no number or a different number is specified
#' @param codonTable a number to assign the type of codon table to be used for the translation. 0 = standard codon table, 2 = vertibrate mitchodrial, . If no number is specified or an option other than those provided above is specified, there will be a warning and the standard codon table will be used
#' @return a character array of animo-acids
#' @examples
#' data('CodonTable0')
#' sqs = sample(1:64, 100, replace=T)
#' codonSequence = CodonTable0[sqs,1] # an arbitrary sequence of codons is genereated based on the standard codon table
#' cat("The sequence of aa based on the standard codon table:\n")
#' cat(CodonTable0[sqs,2])
#' cat("The translated aa sequence:\n")
#' cat(dna2aa(codonSequence, 1, 0))
#'
#' @export
dna2aa <- function(dnaSequence, aaNoLetters = 3, codonTable = 0){
  cTable = arg_check(dnaSequence, 'character', codonTable)
  if(aaNoLetters==1){
    cTable = setNames(cTable[,2],cTable[,1])
  } else {
    cTable = setNames(cTable[,3],cTable[,1])
  }
  return(unname(cTable[dnaSequence]))
}
