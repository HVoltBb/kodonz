#' Translating a codon sequence to into animo-acids
#'
#' This function translate an array of codons into the corresponding amino acid sequence using the codon table of choice.
#'
#' @param dnaSequence a character array of codons
#' @param aaNoLetters a number to assign the format of animo acid. 1 = one-letter abbreviation, 3 = three-letter abbreviation. If it defaults into the 3-letter abbreviation when no number or a different number is specified
#' @param codonTable a number to assign the type of codon table to be used for the translation. 0 = standard codon table, 2 = vertibrate mitchodrial, 3 = yeast mitochondrial, 4 = mold, protozoan, and Coelenterate mitochondrial and Mycoplasma/Spiroplasma, 5 = invertebrate mitochondrial, 6 = Ciliate, Dasycladacean and Hexamita nuclear, 9 = Echinoderm and flatworm mitochondrial, 10 = Euplotid nuclear, 11 = bacterial, Archaeal and plant plastid, 12 = alternative yeast nuclear, 13 = Ascidian mitochondrial, 14 = alternative flatworm mitochondrial, 16 = Chlorophycean mitochondrial, 21 = Trematode mitochondrial, 22 = Scenedesmus obliquus mitochondrial, 23 = Thraustochytrium mitochondrial, 24 = Pterobranchia mitochondrial, 25 = candidate division SR1 and Gracilibacteria, 26 = Pachysolen tannophilus nuclear, 29 = Mesodinium nuclear, 30 = Peritrich nuclear. If no number is specified or an option other than those provided above is specified, there will be a warning and the standard codon table will be used
#' @return a character array of animo-acids
#' @details Codon tables were based on the list compiled by Andrzej Elzanowski and Jim Ostell [1] at NCBI, Bethesda, Maryland, USA. The numbering of codon tables respects the order in Elzanowski and Ostell [1] (except the standard code, which was numbered 0 in this package). There are three addditional codon tables listed in Elzanowski and Ostell [1] (27 Karyorelict nuclear, 28 Condylostoma nuclear and 31 Blastocrithidia nuclear), but were excluded in this package, because these codon tables include codons that can be either intepreted as a AA coding codon or stop codon, and the AA sequence can not be determined based on the DNA sequence alone.
#'The codon table data were available as "CodonTable#", where "#" denote the numeric order of the codon table. For example, "CodonTable21" is the Trematode mitochonrial code table.
#' @references [1] https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG13 [retrieved on 2/4/2019]
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
