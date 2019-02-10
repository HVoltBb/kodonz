#' Arguments checks
#'
#' Performs the usual checks of arguments consistency of codon sequences and codon table to be used
#'
#' @param x a char array of codons
#' @param mode a char array denoting the expected mode of the argument x
#' @param y a number denoting the codon table to be used. 0 = standard codon table, 2 = vertibrate mitchodrial, 3 = yeast mitochondrial, 4 = mold, protozoan, and Coelenterate mitochondrial and Mycoplasma/Spiroplasma, 5 = invertebrate mitochandrial, 6 = Ciliate, Dasycladacean and Hexamita nuclear, 9 = Echinoderm and flatworm mitochondrial, 10 = Euplotid nuclear, 11 = bacterial, Archaeal and plant plastid, 12 = alternative yeast nuclear, 13 = Ascidian mitochondrial, 14 = alternative flatworm mitochondrial, 16 = Chlorophycean mitochondrial, 21 = Trematode mitochondrial, 22 = Scenedesmus obliquus mitochondrial, 23 = Thraustochytrium mitochondrial, 24 = Pterobranchia mitochondrial, 25 = candidate division SR1 and Gracilibacteria, 26 = Pachysolen tannophilus nuclear, 29 = Mesodinium nuclear, 30 = Peritrich nuclear.
#' @return a numeric matrix of the appropriate codon table
#'

arg_check <- function(x, mode, y){
  if(class(x)!=mode) stop(paste0("A ", mode, " object is expected."))

  cTable = CodonTable0
  switch(as.character(y),
         '0' = {},
         '2' = {cTable = CodonTable2},
         '3' = {cTable = CodonTable3},
         '4' = {cTable = CodonTable4},
         '5' = {cTable = CodonTable5},
         '6' = {cTable = CodonTable6},
         '9' = {cTable = CodonTable9},
         '10' = {cTable = CodonTable10},
         '11' = {cTable = CodonTable11},
         '12' = {cTable = CodonTable12},
         '13' = {cTable = CodonTable13},
         '14' = {cTable = CodonTable14},
         '16' = {cTable = CodonTable16},
         '21' = {cTable = CodonTable21},
         '22' = {cTable = CodonTable22},
         '23' = {cTable = CodonTable23},
         '24' = {cTable = CodonTable24},
         '25' = {cTable = CodonTable25},
         '26' = {cTable = CodonTable26},
         '29' = {cTable = CodonTable29},
         '30' = {cTable = CodonTable30},
         {
           warning("Not a valid type of codon table, and the standard codon table will be used.")
         }
  )
  return(cTable)
}
