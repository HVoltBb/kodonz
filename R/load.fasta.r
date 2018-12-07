#' A function to load fasta files
#'
#' Loads the specified fasta file
#'
#' @param file Name of the fasta file.
#' @param ... additional parameters you might want to pass to the read.fasta function.
#'
#' @return A list of objects of class KZsqns.
#' @examples
#' dna1 = load.fasta("mitochondrial01.fasta", as.string = TRUE) # not run
#' @export

load.fasta <- function(file, code = 0, ...){
  if(!file.exists(file)) {stop("\nThe specified file does not exist.\n Please check the file name and try again.\n")}
  temp = seqinr::read.fasta(file, seqtype = "DNA", ...)
  ans = vector("list", length(temp))
  for(i in 1:length(temp)){
    element = list()
    attr(element, "class") <- "KZsqns"
    element$sequence <- temp[[i]]

    if(length(temp[[i]])%%3!=0) stop(paste0("\nThe length of sequence ",i," is not a multiple of 3. Check your sequence and try again.\n "))
    ubases = unique(temp0)
    element$containAmbiguousBases <- FALSE
    if(bases[5] %in% ubases || bases[6] %in% ubases || bases[7] %in% ubases || bases[8] %in% ubases || bases[8] %in% ubases || bases[9] %in% ubases || bases[10] %in% ubases || bases[11] %in% ubases || bases[8] %in% ubases || bases[8] %in% ubases || bases[8] %in% ubases) {warning(paste0("\nSequence", i, " contains ambiguous base(s), which will be removed from further analysis.\n"))
      element$containAmbiguousBases <- TRUE
    }
    bases0 = table(paste0(element$sequence[c(T,F,F)], element$sequence[c(F,T,F)], element$sequence[c(F,F,T)]))[c()] # only keep codons with unambiguous bases
    names(bases0) <- c()
    element$codonComposition <- bases0
    ans[[i]] <- element
  }

  return(ans)
}
