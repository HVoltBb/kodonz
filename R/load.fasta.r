#' A function to load fasta files
#'
#' Loads the specified fasta file
#'
#' @param file Name of the fasta file
#' @param x a character array of DNA sequence. Use "A", "T", "G" and "C" to denote the nucleatide composition. When a fasta file is also specified, the fasta file takes priority, and this character sequence is ignored
#' @param ... additional parameters you might want to pass to the seqinr::read.fasta function
#'
#' @return A list of objects of class KZsqns
#' @examples
#' dna1 = load.fasta("mitochondrial01.fasta", as.string = TRUE) # not run
#' @export

load.fasta <- function(file='n', x,...){
  if(file!='n'){
    if(!file.exists(file)) {stop("\nThe specified file does not exist.\n Please check the file name and try again.\n")}
    temp = seqinr::read.fasta(file, seqtype = "DNA", ...)
  } else if(exists('x')){
    temp = strsplit(x,'')
  } else {
    stop("Either a file or DNA sequence need to be provided.")
  }

  fixx = function(x){
    xx = temp[[i]][3*x+c(-2,-1,0)]
    if(all(xx %in% c("A", "T", "G", "C", "a", "t", "g", "c"))){
      return(toupper(xx))
    }
  }

  for(i in 1:length(temp)){
    ele = temp[[i]]
    if(length(ele)%%3!=0) stop(paste0("\nThe length of sequence ",i," is not a multiple of 3. Check your sequence and try again.\n "))
    ubases = unique(ele)
    if(!all(ubases %in% c("A", "T", "G", "C"))){
      warning(paste0("\nSequence", i, " contains ambiguous base(s) or spacing(s), which will be removed from further analysis.\n"))
      ele = unlist(sapply(1:(length(ele)/3), fixx))
    }
    bases0 = paste0(ele[c(T,F,F)], ele[c(F,T,F)], ele[c(F,F,T)]) # only keep codons with unambiguous bases
    temp[[i]] = bases0
    attr(temp[[i]], "class") <- "KZsqns"
  }

  return(temp)
}
