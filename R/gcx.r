#' GCx of a list of genes
#'
#' Calculates the GC1,2,3 and GC12 contents of a list of genes
#'
#' @param x a list of KZsqns objects.
#' @param pos the position of the codon for which GCx is calculated. 1: GC1; 2:GC2; 3: GC3; 4:GC12
#' @return a matrix of the GCx of the list of genes
#' @examples
#'
#' @export

gcx <- function(x, pos=1){
  if(!is.list(x)){
    cat("Just one perhap very long sequence?\n")
    x = list(x)
  }
  ans = matrix(0, length(x), 1)
  rownames(ans) <- paste0('s_', 1:length(x))
  indx = c(T,T,T)
  adj = 1
  switch(as.character(pos),
         '1' = {
           indx[2:3] = F
           colnames(ans) <- 'GC1'
         },
         '2' = {
           indx[c(1,3)] = F
           colnames(ans) <- 'GC2'
         },
         '3' = {
           indx[1:2] = F
           colnames(ans) <- 'GC3'
         },
         '4' = {
           indx[3] = F
           colnames(ans) <- 'GC12'
           adj = 2
         },
         {
           stop("Not a valid pos value")
         })
  for(i in 1:length(x)){
    temp = table(strsplit(paste0(paste(x[[i]],collapse = ''), 'GGGCCC', sep=''), split = '')[[1]][indx])
    ans[i,1] = (temp["G"] + temp["C"] - 2*adj)/length(x[[i]])/adj
  }
  return(ans)
}
