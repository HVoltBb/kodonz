#' Expected number of ENC given GC3 content
#' 
#' Calculates the expected ENC based on Wright's formula given GC3 percentage
#' @param gc3 a number of GC3 content
#' @return a number of the expected ENC given the GC3 content
#' 
#' @export

enc_exp <- function(gc3=0){
  return(2+gc3+29/(gc3^2+(1-gc3)^2))
}