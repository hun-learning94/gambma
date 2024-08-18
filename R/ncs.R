#' Title
#'
#' @param x
#' @param nk
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
ncs = function(x, nk = 10, lambda = 0.2){
  term = as.character(substitute(x))
  label = paste0("ncs(", term, ")")
  ret = list(term = term, label = label, nk = nk, lambdaVec = lambda)
  class(ret) = "ncsSmoothSpec"
  ret
}
