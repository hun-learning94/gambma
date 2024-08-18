#' Title
#'
#' @param x
#' @param nk
#'
#' @return
#' @export
#'
#' @examples
ncs = function(x, nk = 10){
  term = as.character(substitute(x))
  label = paste0("ncs(", term, ")")
  ret = list(term = term, label = label, nk = nk)
  class(ret) = "ncsSmoothSpec"
  ret
}
