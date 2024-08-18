#' A function to construct a natural cubic spline basis function
#'
#' @description
#' for internal usage 
#'
#' @param x variable
#' @param nk integer; number of knots for each variable
#'
#' @return ncsSmoothSpec object
#' @export
#' @seealso [gambms()]
#' 
#' @examples
ncs = function(x, nk = 10){
  term = as.character(substitute(x))
  label = paste0("ncs(", term, ")")
  ret = list(term = term, label = label, nk = nk)
  class(ret) = "ncsSmoothSpec"
  ret
}
