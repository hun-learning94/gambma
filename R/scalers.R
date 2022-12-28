#' Scale to a user-determined interval
#'
#' @param X A matrix to be scaled column-wise.
#' @param mins A vector of minimums of scaled intervals. Its length should be equal to the number of columns of `X`.
#' @param maxs A vector of maximums of scaled intervals Its length should be equal to the number of columns of `X`.
#'
#' @return A scaled matrix.
#' @export
#'
#' @examples
scaleto01 = function(X, mins=NULL, maxs=NULL){
  if(is.null(mins)) stop("need mins")
  if(is.null(maxs)) stop("need maxs")
  X = sweep(X, 2, mins, "-")
  X = sweep(X, 2, maxs-mins, "/")
  return(X)
}

#' Scale back to the original scale
#'
#' @param X A matrix to be scaled column-wise.
#' @param mins A vector of minimums of scaled intervals. Its length should be equal to the number of columns of `X`.
#' @param maxs A vector of maximums of scaled intervals. Its length should be equal to the number of columns of `X`.
#'
#' @return A scaled matrix.
#' @export
#'
#' @examples
scalefrom01 = function(X, mins=NULL, maxs = NULL){
  if(is.null(mins)) stop("need mins")
  if(is.null(maxs)) stop("need maxs")
  X = sweep(X, 2, maxs-mins, "*")
  X = sweep(X, 2, mins, "+")
  return(X)
}
