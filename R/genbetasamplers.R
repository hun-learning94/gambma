#' Title
#'
#' @param n
#' @param a
#' @param b
#' @param s
#'
#' @return
#' @export
#'
#' @examples
rCH = function(n, a, b, s){
  return(.rCH(n, 100, a, b, s))
}

#' Title
#'
#' @param n
#' @param a
#' @param b
#' @param x
#' @param z
#'
#' @return
#' @export
#'
#' @examples
rGH = function(n, a, b, x, z){
  return(.rGH(n, 100, a, b, x, z))
}

#' Title
#'
#' @param n
#' @param a
#' @param b
#' @param z
#' @param s
#' @param nu
#' @param theta
#'
#' @return
#' @export
#'
#' @examples
rtCCH = function(n, a, b, z, s, nu, theta){
  return(.rtCCH(n, 100, a, b, z, s, nu, theta))
}

#' Title
#'
#' @param n
#' @param a
#' @param b
#' @param z
#' @param w
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
rAH = function(n, a, b, z, w, x, y){
  return(.rAH(n, 100, a, b, z, w, x, y))
}
