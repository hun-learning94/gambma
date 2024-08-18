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
#' @param u
#' @param a
#' @param b
#' @param z
#' @param s
#' @param nu
#' @param theta
#' @param log
#'
#' @return
#' @export
#'
#' @examples
dtCCH = function(u, a, b, z, s, nu, theta, log=F){
  x = 1/theta - 1
  snu = s/nu
  v = u*nu
  res = ((a-1)*log(v) + (b-1)*log1p(-v) -z*log1p(x*v)- snu*v + snu +
           z*log1p(x) - .logPhi1(b, z, a+b, snu, x/(x+1)) - lbeta(a,b))*nu
  if(!log) res = exp(res)
  return(res)
}
