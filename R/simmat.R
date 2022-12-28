#' Title
#'
#' @param flist
#' @param xmin
#' @param xmax
#' @param n
#' @param family
#' @param link
#' @param sig
#' @param intercept
#' @param beta
#'
#' @return
#' @export
#'
#' @examples
simmat = function(flist,
                  xmin = 0, xmax = 1,
                  n = 1000,
                  family = c("bernoulli", "poisson", "gaussian"),
                  link = NULL,
                  sig = NULL,
                  intercept = 0,
                  beta = NULL)
{
  if(family == "bernoulli"){family = "binomial"; glmWeight = 1.0}
  if(!(family %in%  c("binomial", "poisson", "gaussian")))
    stop("unsupported exponential family")

  ## y_scale, glmWeight
  if(family == "binomial"){
    glmWeight = 1
  } else if (family == "poisson"){
    glmWeight = 1.0;
  } else if (family == "gaussian"){
    glmWeight = 1.0;
    if(is.null(sig)) sig = .1
  }

  ## link function
  if(is.null(link)) {
    if (family == "poisson"){link <- "log"}
    if (family == "binomial"){link <- "logit"}
    if (family == "gaussian"){link <- "identity"}
  }

  if(!is.list(flist)) stop("input list of functions")
  p = length(flist)
  Xgrid = seq(xmin, xmax, len = n)

  # include linear predictor
  if(!is.null(beta)){
    plin = length(c(beta))
    Xlin = matrix(runif(n*plin, xmin, xmax), nrow = n, ncol = plin)
    Xlin = sweep(Xlin, 2, colMeans(Xlin), "-")
  }

  # generate design matrix
  X = matrix(runif(n*p, xmin, xmax), nrow = n, ncol = p)

  # generate target response
  eta = rep(0, n) + intercept
  for(i in 1:p){
    f = flist[[i]]
    eta = eta + f(X[,i])
  }
  if(!is.null(beta)){
    eta = eta + Xlin %*% beta
    X = cbind(Xlin, X)
  }

  if(family == "binomial"){
    if(link == "logit") y = rbinom(n, glmWeight, 1/(1+exp(-eta)))
    if(link == "probit") y = rbinom(n, glmWeight, pnorm(eta))
  }
  if(family == "poisson"){
    if(link == "log") y = rpois(n, exp(eta))
  }
  if(family == "gaussian"){
    if(link == "identity") y = rnorm(n, mean = eta, sd = sig)
  }

  res = data.frame(cbind(y, X))
  if(!is.null(beta)){
    colnames(res) = c("y",
                      paste("lin_x", 1:plin, sep=""),
                      paste("x", 1:p, sep="")
    )
  } else {
    colnames(res) = c("y",
                      paste("x", 1:p, sep="")
    )
  }
  return(res)
}

