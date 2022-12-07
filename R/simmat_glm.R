simmat_glm = function(xmin = 0, xmax = 1, 
                      n = 1000,  
                      family = c("poisson", "binomial", "gaussian"),
                      link = c("logit", "probit", "log", "identity"),
                      nbinom = NULL, sig = NULL,
                      A = NULL,
                      intercept = 0,
                      beta = c(-1, 1))
{
  if(family == "gaussian") family = "normal"
  EF_available = c("poisson", "binomial", "normal")
  if(!(family %in% EF_available)) stop("not supported glm family")
  p = length(c(beta))
  
  # generate design matrix
  X = matrix(runif(n*p, xmin, xmax), nrow = n, ncol = p)
  
  # generate target response
  eta = X %*% beta + intercept
  if(!is.null(A)) eta = eta+A
  
  if(family == "binomial"){
    if(link == "logit") y = rbinom(n, nbinom, 1/(1+exp(-eta)))
    if(link == "probit") y = rbinom(n, nbinom, pnorm(eta))
  }
  if(family == "poisson"){
    if(link == "log") y = rpois(n, exp(eta))
  }
  if(family == "normal"){
    if(link == "identity") y = rnorm(n, mean = eta, sd = sig)
  }
  
  res = data.frame(cbind(y, X))
  colnames(res) = c("y", 
                    paste("x", 1:p, sep=""))
  return(res)
}

