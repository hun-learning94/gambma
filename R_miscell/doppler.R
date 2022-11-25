plotf = function(f){
  u = seq(0, 1, by = 0.01)
  plot(u, f(x=u), type="l")
}


## target functions
doppler = function(x, alpha){
  1.5/exp(4*x) * sin((2 * pi * (1 + 2 ^ (( 9 - 4 * alpha ) / 5))) /(x + 2 ^ ((9 - 4 * alpha) / 5)))
}

hun1 = function(x){
  -4*log(x+1) * sin(4*pi*x)
}

hun2 = function(x){
  x = 2*x - 1
  exp(1-x^3)*sin(2.5*pi*x^2) * 0.4
}

hun3 = function(x){
  sin(4*pi*x)*2
}

hun4 = function(x){
  x = 2*x -1
  1.6 * x^2 * (x^3 + 2*exp(-3*x^4 + log(2*x + pi)))
}

hun5 = function(x) doppler(x, 4.5)*1.7
hun6 = function(x) cos(4*pi*x)*2
hun7 = function(x){
  2*(-as.numeric(x>=0 & x<0.5) +
    0.5*as.numeric(x<=1 & x>=0.5) + 0.5*as.numeric(x<=1 & x>=0.75))
}

hun8 = function(x){
  sin(2*pi*x)*2
}

hun9 = function(x){
  exp(-abs(2*x -1))*2
}

hun10 = function(x){
  as.numeric(x > 0.5)* x +
    as.numeric(x <= 0.5) * (1-x)
}

f_list1 = list(f1 = hun8,
               f2 = hun3,
               f3 = hun7,
               f4 = hun10,
               f5 = hun5,
               f6 = hun6,
               f7 = hun7,
               f8 = hun8)
opar = par(no.readonly = T)
par(mfrow = c(2, 4),
    mar = c(2.5, 3, 2, 1), mgp = c(1.5, 0.5, 0), oma = c(0, 0, 0, 0))
for(i in 1:length(f_list1)) plot(f_list1[[i]])
par(opar)


# helper generator for a test matrix
# not CRAD basis yet
# basis expansion should be done in fit level
simmat = function(xmin = 0, xmax = 1, 
                  n = 1000, p, flist, Nk = 20, 
                  family = c("normal", "poisson", "bernoulli", "binomial"),
                  gauss_scale = NA_real_, nbinom = NA_integer_, sig2 = 0.01,
                  corr = NULL)
{
  EF_available = c("normal", "poisson", "bernoulli", "binomial")
  if(!is.list(flist)) stop("input list of functions")
  if(length(flist) < p) stop("number of functions in flist should be larger than p")
  if(Nk >= n) stop("must be Nk < n")
  if(!(family %in% EF_available)) stop("not supported glm family")
  
  # generate design matrix
  if(is.null(corr)){
    X = matrix(runif(n*p, xmin, xmax), nrow = n, ncol = p)
  } else {
    gen.gauss.cop = function(num, Sig){
      Z = matrix(rnorm(num * nrow(Sig)), nrow = nrow(Sig))
      return(pnorm(chol(Sig) %*% Z))
    }
    Sig = toeplitz(c(1, rep(corr, p-1)))
    X = t(gen.gauss.cop(n, Sig))
  }
  
  
  # generate target response
  eta = rep(0, n)
  for(i in 1:p){
    f = flist[[i]]
    eta = eta + f(X[,i]) - mean(f(X[,i])) # to make targets mean centered
  }
  
  if(family == "bernoulli") y = rbinom(n, 1, 1/(1+exp(-eta)))
  if(family == "binomial") y = rbinom(n, nbinom, 1/(1+exp(-eta)))
  if(family == "poisson") y = rpois(n, exp(eta))
  if(family == "normal") y = rnorm(n, eta, sqrt(sig2))
  
  res = data.frame(cbind(y, X))
  colnames(res) = c("y", paste("x", 1:p, sep=""))
  return(res)
}

