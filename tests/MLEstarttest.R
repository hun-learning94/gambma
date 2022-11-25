library(microbenchmark)
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
source("R/simmat.R")
source("R/scale01.R")
source("R_miscell/hunsmooths.R")
Rcpp::sourceCpp("src/gambma_series_diag.cpp")

mletest = function(n, nk, 
                   f_list, xrange,
                   family, link,
                   nbinom = NULL, sig = NULL,
                   prevstart = T,
                   Rglm = F,
                   iter = 10){
  
  glmfamily = switch(family,
                     "binomial" = binomial(logit),
                     "poisson" = poisson(log),
                     "normal" = gaussian(identity))
  
  p = length(f_list)
  A = matrix(0, nrow = n, ncol = 1)
  dat = simmat(f_list, xrange[1], xrange[2], n, family, link, nbinom = nbinom, sig = sig, A)
  y = dat$y
  y = matrix(y, ncol=1)
  X = as.matrix(dat[,-1])
  mins = apply(X, 2, "min")
  maxs = apply(X, 2, "max")
  X_01 = scaleto01(X, mins, maxs)
  
  knots = seq(0, 1, length = nk+2);
  knots = knots[-c(1, length(knots))]
  knotsidx = rep(1:p, each = length(knots))
  knots = rep(knots, p)
  
  ## initial value
  zP = rbinom(length(knots), 1, prob=1/2)
  knotsP = knots[zP == 1]
  knotsidxP = knotsidx[zP == 1]
  tmpP = CRAD_cpp(X_01, knotsP, knotsidxP, T, 0)
  BP = tmpP$X
  
  EtaP = etastart(nbinom, 0, y, family)
  BetaP = Rglmcpp(y, BP, EtaP, F, A, nbinom, 0, family, link)
  EtaP = BP %*% BetaP
  
  if(Rglm){
    for(i in 1:iter){
      zP = rbinom(length(knots), 1, prob=1/2)
      knotsP = knots[zP == 1]
      knotsidxP = knotsidx[zP == 1]
      tmpP = CRAD_cpp(X_01, knotsP, knotsidxP, T, 0)
      BP = tmpP$X
      BetaP = coef(glm.fit(x = BP, y = y, family = glmfamily, offset = A))
      EtaP = BP %*% BetaP
      # if(i %% 10 == 0){
      #   plot(X_01[,1], f_list[[1]](X[,1]), col = "red")
      #   points(X_01[,1], EtaP)
      # }
    }
  } else {
    for(i in 1:iter){
      zP = rbinom(length(knots), 1, prob=1/2)
      knotsP = knots[zP == 1]
      knotsidxP = knotsidx[zP == 1]
      tmpP = CRAD_cpp(X_01, knotsP, knotsidxP, T, 0)
      BP = tmpP$X
      BetaP = Rglmcpp(y, BP, EtaP, prevstart, A, nbinom, prevstart, family, link)
      EtaP = BP %*% BetaP
      # if(i %% 10 == 0){
      #   plot(X_01[,1], f_list[[1]](X[,1]), col = "red")
      #   points(X_01[,1], EtaP)
      # }
    }
  }
  
  return(0)
}

n = 1000; nk = 20; f_list = functions[1]; 
comp = microbenchmark(
  try(mletest(n, nk, f_list, xrange, "poisson", "log", nbinom = 1, sig = 0, prevstart = T, Rglm = F, iter = 50)),
  try(mletest(n, nk, f_list, xrange, "binomial", "logit", nbinom = 1, sig = 0, prevstart = T, Rglm = F, iter = 50)),
  try(mletest(n, nk, f_list, xrange, "normal", "identity", nbinom = 1, sig = 0, prevstart = T, Rglm = F, iter = 50)),
  times = 100
)
comp
autoplot(comp)
comp1= comp

n = 1000; nk = 20; f_list = functions[1]; 
family = "poisson"; link = "log"
comp = microbenchmark(
  try(mletest(n, nk, f_list, xrange, family, link, nbinom = 1, sig = 0, prevstart = T, Rglm = F, iter = 50)),
  try(mletest(n, nk, f_list, xrange, family, link, nbinom = 1, sig = 0, prevstart = F, Rglm = F, iter = 50)),
  try(mletest(n, nk, f_list, xrange, family, link, nbinom = 1, sig = 0, prevstart = F, Rglm = T, iter = 50)),
  times = 100
)
comp
autoplot(comp)
comp2= comp

n = 1000; nk = 20; f_list = functions[1]; 
family = "binomial"; link = "logit"
comp = microbenchmark(
  try(mletest(n, nk, f_list, xrange, family, link, nbinom = 1, sig = 0, prevstart = T, Rglm = F, iter = 50)),
  try(mletest(n, nk, f_list, xrange, family, link, nbinom = 1, sig = 0, prevstart = F, Rglm = F, iter = 50)),
  try(mletest(n, nk, f_list, xrange, family, link, nbinom = 1, sig = 0, prevstart = F, Rglm = T, iter = 50)),
  times = 100
)
comp
autoplot(comp)
comp3= comp












