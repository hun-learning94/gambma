library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
source("R/simmat.R")
source("R/scale01.R")
source("R_miscell/hunsmooths.R")
library(ggplot2)
Rcpp::sourceCpp("src/gambma_series_diag.cpp")

mletest = function(n, nk, 
                   f_list, xrange,
                   family, link,
                   nbinom = NULL, sig = NULL,
                   cpp = F){
  p = length(f_list)
  A = matrix(0, nrow = n, ncol = 1)
  dat = simmat(f_list, xrange[1], xrange[2], n, family, link, nbinom = nbinom, sig = sig, A)
  y = dat$y
  X = as.matrix(dat[,-1])
  mins = apply(X, 2, "min")
  maxs = apply(X, 2, "max")
  X_01 = scaleto01(X, mins, maxs)
  
  knots = seq(0, 1, length = nk+2);
  knots = knots[-c(1, length(knots))]
  knotsidx = rep(1:p, each = length(knots))
  knots = rep(knots, p)
  
  tmp = CRAD_cpp(X_01, knots, knotsidx, T, 0)
  
  if (family == "binomial") {
    if(link == "logit"){
      Rglm = function(y, X, A = NULL, nbinom){
        coef(glm.fit(x = X[,-1], y = y, family=binomial(logit), offset = A))
      }
    }
  }
  if (family == "poisson"){
    if(link == "log"){
      Rglm = function(y, X, A = NULL, nbinom){
        coef(glm.fit(x = X[,-1], y = y, family=poisson(log), offset = A))
      } 
    }
  }
  if (family == "normal"){
    if(link == "identity"){
      Rglm = function(y, X, A = NULL, nbinom){
        solve(crossprod(X), t(X) %*% (y-A))
      } 
    }
  }
  # 
  # # print(c(Rglmcpp(matrix(y, ncol=1), tmp$X, A,nbinom, 1, family, link)))
  # # print(c(Rglm(y, tmp$X, NULL, nbinom)))
  # nearPDres = function(x) as.matrix(Matrix::nearPD(x, doSym=T)$mat)
  # 
  # res1 = Rglmcpp(matrix(y, ncol=1), tmp$X, A, nbinom, 1, family, link, nearPDres)
  # res2 = Rglm(y, tmp$X, A, nbinom)
  # print(sqrt(mean((abs((tmp$X %*% (res1 - res2)))))^2))
  # print(c(res1))
  # print(c(res2))
  # 
  if(cpp){
    res = Rglmcpp(matrix(y, ncol=1), tmp$X, matrix(0, length(y), 1), F,
                  A, nbinom, 1, family, link)
  } else {
    res = Rglm(y, tmp$X, A, nbinom)
  }
  res
}

n = 1000; nk = 20; f_list = functions[1:3]; 
comp = microbenchmark(
  mletest(n, nk, f_list, xrange,
          "normal", "identity", nbinom = 1, sig = 1),
  mletest(n, nk, f_list, xrange,
          "normal", "identity", nbinom = 1, sig = 1, cpp=T),
  mletest(n, nk, f_list, xrange,
          "binomial", "logit", nbinom = 1, sig = 0),
  mletest(n, nk, f_list, xrange,
          "binomial", "logit", nbinom = 1, sig = 0, cpp=T),
  mletest(n, nk, f_list, xrange,
          "poisson", "log", nbinom = 1, sig = 0),
  mletest(n, nk, f_list, xrange,
          "poisson", "log", nbinom = 1, sig = 0, cpp=T),
  times = 200
)
comp
autoplot(comp)

n = 200; nk = 30; f_list = functions[1:3]; 
comp = microbenchmark(
  mletest(n, nk, f_list, xrange,
          "normal", "identity", nbinom = 1, sig = 1),
  mletest(n, nk, f_list, xrange,
          "normal", "identity", nbinom = 1, sig = 1, cpp=T),
  times = 100
)
comp
autoplot(comp)


n = 1000; nk = 20; f_list = functions[1:3]; 
comp = microbenchmark(
  try(mletest(n, nk, f_list, xrange,
            "binomial", "logit", nbinom = 1, sig = 0)),
  try(mletest(n, nk, f_list, xrange,
            "binomial", "logit", nbinom = 1, sig = 0, cpp=T)),
  times = 100
)
comp
autoplot(comp)

n = 200; nk = 30; f_list = functions[1:3]; 
comp = microbenchmark(
  try(mletest(n, nk, f_list, xrange,
              "poisson", "log", nbinom = 1, sig = 0)),
  try(mletest(n, nk, f_list, xrange,
              "poisson", "log", nbinom = 1, sig = 0, cpp=T)),
  times = 100
)
comp
autoplot(comp)
