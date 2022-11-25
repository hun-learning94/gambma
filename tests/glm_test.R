set.seed(2021311165)
library(tidyverse)
library(ggplot2)
source("R/simmat_glm.R")
source("R/glmbma_diag.R")
library(Rcpp)
library(RcppArmadillo)

MCMCiter = 2500
burnin = 500
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
sourceCpp(paste0("src/glmbma_diag.cpp"))

family = "binomial"; 
if(family == "binomial"){
  # N = c(1000, 2000);
  N = 1000
  link = "logit"; 
  nbinom = 1
  Sig = 0;
}else if(family == "poisson"){
  # N = c(150, 200, 300);
  N = 150
  link = "log"; 
  nbinom = 0
  Sig = 0;
}else if(family == "gaussian"){
  N = c(200, 400)
  # N = 200
  link = "identity"; 
  nbinom = 0
  Sig = c(1, 2)
}
intercept = ifelse(family == "poisson", 0, 0)

xrange = c(-1, 1)
n = 500
sig = 1
p = 100
truep = 25
nullp = p - truep
betamax = 0.75
Beta = c(rep(0, nullp), seq(-betamax, betamax, len = p - nullp))
dat = simmat_glm(xrange[1], xrange[2], n, family, link, nbinom = nbinom, sig = sig,
                 beta = Beta, intercept = 1)

res1 = glmbma(y ~ .,
             dat = dat,
             prior = "fixedg",
             g_manual = n,
             family = family, 
             nbinom = 1,
             link = link,
             MCMCiter = MCMCiter,
             burnin = burnin)

res = res1
par(mfrow=c(1,2))
colMeans(res$Z) %>% plot(ylim = c(0,1))
points(res$Z[which.max(res$lpydiagnosis$lpy), ], col = "blue", pch = 3)
abline(v = nullp+0.5)
title(paste0("MAP dim = ", sum(res$Z[which.max(res$lpydiagnosis$lpy), ]), ", wrong = ", sum(abs(res$Z[which.max(res$lpydiagnosis$lpy), ] - as.numeric(Beta != 0))), " out of ", length(Beta)))

colMeans(res$BetaSamp) %>% plot()
abline(v = nullp+0.5)
lines(Beta, col = "red", lty=2)
title(paste0("rmse = ", signif(sqrt(mean((colMeans(res$BetaSamp) - Beta)^2)), 3)))
mtext("fixedg=n", outer  = T, font = 2, cex = 1.5)

res3 = glmbma(y ~ .,
              dat = dat,
              prior = "fixedg",
              g_manual = n/3,
              family = family, 
              nbinom = 1,
              link = link,
              MCMCiter = MCMCiter,
              burnin = burnin)

res = res3
par(mfrow=c(1,2))
colMeans(res$Z) %>% plot(ylim = c(0,1))
points(res$Z[which.max(res$lpydiagnosis$lpy), ], col = "blue", pch = 3)
abline(v = nullp+0.5)
title(paste0("MAP dim = ", sum(res$Z[which.max(res$lpydiagnosis$lpy), ]), ", wrong = ", sum(abs(res$Z[which.max(res$lpydiagnosis$lpy), ] - as.numeric(Beta != 0))), " out of ", length(Beta)))

colMeans(res$BetaSamp) %>% plot()
abline(v = nullp+0.5)
lines(Beta, col = "red", lty=2)
title(paste0("rmse = ", signif(sqrt(mean((colMeans(res$BetaSamp) - Beta)^2)), 3)))
mtext("fixedg=n/3", outer  = T, font = 2, cex = 1.5)

res2 = glmbma(y ~ .,
              dat = dat,
              prior = "betaprime",
              g_manual = n,
              family = family, 
              nbinom = 1,
              link = link,
              MCMCiter = MCMCiter,
              burnin = burnin)

res = res2
par(mfrow=c(1,2))
colMeans(res$Z) %>% plot(ylim = c(0,1))
points(res$Z[which.max(res$lpydiagnosis$lpy), ], col = "blue", pch = 3)
abline(v = nullp+0.5)
title(paste0("MAP dim = ", sum(res$Z[which.max(res$lpydiagnosis$lpy), ]), ", wrong = ", sum(abs(res$Z[which.max(res$lpydiagnosis$lpy), ] - as.numeric(Beta != 0))), " out of ", length(Beta)))

colMeans(res$BetaSamp) %>% plot()
abline(v = nullp+0.5)
lines(Beta, col = "red", lty=2)
title(paste0("rmse = ", signif(sqrt(mean((colMeans(res$BetaSamp) - Beta)^2)), 3)))
mtext("betaprime", outer  = T, font = 2, cex = 1.5)

