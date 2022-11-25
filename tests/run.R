################################################################################
## data generation
################################################################################
set.seed(2021311165)
functionname = "hunsmooths"
source(paste0("R_miscell/", functionname, ".R"))
source("R/simmat.R")

family = "poisson"; 
glmWeight = 1

n = 200
target = c(1,3,4)
f_list = functions[target]
p = length(target)
linadj = Linadj[target]

dat = simmat(f_list, xrange[1], xrange[2], n = n, family = family, glmWeight = glmWeight,
             beta = c(-1,1))


################################################################################
## fit gambma
################################################################################
source("R/gambma_diag.R")
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("src/gambma.cpp")
PRIOR = c("Beta-prime", "ZS-adapted", "Robust", "Intrinsic","constant", "Hyper-g", "Uniform", "Hyper-g/n")

knotConfig = "EVEN" # series, grid, free
maxk = 20
prior = "Robust"

if(knotConfig == "EVEN"){
  Ctrl = list(numMCmodels = 100,
                  enumerate = F,
                  burnIn = 500,
                  mcIter = 500,
                  mcmcIter = 5000)
} else if (knotConfig == "VS"){
  Ctrl = list(burnIn = 500,
                mcmcIter = 2000)
} else if (knotConfig == "FREE"){
  Ctrl = list(nu = 50, 
                  bir_p = 0.4, 
                  dea_p = 0.4,
                  initIter = 200, 
                  burnIn = 500,
                  mcmcIter = 2000,
                  thin = maxk)
}

mf = y~ lin_x1 + lin_x2 +
  ncs(x1, nk=maxk) + ncs(x2, nk = maxk) + ncs(x3, nk = maxk)
# mf = y~ncs(x1, nk=maxk)

tic = Sys.time()
fit = tryCatch(
  gambma(mf, dat, 
         knotConfig = knotConfig, 
         prior = prior, 
         family = family, 
         glmWeight = glmWeight,
         Ctrl = Ctrl),
  error = function(cnd)cnd
)
toc = Sys.time()

summary(fit)
tmp = plot(fit, flist = f_list, show_plot = T, ylim = ylim)
bayesResiduals(fit)

