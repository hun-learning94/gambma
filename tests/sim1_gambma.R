SIMNAME = "sim1"
library(tidyverse)
library(ggplot2)

source("R/gambma_diag.R")
library(Rcpp)
library(RcppArmadillo)
knotConfig = "VS" # series, grid, free
maxk = 30
if(knotConfig == "EVEN"){
  evenCtrl = list(numMCmodels = 100,
                  enumerate = T,
                  burnIn = 500,
                  mcIter = 500,
                  mcmcIter = 2000)
} else if (knotConfig == "VS"){
  vsCtrl = list(burnIn = 500,
                mcmcIter = 2000)
} else if (knotConfig == "FREE"){
  freeCtrl = list(nu = 50, 
                  bir_p = 0.4, 
                  dea_p = 0.4,
                  initIter = 200, 
                  burnIn = 500,
                  mcmcIter = 2000,
                  thin = maxk)
}

Rcpp::sourceCpp("src/gambma.cpp")

## generate data
# source("R/gressani.R")
functionname = "hunsmooths"
source(paste0("R_miscell/", functionname, ".R"))
source("R/simmat.R")

family = "gaussian"; 
if(family == "binomial"){
  N = c(500, 1000, 2000);
  # N = 500
  link = "logit"; 
  glmWeight = 1
  Sig = 0;
}else if(family == "poisson"){
  # N = c(100, 200, 300);
  N=100
  link = "log"; 
  glmWeight = 0
  Sig = 0;
}else if(family == "gaussian"){
  N = c(100, 200, 300);
  # N = 300
  link = "identity"; 
  glmWeight = 0
  # Sig = c(1, 2)
  Sig = 1
}

Nsim_max = 100
Nsim_real = 100
fitted = 0
# TARGET = list(6, 5,4,1)
TARGET = list(6)
np = 200


# PRIOR = c("constant", "Hyper-g", "Uniform", "Hyper-g/n")
PRIOR = c("Beta-prime", "ZS-adapted", "Robust", "Intrinsic","constant", "Hyper-g", "Uniform", "Hyper-g/n")
# PRIOR = c("constant") 

for(target in TARGET){
  f_list = functions[target]
  p = length(target)
  linadj = Linadj[target]
  
  for(n in N){
    for(prior in PRIOR){
      for(sig in Sig){
        fitted = 0
        ## create/read directory
        datdir = paste0("results/", SIMNAME, "/" , family , "-", link, "/n", n)
        if(family == "binomial"){
          datdir = paste0(datdir, "_nbinom", glmWeight)
        } else if (family == "gaussian"){
          datdir = paste0(datdir, "_sig", sig)
        }
        resdir = datdir
        datdir = paste0(datdir, "/f",target ,"/dat")
        resdir = paste0(resdir, "/f",target , "/", prior)
        if(!file.exists(datdir)) dir.create(datdir, recursive = T)
        if(!file.exists(resdir)) dir.create(resdir, recursive = T)
        
        ## begin iteration
        for(i in 1:Nsim_max){
          print(fitted)
          
          
          filename = paste0(resdir, "/dat", i)
          cat(paste0("attempting ", filename, "...\n"))
          if(file.exists(paste0(filename, ".rdata"))){
            # load(paste0(filename, ".rdata"))
            fitted = fitted + 1
            if(fitted >= Nsim_real) break
            next
          }
          if(fitted >= Nsim_real) break
          
          ## 1. read data
          dat = tryCatch(
            read.csv(paste0(datdir, "/dat", i, ".csv")),
            error = function(cnd) cnd
          )
          if(inherits(dat, "error")){
            cat(paste0("Error occured in reading ","dat", i, ".csv"), "\n")
            cat(paste0("Error message: ", dat$message), "\n")
            next
          } else {
            X = as.matrix(dat[,-c(1)])
            mins = apply(X, 2, "min")
            maxs = apply(X, 2, "max")
            X_01 = scaleto01(X, mins, maxs)
            xgrid_01 = seq(0, 1, length = np+2); 
            xgrid_01 = xgrid_01[-c(1, length(xgrid_01))]
            Xgrid_01 = matrix(rep(xgrid_01, p), nrow=np)
            Xgrid = scalefrom01(Xgrid_01, mins, maxs)
          }
          
          ## 2. fit, plot, make PREDtable and save
          tic = Sys.time()
          fit = tryCatch(
            gambma(y ~ ncs(x1, nk = maxk), dat, 
                   knotConfig = knotConfig, 
                   prior = prior, 
                   family = family, 
                   link = link,  
                   glmWeight = glmWeight,
                   vsCtrl = vsCtrl, storeFit = F, g=n),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit, "error")){
            cat(paste0("Error occured in fitting ", knotConfig, " for dat", i), "\n")
            cat(paste0("Error message: ", fit$message), "\n")
            stop("watch yo jet")
            next
          }
          
          # tmp = plot(fit, flist = f_list, show_plot = T, elapsed = elapsed, ylim = ylim)
          # diagnosis(fit)
          # stop("watch yo jet")
          elapsed = difftime(toc, tic, units = c("secs"))
          save(fit, elapsed,
               file = paste0(filename, ".rdata"))
          load(paste0(filename, ".rdata"))
          
          png(paste0(filename, ".png"), width =1000, height = 600)
          tmp = plot(fit, flist = f_list, show_plot = T, elapsed = elapsed, ylim = ylim)
          tmp$knotConfig = paste0(knotConfig)
          dev.off()
          
          
          
          write.csv(tmp, file = paste0(filename, ".csv"), row.names = F)
          
          ## 3. fit until Nsim_real
          fitted = fitted + 1
          if(fitted >= Nsim_real) break
        }
      }
    }
  }
}