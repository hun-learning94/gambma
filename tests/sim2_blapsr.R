SIMNAME = "sim2"
library(tidyverse)
library(ggplot2)

method = "blapsr"
maxk = 30
source("R/gambma_diag.R")
library(blapsr)
source("R_miscell/blapsr/gamlps_plot.R")

## generate data
# source("R/gressani.R")
functionname = "hunsmooths"
source(paste0("R_miscell/", functionname, ".R"))
source("R/simmat.R")

family = "poisson"; 
if(family == "gaussian"){
  N = c(100, 200); glmWeight = 10;
}else if(family == "poisson"){
  N = c(100, 200, 300);
}else if(family == "gaussian"){
  N = c(100, 200, 300);
}else if(family == "bernoulli"){
  glmWeight = 1;
  N = c(1000, 1500, 2000);
}

blapsr_family = family
if(family == "binomial" & glmWeight == 1) blapsr_family = "bernoulli"
blapsr_failed = FALSE

Nsim_max = 501
Nsim_real = 501
fitted = 0
target = c(1,4,3)
f_list = functions[target]
p = length(target)
linadj = Linadj[target]
np = 200

for(n in N){
  fitted = 0
  ## create/read directory
  datdir = paste0("results/", SIMNAME ,"/", family ,"/n", n, "/f", paste0(target, collapse = ""))
  fitted = 0
  resdir = datdir
  datdir = paste0(datdir, "/dat")
  resdir = paste0(resdir, "/", method)
  if(!file.exists(datdir)) dir.create(datdir, recursive = T)
  if(!file.exists(resdir)) dir.create(resdir, recursive = T)
  
  ## begin iteration
  for(i in 1:Nsim_max){
    filename = paste0(resdir, "/dat", i)
    cat(paste0("attempting ", filename, "...\n"))
    # stop("watch yo jet")
    if(file.exists(paste0(filename, ".rdata"))){
      fitted = fitted + 1
      next
    }
    
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
    for(j in 1:1){
      tic = Sys.time()
      fit = tryCatch(
        gamlps(y ~ 1+sm(x1)+sm(x2)+sm(x3),
               data = dat, K = maxk, penorder = 3 - 1, family=blapsr_family),
        error = function(cnd)cnd
      )
      toc = Sys.time()
      if(inherits(fit, "error")){
        cat(paste0("Error occured in fitting ", method, " for dat", i), "\n")
        cat(paste0("Error message: ", fit$message), "\n")
        next
      } else {
        elapsed = difftime(toc, tic, units = c("secs"))
        save(fit, elapsed,
             file = paste0(filename, ".rdata"))
        load(paste0(filename, ".rdata"))
        png(paste0(filename, ".png"), width =1000, height = 600)
        tmp = plot_gamlps(fit, Xgrid, f_list, 0.95, elapsed = elapsed, ylim = ylim, Linadj = linadj)
        tmp$method = paste0(method)
        dev.off()
        
        write.csv(tmp, file = paste0(filename, ".csv"), row.names = F)
        fitted = fitted + 1
        break
      }
    }
    
    ## 3. fit until Nsim_real
    if(fitted >= Nsim_real) break
  }
  
}




