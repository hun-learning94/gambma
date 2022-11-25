SIMNAME = "sim2"
################################################################################
## data setting
################################################################################
set.seed(2021311165)
functionname = "hunsmooths"
source(paste0("R_miscell/", functionname, ".R"))
source("R/simmat.R")

startfrom = 1
family = "poisson"; 
if(family == "binomial"){
  N = c(100, 200); glmWeight = 10;
}else if(family == "poisson"){
  N = c(100, 200, 300);
}else if(family == "gaussian"){
  N = c(100, 200, 300);
}else if(family == "bernoulli"){
  glmWeight = 1;
  N = c(1000, 1500, 2000);
}

TARGET = list(c(1,4,3))
plotWidth = 1000; plotHeight = 400
Nsim_max = 501
Nsim_real = 501
np = 200

################################################################################
## fit gambma
source("R/gambma_diag.R")
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("src/gambma.cpp")

priorCode = function(prior){
  # c("Unit", "Hyper-g", "Uniform", "Hyper-g/n", "Beta-prime", "ZS-adapted", "Robust", "Intrinsic", "constant", "CH")
  switch(prior, 
         "Unit" = "unit", "Hyper-g" = "hyperg", "Uniform" = "uniform", "Hyper-g/n" = "hypergn", "Beta-prime" = "betaprime", 
         "ZS-adapted" = "zsadapted", "Robust" = "robust", "Intrinsic" = "intrinsic", "constant" = "constant", "CH" = "ch",
         NULL)
}
PRIOR = c("Unit", "Uniform","Hyper-g", "Hyper-g/n", "Beta-prime", "ZS-adapted","Robust", "Intrinsic")
prior = "Robust"

knotConfig = "FREE" # EVEN, VS, FREE
maxk = 30
mfGen = function(p){
  mf = paste0("y~")
  for(i in 1:p){mf = paste0(mf, "ncs(x", i, ", nk=maxk)+")}
  mf = substring(mf, 1, nchar(mf)-1)
  return(formula(mf))
}
mf = mfGen(length(TARGET[[1]]))

if(knotConfig == "EVEN"){
  Ctrl = list(numMCmodels = 100,
              enumerate = T,
              burnIn = 500,
              mcIter = 1000,
              mcmcIter = 10000)
} else if (knotConfig == "VS"){
  Ctrl = list(burnIn = 500,
              mcmcIter = 2000)
} else if (knotConfig == "FREE"){
  Ctrl = list(nu = 50, 
              bir_p = 0.4, 
              dea_p = 0.4,
              initIter = 100, 
              burnIn = 500,
              mcmcIter = 2000,
              thin = maxk)
}

################################################################################
for(target in TARGET){
  f_list = functions[target]
  p = length(target)
  linadj = Linadj[target]
  
  ## iterate through sample size ###############################################
  for(n in N){
    datdir = paste0("results/", SIMNAME ,"/", family ,"/n", n, "/f", paste0(target, collapse = ""))
    fitted = 0
    
    for(i in startfrom:Nsim_max){
      fitted = fitted + 1
      if(fitted >= Nsim_real) break
      print(fitted)
      
      ## read in data ##########################################################
      dat = tryCatch(
        read.csv(paste0(datdir, "/dat/", "/dat", i, ".csv")),
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
      
      resdir = paste0(datdir, "/", knotConfig)
      if(!file.exists(resdir)) dir.create(resdir, recursive = T)
      filename = paste0(resdir, "/dat", i)
      cat(paste0("attempting ", filename, "...\n"))
      if(file.exists(paste0(filename, ".rdata"))){ 
        load(paste0(filename, ".rdata"))
        png(paste0(filename, "_trace.png"), width =plotWidth, height = plotHeight)
        plotTrace(fit)
        dev.off()
        next 
      }
      
      ## run 
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
      elapsed = difftime(toc, tic, units = c("secs"))
      if(inherits(fit, "error")){
        cat(paste0("Error occured in fitting ", knotConfig, " for dat", i), "\n")
        cat(paste0("Error message: ", fit$message), "\n")
        # stop("")
        next
      }
      ## save
      plot(fit, flist = f_list, show_plot = T, elapsed = elapsed, ylim = ylim)
      # stop("watch yo jet")
      save(fit, elapsed, file = paste0(filename, ".rdata"))
      load(paste0(filename, ".rdata"))
      png(paste0(filename, ".png"), width =plotWidth, height = plotHeight)
      tmp = plot(fit, flist = f_list, show_plot = T, elapsed = elapsed, ylim = ylim)
      dev.off()
      write.csv(tmp, file = paste0(filename, ".csv"), row.names = F)
    }
  }
}


