SIMNAME = "sim2"
library(tidyverse)
library(ggplot2)

source("R/gambma_diag.R")
source("R/t_col.R")

METHODS = c("Free-knot","Even-knot","VS-knot")
METHODS_col = c("#c7e9b4", "#7fcdbb", "#41b6c4")
methodCode = Vectorize(function(prior){
  switch(prior, 
         "VS-knot"="VS", "Even-knot"="EVEN", "Free-knot" = "FREE",
         "BayesX"="bayesx", "Blapsr"="blapsr", "Mgcv-ps"="mgcvps", "Mgcv-ad"="mgcvad",
         NULL)
}, vectorize.args = "prior")
num_methods = length(METHODS)

################################################################################
## data setting
################################################################################
set.seed(2021311165)
functionname = "hunsmooths"
source(paste0("R_miscell/", functionname, ".R"))
source("R/simmat.R")

family = "poisson"; glmWeight= 1
if(family == "binomial"){
  N = c(100, 200); glmWeight = 10;
}else if(family == "poisson"){
  N = c(100, 200, 300);
}else if(family == "gaussian"){
  N = c(100, 200, 300);
  METHODS = c("Even-knot","VS-knot", "Free-knot")
  METHODS_col = c("#c7e9b4", "#7fcdbb", "#41b6c4")
  num_methods = length(METHODS)
}else if(family == "bernoulli"){
  N = c(1000, 1500, 2000); glmWeight = 1;
}

# TARGET = list(1,2,3)
TARGET = list(c(1,4,3))
plotWidth = 1000; plotHeight = 500
Nsim_max = 500+1
Nsim_plot = 500+1

np = 200
xgrid = seq(xrange[1], xrange[2], length = np+2)
xgrid = xgrid[-c(1, np+2)]


ESS_all = matrix(0, nrow = length(N)*(Nsim_max-1), ncol = length(METHODS))
fuck = 0
for(target in TARGET){
  cnt = 0
  for(n in N){
    cnt = 0
    datdir = paste0("results/", SIMNAME ,"/", family ,"/n", n, "/f", paste0(target, collapse = ""))
    for(i in 1:Nsim_max){
      if(cnt > Nsim_plot) break
      if(i %% 10 == 0) cat(paste(datdir, i, ", cnt", cnt, "\n"))
      if(!all(file.exists(paste0(datdir, "/", methodCode(METHODS), "/dat", i, ".csv")))) next # all fitted idx
      cnt = cnt + 1; 
      
      fuck = fuck+1
      for(ffuck in seq_along(METHODS)){
        load(paste0(datdir, "/", methodCode(METHODS[ffuck]), "/dat", i, ".rdata"))
        ESS_all[fuck, ffuck] = plotTrace(fit)/as.numeric(elapsed)
        stop("watch yo jet")
      }
    }
  }
}

colnames(ESS_all) = METHODS
par(mfrow=c(1,3))
boxplot(ESS_all[1:500,2:3])
boxplot(ESS_all[1:500+500,2:3])
boxplot(ESS_all[1:500+1000,2:3])

write.csv(ESS_all, row.names = F, 
          file = paste0(datdir, "/ESS.csv"))
