SIMNAME = "sim2"

set.seed(2021311165)
library(tidyverse)
library(ggplot2)

## generate data
# source("R/gressani.R")
functionname = "hunsmooths"
source(paste0("R_miscell/", functionname, ".R"))
source("R/simmat.R")
glmWeight = 1
family = "bernoulli"; 
if(family == "binomial"){
  N = c(100, 200); glmWeight = 10;
}else if(family == "poisson"){
  N = c(100, 200, 300);
}else if(family == "gaussian"){
  N = c(100, 200, 300);
}else if(family == "bernoulli"){
  N = c(500, 1000, 1500, 2000); glmWeight = 1;
}

Nsim_max = 600
Nsim_real = 600
# TARGET = list(1,2,3,4)
# TARGET =list(4)
TARGET = list(c(1,4,3))
for(target in TARGET){
  f_list = functions[target]
  p = length(target)
  linadj = Linadj[target]
  
  for(n in N){
    ## create directory
    datdir = paste0("results/", SIMNAME ,"/", family ,"/n", n)
    datdir = paste0(datdir, "/f", paste0(target, collapse = "") ,"/dat")
    if(!file.exists(datdir)) dir.create(datdir, recursive = T)
    
    for(i in 1:Nsim_max){
      dat = simmat(f_list, xrange[1], xrange[2], n, family, glmWeight = glmWeight)
      if(!file.exists(paste0(datdir, "/dat", i, ".csv"))) 
        write.csv(dat, file = paste0(datdir, "/dat", i, ".csv"), row.names = F)
      dat = read.csv(paste0(datdir, "/dat", i, ".csv"))
    }
  }
  
}


