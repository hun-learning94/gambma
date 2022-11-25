SIMNAME = "sim1"
################################################################################
## data setting
################################################################################
set.seed(2021311165)
functionname = "hunsmooths"
source(paste0("R_miscell/", functionname, ".R"))
source("R/simmat.R")

family = "gaussian"; glmWeight= 1
if(family == "binomial"){
  N = c(100, 200); glmWeight = 10;
}else if(family == "poisson"){
  N = c(100, 200, 300);
}else if(family == "gaussian"){
  N = c(100, 200, 300);
  # N = c(200, 300)
}else if(family == "bernoulli"){
  # N = c(500, 1000, 2000); glmWeight = 1;
  N = 500
}

# TARGET = list(1,5,6)
TARGET = list(c(1,2,3))
# TARGET = list(1)
plotWidth = 1000; plotHeight = 500
Nsim_max = 100
Nsim_real = 100
np = 200

################################################################################
## fit gambma
priorCode = function(prior){
  # c("Unit", "Hyper-g", "Uniform", "Hyper-g/n", "Beta-prime", "ZS-adapted", "Robust", "Intrinsic", "constant", "CH")
  switch(prior, 
         "Unit" = "unit", "Hyper-g" = "hyperg", "Uniform" = "uniform", "Hyper-g/n" = "hypergn", "Beta-prime" = "betaprime", 
         "ZS-adapted" = "zsadapted", "Robust" = "robust", "Intrinsic" = "intrinsic", "constant" = "constant", "CH" = "ch",
         NULL)
}
################################################################################
modelSpace = function(FIT, family){
  n = nrow(FIT[[1]]$X)
  nprior = length(FIT)
  nsim = nrow(FIT[[1]]$lpydiagnosis)
  modelSpaceDf = data.frame(
    R2 = double(),
    p = double(),
    prior = character()
  )
  for(i in 1:nprior){
    modelSpaceDf = rbind(modelSpaceDf,
                         data.frame(R2 = FIT[[i]]$lpydiagnosis$r2Qm,
                                    p = rowSums(FIT[[i]]$numknots),
                                    prior = FIT[[i]]$prior))
  }
  if(family!="gaussian") modelSpaceDf$R2 = 1 - exp(-modelSpaceDf$R2/n)
  return(modelSpaceDf)
}
plotModelSpace = function(modelSpaceDf){
  PRIOR = unique(modelSpaceDf$prior)
  nprior = length(PRIOR)
  xlim = range(modelSpaceDf$R2); ylim = range(modelSpaceDf$p)
  par(mfrow=c(2,ceiling(nprior/2)))
  par(mar=c(3, 3, 2, 1), mgp=c(1.9,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)
  # plot(1, xlim=xlim, ylim=ylim, type="n", xlab="R2", ylab="p")
  for(i in seq_along(PRIOR)){
    plotdf = modelSpaceDf[modelSpaceDf$prior == PRIOR[i], ]
    # plot(plotdf$R2, plotdf$p, pch = "+", 
    #      xlim=xlim, ylim=ylim,xlab="R2", ylab="p",
    #      main = PRIOR[i])
    smoothScatter(plotdf$p~plotdf$R2, 
                  xlim=xlim, ylim=ylim,xlab="R2", ylab="p",main = PRIOR[i])
  }
}
# plotModelSpace(tmp)
################################################################################
PRIOR =  c("Unit","Uniform", "Hyper-g", "Hyper-g/n", "Beta-prime", "ZS-adapted","Robust", "Intrinsic")
maxk = 30


for(target in TARGET){
  f_list = functions[target]
  p = length(target)
  linadj = Linadj[target]
  
  ## iterate through sample size ###############################################
  for(n in N){
    datdir = paste0("results/", SIMNAME ,"/", family ,"/n", n, "/f", paste0(target, collapse = ""))
    fitted = 0
    
    for(i in 1:Nsim_max){
      fitted = fitted + 1
      if(fitted >= Nsim_real) break
      print(fitted)
      
      FIT = list()            
      # stop("watch yo jet")
      ## iterate through prior #################################################
      for(pr in seq_along(PRIOR)){
        prior = PRIOR[pr]
        resdir = paste0(datdir, "/", priorCode(prior))
        if(!file.exists(resdir)) dir.create(resdir, recursive = T)
        filename = paste0(resdir, "/dat", i)
        cat(paste0("reading in ", paste0(filename, ".rdata") , "...\n"))
        
        ## save
        # tryCatch(expr = load(paste0(filename, ".rdata")),
        #          error = function(e) print(paste0("cannot find ", paste0(filename, ".rdata"))),
        #          finally = break)
        load(paste0(filename, ".rdata"))
        FIT[[pr]] = fit
      }
      png(paste0(datdir, "/modelDiag", i,".png"), width =800, height = 400)
      plotModelSpace(modelSpace(FIT, family))
      dev.off()
    }
  }
}


