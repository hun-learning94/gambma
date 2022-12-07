SIMNAME = "sim1"
library(tidyverse)
library(ggplot2)

source("R/gambma_diag.R")
source("R/t_col.R")

METHODS = c("uniform", "hyperg", "benchmark", "zsadapted", "betaprime", "fixedg")
METHODS_plot = c("Uniform", "Hyperg", "Benchmark", "ZSadapted", "Betaprime","Fixed-g")
METHODS_col = c(t_col("palevioletred1"),"red",   t_col("steelblue1"), t_col("turquoise3"),  "blue","olivedrab")
num_methods = length(METHODS)

library(Rcpp)
library(RcppArmadillo)

## generate data
functionname = "hunsmooths"
source(paste0("R_miscell/", functionname, ".R"))
source("R/simmat.R")

family = "binomial"; 
if(family == "binomial"){
  N = c(500, 1000, 2000);
  # N = 500
  link = "logit"; 
  nbinom = 1
  Sig = 0;
}else if(family == "poisson"){
  N = c(100, 200, 300);
  link = "log"; 
  nbinom = 0
  Sig = 0;
}else if(family == "gaussian"){
  N = c(100, 200, 300);
  link = "identity"; 
  nbinom = 0
  # Sig = c(1, 2)
  Sig=1
}

Nsim_max = 600
Nsim_plot = 500
TARGET = c(1,2,3,4,5,6)
np = 200
xgrid = seq(xrange[1], xrange[2], length = np+2)
xgrid = xgrid[-c(1, np+2)]
PREDTABLE_all = data.frame(
  f = integer(),
  grid = double(),
  mean = double(),
  lb = double(),
  ub = double(),
  true = double(),
  cover = integer(),
  method = character(),
  id = integer(),
  n = double()
)
MSE_all = data.frame(
  f = integer(), 
  method = character(), 
  rmse= double(), 
  sim_iter = integer(),
  n = double()
)
COVER_all = data.frame(
  f = integer(), 
  method =character(), 
  cover =  double(),
  n = double()
)

for(target in TARGET){
  
  cnt = 0
  f_list = functions[target]
  p = length(target)
  linadj = Linadj[target]
  
  # EXISTS= list()
  # exist = 1
  for(n in N){
    for(sig in Sig){
      cnt = 0
      ## create/read directory
      datdir = paste0("results/", SIMNAME, "/" , family , "-", link, "/n", n)
      if(family == "binomial"){
        datdir = paste0(datdir, "_nbinom", nbinom)
      } else if (family == "gaussian"){
        datdir = paste0(datdir, "_sig", sig)
      }
      resdir = datdir
      datdir = paste0(datdir, "/f",target ,"/dat")
      resdir = paste0(resdir, "/f",target)
      
      PREDTABLE = data.frame(
        f = integer(),
        grid = double(),
        mean = double(),
        lb = double(),
        ub = double(),
        true = double(),
        cover = integer(),
        method = character(),
        id = integer()
      )
      MSE = data.frame(
        f = integer(), 
        method = character(), 
        rmse= double(), 
        sim_iter = integer()
      )
      COVER = data.frame(
        f = rep(rep(1:p, each = np), num_methods),
        method = rep(METHODS, each = np*p),
        cover = integer(np*p*num_methods)
      )
      
      for(i in 1:Nsim_max){
        if(cnt > Nsim_plot) break
        if(i %% 50 == 0) cat(paste(resdir, i, ", cnt", cnt, "\n"))
        # EXIST[i, ] = file.exists(paste0(resdir, "/", METHODS, "/dat", i, ".csv"))
        # if(!all(EXIST[i, ])) next # all fitted idx
        if(!all(file.exists(paste0(resdir, "/", METHODS, "/dat", i, ".csv")))) next # all fitted idx
        cnt = cnt + 1; 
        
        PREDtable = data.frame(
          f = integer(),grid = double(),mean = double(),
          lb = double(),ub = double(),
          true = double(),cover = integer(),
          method = character(),id = integer()
        )
        for(method in METHODS){
          predtable = read.csv(paste0(resdir, "/", method, "/dat", i, ".csv"))
          predtable$method = method
          PREDtable = rbind(PREDtable, predtable)
        }
        PREDtable$f = target
        PREDtable$id = i
        PREDtable$n = n
        PREDTABLE = rbind(PREDTABLE, PREDtable)
        
        ## MSE
        MSE_tmp = PREDtable %>% 
          mutate(sq_diff = (true-mean)^2) %>% 
          select(f, method, sq_diff) %>% 
          group_by(f, method) %>% summarise(logrmse = log(sqrt(mean(sq_diff))), .groups = "drop")
        MSE_tmp$sim_iter = cnt
        MSE_tmp$f = target
        MSE_tmp$n = n
        MSE = rbind(MSE, MSE_tmp)
        
        ## COVER
        COVER$cover = COVER$cover + PREDtable$cover
      }
      # EXISTS[[exist]] = na.omit(EXIST)
      # exist = exist+1
      COVER$cover = COVER$cover / (cnt-1)
      COVER$x = xgrid
    }
    # PREDTABLE$n = n
    COVER$n = n
    COVER$f = target
    # MSE$n = n
    PREDTABLE_all = rbind(PREDTABLE_all, PREDTABLE)
    COVER_all = rbind(COVER_all, COVER)
    MSE_all = rbind(MSE_all, MSE)
  }
}


## save results
if(family!= "gaussian"){
  PATH_PREDTABLE_all = paste0("results/", SIMNAME, "/" , family , "-", link,  "/f", paste0(TARGET, collapse = ""), "_PREDTABLE.rds")
  PATH_COVER = paste0("results/", SIMNAME, "/" , family , "-", link, "/f", paste0(TARGET, collapse = ""), "_COVER.csv")
  PATH_MSE = paste0("results/", SIMNAME, "/" , family , "-", link, "/f", paste0(TARGET, collapse = ""), "_MSE.csv")
  PATH_numfailed = paste0("results/", SIMNAME, "/" , family , "-", link, "/f", paste0(TARGET, collapse = ""), "_numfailed.csv")
  
  saveRDS(PREDTABLE_all, file = PATH_PREDTABLE_all)
  write.csv(COVER_all, file = PATH_COVER, row.names=F)
  write.csv(MSE_all, file = PATH_MSE, row.names=F)
  # write.csv(numfailed, file = PATH_numfailed)
  
  PREDTABLE_all = readRDS(file = PATH_PREDTABLE_all)
  COVER_all = read.csv(PATH_COVER)
  MSE_all = read.csv(PATH_MSE)
  # numfailed = read.csv(PATH_numfailed)
  
  PLOTNAME = paste0("results/", SIMNAME, "/" , family , "-", link, "/", family , "-", link,"_f", paste0(TARGET, collapse = ""), ".pdf")
  DIAGNAME = paste0("results/", SIMNAME, "/" , family , "-", link, "/", family , "-", link,"_f", paste0(TARGET, collapse = ""), "_fit.pdf")
  dropbox = "C:/Users/USER/Dropbox/Kang/GAM3/plots/"
  PLOTNAMEdropbox = paste0(dropbox, family , "-", link,"_f", paste0(TARGET, collapse = ""), ".pdf")
  DIAGNAMEdropbox = paste0(dropbox, family , "-", link,"_f", paste0(TARGET, collapse = ""), "_fit.pdf")
} else {
  PATH_PREDTABLE_all = paste0("results/", SIMNAME, "/" , family , "-", link,  "/f", paste0(TARGET, collapse = ""),"_sig", Sig, "_PREDTABLE.rds")
  PATH_COVER = paste0("results/", SIMNAME, "/" , family , "-", link, "/f", paste0(TARGET, collapse = ""),"_sig", Sig, "_COVER.csv")
  PATH_MSE = paste0("results/", SIMNAME, "/" , family , "-", link, "/f", paste0(TARGET, collapse = ""),"_sig", Sig, "_MSE.csv")
  PATH_numfailed = paste0("results/", SIMNAME, "/" , family , "-", link, "/f", paste0(TARGET, collapse = ""), "_sig", Sig,"_numfailed.csv")
  
  # saveRDS(PREDTABLE_all, file = PATH_PREDTABLE_all)
  # write.csv(COVER_all, file = PATH_COVER, row.names=F)
  # write.csv(MSE_all, file = PATH_MSE, row.names=F)
  # write.csv(numfailed, file = PATH_numfailed)
  
  PREDTABLE_all = readRDS(file = PATH_PREDTABLE_all)
  COVER_all = read.csv(PATH_COVER)
  MSE_all = read.csv(PATH_MSE)
  # numfailed = read.csv(PATH_numfailed)
  
  PLOTNAME = paste0("results/", SIMNAME, "/" , family , "-", link, "/", family , "-", link,"_f", paste0(TARGET, collapse = ""), "_sig", Sig,".pdf")
  DIAGNAME = paste0("results/", SIMNAME, "/" , family , "-", link, "/", family , "-", link,"_f", paste0(TARGET, collapse = ""), "_sig", Sig,"_fit.pdf")
  dropbox = "C:/Users/USER/Dropbox/Kang/GAM3/plots/"
  PLOTNAMEdropbox = paste0(dropbox, family , "-", link,"_f", paste0(TARGET, collapse = ""), "_sig", Sig,".pdf")
  DIAGNSMEdropbox = paste0(dropbox, family , "-", link,"_f", paste0(TARGET, collapse = ""), "_sig", Sig,"_fit.pdf")
}

## diagnose result
# pdf(DIAGNAMEdropbox, width = 3*length(comp_methods), height = 2.2*length(TARGET))
pdf(DIAGNAME, width = 3*length(comp_methods), height = 2.2*length(TARGET))
comp_methods = c("hyperg", "betaprime", "fixedg")
plotlabels = function(x){
  switch(x, hyperg = "Hyperg", betaprime = "Betaprime", fixedg = "Fixed g = n")
}
target = 1; n = N[1]; maxlines = 500;
par(mar=c(1.5, 2, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,1,0), cex.axis = 1.0, cex.main= 1)
par(mfrow=c(length(TARGET),length(comp_methods)))
for(target in TARGET){
  for(method in comp_methods){
    plot(xgrid, functions[[target]](xgrid), type= "n", ylab = "", xlab = "", ylim = c(ylim[1]-1, ylim[2]+1),
         main = plotlabels(method))
    plotdf = PREDTABLE_all[PREDTABLE_all$n == n & PREDTABLE_all$method == method & PREDTABLE_all$f == target, ]
    meancurve = rep(0, length(xgrid))
    medcurve = matrix(NA, nrow = max(plotdf$id), ncol = length(xgrid))
    fitted = 0
    for(i in 1:max(plotdf$id)){
      if(sum(plotdf$id == i) == 0 ) next
      meancurve = meancurve + plotdf[(plotdf$id == i), "mean"]
      medcurve[i, ] = plotdf[(plotdf$id == i), "mean"]
      fitted = fitted + 1
      if(fitted < maxlines) {
        lines(xgrid, plotdf[(plotdf$id == i), "mean"],
              lty = 1, col = t_col("grey"), lwd = 0.25)
      } 
    }
    lines(xgrid, functions[[target]](xgrid), col = "black", lwd = 0.5)
    # lines(xgrid, meancurve / fitted, col = "blue", lty = 2, lwd = 1.2)
    lines(xgrid, apply(medcurve, 2, function(x) median(x, na.rm=T)), col = "blue", lty=2, lwd=1.2)
  } 
}
dev.off()

# ## diagnose result (PPT)
# comp_methods = c("hyperg", "betaprime", "fixedg")
# plotlabels = function(x){
#   switch(x, hyperg = "p(g), Hyperg", betaprime = "p(g;n), Betaprime", fixedg = "g = n")
# }
# n = N[1]; maxlines = 200; 
# for(targett in c(1,5,4)){
#   pdf(paste0("C:/Users/USER/Naver MYBOX/hun_workstation/Papers/BVSGAM/docs/reports/DM0827/plots/", family , "-", link,"_f", targett, "_fit.pdf"), 
#       width = 2.5*length(comp_methods), height = 3)
#   par(mar=c(1.5, 2, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)
#   par(mfrow=c(1,length(comp_methods)))
#   for(target in targett){
#     for(method in comp_methods){
#       plot(xgrid, functions[[target]](xgrid), type= "n", ylab = "", xlab = "", ylim = c(ylim[1]-1, ylim[2]+1),
#            main = plotlabels(method))
#       plotdf = PREDTABLE_all[PREDTABLE_all$n == n & PREDTABLE_all$method == method & PREDTABLE_all$f == target, ]
#       meancurve = rep(0, length(xgrid))
#       medcurve = matrix(NA, nrow = max(plotdf$id), ncol = length(xgrid))
#       fitted = 0
#       for(i in 1:max(plotdf$id)){
#         if(sum(plotdf$id == i) == 0 ) next
#         meancurve = meancurve + plotdf[(plotdf$id == i), "mean"]
#         medcurve[i, ] = plotdf[(plotdf$id == i), "mean"]
#         fitted = fitted + 1
#         if(fitted < maxlines) {
#           lines(xgrid, plotdf[(plotdf$id == i), "mean"],
#                 lty = 1, col = t_col("grey"), lwd = 0.5)
#         } 
#       }
#       lines(xgrid, functions[[target]](xgrid), col = "black", lwd = 0.5)
#       lines(xgrid, meancurve / fitted, col = "blue", lty = 2, lwd = 0.5)
#       # lines(xgrid, apply(medcurve, 2, function(x) median(x, na.rm=T)), col = "blue", lty=2, lwd=0.5)
#     } 
#   }
#   dev.off()
# }

## plot results
pdf(PLOTNAMEdropbox, width = 10*1.1, height = length(TARGET)*1.9)
par(mar=c(1.5, 2, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,1,0), cex.axis = 1.0, cex.main= 1)
layout(matrix(rep(1:(p*(1+length(N))), rep(c(length(N)-1,rep(1, length(N))), p)), ncol=(length(N)+2), byrow=T))
p = length(TARGET)
plotdf = MSE_all
plotdf$method = paste0(plotdf$method, plotdf$n)
plotat = c(seq_along(METHODS), 
           seq_along(METHODS) + length(METHODS) + 1, 
           seq_along(METHODS) + length(METHODS)*2 + 2)
ylim2= c(0.45, 0.999)
ytick = c(0.5, 0.8, 0.9, 0.95, 0.99)
coverf = function(x) 1/(1-x/1.1)-1
plotdf2 = COVER_all %>% mutate(cover = coverf(cover))
plotdf2 = plotdf2 %>% pivot_wider(names_from = method, values_from = cover)
plotdf2$x = round(plotdf2$x, 2)
gridpoints = round(seq(-0.95, 0.95, length = 11), 2)
for(pp in TARGET){
  # plot(xgrid, f_list[[pp]](xgrid), type= "l", ylab = "", xlab = "", ylim = ylim, main = paste0("true f", pp))
  plotdf %>% 
    pivot_wider(names_from = method, values_from = logrmse) %>% ungroup %>% filter(f==pp) %>%
    select(any_of(c(outer(METHODS, N, paste0)))) %>% 
    boxplot(outline = F, col = METHODS_col, border = "black", main = paste0("logrmse"),
            xaxt = "n",
            at = plotat)
  axis(1, at = c(0, 1+length(METHODS), 2+2*length(METHODS))+4, labels = paste0("n=", N))
  abline(v = c(0, 1+length(METHODS))+(length(METHODS)+1), lty=2, col = "grey")
  legend("topright", col = METHODS_col, inset = c(-0.006, -0.01), bty="n",
         lty=1, lwd = 3, legend = METHODS_plot, ncol=2, cex = 0.8)
  
  for(nn in N){
    plotdf2 %>% filter(x %in% gridpoints) %>% filter(f == pp, n == nn) %>% distinct(x, .keep_all = T) %>% 
      select(any_of(METHODS)) %>% matplot(type="b", pch = 1, cex = 0.7, col = (METHODS_col), lty=1, lwd=1.2, 
                                          ylim = coverf(ylim2), x = gridpoints, main = paste0("n=", nn), xlab="", yaxt="n")
    axis(side=2, at = coverf(ytick), labels = ytick)
    abline(h = coverf(0.95), lty=2, col = "gray40")
  }
  # legend("bottomright", col = METHODS_col, lty=1, pch = 1, lwd = 1.2, legend = METHODS_plot, ncol=4)
}
dev.off()












