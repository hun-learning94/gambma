SIMNAME = "sim2"
library(tidyverse)
library(ggplot2)

source("R/gambma_diag.R")
source("R/t_col.R")

METHODS = c("series", "grid", "free", "bayesx", "blapsr", "mgcvps", "mgcvad")
METHODS_plot = c("Even-knot", "VS-knot", "Free-knot", "BayesX", "Blapsr","Mgcv-ps", "Mgcv-ad")
METHODS_col = c(t_col("steelblue1"), "blue", t_col("turquoise3"), "olivedrab",t_col("lightgreen"), "red", t_col("palevioletred1"))
# METHODS = c("series", "grid")
num_methods = length(METHODS)

library(Rcpp)
library(RcppArmadillo)
library(R2BayesX)
source("R_miscell/r2bayesX/R2BayesX_plot.R")
library(mgcv)
source("R_miscell/mgcv/mgcv_plot.R")
library(blapsr)
source("R_miscell/blapsr/gamlps_plot.R")

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
  Sig = 1
}

if(family == "gaussian"){
  METHODS = c("series", "grid", "free", "bayesx", "mgcvps", "mgcvad")
  METHODS_plot = c("Even-knot", "VS-knot", "Free-knot", "BayesX", "Mgcv-ps", "Mgcv-ad")
  METHODS_col = c(t_col("steelblue1"), "blue", t_col("turquoise3"), "olivedrab", "red", t_col("palevioletred1"))
  num_methods = length(METHODS)
}

Nsim_max = 1500
Nsim_plot = 500
cnt = 0
TARGET= target = c(1,2,3,6)
f_list = functions[target]
p = length(target)
linadj = Linadj[target]
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
EXISTS= list()
exist = 1

for(n in N){
  for(sig in Sig){
    EXIST = matrix(NA, nrow = Nsim_max, ncol = length(METHODS))
    cnt = 1
    ## create/read directory
    datdir = paste0("results/", SIMNAME, "/" , family , "-", link, "/n", n)
    if(family == "binomial"){
      datdir = paste0(datdir, "_nbinom", nbinom)
    } else if (family == "gaussian"){
      datdir = paste0(datdir, "_sig", sig)
    }
    resdir = datdir
    datdir = paste0(datdir, "/dat")
    
    
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
      EXIST[i, ] = file.exists(paste0(resdir, "/", METHODS, "/dat", i, ".csv"))
      if(!all(EXIST[i, ])) next # all fitted idx
      cnt = cnt + 1; 
      
      PREDtable = data.frame(
        f = integer(),grid = double(),mean = double(),
        lb = double(),ub = double(),
        true = double(),cover = integer(),
        method = character(),id = integer()
      )
      for(method in METHODS){
        predtable = read.csv(paste0(resdir, "/", method, "/dat", i, ".csv"))
        PREDtable = rbind(PREDtable, predtable)
      }
      PREDtable$id = i
      PREDTABLE = rbind(PREDTABLE, PREDtable)
      
      ## MSE
      MSE_tmp = PREDtable %>% 
        mutate(sq_diff = (true-mean)^2) %>% 
        select(f, method, sq_diff) %>% 
        group_by(f, method) %>% summarise(logrmse = log(sqrt(mean(sq_diff))), .groups = "drop")
      MSE_tmp$sim_iter = cnt
      MSE = rbind(MSE, MSE_tmp)
      
      ## COVER
      COVER$cover = COVER$cover + PREDtable$cover
    }
    EXISTS[[exist]] = na.omit(EXIST)
    exist = exist+1
    COVER$cover = COVER$cover / (cnt-1)
    COVER$x = xgrid
  }
  PREDTABLE$n = n
  COVER$n = n
  MSE$n = n
  PREDTABLE_all = rbind(PREDTABLE_all, PREDTABLE)
  COVER_all = rbind(COVER_all, COVER)
  MSE_all = rbind(MSE_all, MSE)
}

## number of failed results
numfailed = t(sapply(EXISTS, function(x) colSums(!x)))
numtried = sapply(EXISTS, function(x) nrow(x))
colnames(numfailed) = METHODS
rownames(numfailed) = N
numfailed = cbind(numfailed, numtried)

## save results
if(family!= "gaussian"){
  # saveRDS(PREDTABLE_all, file = paste0("results/", SIMNAME, "/" , family , "-", link, "/PREDTABLE.rds"))
  # write.csv(COVER_all, file = paste0("results/", SIMNAME, "/" , family , "-", link, "/COVER.csv"),
  #           row.names=F)
  # write.csv(MSE_all, file = paste0("results/", SIMNAME, "/" , family , "-", link, "/MSE.csv"),
  #           row.names=F)
  # write.csv(numfailed, file = paste0("results/", SIMNAME, "/" , family , "-", link, "/numfailed.csv"))
  
  PREDTABLE_all = readRDS(file = paste0("results/", SIMNAME, "/" , family , "-", link, "/PREDTABLE.rds"))
  COVER_all = read.csv(paste0("results/", SIMNAME, "/" , family , "-", link, "/COVER.csv"))
  MSE_all = read.csv(paste0("results/", SIMNAME, "/" , family , "-", link, "/MSE.csv"))
  
  PLOTNAME = paste0("results/", SIMNAME, "/" , family , "-", link, "/", family,"-",link, ".pdf")
  DIAGNAME = paste0("results/", SIMNAME, "/" , family , "-", link, "/", family,"-",link, "_fit.pdf")
  dropbox = "C:/Users/USER/Dropbox/Kang/GAM3/plots/"
  PLOTNAMEdropbox = paste0(dropbox, family,"-",link, ".pdf")
  DIAGNAMEdropbox = paste0(dropbox, family,"-",link, "_fit.pdf")
  
} else {
  # saveRDS(PREDTABLE_all, file = paste0("results/", SIMNAME, "/" , family , "-", link, "/sig", Sig, "_PREDTABLE.rds"))
  # write.csv(COVER_all, file = paste0("results/", SIMNAME, "/" , family , "-", link, "/sig", Sig, "_COVER.csv"),
  #           row.names=F)
  # write.csv(MSE_all, file = paste0("results/", SIMNAME, "/" , family , "-", link,"/sig", Sig,  "_MSE.csv"),
  #           row.names=F)
  # write.csv(numfailed, file = paste0("results/", SIMNAME, "/" , family , "-", link,"/sig", Sig, "_numfailed.csv"))
  
  PREDTABLE_all = readRDS(file = paste0("results/", SIMNAME, "/" , family , "-", link, "/sig", Sig, "_PREDTABLE.rds"))
  COVER_all = read.csv(paste0("results/", SIMNAME, "/" , family , "-", link, "/sig", Sig, "_COVER.csv"))
  MSE_all = read.csv(paste0("results/", SIMNAME, "/" , family , "-", link, "/sig", Sig, "_MSE.csv"))
  
  PLOTNAME = paste0("results/", SIMNAME, "/" , family , "-", link, "/", family,"-",link,"_sig", Sig, ".pdf")
  DIAGNAME = paste0("results/", SIMNAME, "/" , family , "-", link, "/", family,"-",link,"_sig", Sig, "_fit.pdf")
  dropbox = "C:/Users/USER/Dropbox/Kang/GAM3/plots/"
  PLOTNAMEdropbox = paste0(dropbox, family,"-",link,"_sig", Sig, ".pdf")
  DIAGNSMEdropbox = paste0(dropbox, family,"-",link,"_sig", Sig, "_fit.pdf")
}

## diagnose result
comp_methods = c("mgcvad","mgcvps","grid", "bayesx");
comp_methods_name = c("Mgcv-ad","Mgcv-ps", "VS-knot", "BayesX");
pdf(DIAGNAMEdropbox, width = 3*length(comp_methods) /1.2, height = 2.2*length(TARGET)/1.2)
TARGET = c(1,6,2,3)
plotlabels = function(x){
  switch(x, mgcvad = "Mgcv-ad", mgcvps = "Mgcv-ps", grid = "VS-knot", bayesx = "BayesX")
}
plotfuncs = function(x){
  switch(as.character(x), "1"=1, "6"=4, "2"=2,"3"=3)
}

n = N[1]; maxlines = 500; fuck2 = fuck = 0
par(mar=c(1.5, 2, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,1,0), cex.axis = 1.0, cex.main= 1)
par(mfrow=c(length(TARGET),length(comp_methods)))
for(target in TARGET){
  fuck = fuck + 1
  fuck2=0
  for(method in comp_methods){
    fuck2 = fuck2+1
    plot(xgrid, functions[[target]](xgrid), type= "n", ylab = "", xlab = "", ylim = c(ylim[1]-1, ylim[2]+1),
         main = plotlabels(method))
    plotdf = PREDTABLE_all[PREDTABLE_all$n == n & PREDTABLE_all$method == method & PREDTABLE_all$f == plotfuncs(target), ]
    meancurve = rep(0, length(xgrid))
    fitted = 0
    for(i in 1:max(plotdf$id)){
      if(sum(plotdf$id == i) == 0 ) next
      meancurve = meancurve + plotdf[(plotdf$id == i), "mean"]
      fitted = fitted + 1
      if(fitted < maxlines) {
        lines(xgrid, plotdf[(plotdf$id == i), "mean"],
              lty = 1, col = t_col("grey"), lwd = 0.5)
      } 
    }
    lines(xgrid, functions[[target]](xgrid), col = "black", lwd = 0.5)
    lines(xgrid, meancurve / fitted, col = "blue", lty = 2, lwd = 1.2)
    # lines(xgrid, apply(medcurve, 2, function(x) median(x, na.rm=T)), col = "blue", lty=2, lwd=1.2)
  } 
}
dev.off()

# ## diagnose result (PPT)
# comp_methods = c("grid", "mgcvps", "bayesx");
# comp_methods_name = c("VS-knot", "Mgcv-ps", "BayesX");
# n = N[1]; maxlines = 200; fuck2 = fuck = 0
# for(target in TARGET){
#   pdf(paste0("C:/Users/USER/Naver MYBOX/hun_workstation/Papers/BVSGAM/docs/reports/DM0827/plots/", family , "-", link,"_f", target, "_sim2_fit.pdf"), 
#       width = 2.5*length(comp_methods), height = 3)
#   par(mar=c(1.5, 2, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,1,0), cex.axis = 1.0, cex.main= 1)
#   par(mfrow=c(1,length(comp_methods)))
#   fuck = fuck + 1
#   fuck2=0
#   for(method in comp_methods){
#     fuck2 = fuck2+1
#     plot(xgrid, functions[[target]](xgrid), type= "n", ylab = "", xlab = "", ylim = c(ylim[1]-1, ylim[2]+1),
#          main = comp_methods_name[fuck2])
#     plotdf = PREDTABLE_all[PREDTABLE_all$n == n & PREDTABLE_all$method == method & PREDTABLE_all$f == fuck, ]
#     meancurve = rep(0, length(xgrid))
#     fitted = 0
#     for(i in 1:max(plotdf$id)){
#       if(sum(plotdf$id == i) == 0 ) next
#       meancurve = meancurve + plotdf[(plotdf$id == i), "mean"]
#       fitted = fitted + 1
#       if(fitted < maxlines) {
#         lines(xgrid, plotdf[(plotdf$id == i), "mean"],
#               lty = 1, col = t_col("grey"), lwd = 0.5)
#       } 
#     }
#     lines(xgrid, functions[[target]](xgrid), col = "black", lwd = 0.5)
#     lines(xgrid, meancurve / fitted, col = "blue", lty = 2, lwd = 0.5)
#   } 
#   dev.off()
# }

## plot results
## 1. MSE and COVER
pdf(PLOTNAMEdropbox, width = 10*1.1, height = length(TARGET)*1.5)
TARGET = c(1,6,2,3)
plotfuncs = function(x){
  switch(as.character(x), "1"=1, "2"=4, "3"=2,"4"=3)
}
plotlabels = function(x){
  switch(as.character(x), "1"=1, "2"=4, "3"=5,"4"=6)
}
show_legend = c(F,T,F,F)
par(mar=c(1.5, 2, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,1,0), cex.axis = 1.0, cex.main= 1)
layout(matrix(rep(1:(p*(1+length(N))), rep(c(length(N)-1,rep(1, length(N))), p)), ncol=(length(N)+2), byrow=T))
plotdf = MSE_all
plotdf$method = paste0(plotdf$method, plotdf$n)
plotat = c(seq_along(METHODS), 
           seq_along(METHODS) + length(METHODS) + 1, 
           seq_along(METHODS) + length(METHODS)*2 + 2)
# ylim= c(0.48, 0.99)
ylim2= c(0.45, 0.999)
ytick = c(0.5, 0.8, 0.9, 0.95, 0.98)
coverf = function(x) 1/(1-x/1.1)-1
plotdf2 = COVER_all %>% mutate(cover = coverf(cover))
plotdf2 = plotdf2 %>% pivot_wider(names_from = method, values_from = cover)
plotdf2$x = round(plotdf2$x, 2)
gridpoints = round(seq(-0.95, 0.95, length = 11), 2)
for(pp in 1:p){
  # plot(xgrid, f_list[[pp]](xgrid), type= "l", ylab = "", xlab = "", ylim = ylim, main = paste0("true f", pp))
  plotdf %>% 
    pivot_wider(names_from = method, values_from = logrmse) %>% ungroup %>% filter(f==plotfuncs(pp)) %>%
    select(any_of(c(outer(METHODS, N, paste0)))) %>% 
    boxplot(outline = F, col = METHODS_col, border = "black", main = paste0("f", plotlabels(pp), " logrmse"),
            xaxt = "n",
            at = plotat)
  axis(1, at = c(0, 1+length(METHODS), 2+2*length(METHODS))+4, labels = paste0("n=", N))
  abline(v = c(0, 1+length(METHODS))+(length(METHODS)+1), lty=2, col = "grey")
  if(show_legend[pp]) legend("topright", col = METHODS_col, inset = c(-0.015, -0.01), bty="n",
                             lty=1, lwd = 3, legend = METHODS_plot, ncol=2, cex = 0.8)
  
  for(nn in N){
    plotdf2 %>% filter(x %in% gridpoints) %>% filter(f == plotfuncs(pp), n == nn) %>% distinct(x, .keep_all = T) %>% 
      select(any_of(METHODS)) %>% matplot(type="b", pch = 1, cex = 0.7, col = (METHODS_col), lty=1, lwd=1.2, 
                                          ylim = coverf(ylim2), x = gridpoints, main = paste0("n=", nn), xlab="", yaxt="n")
    axis(side=2, at = coverf(ytick), labels = ytick)
    abline(h = coverf(0.95), lty=2, col = "gray40")
  }
  # legend("bottomright", col = METHODS_col, lty=1, pch = 1, lwd = 1.2, legend = METHODS_plot, ncol=4)
}
dev.off()

# ## 3. COVER: table
# plotdf = COVER_all
# plotdf = plotdf %>% pivot_wider(names_from = method, values_from = cover)
# plotdf$x = round(plotdf$x, 2)
# 
# gridpoints = c(-0.95, -0.7, -0.5, -0.2, 0.0, 0.2, 0.5, 0.7, 0.95)
# TMP = data.frame()
# for(pp in 1:p){
#   for(nn in 500){
#     tmp = plotdf %>% filter(x %in% gridpoints) %>% filter(f == pp, n == nn) %>% 
#       distinct(x, .keep_all = T) %>% 
#       select(any_of(METHODS)) %>% t() %>% as.data.frame()
#     colnames(tmp) = gridpoints
#     tmp$method = METHODS
#     tmp$f = pp; #tmp$n = nn
#     TMP = rbind(TMP, tmp)
#   }
# }
# rownames(TMP) = NULL
# TMP = TMP %>% relocate(f, method)
# TMP
# write.csv(TMP, file = paste0("results/", SIMNAME, "/" , family , "-", link, "/COVER_df.csv"),
#           row.names=F)

# ## 3. COVER: table
# plotdf = COVER_all
# plotdf = plotdf %>% pivot_wider(names_from = method, values_from = cover)
# plotdf$x = round(plotdf$x, 2)
# 
# TMP = matrix(0, ncol= 2+length(METHODS), nrow = p*length(N))
# i = 1
# for(pp in 1:p){
#   for(nn in N){
#     tmp = plotdf %>% filter(f==pp, n == nn) %>% select(any_of(METHODS)) %>% colMeans %>% round(2)
#     TMP[i, 3:(length(METHODS)+2)] = tmp
#     TMP[i, 1] = pp; TMP[i,2] = nn
#     i = i+1
#   }
# }
# colnames(TMP) = c("f", "n", METHODS)
# write.csv(as.data.frame(TMP),
#           file = paste0("results/", SIMNAME, "/" , family , "-", link, "/COVER_df2.csv"),
#           row.names=F)






