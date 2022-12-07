SIMNAME = "sim2"
library(tidyverse)
library(ggplot2)

source("R/gambma_diag.R")
source("R/t_col.R")

METHODS = c("Blapsr","BayesX", "Mgcv-ps", "Mgcv-ad", "Even-knot","VS-knot", "Free-knot")
METHODS_col = c("#fed976", "#fd8d3c",  "#fa9fb5", "#dd3497","#c7e9b4", "#7fcdbb", "#41b6c4")
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

family = "bernoulli"; glmWeight= 1
if(family == "binomial"){
  N = c(100, 200); glmWeight = 10;
  N = c(100)
}else if(family == "poisson"){
  N = c(100, 200, 300);
}else if(family == "gaussian"){
  N = c(100, 200, 300);
  METHODS = c("BayesX", "Mgcv-ps", "Mgcv-ad", "Even-knot","VS-knot", "Free-knot")
  METHODS_col = c("#fd8d3c",  "#fa9fb5", "#dd3497","#c7e9b4", "#7fcdbb", "#41b6c4")
  num_methods = length(METHODS)
}else if(family == "bernoulli"){
  N = c(1000, 1500, 2000); glmWeight = 1;
}

# TARGET = list(1,2,3)
TARGET = list(c(1,4,3))
plotWidth = 1000; plotHeight = 500
Nsim_max = 500
Nsim_plot = 500+1

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

EXISTS_ALL = matrix(0, nrow = length(N), ncol = num_methods); nnn = 1
colnames(EXISTS_ALL) = METHODS

for(target in TARGET){
  
  cnt = 0
  f_list = functions[target]
  p = length(target)
  linadj = Linadj[target]
  
  for(n in N){
    
    cnt = 0
    datdir = paste0("results/", SIMNAME ,"/", family ,"/n", n, "/f", paste0(target, collapse = ""))
    
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
    EXISTS = matrix(0, nrow = Nsim_max, ncol = num_methods)
    for(i in 1:Nsim_max){
      if(cnt > Nsim_plot) break
      if(i %% 50 == 0) cat(paste(datdir, i, ", cnt", cnt, "\n"))
      EXISTS[i,] = file.exists(paste0(datdir, "/", methodCode(METHODS), "/dat", i, ".csv"))
      if(!all(EXISTS[i,])) next # all fitted idx
      cnt = cnt + 1; 
      
      PREDtable = data.frame(
        f = integer(),grid = double(),mean = double(),
        lb = double(),ub = double(),
        true = double(),cover = integer(),
        method = character(),id = integer()
      )
      
      for(method in METHODS){
        predtable = read.csv(paste0(datdir, "/", methodCode(method), "/dat", i, ".csv"))
        predtable$method = method
        PREDtable = rbind(PREDtable, predtable)
      }
      PREDtable$id = i
      PREDtable$n = n
      PREDTABLE = rbind(PREDTABLE, PREDtable)
      
      ## MSE
      MSE_tmp = PREDtable %>% 
        mutate(sq_diff = (true-mean)^2) %>% 
        select(f, method, sq_diff) %>% 
        group_by(f, method) %>% summarise(logrmse = log(sqrt(mean(sq_diff))), .groups = "drop")
      MSE_tmp$sim_iter = cnt
      MSE_tmp$n = n
      MSE = rbind(MSE, MSE_tmp)
      
      ## COVER
      COVER$cover = COVER$cover + PREDtable$cover
    }
    COVER$cover = COVER$cover / (cnt-1)
    COVER$x = xgrid
    COVER$n = n
    PREDTABLE_all = rbind(PREDTABLE_all, PREDTABLE)
    COVER_all = rbind(COVER_all, COVER)
    MSE_all = rbind(MSE_all, MSE)
    
    EXISTS_ALL[nnn, ] = colSums(EXISTS)
    nnn = nnn+1
  }
  
  ## save results
  datdir = paste0("results/", SIMNAME ,"/", family ,"/n", tail(N,1), "/f", paste0(target, collapse = ""))
  PATH_PREDTABLE_all = paste0(datdir, "/PREDTABLE.rds")
  PATH_COVER = paste0(datdir, "/COVER.csv")
  PATH_MSE = paste0(datdir, "/MSE.csv")
  PATH_EXISTS = paste0(datdir, "/EXISTS.csv")
  
  saveRDS(PREDTABLE_all, file = PATH_PREDTABLE_all)
  write.csv(COVER_all, file = PATH_COVER, row.names=F)
  write.csv(MSE_all, file = PATH_MSE, row.names=F)
  write.csv(EXISTS_ALL, file = PATH_EXISTS, row.names=F)
}

target = TARGET[[1]]
datdir = paste0("results/", SIMNAME ,"/", family ,"/n", tail(N,1), "/f", paste0(target, collapse = ""))
PATH_PREDTABLE_all = paste0(datdir, "/PREDTABLE.rds")
PATH_COVER = paste0(datdir, "/COVER.csv")
PATH_MSE = paste0(datdir, "/MSE.csv")
PREDTABLE_all = readRDS(file = PATH_PREDTABLE_all)
COVER_all = read.csv(PATH_COVER)
MSE_all = read.csv(PATH_MSE)
# EXISTS = read.csv(PATH_EXISTS)


## diagnose result
p = length(TARGET[[1]])
if(family== "gaussian")mfrow=c(2,3) else mfrow=c(2,4)
# if(family == "gaussian") mfrow=c(2*p,3)
for(j in seq_along(N)){
  target = 1; n = N[j]; maxlines = 50;
  fuckTARGET = seq_along(TARGET[[1]])
  for(target in seq_along(TARGET[[1]])){
    DIAGNAME = paste0(datdir, "/",SIMNAME, "_", family,  "_diagPlotsN", N[j], "f", target,".pdf")
    pdf(DIAGNAME, width = 2*mfrow[2], height = 1.5*mfrow[1])
    par(mar=c(1.5, 2, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), 
        cex.axis = 1., cex.main= 1.1)
    par(mfrow=mfrow)
    for(method in METHODS){
      plot(xgrid, functions[[(TARGET[[1]])[target]]](xgrid), type= "n",
           ylab = "", xlab = "", ylim = c(ylim[1]-1, ylim[2]+1),
           main = method)
      plotdf = PREDTABLE_all[PREDTABLE_all$n == n &
                               PREDTABLE_all$method == method &
                               PREDTABLE_all$f == target, ]
      meancurve = rep(0, length(xgrid))
      medcurve = matrix(NA, nrow = max(plotdf$id), ncol = length(xgrid))
      fitted = 0
      for(i in sample(unique(plotdf$id), maxlines)){
        if(sum(plotdf$id == i) == 0 ) next
        meancurve = meancurve + plotdf[(plotdf$id == i), "mean"]
        medcurve[i, ] = plotdf[(plotdf$id == i), "mean"]
        fitted = fitted + 1
        if(fitted < maxlines) {
          lines(xgrid[seq(1, 200, by=2)],
                plotdf[(plotdf$id == i), "mean"][seq(1, 200, by=2)],
                lty = 1, col = t_col("grey"), lwd = 0.8)
        }
      }
      lines(xgrid, functions[[(TARGET[[1]])[target]]](xgrid), col = "blue", lwd = 1.2, lty=5)
    }
    if(family!="gaussian")plot.new()
    dev.off()
  }
}


## plot results
fucklayout = function(p, n){
  tmp = matrix(
    rep(c(rep(1, n), 1:n+1), p), 
    nrow= p, byrow=T)
  for(i in 1:p){
    tmp[i,] = tmp[i,] + (length(tmp[i,])-n+1)*(i-1)
  }
  tmp
}

fuckfuck = fucklayout(1, length(N))
plotdf = MSE_all
plotdf$method = paste0(plotdf$method, plotdf$n)
plotat = function(METHODS, N){
  out= c()
  for(i in seq_along(N)){
    out = c(out,
            seq_along(METHODS) + (i-1)*length(METHODS) + (i-1))
  }
  out
}
plotat2 = function(METHODS, N){
  out= c()
  for(i in seq_along(N)){
    out = c(out, 
            (length(METHODS)+1)*(i-1))
  }
  out = out + length(METHODS)/2
  if(length(METHODS)%%2 == 0) out = out+0.5
  out
}
plotat3 = function(METHODS, N){
  out=c()
  for(i in seq_along(N)){
    out = c(out, 
            (length(METHODS)+1)*(i-1))
  }
  out
}
ylim2 = c(0.74, 1.01)
ytick = c(0.75, 0.8, 0.85, 0.9, 0.95, 0.99)
coverf = function(x){
  x
}
plotdf2 = COVER_all %>% mutate(cover = coverf(cover))
plotdf2 = plotdf2 %>% pivot_wider(names_from = method, values_from = cover)
plotdf2$x = round(plotdf2$x, 2)
gridpoints = round(seq(-0.95, 0.95, length = 11), 2)
for(pp in seq_along(TARGET[[1]])){
  PLOTNAME = paste0(datdir, "/",SIMNAME, "_", family, "_f", pp, "_summary.pdf")
  pdf(PLOTNAME, width = ncol(fuckfuck)*3, height = nrow(fuckfuck)*4)
  par(mar=c(3.5, 4, 3, 1), mgp=c(2,1.5,0), oma = c(0,0,0,0),
      cex.axis = 2.4, cex.main =2.5)
  layout(fuckfuck)
  
  
  plotdf %>% 
    pivot_wider(names_from = method, values_from = logrmse) %>% ungroup %>% filter(f==pp) %>%
    select(any_of(c(outer(METHODS, N, paste0)))) %>% 
    boxplot(outline = F, col = METHODS_col, border = "black", main = paste0("logRMSE"),
            xaxt = "n",
            at = plotat(METHODS, N))
  axis(1, at = plotat2(METHODS, N), labels = paste0("n=", N))
  abline(v = plotat3(METHODS, N), lty=2, col = "grey")
  # if(pp!=3) legend("topright", col = METHODS_col, bty="n", lwd = 5, legend = METHODS, ncol=2, cex = 1.5)
  
  for(nn in N){
    plotdf2 %>% filter(x %in% gridpoints) %>% filter(f == pp, n == nn) %>% distinct(x, .keep_all = T) %>% 
      select(any_of(METHODS)) %>% matplot(type="b", pch = 1, cex = 1, col = (METHODS_col), lty=1, lwd=2, 
                                          ylim = coverf(ylim2), x = gridpoints, main = paste0("n=", nn), xlab="", yaxt="n")
    axis(side=2, at = coverf(ytick), labels = ytick)
    abline(h = coverf(0.95), lty=2, col = "gray40")
  }
  if(pp==3) legend("bottomright", col = METHODS_col, bty="n", lwd = 5, legend = METHODS, ncol=1, cex = 2)
  
  dev.off()
}












