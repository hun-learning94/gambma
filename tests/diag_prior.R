set.seed(2021311165)
library(tidyverse)
library(ggplot2)

## generate data
# source("R/gressani.R")
functionname = "hunsmooths"
source(paste0("R_miscell/", functionname, ".R"))
source("R/simmat.R")
source("R/gambma_diag.R")
library(Rcpp)
library(RcppArmadillo)

MCMCiter = 2000
burnin = 500
method = "series"
enumerate_MC = T

sourceCpp(paste0("src/gambma_", method, "_diag.cpp"))
# sourceCpp(paste0("src/gambma_", "grid", "_diag.cpp"))
# sourceCpp(paste0("src/gambma_", "free", "_diag.cpp"))
maxk = 20; np = 200; degree = 3;

family = "binomial"; 
if(family == "binomial"){
  N = c(500, 1000, 2000);
  # N = 500
  link = "logit"; 
  nbinom = 1
  Sig = 0;
}else if(family == "poisson"){
  # N = c(100, 150);
  N = 500
  link = "log"; 
  nbinom = 0
  Sig = 0;
}else if(family == "gaussian"){
  # N = c(100, 200)
  N = 300
  link = "identity"; 
  nbinom = 0
  # Sig = c(1, 2)
  Sig = 1
}
intercept = ifelse(family == "poisson", 0, 0)

np = 200

itermax = 200
iter_mse = 200
startmid = 0

# library(mgcv)
# source("R_miscell/mgcv/mgcv_plot.R")
# library(blapsr)
# source("R_miscell/blapsr/gamlps_plot.R")
# library(R2BayesX)
# source("R_miscell/r2bayesX/R2BayesX_plot.R")

###############################################################################
TARGET= list(c(1))
for(target in TARGET){
  MSE_PLOTS = list()
  COVER_PLOTS = list()
  SECD_PLOTS = list()
  FIRD_PLOTS = list()
  
  p = length(target)
  f_list = functions[target]
  linadj = Linadj[target]
  
  NS = T
  lambda = 0
  
  bdmargin = 0
  
  method1 = "fixedg"
  method2 = "betaprime"
  method3 = "zsadapted"
  method4 = "benchmark"
  method5 = "hypergn"
  method6 = "hyperg"
  method7 = "uniform"
  METHODS = c(method1, method2, method3, method4, method5, method6, method7)
  METHODS_col = seq_along(METHODS) + 1
  num_methods = length(METHODS)
  
  # iterate over N and Sig
  plotidx = 1
  for(n in N){
    g_manual = n/4
    for(sig in Sig){
      resultdir = paste0("results/prioreffect_", family, "_NS", as.numeric(NS), "_bdm", bdmargin, "_poi", lambda)
      if(!file.exists(resultdir)) dir.create(resultdir)
      resultdir = paste0(resultdir, "/", 
                         functionname , 
                         paste(target, collapse = ""),
                         "_maxk", maxk,
                         "_n", n,
                         "_nb", nbinom, 
                         "_sig", sig)
      if(intercept != 0) resultdir = paste0(resultdir, "_int", intercept)
      if(!file.exists(resultdir)) dir.create(resultdir)
      resultdir_sum = paste0(resultdir, "/summary")
      if(!file.exists(resultdir_sum)) dir.create(resultdir_sum)
      xgrid = seq(xrange[1], xrange[2], length = np)
      ###############################################################################
      
      i_done = 0
      for(i in 1:(itermax*2)){
        cat(paste0("n", n, " i", i, "\n"))
        stop("watch yo jet")
        if(i_done > itermax) break
        filename = paste0(resultdir, "/sim", i, ".Rdata")
        if(!file.exists(filename)){
          dat = simmat(f_list, xrange[1], xrange[2], n, family, link, nbinom = nbinom, sig = sig,
                       intercept = intercept)
          X = as.matrix(dat[,-1])
          mins = apply(X, 2, "min")
          maxs = apply(X, 2, "max")
          X_01 = scaleto01(X, mins, maxs)
          xgrid_01 = seq(0, 1, length = np+2); 
          xgrid_01 = xgrid_01[-c(1, length(xgrid_01))]
          Xgrid_01 = matrix(rep(xgrid_01, p), nrow=np)
          Xgrid = scalefrom01(Xgrid_01, mins, maxs)
          if(i < startmid) next
          
          # Xgrid = matrix(rep(seq(xrange[1], xrange[2], len = np), p),ncol=p)
          # plot(dat$x1, dat$y)
          
          ###############################################################################
          method1 = "fixedg"
          tic = Sys.time()
          fit1 = tryCatch(
            gambma(y ~ sm(x1), dat, 
                   maxk = maxk, 
                   knot_config = "series", prior = method1, family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = n/4,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   enumerate = T),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit1, "error")){
            cat(paste0("Error occured in fitting ", method1), "\n")
            cat(paste0("Error message: ", fit1$message), "\n")
            next
          }
          elapsed1 = difftime(toc, tic, units = c("secs"))
          tmp1 = plot(fit1, flist = f_list, show_plot = T, elapsed = elapsed1, ylim = ylim)
          tmp1$method = paste0(method1)
          
          ###############################################################################
          method2 = "betaprime"
          tic = Sys.time()
          fit2 = tryCatch(
            gambma(y ~ sm(x1), dat, 
                   maxk = maxk, 
                   knot_config = "series", prior = method2, family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = g_manual,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   enumerate = T),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit2, "error")){
            cat(paste0("Error occured in fitting ", method2), "\n")
            cat(paste0("Error message: ", fit2$message), "\n")
            next
          }
          elapsed2 = difftime(toc, tic, units = c("secs"))
          tmp2 = plot(fit2, flist = f_list, show_plot = T, elapsed = elapsed2, ylim = ylim)
          tmp2$method = paste0(method2)
          ###############################################################################
          method3 = "zsadapted"
          tic = Sys.time()
          fit3 = tryCatch(
            gambma(y ~ sm(x1), dat, 
                   maxk = maxk, 
                   knot_config = "series", prior = method3, family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = g_manual,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   enumerate_MC = enumerate_MC,
                   NS = NS,
                   bdmargin = bdmargin),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit3, "error")){
            cat(paste0("Error occured in fitting ", method3), "\n")
            cat(paste0("Error message: ", fit3$message), "\n")
            next
          }
          elapsed3 = difftime(toc, tic, units = c("secs"))
          tmp3 = plot(fit3, flist = f_list, show_plot = T, elapsed = elapsed3, ylim = ylim)
          tmp3$method = paste0(method3)
          ###############################################################################
          method4 = "benchmark"
          tic = Sys.time()
          fit4 = tryCatch(
            gambma(y ~ sm(x1), dat, 
                   maxk = maxk, 
                   knot_config = "series", prior = method4, family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = g_manual,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   enumerate_MC = enumerate_MC,
                   NS = NS,
                   bdmargin = bdmargin),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit4, "error")){
            cat(paste0("Error occured in fitting ", method4), "\n")
            cat(paste0("Error message: ", fit4$message), "\n")
            next
          }
          elapsed4 = difftime(toc, tic, units = c("secs"))
          tmp4 = plot(fit4, flist = f_list, show_plot = T, elapsed = elapsed4, ylim = ylim)
          tmp4$method = paste0(method4)
          ###############################################################################
          method5 = "hypergn"
          tic = Sys.time()
          fit5 = tryCatch(
            gambma(y ~ sm(x1), dat, 
                   maxk = maxk, 
                   knot_config = "series", prior = method5, family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = g_manual,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   enumerate_MC = enumerate_MC,
                   NS = NS,
                   bdmargin = bdmargin),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit5, "error")){
            cat(paste0("Error occured in fitting ", method5), "\n")
            cat(paste0("Error message: ", fit5$message), "\n")
            next
          }
          elapsed5 = difftime(toc, tic, units = c("secs"))
          tmp5 = plot(fit5, flist = f_list, show_plot = T, elapsed = elapsed5, ylim = ylim)
          tmp5$method = paste0(method5)
          
          ###############################################################################
          method6 = "hyperg"
          tic = Sys.time()
          fit6 = tryCatch(
            gambma(y ~ sm(x1), dat, 
                   maxk = maxk, 
                   knot_config = "series", prior = method6, family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = g_manual,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   enumerate_MC = enumerate_MC,
                   NS = NS,
                   bdmargin = bdmargin),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit6, "error")){
            cat(paste0("Error occured in fitting ", method6), "\n")
            cat(paste0("Error message: ", fit6$message), "\n")
            next
          }
          elapsed6 = difftime(toc, tic, units = c("secs"))
          tmp6 = plot(fit6, flist = f_list, show_plot = T, elapsed = elapsed6, ylim = ylim)
          tmp6$method = paste0(method6)
          
          ###############################################################################
          method7 = "uniform"
          tic = Sys.time()
          fit7 = tryCatch(
            gambma(y ~ sm(x1), dat, 
                   maxk = maxk, 
                   knot_config = "series", prior = method7, family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = g_manual,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   enumerate_MC = enumerate_MC,
                   NS = NS,
                   bdmargin = bdmargin),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit7, "error")){
            cat(paste0("Error occured in fitting ", method7), "\n")
            cat(paste0("Error message: ", fit7$message), "\n")
            next
          }
          elapsed7 = difftime(toc, tic, units = c("secs"))
          tmp7 = plot(fit7, flist = f_list, show_plot = T, elapsed = elapsed7, ylim = ylim)
          tmp7$method = paste0(method7)
          
          ###############################################################################
          i_done = i_done + 1
          PREDtable = rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7)
          save(fit1, elapsed1,
               fit2, elapsed2,
               fit3, elapsed3,
               fit4, elapsed4,
               fit5, elapsed5,
               fit6, elapsed6,
               fit7, elapsed7,
               PREDtable, file = filename)
        } else {
          load(filename)
          Xgrid = fit1$Xgrid
          i_done = i_done + 1
        }
        
        print(i_done)
        png(paste0(resultdir, "/sim", i, "_", method1, ".png"), width = 500, height = 500)
        tmp1 = plot(fit1, flist = f_list, elapsed = elapsed1, ylim = ylim)
        dev.off()
        tmp1$method = method1
        png(paste0(resultdir, "/sim", i, "_", method1, "_diagnosis.png"), width = 1000, height = 500)
        diagnosis(fit1)
        dev.off()

        png(paste0(resultdir, "/sim", i, "_", method2, ".png"), width = 500, height = 500)
        tmp2 = plot(fit2, flist = f_list, elapsed = elapsed2, ylim = ylim)
        dev.off()
        tmp2$method = method2
        png(paste0(resultdir, "/sim", i, "_", method2, "_diagnosis.png"), width = 1000, height = 500)
        diagnosis(fit2)
        dev.off()

        png(paste0(resultdir, "/sim", i, "_", method3, ".png"), width = 500, height = 500)
        tmp3 = plot(fit3, flist = f_list, elapsed = elapsed3, ylim = ylim)
        dev.off()
        tmp3$method = method3
        png(paste0(resultdir, "/sim", i, "_", method3, "_diagnosis.png"), width = 1000, height = 500)
        diagnosis(fit3)
        dev.off()
        
        png(paste0(resultdir, "/sim", i, "_", method4, ".png"), width = 500, height = 500)
        tmp4 = plot(fit4, flist = f_list, elapsed = elapsed4, ylim = ylim)
        dev.off()
        tmp4$method = method4
        png(paste0(resultdir, "/sim", i, "_", method4, "_diagnosis.png"), width = 1000, height = 500)
        diagnosis(fit4)
        dev.off()
        
        png(paste0(resultdir, "/sim", i, "_", method5, ".png"), width = 500, height = 500)
        tmp5 = plot(fit5, flist = f_list, elapsed = elapsed5, ylim = ylim)
        dev.off()
        tmp5$method = method5
        png(paste0(resultdir, "/sim", i, "_", method5, "_diagnosis.png"), width = 1000, height = 500)
        diagnosis(fit5)
        dev.off()
        
        png(paste0(resultdir, "/sim", i, "_", method6, ".png"), width = 500, height = 500)
        tmp6 = plot(fit6, flist = f_list, elapsed = elapsed6, ylim = ylim)
        dev.off()
        tmp6$method = method6
        png(paste0(resultdir, "/sim", i, "_", method6, "_diagnosis.png"), width = 1000, height = 500)
        diagnosis(fit6)
        dev.off()
        
        png(paste0(resultdir, "/sim", i, "_", method7, ".png"), width = 500, height = 500)
        tmp7 = plot(fit6, flist = f_list, elapsed = elapsed6, ylim = ylim)
        dev.off()
        tmp7$method = method7
        png(paste0(resultdir, "/sim", i, "_", method7, "_diagnosis.png"), width = 1000, height = 500)
        diagnosis(fit7)
        dev.off()

        PREDtable = rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7)
        save(fit1, elapsed1,
             fit2, elapsed2,
             fit3, elapsed3,
             fit4, elapsed4,
             fit5, elapsed5,
             fit6, elapsed6,
             fit7, elapsed7,
             PREDtable, file = filename)
        rm(fit1, fit2, fit3, fit4, fit5, fit6, fit7)
        
        ##############################################################################
        ## summary: Coverage, MSE
        if(i_done %% iter_mse == 0 & startmid ==0){
          
          # write.table(numfailed, file = paste0(resultdir_sum, "/numfailed.csv"), col.names = F)
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
          
          cnt = 0
          for(ii in 1:i){
            filename = paste0(resultdir, "/sim", ii, ".Rdata")
            if(!file.exists(filename)) next
            load(filename)
            PREDtable = PREDtable %>% filter(method %in% METHODS)
            cnt = cnt + 1
            print(cnt)
            PREDtable$id = cnt
            PREDTABLE = rbind(PREDTABLE, PREDtable)
            ## 1. MSE
            MSE_tmp = PREDtable %>% 
              mutate(sq_diff = (true-mean)^2) %>% 
              select(f, method, sq_diff) %>% 
              group_by(f, method) %>% summarise(logrmse = log(sqrt(mean(sq_diff))))
            MSE_tmp$sim_iter = cnt
            MSE = rbind(MSE, MSE_tmp)
            ## 2. coverage
            COVER$cover = COVER$cover + PREDtable$cover
          }
          COVER$cover = COVER$cover / cnt
          COVER$x = xgrid
          
          # stop("watch yo jet")
          
          library(plot.matrix)
          plotcoverage = function(COVER, METHODS, pp, nth){
            tmptable = COVER %>% pivot_wider(names_from = method, values_from = cover) %>% filter(f== pp) %>%
              select(METHODS) %>% filter(row_number() %% nth == 0) 
            for(i in 1:nrow(tmptable)){
              idx = which.min(abs(tmptable[i,] - 0.95))
              tmptable[i, idx] = 1
              tmptable[i, -idx] = 0
            }
            plot((as.matrix(tmptable)), col = c("gray88", "deepskyblue3"), key= NULL,
                 xlab="", ylab="", main="Closest 95% coverage at each quantile", border = "white",
                 axis.row = NULL)
          }
          
          png(paste0(resultdir_sum, "/summary_",cnt,  ".png"), width = 1500, height = 400*p)
          layout(matrix(rep(1:(p*3), rep(c(1,1,1), p)), ncol=3, byrow=T))
          par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.2)
          for(pp in 1:p){
            plot(xgrid, f_list[[pp]](xgrid), type= "l", ylab = "", xlab = "", ylim = ylim, main = paste0("true f", pp))
            MSE %>% pivot_wider(names_from = method, values_from = logrmse) %>% ungroup %>% filter(f==pp) %>%
              select(METHODS) %>% boxplot(outline = F, col=METHODS_col, border = "black", main = paste0("f", pp, " log10rmse"))
            legend("bottomleft", legend  = METHODS, col = METHODS_col, lty=1)
            COVER %>% pivot_wider(names_from = method, values_from = cover) %>% filter(f== pp) %>%
              select(METHODS) %>% matplot(type="l", col = METHODS_col, lty=1, lwd=2, ylim = c(0.5, 1.02), x = xgrid, main = paste0("f", pp, " 95% coverage"), xlab="")
            abline(h = 0.95, lty=2, col = "gray40")
            legend("bottomleft", legend  = METHODS, col = METHODS_col, lty=1)
            # plotcoverage(COVER, METHODS, pp, 5)
          }
          dev.off()
          
          
          for(whichmethod in METHODS){
            png(paste0(resultdir_sum, "/summary_",whichmethod,  ".png"), width = 1200, height = 600)
            par(mfrow=c(1,1))
            par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
            for(p in seq_along(target)){
              crit = PREDTABLE$f == p & PREDTABLE$method == whichmethod
              plot(PREDTABLE$grid[PREDTABLE$id == 1 & crit], PREDTABLE$true[PREDTABLE$id == 1 & crit],
                   type="n", col = "red", ylim = ylim, xlab = "", ylab = paste0("f", p),
                   main = whichmethod)
              for(i in 1:cnt){
                lines(PREDTABLE$grid[PREDTABLE$id == i & crit], PREDTABLE$mean[PREDTABLE$id == i & crit],
                      col = "grey", lty=2)
              }
              lines(PREDTABLE$grid[PREDTABLE$id == 1 & crit], apply(matrix(PREDTABLE$mean[crit], byrow = T, ncol = length(xgrid)), 2, mean),
                    col = "blue", lty =2)
              lines(PREDTABLE$grid[PREDTABLE$id == 1 & crit], PREDTABLE$true[PREDTABLE$id == 1 & crit],
                    col = "red")
            }
            dev.off()
          }
          
          for(whichmethod in METHODS){
            png(paste0(resultdir_sum, "/worstfitseach_",whichmethod,  ".png"), width = 1200, height = 600)
            par(mfrow=c(1,1))
            par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
            for(p in seq_along(target)){
              tmp = MSE %>% filter(f == p & method == whichmethod) %>% group_by(sim_iter) %>% summarise(logrmse = sum(logrmse)) %>% arrange(desc(logrmse))
              worstfits_id = head(tmp, nrow(tmp) %/% 5)$sim_iter[1:5]
              tmpPREDTABLE = PREDTABLE %>% filter(f == p & method == whichmethod & id %in% worstfits_id)
              plot(tmpPREDTABLE$grid[tmpPREDTABLE$id == worstfits_id[1]], tmpPREDTABLE$true[tmpPREDTABLE$id == worstfits_id[1]],
                   type="n", col = "red", ylim = ylim, xlab = "", ylab = paste0("f", p),
                   main = whichmethod)
              for(i in seq_along(worstfits_id)){
                lines(tmpPREDTABLE$grid[tmpPREDTABLE$id == worstfits_id[i]], tmpPREDTABLE$mean[tmpPREDTABLE$id == worstfits_id[i]],
                      col = "grey", lty=2)
              }
              lines(tmpPREDTABLE$grid[tmpPREDTABLE$id == worstfits_id[i]], tmpPREDTABLE$true[tmpPREDTABLE$id == worstfits_id[i]],
                    col = "red")
            }
            dev.off()
          }
          
          # for(whichmethod in METHODS){
          #   png(paste0(resultdir_sum, "/bestfitseach_",whichmethod,  ".png"), width = 800, height = 600)
          #   par(mfrow=c(1,length(target)))
          #   for(p in seq_along(target)){
          #     tmp = MSE %>% filter(f == p & method == whichmethod) %>% group_by(sim_iter) %>% summarise(logrmse = sum(logrmse)) %>% arrange((logrmse))
          #     worstfits_id = head(tmp, nrow(tmp) %/% 10)$sim_iter[1:10]
          #     tmpPREDTABLE = PREDTABLE %>% filter(f == p & method == whichmethod & id %in% worstfits_id)
          #     plot(tmpPREDTABLE$grid[tmpPREDTABLE$id == worstfits_id[1]], tmpPREDTABLE$true[tmpPREDTABLE$id == worstfits_id[1]],
          #          type="n", col = "red", ylim = ylim, xlab = "", ylab = paste0("f", p),
          #          main = whichmethod)
          #     for(i in seq_along(worstfits_id)){
          #       lines(tmpPREDTABLE$grid[tmpPREDTABLE$id == worstfits_id[i]], tmpPREDTABLE$mean[tmpPREDTABLE$id == worstfits_id[i]],
          #             col = "grey", lty=2)
          #     }
          #     lines(tmpPREDTABLE$grid[tmpPREDTABLE$id == worstfits_id[i]], tmpPREDTABLE$true[tmpPREDTABLE$id == worstfits_id[i]],
          #           col = "red")
          #   }
          #   dev.off()
          # }
          
          plotidx = plotidx + 1
        }
        
      }
    }
  }
}



