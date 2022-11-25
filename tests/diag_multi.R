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
sourceCpp(paste0("src/gambma_", "grid", "_diag.cpp"))
# sourceCpp(paste0("src/gambma_", "free", "_diag.cpp"))
maxk = 20; np = 200; degree = 3;

family = "binomial"; 
if(family == "binomial"){
  N = c(500, 800, 1000);
  # N = 500
  link = "logit"; 
  nbinom = 1
  Sig = 0;
}else if(family == "poisson"){
  # N = c(150, 200, 300);
  N = 200
  link = "log"; 
  nbinom = 0
  Sig = 0;
}else if(family == "gaussian"){
  N = c(200, 400)
  # N = 200
  link = "identity"; 
  nbinom = 0
  Sig = c(2, 1)
}
intercept = ifelse(family == "poisson", -0.5, 0)

np = 200

itermax = 200
iter_mse = 200
startmid = 0

library(mgcv)
source("R_miscell/mgcv/mgcv_plot.R")
library(blapsr)
source("R_miscell/blapsr/gamlps_plot.R")
library(R2BayesX)
source("R_miscell/r2bayesX/R2BayesX_plot.R")

###############################################################################
# TARGET= list(c(1), c(2), c(3), c(4), c(5), c(6))
TARGET = list(c(1, 4))
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
  
  method1 = "mgcv"
  method2 = "mgcvad"
  method3 = "series"
  method4 = "grid"
  method5 = "free"
  method6 = "bayesx"
  method7 = "blapsr"
  METHODS = c(method1, method2, method3, method4, method5, method6, method7)
  METHODS_col = c("palevioletred1", "deeppink1", "steelblue1", "turquoise3", "slateblue1", "limegreen", "olivedrab")
  num_methods = length(METHODS)
  
  # iterate over N and Sig
  plotidx = 1
  for(n in N){
    g_manual = n
    for(sig in Sig){
      resultdir = paste0("results/test", family, "_NS", as.numeric(NS), "_bdm", bdmargin, "_poi", lambda)
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
        print(n)
        stop("watch yo jet")
        if(i_done > itermax) break
        filename = paste0(resultdir, "/sim", i, ".Rdata")
        if(!file.exists(filename)){
          dat = simmat(f_list, xrange[1], xrange[2], n, family, link, nbinom = nbinom, sig = sig)
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
          method7 = "blapsr"
          blapsr_family = family
          if(family == "binomial" & nbinom == 1) blapsr_family = "bernoulli"
          tic = Sys.time()
          fit7 = tryCatch(
            gamlps(y ~ sm(x1) + sm(x2),
                   data = dat, K = maxk, penorder = degree - 1, family=blapsr_family),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit7, "error")){
            # numfailed[method3] = numfailed[method3] + 1
            cat(paste0("Error occured in fitting ", method7), "\n")
            cat(paste0("Error message: ", fit7$message), "\n")
            rm(fit7)
            next
          } else {
            elapsed7 = difftime(toc, tic, units = c("secs"))
            tmp7 = plot_gamlps(fit7, Xgrid, f_list, 0.95, elapsed = elapsed7, ylim = ylim, Linadj = linadj)
            tmp7$method = paste0(method7)
            print(paste0(method7, " done"))
          }
          ###############################################################################
          method6 = "bayesx"
          tic = Sys.time()
          fit6 = tryCatch(
            bayesx(y ~ sx(x1, bs = "ps", degree = degree, knots = maxk)+
                     sx(x2, bs = "ps", degree = degree, knots = maxk),
                   method = "MCMC", data = dat, family = family),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          tmp = tryCatch(plot(fit6), warning = function(cnd) cnd, error = function(cnd) cnd)
          if(inherits(fit6, "error") || inherits(tmp, "warning")){
            # numfailed[method4] = numfailed[method4] + 1
            cat(paste0("Error occured in fitting ", method6), "\n")
            if(inherits(fit6, "error")) cat(paste0("Error message: ", fit6$message), "\n")
            if(inherits(tmp, "warning")) cat(paste0("Warning message: ", tmp$message), "\n")
            rm(fit6, tmp)
            next
          } else {
            elapsed6 = difftime(toc, tic, units = c("secs"))
            tmp6 = plot_bayesx(fit6, Xgrid, f_list, 0.95, maxk, elapsed = elapsed6, ylim = ylim)
            tmp6$method = paste0(method6)
            print(paste0(method6, " done"))
          }
          
          ###############################################################################
          method1 = "mgcv"
          tic = Sys.time()
          fit1 = tryCatch(
            gam(y~ s(x1, bs = "ps", k = maxk)+
                  s(x2, bs = "ps", k = maxk),
                data = dat, family = family),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit1, "error")){
            cat(paste0("Error occured in fitting ", method1), "\n")
            cat(paste0("Error message: ", fit1$message), "\n")
            rm(fit1)
            next
          } else {
            elapsed1 = difftime(toc, tic, units = c("secs"))
            tmp1 = plot_mgcv(fit1, Xgrid, f_list, 0.95, maxk, elapsed = elapsed1, ylim = ylim,
                             Linadj = linadj)
            
            tmp1$method = paste0(method1)
            print(paste0(method1, " done"))
          } 
          ###############################################################################
          method2 = "mgcvad"
          # dat = simmat(f_list, xrange[1], xrange[2], n, family, link, nbinom = nbinom, sig = sig)
          tic = Sys.time()
          fit2 = tryCatch(
            gam(y~ s(x1, bs = "ad", k = maxk)+
                  s(x2, bs = "ad", k = maxk),
                data = dat, family = family),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit2, "error")){
            cat(paste0("Error occured in fitting ", method2), "\n")
            cat(paste0("Error message: ", fit2$message), "\n")
            rm(fit2)
            next
          } else {
            elapsed2 = difftime(toc, tic, units = c("secs"))
            tmp2 = plot_mgcv(fit2, Xgrid, f_list, 0.95, maxk, elapsed = elapsed2, ylim = ylim,
                             Linadj = linadj)
            tmp2$method = paste0(method2)
            print(paste0(method2, " done"))
          } 
          ###############################################################################
          method3 = "series"
          dat = simmat(f_list, xrange[1], xrange[2], n, family, link, nbinom = nbinom, sig = sig)
          tic = Sys.time()
          fit3 = tryCatch(
            gambma(y ~ sm(x1)+sm(x2), 
                   dat, 
                   maxk = maxk, 
                   knot_config = method, prior = "betaprime", family = family, link = link,  
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
            rm(fit3)
            next
          }
          elapsed3 = difftime(toc, tic, units = c("secs"))
          tmp3 = plot(fit3, flist = f_list, show_plot = T, elapsed = elapsed3, ylim = ylim)
          
          tmp3$method = paste0(method3)
          
          ###############################################################################
          method4 = "grid"
          tic = Sys.time()
          fit4 = tryCatch(
            gambma(y ~ sm(x1) + sm(x2), dat, maxk = maxk, 
                   knot_config = "grid", prior = "betaprime", family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = g_manual,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   NS= NS,
                   bdmargin = bdmargin),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit4, "error")){
            cat(paste0("Error occured in fitting ", method4), "\n")
            cat(paste0("Error message: ", fit4$message), "\n")
            rm(fit4)
            next
          }
          elapsed4 = difftime(toc, tic, units = c("secs"))
          tmp4 = plot(fit4, flist = f_list, show_plot = T, elapsed = elapsed4, ylim = ylim)
          tmp4$method = paste0(method4)
          
          ###############################################################################
          method5 = "free"
          tic = Sys.time()
          fit5 = tryCatch(
            gambma(y ~ sm(x1) + sm(x2), dat, maxk = maxk, 
                   knot_config = "free", prior = "betaprime", family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = g_manual,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   NS= NS,
                   bdmargin = bdmargin),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit5, "error")){
            cat(paste0("Error occured in fitting ", method5), "\n")
            cat(paste0("Error message: ", fit5$message), "\n")
            rm(fit5)
            next
          }
          elapsed5 = difftime(toc, tic, units = c("secs"))
          tmp5 = plot(fit5, flist = f_list, show_plot = T, elapsed = elapsed5, ylim = ylim)
          tmp5$method = paste0(method5)
          
          ###############################################################################
          PREDtable = rbind(tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7)
          save(fit1, elapsed1,
               fit2, elapsed2,
               fit3, elapsed3,
               fit4, elapsed4,
               fit5, elapsed5,
               fit6, elapsed6,
               fit7, elapsed7,
               PREDtable, file = filename)
          i_done = i_done + 1
          
        } else {
          load(filename)
          Xgrid = fit3$Xgrid
          dat = fit3$data
          ###############################################################################
          method5 = "free"
          tic = Sys.time()
          fit5 = tryCatch(
            gambma(y ~., dat, maxk = maxk, 
                   knot_config = "free", prior = "betaprime", family = family, link = link,  
                   nbinom = nbinom,
                   g_manual = g_manual,
                   lambda = lambda,
                   MCMCiter = MCMCiter,
                   burnin = burnin,
                   NS= NS,
                   bdmargin = bdmargin),
            error = function(cnd)cnd
          )
          toc = Sys.time()
          if(inherits(fit5, "error")){
            cat(paste0("Error occured in fitting ", method5), "\n")
            cat(paste0("Error message: ", fit5$message), "\n")
            rm(fit5)
            next
          }
          elapsed5 = difftime(toc, tic, units = c("secs"))
          tmp5 = plot(fit5, flist = f_list, show_plot = T, elapsed = elapsed5, ylim = ylim)
          tmp5$method = paste0(method5)
          ##############################################################################
          if(!exists("fit7")){
            if(i < 0) next
            method7 = "blapsr"
            blapsr_family = family
            if(family == "binomial" & nbinom == 1) blapsr_family = "bernoulli"
            tic = Sys.time()
            fit7 = tryCatch(
              gamlps(y ~ 1+sm(x1)+sm(x2)+sm(x3)+sm(x4),
                     data = dat, K = maxk, penorder = degree - 1, family=blapsr_family),
              error = function(cnd)cnd
            )
            toc = Sys.time()
            if(inherits(fit7, "error")){
              cat(paste0("Error occured in fitting ", method7), "\n")
              cat(paste0("Error message: ", fit7$message), "\n")
              rm(fit7)
              next
            } else {
              elapsed7 = difftime(toc, tic, units = c("secs"))
              tmp7 = plot_gamlps(fit7, Xgrid, f_list, 0.95, elapsed = elapsed7, ylim = ylim, Linadj = linadj)
              tmp7$method = paste0(method7)
              print(paste0(method7, " done"))
            }
          }

          ##############################################################################
          if(!exists("fit6")){
            method6 = "bayesx"
            tic = Sys.time()
            fit6 = tryCatch(
              bayesx(y ~ sx(x1, bs = "ps", degree = degree, knots = maxk)+
                       sx(x2, bs = "ps", degree = degree, knots = maxk)+
                       sx(x3, bs = "ps", degree = degree, knots = maxk)+
                       sx(x4, bs = "ps", degree = degree, knots = maxk),
                     method = "MCMC", data = dat, family = family),
              error = function(cnd)cnd
            )
            toc = Sys.time()
            tmp = tryCatch(plot(fit6), warning = function(cnd) cnd, error = function(cnd) cnd)
            if(inherits(fit6, "error") || inherits(tmp, "warning")){
              # numfailed[method4] = numfailed[method4] + 1
              cat(paste0("Error occured in fitting ", method6), "\n")
              if(inherits(fit6, "error")) cat(paste0("Error message: ", fit6$message), "\n")
              if(inherits(tmp, "warning")) cat(paste0("Warning message: ", tmp$message), "\n")
              rm(fit6, tmp)
              next
            } else {
              elapsed6 = difftime(toc, tic, units = c("secs"))
              tmp6 = plot_bayesx(fit6, Xgrid, f_list, 0.95, maxk, elapsed = elapsed6, ylim = ylim)
              tmp6$method = paste0(method6)
              print(paste0(method6, " done"))
            }
          }
          
          i_done = i_done + 1
        }
        print(i_done)
        
        if(i > 0){
          png(paste0(resultdir, "/sim", i, "_", method1, ".png"), width =1200, height = 600)
          tmp1 = plot_mgcv(fit1, Xgrid, f_list, 0.95, maxk, elapsed = elapsed1, ylim = ylim,
                           Linadj = linadj)
          dev.off()
          tmp1$method = method1
          
          png(paste0(resultdir, "/sim", i, "_", method2, ".png"), width =1200, height = 600)
          tmp2 = plot_mgcv(fit2, Xgrid, f_list, 0.95, maxk, elapsed = elapsed2, ylim = ylim,
                           Linadj = linadj)
          dev.off()
          tmp2$method = method2
          
          png(paste0(resultdir, "/sim", i, "_", method3, ".png"), width = 1200, height = 600)
          tmp3 = plot(fit3, flist = f_list, elapsed = elapsed3, ylim = ylim)
          dev.off()
          tmp3$method = method3
          
          png(paste0(resultdir, "/sim", i, "_", method4, ".png"), width = 1200, height = 600)
          tmp4 = plot(fit4, flist = f_list, elapsed = elapsed4, ylim = ylim)
          dev.off()
          tmp4$method = method4
          
          png(paste0(resultdir, "/sim", i, "_", method5, ".png"), width = 1200, height = 600)
          tmp5 = plot(fit5, flist = f_list, elapsed = elapsed5, ylim = ylim)
          dev.off()
          tmp5$method = method5
          
          png(paste0(resultdir, "/sim", i, "_", method6, ".png"), width = 1200, height = 600)
          tmp6 = plot_bayesx(fit6, Xgrid, f_list, 0.95, maxk, elapsed = elapsed6, ylim = ylim)
          dev.off()
          tmp6$method = method6
          
          png(paste0(resultdir, "/sim", i, "_", method7, ".png"), width = 1200, height = 600)
          tmp7 = plot_gamlps(fit7, Xgrid, f_list, 0.95, elapsed = elapsed7, ylim = ylim, Linadj = linadj)
          dev.off()
          tmp7$method = method7
          
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
        }
        
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
          
          png(paste0(resultdir_sum, "/summary_",cnt,  ".png"), width = 1000, height = 250*p)
          layout(matrix(rep(1:(p*3), rep(c(1,2,2), p)), ncol=5, byrow=T))
          par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.2)
          for(pp in 1:p){
            plot(xgrid, f_list[[pp]](xgrid), type= "l", ylab = "", xlab = "", ylim = ylim, main = paste0("true f", pp))
            MSE %>% pivot_wider(names_from = method, values_from = logrmse) %>% ungroup %>% filter(f==pp) %>%
              select(METHODS) %>% boxplot(outline = F, col=METHODS_col, border = "black", main = paste0("f", pp, " logrmse"))
            COVER %>% pivot_wider(names_from = method, values_from = cover) %>% filter(f== pp) %>%
              select(METHODS) %>% matplot(type="l", col = METHODS_col, lty=1, lwd=2, ylim = c(0.5, 1.02), x = xgrid, main = paste0("f", pp, " 95% coverage"), xlab="")
            abline(h = 0.95, lty=2, col = "gray40")
          }
          dev.off()

          
          for(whichmethod in METHODS){
            png(paste0(resultdir_sum, "/summary_",whichmethod,  ".png"), width = 1200, height = 600)
            par(mfrow=c(2,length(target)/2))
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
            par(mfrow=c(2,length(target)/2))
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



