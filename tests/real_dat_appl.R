library(tidyverse)
source("R/gambma_diag.R")
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("src/gambma.cpp")
source("R_miscell/mgcv/mgcv_plot.R")
library("PerformanceAnalytics")

resdir = "results/realData"
if(!file.exists(resdir)) dir.create(resdir, recursive = T)

# if(knotConfig == "EVEN"){
#   Ctrl = list(numMCmodels = 100,
#               enumerate = T,
#               burnIn = 500,
#               mcIter = 500,
#               mcmcIter = 2000)
# } else if (knotConfig == "VS"){
#   Ctrl = list(burnIn = 500,
#               mcmcIter = 2000)
# } else if (knotConfig == "FREE"){
#   Ctrl = list(nu = 50, 
#               bir_p = 0.4, 
#               dea_p = 0.4,
#               initIter = 200, 
#               burnIn = 500,
#               mcmcIter = 2000,
#               thin = maxk)
# }

##########################################################################################
## NORMAL: Boston 
## https://www.cs.toronto.edu/~delve/data/boston/bostonDetail.html
## http://lib.stat.cmu.edu/datasets/boston
##########################################################################################
library(MASS)
?MASS::Boston
data(Boston)
Boston_dat = Boston
head(Boston_dat)

## check na

## data eda
hist(log(Boston_dat$medv), breaks=30)
Boston_dat %>% dplyr::select(-medv) %>% chart.Correlation(histogram=T)
# Boston_dat %>% 
#   dplyr::select(c(crim, indus, nox, age, dis, rad, tax)) %>% # highly correlated
#   chart.Correlation(histogram=T)
# 
# ## dim red. by PCA
# library(factoextra)
# Boston_dat.pca = Boston_dat %>% dplyr::select(c(indus, nox, crim, tax))
# res.pca = prcomp(Boston_dat.pca, scale=T)
# summary(res.pca)
# fviz_eig(res.pca)
# fviz_pca_var(res.pca,
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )
# fviz_pca_biplot(res.pca, repel = TRUE,
#                 col.var = "#2E9FDF", # Variables color
#                 col.ind = "#696969"  # Individuals color
# )
# PCs = get_pca_ind(res.pca)$coord
# PC1 = PCs[,1]
# par(mfrow=c(1,4))
# plot(PC1, Boston_dat.pca$indus)
# plot(PC1, Boston_dat.pca$nox)
# plot(PC1, Boston_dat.pca$crim)
# plot(PC1, Boston_dat.pca$tax)
# Boston_dat$pc1 = PC1

apply(Boston_dat, 2, function(x) length(unique(x)))

maxk=15
mf = log(medv) ~ chas +
  ncs(crim, nk = maxk, lambda=lambda) +
  ncs(zn, nk = maxk, lambda=lambda) +
  ncs(indus, nk = maxk, lambda=lambda) +
  ncs(nox, nk = maxk, lambda=lambda) +
  ncs(rm, nk = maxk, lambda=lambda) +
  ncs(age, nk = maxk, lambda=lambda) +
  ncs(dis, nk = maxk, lambda=lambda) +
  ncs(rad, nk = maxk, lambda=lambda) +
  ncs(tax, nk = maxk, lambda=lambda) +
  ncs(ptratio, nk = maxk, lambda=lambda) +
  ncs(black, nk = maxk, lambda=lambda) +
  ncs(lstat, nk = maxk, lambda=lambda)

# tic = Sys.time()
# fitEV = tryCatch( 
#   gambma(mf, Boston_dat, 
#          knotConfig = "EVEN", 
#          prior = "Robust", 
#          family = "gaussian", 
#          Ctrl = list(enumerate = F,
#                      burnIn = 1000,
#                      mcmcIter = 3000),
#          storeFitted = T),
#   error = function(cnd)cnd
# )
# toc = Sys.time()
# elapsedEV = difftime(toc, tic, units = c("secs"))
# summary(fitEV)
# plot(fitEV, show_title = F)
# plotnumknot(fitEV)
# bayesResiduals(fitEV)

# lambda = 0
# tic = Sys.time()
# fitVS0 =  tryCatch( 
#   gambma(mf, Boston_dat, 
#          knotConfig = "VS", 
#          prior = "Robust", 
#          family = "gaussian", 
#          Ctrl = list(burnIn = 1000,
#                      mcmcIter = 4000),
#          storeFitted = T),
#   error = function(cnd)cnd
# )
# toc = Sys.time()
# elapsedVS0 = difftime(toc, tic, units = c("secs"))
# summary(fitVS0)
# plot(fitVS0, show_title = F)
# plotnumknot(fitVS0)
# bayesResiduals(fitVS0)


# lambda = 1
# tic = Sys.time()
# fitVS1 =  tryCatch( 
#   gambma(mf, Boston_dat, 
#          knotConfig = "VS", 
#          prior = "Robust", 
#          family = "gaussian", 
#          Ctrl = list(burnIn = 1000,
#                      mcmcIter = 4000),
#          storeFitted = T),
#   error = function(cnd)cnd
# )
# toc = Sys.time()
# elapsedVS1 = difftime(toc, tic, units = c("secs"))
# summary(fitVS1)
# plot(fitVS1, show_title = F)
# plotnumknot(fitVS1)
# bayesResiduals(fitVS1)
# 
# 
# save(fitVS1, elapsedVS1,
#      file = paste0(resdir, "/BostonResults.rdata"))
load(paste0(resdir, "/BostonResults.rdata"))

fit = fitVS1
sink(file = paste0(resdir, "/BostonSummary.txt"))
timestamp()
summary(fit)
sink()
# file.show(paste0(resdir, "/BostonSummary.txt"))

pdf(paste0(resdir, "/BostonFit.pdf"), width = 6, height = 4)
par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
plot(fit, show_title = F, n_row=3)
dev.off()

pdf(paste0(resdir, "/BostonNumKnot.pdf"), width = 10, height = 3.5)
par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
plotnumknot(fit)
dev.off()

##########################################################################################
## BINOMIAL 
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2245318/
## https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database
##########################################################################################
library(tidyverse)
library("PerformanceAnalytics")
library("mlbench")
data("PimaIndiansDiabetes2")
Pima_dat = PimaIndiansDiabetes2 %>% dplyr::select(-insulin)
Pima_dat = na.omit(Pima_dat)
Pima_dat$diabetes = ifelse(Pima_dat$diabetes=="pos", 1, 0)
head(Pima_dat)

## eda
Pima_dat %>% dplyr::select(-diabetes) %>% 
  chart.Correlation(histogram=T)

# ## pc..?
# library(factoextra)
# Pima_dat.pca = Pima_dat %>% dplyr::select(c(mass, triceps))
# res.pca = prcomp(Pima_dat.pca, scale=T)
# summary(res.pca)
# fviz_eig(res.pca)
# fviz_pca_var(res.pca,
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )
# fviz_pca_biplot(res.pca, repel = TRUE,
#                 col.var = "#2E9FDF", # Variables color
#                 col.ind = "#696969"  # Individuals color
# )
# PCs = get_pca_ind(res.pca)$coord
# PC1 = PCs[,1]
# plot(PC1, Pima_dat.pca$triceps)
# Pima_dat$pc1 = -PC1
# 
# Pima_dat %>% dplyr::select(-c(diabetes, triceps, mass)) %>% chart.Correlation(histogram=T)
# Pima_dat %>% dplyr::select(-c(diabetes, triceps, mass)) %>% apply(.,2,function(x) length(unique(x)))

maxk=15
mf = diabetes ~ 
  ncs(pregnant, nk = maxk, lambda = lambda) +
  ncs(glucose, nk = maxk, lambda = lambda) +
  ncs(pressure, nk = maxk, lambda = lambda) +
  ncs(triceps, nk = maxk, lambda = lambda) +
  ncs(mass, nk = maxk, lambda = lambda) +
  ncs(pedigree, nk = maxk, lambda = lambda) +
  ncs(age, nk = maxk, lambda = lambda)

# tic = Sys.time()
# fitEV = tryCatch(
#   gambma(mf, Pima_dat,
#          knotConfig = "EVEN",
#          prior = "Robust",
#          family = "bernoulli",
#          Ctrl = list(enumerate = F,
#                      burnIn = 1000,
#                      mcmcIter = 3000),
#          storeFitted = T),
#   error = function(cnd)cnd
# )
# toc = Sys.time()
# elapsedEV = difftime(toc, tic, units = c("secs"))
# summary(fitEV)
# plot(fitEV, show_title = F, ylim = c(-5, 5))
# plotnumknot(fitEV)
# bayesResiduals(fitEV)

# lambda = 0
# tic = Sys.time()
# fitVS0 =  tryCatch( 
#   gambma(mf, Pima_dat, 
#          knotConfig = "VS", 
#          prior = "Robust", 
#          family = "bernoulli", 
#          Ctrl = list(burnIn = 1000,
#                      mcmcIter = 4000),
#          storeFitted = T),
#   error = function(cnd)cnd
# )
# toc = Sys.time()
# elapsedVS0 = difftime(toc, tic, units = c("secs"))
# summary(fitVS0)
# plot(fitVS0,show_title = F, ylim = c(-3,3))
# plotnumknot(fitVS0)
# bayesResiduals(fitVS0)
# 
# lambda = 1
# tic = Sys.time()
# fitVS1 =  tryCatch( 
#   gambma(mf, Pima_dat, 
#          knotConfig = "VS", 
#          prior = "Robust", 
#          family = "bernoulli", 
#          Ctrl = list(burnIn = 1000,
#                      mcmcIter = 4000),
#          storeFitted = T),
#   error = function(cnd)cnd
# )
# toc = Sys.time()
# elapsedVS1 = difftime(toc, tic, units = c("secs"))
# summary(fitVS1)
# plot(fitVS1,show_title = F)
# plotnumknot(fitVS1)
# bayesResiduals(fitVS1)

# pdf(paste0(resdir, "/PimaFit.pdf"), width = 7.5, height = 4)
# par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
# plot(fitVS,show_title = F, ylim=c(-6, 6))
# dev.off()
# 
# pdf(paste0(resdir, "/PimaNumKnot.pdf"), width = 7.5, height = 4)
# par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
# plotnumknot(fitVS)
# dev.off()

# save(fitVS1, 
#      file = paste0(resdir, "/PimaResults.rdata"))
load(paste0(resdir, "/PimaResults.rdata"))

fit = fitVS1
sink(file = paste0(resdir, "/PimaSummary.txt"))
timestamp()
summary(fit)
sink()
# file.show(paste0(resdir, "/BostonSummary.txt"))

pdf(paste0(resdir, "/PimaFit.pdf"), width = 6, height = 2.7)
par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
plot(fit, show_title = F)
dev.off()

pdf(paste0(resdir, "/PimaNumKnot.pdf"), width = 10, height = 3.5)
par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
plotnumknot(fit)
dev.off()


