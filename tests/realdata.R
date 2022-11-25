source("R/gambma_diag.R")
Rcpp::sourceCpp("src/gambma_series_diag.cpp")
Rcpp::sourceCpp("src/gambma_grid_diag.cpp")
Rcpp::sourceCpp("src/gambma_free_diag.cpp")
source("R_miscell/mgcv/mgcv_plot.R")
MCMCiter = 1000; burnin = 200
maxk = 20
##########################################################################################
## poisson
##########################################################################################
library(blapsr)
data("medicaid")
dat = medicaid
formular = numvisits ~ children + race + maritalstat + 
  sm(age) + sm(income1000) + sm(access) + sm(pc1times1000)

if(!inherits(formular, "formula")) stop("Incorrect model formula")
mf = stats::model.frame(formular, dat = dat)
y = matrix(stats::model.extract(mf, "response"), ncol = 1)
X = stats::model.matrix(mf, dat = dat)
if(is.null(dim(X))) X = matrix(X, ncol = 1)
smterms = grepl("sm(", colnames(X), fixed = T)
if(!any(smterms)) stop("no smooth term")
Xlin = X[, !smterms, drop=F]
X = X[, smterms, drop= F]
n = nrow(y); p = ncol(X); plin = ncol(Xlin)

uniqcnts = apply(X, 2, function(x) length(unique(x)))
maxk_vec = rep(maxk, p)
maxk_vec[uniqcnts - maxk < 7] = uniqcnts[uniqcnts - maxk < 7] - 7

maxk = 20
tic = Sys.time()
fit1 = tryCatch(
  gambma(formular, dat, 
         maxk = maxk, 
         maxk_vec = c(maxk, 15, 11, maxk),
         knot_config = "series", prior = "betaprime", family = "poisson", link = "log",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 30000,
         burnin = 10000,
         enumerate = F,
         NS = T,
         bdmargin = 0),
  error = function(cnd)cnd
)
toc = Sys.time()
plot(fit1, elapsed = difftime(toc, tic, units = c("secs")))
fit1$maxk
series_MH_maxk20 = fit1
# direct enumeration numknots marginal for each p

## MH numknots marignal

par(mfrow=c(2,2))
for(i in 1:4){
  fit1$numknots[,i] %>% table() %>% plot()
}

# png(paste0("results/", maxk, "poisson_series.png"), width =800, height = 500)
# plot(fit1, elapsed = difftime(toc, tic, units = c("secs")))
# dev.off()

tic = Sys.time()
fit2 = tryCatch(
  gambma(formular, dat, 
         maxk = maxk, 
         # maxk_vec = c(maxk, 15, 11, maxk),
         knot_config = "grid", prior = "betaprime", family = "poisson", link = "log",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 11000,
         burnin = 1000,
         NS = T,
         bdmargin = 0),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit2)
# png(paste0("results/", maxk, "poisson_grid.png"), width =800, height = 500)
plot(fit2, elapsed = difftime(toc, tic, units = c("secs")))
# dev.off()

tic = Sys.time()
fit3 = tryCatch(
  gambma(formular, dat, 
         maxk = maxk, 
         # maxk_vec = c(maxk, 15, 11, maxk),
         knot_config = "free", prior = "betaprime", family = "poisson", link = "log",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 2000,
         burnin = 500),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit3)
# png(paste0("results/", maxk, "poisson_grid.png"), width =800, height = 500)
plot(fit3, elapsed = difftime(toc, tic, units = c("secs")))
plotnumknot(fit3)

tic = Sys.time()
fit4 = tryCatch( 
  gam(numvisits ~ children + race + maritalstat + 
        s(age, bs = "ps", k = 10)+ 
        s(income1000, bs = "ps", k = 10)+
        s(access, bs = "ps", k = 10)+
        s(pc1times1000, bs = "ps", k = 10),
      data = medicaid, family =  "poisson"),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit4)
png(paste0("results/", maxk, "poisson_mgcvps.png"), width =800, height = 500)
plot_mgcv(fit4, Xgrid = fit2$Xgrid, maxk = maxk, elapsed = difftime(toc, tic, units = c("secs")))
dev.off()

tic = Sys.time()
fit5 = tryCatch( 
  gam(numvisits ~ children + race + maritalstat + 
        s(age, bs = "ad", k = 20)+ 
        s(income1000, bs = "ad", k = 16)+
        s(access, bs = "ad", k = 10)+
        s(pc1times1000, bs = "ad", k = 20),
      data = medicaid, family =  "poisson"),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit5)
png(paste0("results/", maxk, "poisson_mgcvad.png"), width =800, height = 500)
plot_mgcv(fit5, Xgrid = fit2$Xgrid, maxk = maxk, elapsed = difftime(toc, tic, units = c("secs")))
dev.off()

save(fit1, fit2, fit4, fit5, file = paste0("results/", maxk, "poisson_result.rdata"))
load(paste0("results/", maxk, "poisson_result.rdata"))

tic = Sys.time()
fit6 = tryCatch(
  bayesx(numvisits ~ children + race + maritalstat +
           sx(age, bs = "ps", degree = degree, knots = maxk)+
           sx(income1000, bs = "ps", degree = degree, knots = maxk)+
           sx(access, bs = "ps", degree = degree, knots = maxk)+
           sx(pc1times1000, bs = "ps", degree = degree, knots = maxk),
         method = "MCMC", data = medicaid, family = "poisson"),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit6)
plot(fit6)
png(paste0("results/", maxk, "poisson_bayesx.png"), width =800, height = 500)
plot_mgcv(fit5, Xgrid = fit1$Xgrid, maxk = maxk, elapsed = difftime(toc, tic, units = c("secs")))
dev.off()

##########################################################################################
## normal
##########################################################################################
maxk=20

dat = MASS::Boston
# formular = medv ~ chas + rad+ black + zn +age+
#   sm(crim) + sm(indus) + sm(nox) + sm(rm)+ sm(dis)+ sm(tax)+ sm(ptratio)+ sm(lstat)
formular = medv ~ chas + rad+ black + zn + age+ptratio+
   sm(log(crim)) + sm(indus) + sm(nox) + sm(rm) + sm(log(dis)) + sm(tax) + sm(lstat)

if(!inherits(formular, "formula")) stop("Incorrect model formula")
mf = stats::model.frame(formular, dat = dat)
y = matrix(stats::model.extract(mf, "response"), ncol = 1)
X = stats::model.matrix(mf, dat = dat)
if(is.null(dim(X))) X = matrix(X, ncol = 1)
smterms = grepl("sm(", colnames(X), fixed = T)
if(!any(smterms)) stop("no smooth term")
Xlin = X[, !smterms, drop=F]
X = X[, smterms, drop= F]
n = nrow(y); p = ncol(X); plin = ncol(Xlin)

uniqcnts = apply(X, 2, function(x) length(unique(x)))
uniqcnts

tic = Sys.time()
fit = tryCatch(
  gambma(formular, dat,
         maxk =20,
         knot_config = "series", prior = "betaprime", family = "gaussian", link = "identity",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 30000,
         burnin = 10000,
         enumerate = F),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit)
plot(fit, elapsed = difftime(toc, tic, units = c("secs")))
plotnumknot(fit)
fit0 = fit

tic = Sys.time()
fit = tryCatch(
  gambma(formular, dat, 
         maxk =20, 
         knot_config = "grid", prior = "betaprime", family = "gaussian", link = "identity",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 3000,
         burnin = 1000),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit)
plot(fit, elapsed = difftime(toc, tic, units = c("secs")))
plotnumknot(fit)
fit1 = fit


library(mgcv)
maxk = 10
fit4 = tryCatch( 
  gam(medv ~ chas + rad + black + zn+
        s(crim, bs = "ps", k = maxk)+ 
        s(indus, bs = "ps", k = maxk)+ 
        s(nox, bs = "ps", k = maxk)+ 
        s(rm, bs = "ps", k = maxk)+ 
        s(dis, bs = "ps", k = maxk)+
        s(tax, bs = "ps", k = maxk)+ 
        s(ptratio, bs = "ps", k = maxk)+
        s(lstat, bs = "ps", k = maxk),
      data = MASS::Boston, family =  "gaussian"),
  error = function(cnd)cnd
)
summary(fit4)
plot_mgcv(fit4, Xgrid = fit1$Xgrid, maxk = maxk)











##########################################################################################\
library(gamair); data(mpg); head(mpg)
dat = mpg

formular = city.mpg ~ fuel + style + drive + 
  sm(weight) + sm(hp)

if(!inherits(formular, "formula")) stop("Incorrect model formula")
mf = stats::model.frame(formular, dat = dat)
y = matrix(stats::model.extract(mf, "response"), ncol = 1)
X = stats::model.matrix(mf, dat = dat)
if(is.null(dim(X))) X = matrix(X, ncol = 1)
smterms = grepl("sm(", colnames(X), fixed = T)
if(!any(smterms)) stop("no smooth term")
Xlin = X[, !smterms, drop=F]
X = X[, smterms, drop= F]
n = nrow(y); p = ncol(X); plin = ncol(Xlin)

uniqcnts = apply(X, 2, function(x) length(unique(x)))
uniqcnts

maxk = 20
tic = Sys.time()
fit4 = tryCatch( 
  gam(city.mpg ~ fuel + style + drive + 
        s(weight, bs = "ps", k = maxk)+ 
        s(hp, bs = "ps", k = maxk),
      data = dat, family =  "gaussian"),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit4)
plot_mgcv(fit4, Xgrid = fit1$Xgrid, maxk = maxk, elapsed = difftime(toc, tic, units = c("secs")))


tic = Sys.time()
fit1 = tryCatch(
  gambma(formular, dat, 
         maxk =20, 
         knot_config = "series", prior = "betaprime", family = "gaussian", link = "identity",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 30000,
         burnin = 10000,
         enumerate = T),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit1)
# png(paste0("results/", maxk, "normal_series.png"), width =800, height = 500)
plot(fit1, elapsed = difftime(toc, tic, units = c("secs")))
plotnumknot(fit1)

tic = Sys.time()
fit2 = tryCatch(
  gambma(formular, dat, 
         maxk = maxk, 
         knot_config = "grid", prior = "betaprime", family = "gaussian", link = "identity",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 3000,
         burnin = 1000,
         NS = T,
         bdmargin = 0),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit2)
# png(paste0("results/", maxk, "normal_grid.png"), width =800, height = 500)
plot(fit2, elapsed = difftime(toc, tic, units = c("secs")))
# dev.off()
plotnumknot(fit2)

tic = Sys.time()
fit3 = tryCatch(
  gambma(formular, dat, 
         maxk = 20, 
         knot_config = "free", prior = "betaprime", family = "gaussian", link = "identity",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 3000,
         burnin = 1000),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit3)
plot(fit3, elapsed = difftime(toc, tic, units = c("secs")))
plotnumknot(fit3)

##########################################################################################
## normal
##########################################################################################
dat = read.table("plasma.txt")
colnames(dat) = c("age", "sex", "smokestat", "bmi", "vituse", "calories", "fat", "fiber", "alchohol", "cholesterol", "betadiet", "retdiet", "betaplasma", "retplasma")
dat$log_betaplasma = log(dat$betaplasma+0.01)
dat$log_cholesterol = log(dat$cholesterol)
head(dat)
formular = log_betaplasma ~ bmi + betadiet + sex + smokestat + fiber + fat +
  sm(age) + sm(log_cholesterol)

tic = Sys.time()
fit1 = tryCatch(
  gambma(formular, dat, 
         maxk =20, 
         knot_config = "series", prior = "betaprime", family = "gaussian", link = "identity",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 30000,
         burnin = 10000,
         enumerate = T,
         NS = T,
         bdmargin = 0),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit1)
# png(paste0("results/", maxk, "normal_series.png"), width =800, height = 500)
plot(fit1, elapsed = difftime(toc, tic, units = c("secs")))
plotnumknot(fit1)
# dev.off()

tic = Sys.time()
fit2 = tryCatch(
  gambma(formular, dat, 
         maxk = maxk, 
         knot_config = "grid", prior = "betaprime", family = "gaussian", link = "identity",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 3000,
         burnin = 1000,
         NS = T,
         bdmargin = 0),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit2)
# png(paste0("results/", maxk, "normal_grid.png"), width =800, height = 500)
plot(fit2, elapsed = difftime(toc, tic, units = c("secs")))
# dev.off()
plotnumknot(fit2)

tic = Sys.time()
fit3 = tryCatch(
  gambma(formular, dat, 
         maxk = 20, 
         knot_config = "free", prior = "betaprime", family = "gaussian", link = "identity",
         nbinom = nbinom,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 2000,
         burnin = 1000),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit3)
plot(fit3, elapsed = difftime(toc, tic, units = c("secs")))
plotnumknot(fit3)

tic = Sys.time()
fit4 = tryCatch( 
  gam(log_betaplasma ~ bmi + betadiet + sex + smokestat + fiber + fat +
        s(age, bs = "ps", k = maxk)+ 
        s(log_cholesterol, bs = "ps", k = maxk),
      data = dat, family =  "gaussian"),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit4)
png(paste0("results/", maxk, "normal_mgcvps.png"), width =800, height = 500)
plot_mgcv(fit4, Xgrid = fit1$Xgrid, maxk = maxk, elapsed = difftime(toc, tic, units = c("secs")))
dev.off()

tic = Sys.time()
fit5 = tryCatch( 
  gam(log_betaplasma ~ bmi + betadiet + sex + smokestat + fiber + fat +
        s(age, bs = "ad", k = maxk)+ 
        s(log_cholesterol, bs = "ad", k = maxk),
      data = dat, family =  "gaussian"),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit5)
png(paste0("results/", maxk, "normal_mgcvad.png"), width =800, height = 500)
plot_mgcv(fit5, Xgrid = fit1$Xgrid, maxk = maxk, elapsed = difftime(toc, tic, units = c("secs")))
dev.off()

save(fit1, fit2, fit4, fit5, file = paste0("results/", maxk, "normal_result.rdata"))
load(paste0("results/", maxk, "normal_result.rdata"))

tic = Sys.time()
fit6 = tryCatch(
  bayesx(log_betaplasma ~ bmi + betadiet + sex + smokestat + fiber + fat +
           sx(age, bs = "ps", degree = degree, knots = maxk)+
           sx(log_cholesterol, bs = "ps", degree = degree, knots = maxk),
         method = "MCMC", data = dat, family = "gaussian"),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit6)
plot(fit6)
png(paste0("results/", maxk, "normal_bayesx.png"), width =800, height = 500)
plot_bayesx(fit6, Xgrid = fit1$Xgrid, maxk = maxk, elapsed = difftime(toc, tic, units = c("secs")))
dev.off()

##########################################################################################
## binomial
##########################################################################################
library("mlbench")
data("PimaIndiansDiabetes2")
dat = PimaIndiansDiabetes2
dat = dat[, -5] # no insulin
dat = na.omit(dat)
head(dat)

formular = diabetes ~ sm(pregnant) + sm(glucose) + sm(pressure) + sm(triceps) + sm(mass) + sm(pedigree) + sm(age)

maxk = 20
mf = stats::model.frame(formular, dat = dat)
y = matrix(stats::model.extract(mf, "response"), ncol = 1)
X = stats::model.matrix(mf, dat = dat)
if(is.null(dim(X))) X = matrix(X, ncol = 1)
smterms = grepl("sm(", colnames(X), fixed = T)
if(!any(smterms)) stop("no smooth term")
Xlin = X[, !smterms, drop=F]
X = X[, smterms, drop= F]
n = nrow(y); p = ncol(X); plin = ncol(Xlin)
uniqcnts = apply(X, 2, function(x) length(unique(x)))
maxk_vec = rep(maxk, p)
maxk_vec[uniqcnts - maxk < 5] = uniqcnts[uniqcnts - maxk < 5] - 5

tic = Sys.time()
fit1 = tryCatch(
  gambma(formular, dat, 
         maxk = 20, 
         knot_config = "series", prior = "betaprime", family = "binomial", link = "logit",
         nbinom = 1,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 30000,
         burnin = 10000,
         NS = T,
         bdmargin = 0),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit1)
# png(paste0("results/", maxk, "binomial_series.png"), width =800, height = 500)
plot(fit1, elapsed = difftime(toc, tic, units = c("secs")))

pdf("PIMA1.pdf", width = 8, height = 4)
plot(fit1)
dev.off()
plotnumknot(fit1)
# dev.off()

tic = Sys.time()
fit2 = tryCatch(
  gambma(formular, dat, 
         maxk = 20, 
         knot_config = "grid", prior = "betaprime", family = "binomial", link = "logit",
         nbinom = 1,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 3000,
         burnin = 1000,
         NS = T,
         bdmargin = 0),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit2)
# png(paste0("results/", maxk, "binomial_grid.png"), width =800, height = 500)
plot(fit2, elapsed = difftime(toc, tic, units = c("secs")))
plotnumknot(fit2)
# dev.off()
pdf("PIMA2.pdf", width = 8, height = 4)
plot(fit2)
dev.off()


tic = Sys.time()
fit3 = tryCatch(
  gambma(formular, dat, 
         maxk = 20, 
         knot_config = "free", prior = "betaprime", family = "binomial", link = "logit",
         nbinom = 1,
         g_manual = nrow(dat),
         lambda = 0,
         MCMCiter = 3000,
         burnin = 1000,
         NS = T,
         bdmargin = 0),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit3)
# png(paste0("results/", maxk, "binomial_grid.png"), width =800, height = 500)
plot(fit3, elapsed = difftime(toc, tic, units = c("secs")))
plotnumknot(fit3)
# dev.off()


tic = Sys.time()
fit4 = tryCatch(
  gam(diabetes ~
        s(pregnant, bs = "ps", k = maxk)+
        s(glucose, bs = "ps", k = maxk)+
      s(pressure, bs = "ps", k = maxk)+
      s(triceps, bs = "ps", k = maxk)+
      s(insulin, bs = "ps", k = maxk)+
      s(mass, bs = "ps", k = maxk)+
        s(pedigree, bs = "ps", k = maxk)+
        s(age, bs = "ps", k = maxk),
      data = dat, family =  "binomial"),
  error = function(cnd)cnd
)
toc = Sys.time()
png(paste0("results/", maxk, "binomial_mgcvps.png"), width =800, height = 500)
plot_mgcv(fit4, Xgrid = fit1$Xgrid, maxk = maxk, elapsed = difftime(toc, tic, units = c("secs")))
dev.off()

tic = Sys.time()
fit5 = tryCatch(
  gam(diabetes ~
        s(pregnant, bs = "ad", k = maxk)+
        s(glucose, bs = "ad", k = maxk)+
        s(pressure, bs = "ad", k = maxk)+
        s(triceps, bs = "ad", k = maxk)+
        s(insulin, bs = "ad", k = maxk)+
        s(mass, bs = "ad", k = maxk)+
        s(pedigree, bs = "ad", k = maxk)+
        s(age, bs = "ad", k = maxk),
      data = dat, family =  "binomial"),
  error = function(cnd)cnd
)
toc = Sys.time()
png(paste0("results/", maxk, "binomial_mgcvad.png"), width =800, height = 500)
plot_mgcv(fit5, Xgrid = fit1$Xgrid, maxk = maxk, elapsed = difftime(toc, tic, units = c("secs")))
dev.off()

save(fit1, fit2, fit4, fit5, file = paste0("results/", maxk, "binomial_result.rdata"))
load(paste0("results/", maxk, "binomial_result.rdata"))


tic = Sys.time()
fit6 = tryCatch(
  bayesx(diabetes ~ 
           sx(pregnant, bs = "ps", degree = degree, knots = maxk)+
           sx(glucose, bs = "ps", degree = degree, knots = maxk)+
           sx(pressure, bs = "ps", degree = degree, knots = maxk)+
           sx(triceps, bs = "ps", degree = degree, knots = maxk)+
           sx(insulin, bs = "ps", degree = degree, knots = maxk)+
           sx(mass, bs = "ps", degree = degree, knots = maxk)+
           sx(pedigree, bs = "ps", degree = degree, knots = maxk)+
           sx(age, bs = "ps", degree = degree, knots = maxk),
         method = "MCMC", data = dat, family = "binomial"),
  error = function(cnd)cnd
)
toc = Sys.time()
summary(fit6)
plot(fit6)
