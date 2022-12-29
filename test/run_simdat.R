library(gambms)
################################################################################
## data generation
################################################################################
set.seed(2021311165)
f_list = list(f1 = function(x) 0.5 * (2*x^5 + 3*x^2 + cos(3*pi*x) - 1),
              f2 = function(x) x,
              f3 = function(x) 0.75*(0.0035 * (x*3 + 1.5)^3 + (x > -0.5 & x < 0.85) *
                                       0.07 *sin(1.7*pi*(x*3 + 1.5)^2 / 3.2)*(x*3 -2.5)^2 * exp(x*3 + 1.5)))
n = 300
dat = simmat(f_list, -1, 1, n = n, family = "poisson")


################################################################################
## fit gambms
################################################################################
maxk = 20
mf = y~ncs(x1, nk=maxk)+ ncs(x2, nk = maxk)  + ncs(x3, nk = maxk)

fit_sim = tryCatch(
  gambms(mf, dat,
         knotConfig = "FREE",
         prior = "Robust",
         family = "poisson",
         printIter=500,
         freeCtrl=list(mcmcIter = 2000, thin = maxk),
         storeFitted = T),
  error = function(cnd)cnd
)

summary(fit_sim)
plot(fit_sim)
plotnumknot(fit_sim)
plotresiduals(fit_sim)

pdf("test/fit_sim_plot.pdf", width = 9, height=3)
plot(fit_sim)
dev.off()

