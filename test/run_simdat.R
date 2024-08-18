################################################################################
## data generation
################################################################################
set.seed(1)
f_list = list(f1 = function(x) 0.5 * (2*x^5 + 3*x^2 + cos(3*pi*x) - 1),
              f2 = function(x) 0.75*(0.0035 * (x*3 + 1.5)^3 + (x > -0.5 & x < 0.85) * 
                                       0.07 *sin(1.7*pi*(x*3 + 1.5)^2 / 3.2)*
                                       (x*3 -2.5)^2 * exp(x*3 + 1.5)),
              f3 = function(x) x,
              f4 = function(x) x*0)
n = 200
dat = simmat(f_list, -1, 1, n = n, family = "poisson")


################################################################################
## fit gambms
################################################################################
mf = y~ncs(x1, nk = 20)+ 
  ncs(x2, nk = 20)  + 
  ncs(x3, nk = 20) + 
  ncs(x4, nk = 20)

fit_sim = tryCatch(
  gambms(mf, dat,
         knotConfig = "FREE",
         prior = "Intrinsic",
         family = "poisson"),
  error = function(cnd)cnd
)

plot(fit_sim)
plotnumknot(fit_sim)
plotresiduals(fit_sim)

# save(fit_sim, file = "test//fit_sim.rdata")


