library(gambms)
##########################################################################################
## NORMAL: Boston
## https://www.cs.toronto.edu/~delve/data/boston/bostonDetail.html
## http://lib.stat.cmu.edu/datasets/boston
## deprecated
##########################################################################################
data("Boston")

maxk=15; lambda = .1
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

fit_Boston =  gambms(mf, Boston,
                 knotConfig = "EVEN",
                 prior = "Intrinsic",
                 family = "gaussian")
summary(fit_Boston)
plot(fit_Boston)
plotnumknot(fit_Boston)
plotresiduals(fit_Boston)
save(fit_Boston, file = "test//fit_Boston.rdata")

##########################################################################################
## BINOMIAL: Pima
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2245318/
## https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database
##########################################################################################
data("Pima")

maxk=10; lambda = .1
mf = diabetes ~
  ncs(pregnant, nk = maxk, lambda = lambda) +
  ncs(glucose, nk = maxk, lambda = lambda) +
  ncs(pressure, nk = maxk, lambda = lambda) +
  ncs(triceps, nk = 5, lambda = lambda) +
  ncs(mass, nk = maxk, lambda = lambda) +
  ncs(pedigree, nk = maxk, lambda = lambda) +
  ncs(age, nk = 5, lambda = lambda)

fit_Pima =  gambms(mf, Pima,
                   knotConfig = "VS",
                   prior = "Intrinsic",
                   family = "bernoulli",
                   Ctrl = list(mcmcIter = 5000))

summary(fit_Pima)
plot(fit_Pima)
plotnumknot(fit_Pima)
plotresiduals(fit_Pima)

save(fit_Pima, "test//fit_Pima.rdata")

