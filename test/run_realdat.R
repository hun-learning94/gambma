library(gambms)
##########################################################################################
## BINOMIAL: Pima
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2245318/
## https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database
##########################################################################################
data("Pima")

maxk = 20;
mf = diabetes ~
  ncs(pregnant, nk = maxk) +
  ncs(glucose, nk = maxk) +
  ncs(pressure, nk = max) +
  ncs(triceps, nk = maxk) +
  ncs(mass, nk = maxk) +
  ncs(pedigree, nk = maxk) +
  ncs(age, nk = maxk)

fit_Pima =  gambms(mf, Pima,
                   knotConfig = "VS",
                   prior = "Intrinsic",
                   family = "bernoulli",
                   Ctrl = list(mcmcIter = 2000))

summary(fit_Pima)
plot(fit_Pima)
plotnumknot(fit_Pima)
plotresiduals(fit_Pima)

# save(fit_Pima, "test//fit_Pima.rdata")

