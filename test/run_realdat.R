library(gambms)
##########################################################################################
## BINOMIAL: Pima
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2245318/
## https://www.kaggle.com/datasets/uciml/pima-indians-diabetes-database
##########################################################################################
data("Pima")

mf = diabetes ~
  ncs(pregnant, nk = 20) +
  ncs(glucose, nk = 20) +
  ncs(pressure, nk = 20) +
  ncs(triceps, nk = 20) +
  ncs(mass, nk = 20) +
  ncs(pedigree, nk = 20) +
  ncs(age, nk = 20)

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

