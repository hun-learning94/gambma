## Constructor for a class "gambvs"
## i.e. a model fit function for gambvs
library(BAS)
glmbma = function(formula, dat, 
                  prior = c("betaprime", "fixedg", "hyperg", "zsadapted", "benchmark"),
                  family = c("poisson", "binomial", "normal"),
                  link = c("log", "logit", "probit", "identity"),
                  nbinom = NULL, y_scale = NULL,
                  g_manual = NULL,
                  u_sigma = 0.1,
                  MCMCiter = 2000, 
                  burnin = 500)
{
  if(family == "gaussian") family = "normal"
  
  ############################################################
  ############################################################
  ## model design matrix X, and response y
  if(!inherits(formula, "formula")) stop("Incorrect model formula")
  mf = stats::model.frame(formula, dat = dat)
  y = matrix(stats::model.extract(mf, "response"), ncol = 1)
  X = stats::model.matrix(mf, dat = dat)
  if(is.null(dim(X))) X = matrix(X, ncol = 1)
  n = nrow(y); p = ncol(X);

  ## A, offset matrix
  A = stats::model.offset(mf)
  if(is.null(A)){
    A = matrix(rep(0, nrow(X)), ncol = 1)
  } else {
    A = matrix(A, ncol= 1)
  }
  
  ############################################################
  ############################################################
  ## check for family, link
  if(!(family %in%  c("binomial", "poisson", "normal"))) 
    stop("unsupported exponential family")
  
  ## y_scale, nbinom
  if(family == "binomial"){
    if(is.null(nbinom)) stop("no nbinom")
    y_scale = 1.0;
  } else if (family == "poisson"){
    nbinom = y_scale = 1.0;
  } else if (family == "normal"){
    nbinom = y_scale = 1.0;
  }  
  
  ## link function
  if(is.null(link)) {
    if (family == "poisson"){link <- "log"}
    if (family == "binomial"){link <- "logit"}
    if (family == "normal"){link <- "identity"}
  } else {
    if (family == "binomial"){
      if (!(link %in% c("probit", "logit"))) { stop(paste0("Invalid link (", link , ") for ", family)) }
    }
    if (family == "poisson"){
      if (!(link %in% c("log"))){ stop(paste0("Invalid link (", link , ") for ", family)) }
    }
  }
  
  ## code prior
  priorcode = function(prior){
    switch(prior, 
           hyperg = 1L, uniform = 2L, benchmark = 3L, betaprime = 4L, zsadapted = 5L, hypergn = 6L, fixedg = 7L,
           0L)
  }
  
  ## check prior validity
  # if(family %in% c("binomial", "poisson")){
  #   if(!(prior %in% c("fixedg", "betaprime", "hyperg", "zsadapted", "benchmark", "hypergn", "uniform"))){
  #     stop(paste0("Invalid prior (", prior ,") for ", family))
  #   }
  # }
  
  ############################################################
  ############################################################
  ## Define backup Rglm for each Exponential family
  if (family == "binomial") {
    if(link == "logit"){
      Rglm = function(y, X, A, nbinom, etastart){
        # res = matrix(glm(cbind(y, nbinom -y) ~ X - 1, family = binomial(link = "logit"), offset = A)$coefficients, ncol=1)
        res = tryCatch(matrix(glm.fit(X, cbind(y, nbinom - y), etastart = etastart, offset = A, family = binomial(link = "logit"))$coefficients, ncol=1),
                 warning = function(cnd) cnd)
        if(inherits(res, "warnings")) stop("MLE error")
        res[is.na(res)] = 0.0
        return(res)
      }
    }
    if(link == "probit"){
      Rglm = function(y, X, A, nbinom, etastart){
        # res = matrix(glm(cbind(y, nbinom -y) ~ X - 1, family = binomial(link = "probit"), offset = A)$coefficients, ncol=1)
        res = tryCatch(matrix(glm.fit(X, cbind(y, nbinom - y), etastart = etastart, offset = A, family = binomial(link = "probit"))$coefficients, ncol=1),
                       warning = function(cnd) cnd)
        if(inherits(res, "warnings")) stop("MLE error")
        res[is.na(res)] = 0.0
        return(res)
      }
    }
  }
  if (family == "poisson"){
    if(link == "log"){
      Rglm = function(y, X, A, nbinom, etastart){
        # res = matrix(glm(y ~ X - 1, family = "poisson", offset = A)$coefficients, ncol=1)  
        res = tryCatch(matrix(glm.fit(X, y, etastart = etastart, offset = A, family = poisson())$coefficients, ncol=1),
                       warning = function(cnd) cnd)
        if(inherits(res, "warnings")) stop("MLE error")
        res[is.na(res)] = 0.0
        return(res)
      } 
    }
  }
  if (family == "normal"){
    if(link == "identity"){
      Rglm = function(y, X, A, nbinom, etastart){
        # res = matrix(lm(y ~ X - 1, offset = A)$coefficients, ncol=1)  
        res = tryCatch(matrix(glm.fit(X, y, etastart = etastart, offset = A, family = gaussian())$coefficients, ncol=1),
                       warning = function(cnd) cnd)
        res[is.na(res)] = 0.0
        return(res)
      } 
    }
  }
  
  ## Define nearPD for use in RcppArmadillo
  nearPDres = function(x) as.matrix(Matrix::nearPD(x, doSym=T)$mat)
  
  ############################################################
  if(is.null(g_manual)) g_manual = n
  ## RUN
  res = glmbma_cpp(nbinom, 
                   y_scale,
                   y, 
                   X, 
                   A, 
                   family, 
                   priorcode(prior), 
                   link,
                   u_sigma,
                   g_manual,
                   MCMCiter,
                   Rglm, 
                   nearPDres)
  
  ############################################################
  ############################################################
  ## results
  
  out = list(
    formula = formula,
    family = family,
    data = dat,
    X = X,
    link = link,
    N = n,
    p = p,
    prior = prior
  )
  
  out$lpydiagnosis = data.frame(
    lpy = c(res$lpy)[-c(1:burnin)],
    Qm = c(res$Qm)[-c(1:burnin)],
    comp1 = c(res$comp1)[-c(1:burnin)],
    comp2 = c(res$comp2)[-c(1:burnin)],
    comp3 = c(res$comp3)[-c(1:burnin)]
  )
  
  out$Z = res$Z[-c(1:burnin),-1]
  out$InterceptSamp = res$InterceptSamp[-c(1:burnin)]
  out$g_sampled = c(1 / res$U[-c(1:burnin)] - 1)
  
  out$BetaSamp = res$BetaSamp[-c(1:burnin),]
  # BetaSamp = matrix(0, nrow = length(out$BetaSamp), ncol = ncol(out$Z))
  # for(i in 1:length(out$BetaSamp)){
  #   BetaSamp[i, out$Z[i,] > 0] = out$BetaSamp[[i]]
  # }
  # out$BetaSamp = BetaSamp
  
  if(family == "normal") out$sigma2 = res$SIGMA2
  
  attr(out, "class") = "glmbma" # S3 object
  return(out)
}
