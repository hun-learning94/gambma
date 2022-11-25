## Constructor for a class "gambvs"
## i.e. a model fit function for gambvs
source("R/scale01.R")
sm = function(x) x
ncs = function(x, nk = 10, lambda = 0){
  term = as.character(substitute(x))
  label = paste0("ncs(", term, ")")
  ret = list(term = term, label = label, nk = nk, lambdaVec = lambda)
  class(ret) = "ncsSmoothSpec"
  ret
}

gambma = function(fm, dat, 
                  knotConfig = c("EVEN", "VS", "FREE"),
                  prior = c("Unit", "Hyper-g", "Uniform", "Hyper-g/n", "Beta-prime", "ZS-adapted",
                            "Robust", "Intrinsic", "constant", "CH"),
                  family = c("poisson", "gaussian", "bernoulli"),
                  link = NULL,
                  glmWeight = NULL, 
                  # params related to manual priors
                  g_manual = NA_real_, 
                  a_manual = NA_real_, 
                  b_manual = NA_real_, 
                  s_manual = NA_real_,
                  # params for even, vs, free knots
                  Ctrl = list(),
                  # evenCtrl = list(numMCmodels = 100,
                  #                 enumerate = FALSE,
                  #                 burnIn = 500,
                  #                 mcIter = 500,
                  #                 mcmcIter = 2000),
                  # vsCtrl = list(burnIn = 500,
                  #               mcmcIter = 2000),
                  # freeCtrl = list(nu = 50, 
                  #                 bir_p = 0.4, 
                  #                 dea_p = 0.4,
                  #                 initIter = 200, 
                  #                 burnIn = 500,
                  #                 mcmcIter = 2000,
                  #                 thin = NULL),
                  # params for progress monitoring
                  printIter = 200,
                  # grid length for prediction
                  np = 200,
                  # to save the fitted values
                  storeFitted = F)
{
  ############################################################
  ############################################################
  ## model design matrix X, and response y
  if(!inherits(fm, "formula")) stop("Incorrect model fm")
  tf = terms.formula(fm, specials = "ncs")
  smoothTermsIdx = attr(tf, "specials")$ncs
  if(is.null(smoothTermsIdx)) stop("no smooth term")
  allTerms = rownames(attr(tf, "factors"))
  smoothTerms = allTerms[smoothTermsIdx]
  smoothTermsSpec = list()
  for(i in seq_along(smoothTerms)){ smoothTermsSpec[[i]] = try(eval(parse(text = smoothTerms[i]))) }
  maxkVec = sapply(smoothTermsSpec, function(x) x$nk)
  lambdaVec = sapply(smoothTermsSpec, function(x) x$lambdaVec)
  smoothTermsName = sapply(smoothTermsSpec, function(x) x$term)
  fakeFormula = ""; j=1
  for(i in seq_along(allTerms)){
    if(i == 1){ fakeFormula = paste0(fakeFormula, allTerms[i])
    } else { 
      if(i %in% smoothTermsIdx) {fakeFormula = paste0(fakeFormula, smoothTermsSpec[[j]]$term); j=j+1}
      else fakeFormula = paste0(fakeFormula, allTerms[i])
    }
    if(i == 1){ fakeFormula = paste0(fakeFormula, "~")
    } else if(i < length(allTerms)) { fakeFormula = paste0(fakeFormula, "+") }
    
  }
  fakeFormula = formula(fakeFormula)
  
  mf = stats::model.frame(fakeFormula, dat = dat)
  y = stats::model.extract(mf, "response")
  X = stats::model.matrix(mf, dat = dat)
  if(is.null(dim(X))) X = matrix(X, ncol = 1)
  smterms = (colnames(X) %in% smoothTermsName)
  XLin = X[, !smterms, drop=F]
  X = X[, smterms, drop= F]
  n = length(y); p = ncol(X); pLin = ncol(XLin)
  
  ## X should be normalized to lie between 0 and 1, inclusive
  mins = apply(X, 2, "min")
  maxs = apply(X, 2, "max")
  X_01 = scaleto01(X, mins, maxs)
  xgrid_01 = seq(0, 1, length = np+2); 
  xgrid_01 = xgrid_01[-c(1, length(xgrid_01))]
  Xgrid_01 = matrix(rep(xgrid_01, p), nrow=np)
  Xgrid = scalefrom01(Xgrid_01, mins, maxs)
  
  ## orthogonalize XLin (column-wise centering and rescaling)
  XLinOrg = XLin
  if(ncol(XLin)>1){
    if(knotConfig == "FREE") stop("currently free-knot does not support linear covariates")
    XLinSub = XLin[,-1, drop=F]
    XLinSubMn = colMeans(XLinSub)
    XLinSubSd = apply(XLinSub, 2, sd)
    XLinSubStd = sweep(XLinSub, 2, XLinSubMn)
    XLinSubStd = sweep(XLinSubStd, 2, XLinSubSd, FUN="/")
    XLin = cbind(1, XLinSubStd)
  }
  
  ## A, offset matrix
  A = stats::model.offset(mf)
  if(is.null(A)) A = rep(0, nrow(X))
  
  ############################################################
  ############################################################
  ## check for family, link
  if(family == "bernoulli"){family = "binomial"; glmWeight = 1.0}
  if(!(family %in%  c("binomial", "poisson", "gaussian"))) 
    stop("unsupported exponential family")
  
  ## y_scale, glmWeight
  if(family == "binomial"){
    if(is.null(glmWeight)) stop("no glmWeight")
  } else if (family == "poisson"){
    glmWeight = 1.0;
  } else if (family == "gaussian"){
    glmWeight = 1.0;
  }  
  
  ## link function
  if(is.null(link)) {
    if (family == "poisson"){link <- "log"}
    if (family == "binomial"){link <- "logit"}
    if (family == "gaussian"){link <- "identity"}
  } else {
    if (family == "binomial"){
      if (!(link %in% c("probit", "logit"))) { stop(paste0("Invalid link (", link , ") for ", family)) }
    }
    if (family == "poisson"){
      if (!(link %in% c("log"))){ stop(paste0("Invalid link (", link , ") for ", family)) }
    }
  }
  
  ## code familyLink
  getFamilyLinkCode = function(family, link){
    familyLink = NULL
    if (family == "poisson") {
      familyLink = 31L
    } else if(family == "binomial"){
      if(link == "logit") familyLink = 21L
      if(link == "probit") familyLink = 22L
    } else if (family == "gaussian"){
      familyLink= 11L
    }
    return(familyLink)
  }
  familyLinkCode = getFamilyLinkCode(family, link)
  if(is.null(familyLinkCode)) stop(paste0("Unsupported family - link : ", family, " - ", link))
  
  ## code prior
  getPriorCode = function(prior){
    # c("Unit", "Hyper-g", "Uniform", "Hyper-g/n", "Beta-prime", "ZS-adapted", "Robust", "Intrinsic", "constant", "CH")
    switch(prior, 
           "Unit" = 0L, "Hyper-g" = 1L, "Uniform" = 2L, "Hyper-g/n" = 3L, "Beta-prime" = 4L, 
           "ZS-adapted" = 5L, "Robust" = 6L, "Intrinsic" = 7L, "constant" = 8L, "CH" = 9L,
           NULL)
  }
  priorCode = getPriorCode(prior)
  if(is.null(priorCode)) stop(paste0("Unsupported prior : ", prior))
  if(priorCode == 8 & is.na(g_manual)) stop("input g for constant g prior")
  if(priorCode == 9 & any(is.na(c(a_manual, b_manual, s_manual)))) stop("input a, b, s for CH g prior")
  
  ############################################################
  ############################################################
  ## Define backup Rglm for each Exponential family
  if (family == "binomial") {
    if(link == "logit"){
      Rglm = function(y, X, A, glmWeight, etastart){
        # res = matrix(glm(cbind(y, glmWeight -y) ~ X - 1, family = binomial(link = "logit"), offset = A)$coefficients, ncol=1)
        res = tryCatch(matrix(glm.fit(X, cbind(y, glmWeight - y), etastart = etastart, offset = A, family = binomial(link = "logit"))$coefficients, ncol=1),
                 warning = function(cnd) cnd)
        if(inherits(res, "warnings")) stop("MLE error")
        res[is.na(res)] = 0.0
        return(res)
      }
    }
    if(link == "probit"){
      Rglm = function(y, X, A, glmWeight, etastart){
        # res = matrix(glm(cbind(y, glmWeight -y) ~ X - 1, family = binomial(link = "probit"), offset = A)$coefficients, ncol=1)
        res = tryCatch(matrix(glm.fit(X, cbind(y, glmWeight - y), etastart = etastart, offset = A, family = binomial(link = "probit"))$coefficients, ncol=1),
                       warning = function(cnd) cnd)
        if(inherits(res, "warnings")) stop("MLE error")
        res[is.na(res)] = 0.0
        return(res)
      }
    }
  }
  if (family == "poisson"){
    if(link == "log"){
      Rglm = function(y, X, A, glmWeight, etastart){
        # res = matrix(glm(y ~ X - 1, family = "poisson", offset = A)$coefficients, ncol=1)  
        res = tryCatch(matrix(glm.fit(X, y, etastart = etastart, offset = A, family = poisson())$coefficients, ncol=1),
                       warning = function(cnd) cnd)
        if(inherits(res, "warnings")) stop("MLE error")
        res[is.na(res)] = 0.0
        return(res)
      } 
    }
  }
  if (family == "gaussian"){
    if(link == "identity"){
      Rglm = function(y, X, A, glmWeight, etastart){
        # res = matrix(lm(y ~ X - 1, offset = A)$coefficients, ncol=1)  
        # res = tryCatch(matrix(glm.fit(X, y, etastart = etastart, offset = A, family = gaussian())$coefficients, ncol=1),
        #                warning = function(cnd) cnd)
        res = tryCatch(matrix(lm.fit(X, y)$coefficients, ncol=1),
                       warning = function(cnd) cnd)
        res[is.na(res)] = 0.0
        return(res)
      } 
    }
  }
  
  ############################################################
  ############################################################
  ## Define nearPD for use in RcppArmadillo
  nearPDres = function(x) as.matrix(Matrix::nearPD(x, doSym=T)$mat)
  
  ## adjust maxk
  if(is.null(maxkVec)){
    uniqCnts = apply(X, 2, function(x) length(unique(x)))
    maxkVec = rep(maxk, p)
    if(any(uniqCnts < maxkVec)){
      defects = colnames(X)[uniqCnts < maxkVec]
      defectsUniqCnts = uniqCnts[uniqCnts < maxkVec]
      defectsMaxkVec = maxkVec[uniqCnts < maxkVec]
      for(i in seq_along(defects)){
        warning(paste0("possible issue with variable ", defects[i], 
                       " as the number of unique counts ", defectsUniqCnts[i],
                       " is less than maxk ",  defectsMaxkVec[i]))
      }
      warning("readjusting to be maxk < unique counts")
      maxkVec[uniqCnts < maxkVec] = uniqCnts[uniqCnts < maxkVec] - 1
    }
    # maxkVec[uniqCnts - maxk < 5] = uniqCnts[uniqCnts - maxk < 5] - 5
  }
  
  
  ############################################################
  ############################################################
  ## RUN
  # if(CRADNS) maxk = maxk - 2
  if(knotConfig == "EVEN"){
    evenCtrl = Ctrl
    if(is.null(evenCtrl$numMCmodels)) evenCtrl$numMCmodels = 100L
    if(is.null(evenCtrl$enumerate)) evenCtrl$enumerate = F
    if(is.null(evenCtrl$burnIn)) evenCtrl$burnIn = 500L
    if(is.null(evenCtrl$mcIter)) evenCtrl$mcIter = 2000L
    if(is.null(evenCtrl$mcmcIter)) evenCtrl$mcmcIter = 2000L
    res = gambmaEVEN(y, 
                     glmWeight, 
                     X_01, 
                     Xgrid_01,
                     XLin,
                     A, 
                     maxkVec, 
                     lambdaVec,
                     familyLinkCode,
                     priorCode,
                     a_manual,
                     b_manual,
                     s_manual,
                     g_manual,
                     evenCtrl$enumerate,
                     evenCtrl$numMCmodels,
                     evenCtrl$mcIter,
                     evenCtrl$mcmcIter + evenCtrl$burnIn,
                     Rglm, nearPDres, 
                     storeFitted, printIter)
  } else if(knotConfig == "FREE"){
    freeCtrl = Ctrl
    if(is.null(freeCtrl$nu)) freeCtrl$nu = 100
    if(is.null(freeCtrl$bir_p)) freeCtrl$bir_p = 0.4
    if(is.null(freeCtrl$dea_p)) freeCtrl$dea_p = 0.4
    if(is.null(freeCtrl$initIter)) freeCtrl$initIter = 200L
    if(is.null(freeCtrl$burnIn)) freeCtrl$burnIn = 500L
    if(is.null(freeCtrl$mcmcIter)) freeCtrl$mcmcIter = 2000L
    if(is.null(freeCtrl$thin)) freeCtrl$thin = maxk
    res = gambmaFREE(y, 
                     glmWeight, 
                     X_01, 
                     Xgrid_01,
                     XLin,
                     A, 
                     maxkVec, 
                     lambdaVec,
                     familyLinkCode,
                     priorCode,
                     a_manual,
                     b_manual,
                     s_manual,
                     g_manual,
                     freeCtrl$initIter,
                     freeCtrl$mcmcIter + freeCtrl$burnIn,
                     freeCtrl$thin,
                     freeCtrl$bir_p, freeCtrl$dea_p, freeCtrl$nu, 
                     Rglm, nearPDres, storeFitted, printIter)
  } else if(knotConfig == "VS"){
    vsCtrl = Ctrl
    if(is.null(vsCtrl$burnIn)) vsCtrl$burnIn = 500L
    if(is.null(vsCtrl$mcmcIter)) vsCtrl$mcmcIter = 2000L
    res = gambmaVS(y, 
                   glmWeight, 
                   X_01, 
                   Xgrid_01,
                   XLin,
                   A, 
                   maxkVec, 
                   lambdaVec,
                   familyLinkCode,
                   priorCode,
                   a_manual,
                   b_manual,
                   s_manual,
                   g_manual,
                   vsCtrl$mcmcIter + vsCtrl$burnIn,
                   Rglm, nearPDres, F, storeFitted, printIter)
  } else {stop("invalid knotConfig")}
  print("cpp done")
  ############################################################
  ############################################################
  ## results
  out = list(
    fm = fm,
    family = family,
    data = dat,
    y = y,
    X = X,
    XLin = XLin,
    XLinOrg = XLinOrg,
    link = link,
    N = n,
    p = p,
    pLin = pLin,
    maxk = maxkVec,
    lambdaVec = lambdaVec,
    knotConfig = knotConfig,
    prior = prior,
    Xgrid = Xgrid,
    np = np
  )
  
  if(knotConfig == "EVEN"){
    if(res$enumerate){  # DEMC
      out$enumerate = T
      out$mcIter = evenCtrl$mcIter
      out$knotnumsAll = res$KNOTNUMS
      out$knotLocs = out$MCKNOTLOCS
      out$knotLocsIdx = out$MCKNOTLOCSIDX
      getIdx = 1:length(out$knotnumsAll)
    } else { # RWMH
      getIdx = (evenCtrl$burnIn+1):evenCtrl$mcmcIter
      out$mcmcIter = evenCtrl$mcmcIter
      out$numknots = res$KNOTNUMS
      out$acceptProb = res$AcceptProp
    }
  } else if (knotConfig == "VS"){
    getIdx = (vsCtrl$burnIn+1):vsCtrl$mcmcIter
    out$mcmcIter = vsCtrl$mcmcIter
    out$knots = res$KNOTS[getIdx,]
    tmp = list()
    start = out$pLin + 2
    for(i in 1:out$p){
      end = start + out$maxk[i]-1
      tmp[[i]] = out$knots[, start:end]
      start = end + 2
    }
    out$numknots = sapply(tmp, function(x) rowSums(x))
  } else if (knotConfig == "FREE"){
    getIdx = (freeCtrl$burnIn+1):(freeCtrl$mcmcIter+freeCtrl$burnIn)
    out$mcmcIter = freeCtrl$mcmcIter
    out$knotLocs = res$KNOTS[getIdx]
    out$knotLocsIdx = res$KNOTSIDX[getIdx]
    out$acceptProb = res$AcceptProp
    
    ppp = matrix(0, nrow = length(out$knotLocs), ncol = out$p)
    for(i in seq_along(out$knotLocs)){
      for(j in 1:out$p){
        ppp[i,j] = sum((out$knotLocsIdx[[i]] == j) & (out$knotLocs[[i]] >0))
      }
    }
    out$numknots = ppp
    
    # out$fuckyou = res$fuckyou
  }
  
  out$lpydiagnosis = data.frame(
    lpy = c(res$lpy),
    r2Qm = c(res$r2Qm),
    comp1 = c(res$comp1),
    comp2 = c(res$comp2),
    comp3 = c(res$comp3)
  )[getIdx, ]
  
  out$linears = res$PREDLINEARS[getIdx, ]
  out$predicted = res$PREDSMOOTHS[getIdx, ]
  out$fitted = res$FITTEDSMOOTHS[getIdx, ]
  out$g_sampled = res$G[getIdx, ]
  if(family == "gaussian") out$sigma2 = 1/res$PHI[getIdx, ] 
  
  # adjust for centering of linear coefficients
  if(ncol(XLin)>1){
    linbetasidx = 2:ncol(XLin)
    out$linears[,1] = out$linears[,1] + out$linears[,linbetasidx, drop=F] %*% (XLinSubMn / XLinSubSd)
    out$linears[,-1] = sweep(out$linears[,-1, drop=F], 2, XLinSubSd, FUN="/")
  }
  
  attr(out, "class") = "gambma" # S3 object
  return(out)
}

###### print method ######
print.gambma = function(x){
  cat("fm: \n")
  print(x$fm)
  cat("\n")
  cat("Family:                                          ", x$family, "\n")
  cat("Link function:                                   ", x$link, "\n")
  cat("g-prior:                                         ", x$prior, "\n")
  cat("------------------------------------------------------------ \n")
  cat("Sample size:                                     ", x$N, "\n")
  cat("Number of smooth covariates:                     ", x$p, "\n")
  cat("Number of linear terms (including intercept):    ", x$pLin, "\n")
  cat("Knot configuration:                              ", x$knotConfig, "\n")
  cat("Maximum number of knots:                         ", x$maxk, "\n")
  # cat("Natural spLine:                                  ", x$NS, "\n")
  # if(x$NS){
  #   cat("  Natural spLine boundary at:                    ", x$bdmargin, "\n")
  # } 
  cat("------------------------------------------------------------ \n")
}

summary.gambma = function(x){
  if(is.null(dim(x$X))) x$X = matrix(x$X, ncol = 1)
  if(is.null(dim(x$linears))) x$linears = matrix(x$linears, ncol = 1)
  mean_post = colMeans(x$linears)
  median_post = apply(x$linears,2, median)
  sd_post = apply(x$linears, 2, sd)
  z_score = mean_post / sd_post
  lower95 = apply(x$linears, 2, function(x) quantile(x, 0.025))
  upper95 = apply(x$linears, 2, function(x) quantile(x, 1-0.025))
  linear_coeff = data.frame(
    mean = mean_post,
    median = median_post,
    sd = sd_post,
    z_score = z_score,
    lower95 = lower95,
    upper95 = upper95
  )
  rownames(linear_coeff) = colnames(x$XLinOrg)
  
  g_samp = data.frame(
    mean = mean(fit$g_sampled),
    median = median(fit$g_sampled),
    sd = sd(fit$g_sampled),
    lower95 = quantile(fit$g_sampled, 0.025),
    upper95 = quantile(fit$g_sampled, 1-0.025)
  )
  rownames(g_samp) = "g"
  
  if(x$family == "gaussian"){
    sigsamp = sqrt(fit$sigma2)
    sig_samp = data.frame(
      mean = mean(sigsamp),
      median = median(sigsamp),
      sd = sd(sigsamp),
      lower95 = quantile(sigsamp, 0.025),
      upper95 = quantile(sigsamp, 1-0.025)
    )
    rownames(sig_samp) = "sig"
  }
  
  cat("fm: \n")
  print(x$fm)
  cat("\n")
  cat("Family:                                          ", x$family, "\n")
  cat("Link function:                                   ", x$link, "\n")
  cat("g-prior:                                         ", x$prior, "\n")
  cat("------------------------------------------------------------ \n")
  cat("Sample size:                                     ", x$N, "\n")
  cat("Number of smooth covariates:                     ", x$p, "\n")
  cat("Number of linear terms (including intercept):    ", x$pLin, "\n")
  cat("Knot configuration:                              ", x$knotConfig, "\n")
  cat("Maximum number of knots:                         ", x$maxk, "\n")
  cat("Prior mean on number of knots:                   ", x$lambdaVec, "\n")
  # cat("Natural spLine:                                  ", x$NS, "\n")
  # if(x$NS){
  # cat("  Natural spLine boundary at:                    ", x$bdmargin, "\n")
  # } 
  cat("------------------------------------------------------------- \n")
  cat("Marginal probability that sm is linear: \n")
  tmp = colMeans(x$numknots == 0)
  # names(tmp) = paste0("sm(",1:x$p, ")")
  names(tmp) = colnames(x$X)
  print(signif(tmp, 2))
  cat("------------------------------------------------------------- \n")
  cat("Linear coefficients:\n")
  print.table(format(round(as.matrix(linear_coeff), 4), nsmall=4), right=T)
  cat("------------------------------------------------------------ \n")
  cat("Sampled g values:\n")
  print.table(format(round(as.matrix(g_samp), 4), nsmall=4), right=T)
  cat("------------------------------------------------------------- \n")
  if(x$family == "gaussian"){
    cat("Sampled sigma values:\n")
    print.table(format(round(as.matrix(sig_samp), 4), nsmall=4), right=T)
    cat("------------------------------------------------------------- \n")
  }
}

###### plot method ######
plot.gambma = function(x, alpha = 0.025, flist = NULL, elapsed = NULL, show_plot = T, ylim= NULL, 
                       show_title = T, Linadj = NULL, n_row = NULL)
{
  old_pars <- par(no.readonly = TRUE)
  on.exit(par(old_pars), add = TRUE)
  
  p = x$p
  pLin = x$pLin
  family = x$family
  knotConfig = x$knotConfig
  maxk = x$maxk
  Xgrid = x$Xgrid; np = x$np
  
  ## PREDtable
  PREDtable = data.frame(
    f = integer(),
    grid = double(),
    mean = double(),
    lb = double(),
    ub = double()
  )
  if(!is.null(flist)){
    PREDtable$true = double()
    PREDtable$cover = integer()
  }
  if(is.null(x$enumerate)) x$enumerate = T
  if(is.null(x$whichmethod)) x$whichmethod = "MH"
  idx = 1:np
  
  if(is.null(colnames(x$X))) x$X = paste0("sm(x", 1:p, ")")
  substr_xnames = function(x, n){
    substr(x, n, nchar(x) - 1)
  }
  xnames = colnames(x$X)
  
  for(j in 1:p){
    xdata = x$data[, j+pLin]
    xgrid = Xgrid[,j]
    Fits = x$predicted[, idx]
    Fits_mean = colMeans(Fits) 
    Fits_ub = apply(Fits, 2, function(x) quantile(x, 1 - alpha)) 
    Fits_lb = apply(Fits, 2, function(x) quantile(x, alpha))
    predtable = data.frame(f = j, grid = xgrid, mean = Fits_mean, lb = Fits_lb, ub = Fits_ub)
    if(!is.null(flist)){
      ytrue = flist[[j]](xgrid); ytrue = ytrue - mean(flist[[j]](xdata))
      cover = ifelse((ytrue > Fits_lb) & (ytrue < Fits_ub) , 1, 0)
      predtable$true = ytrue
      predtable$cover = cover
      if(!is.null(Linadj)){
        if(Linadj[[j]] & sum(cover)/length(xgrid) < 0.3){ predtable$cover = 0 }
      }
    }
    PREDtable = rbind(PREDtable, predtable)
    idx = idx + np
  }
  
  ## plot
  if(show_plot == T){
    opar = par(no.readonly = T)
    if(p < 4){
      plotarray = c(1, (p))
    } else {
      if(is.null(n_row)) n_row = 2;
      plotarray = c(n_row, ceiling((p)/n_row))
    }
    par(mfrow=plotarray)
    # par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
    if(show_title == F) par(oma = c(0,0,0,0))
    
    idx = 1:np
    for(j in 1:p){
      predtable = PREDtable[PREDtable$f ==j ,]
      xgrid = predtable$grid
      Fits_mean = predtable$mean
      Fits_ub = predtable$ub
      Fits_lb = predtable$lb
      ytrue = predtable$true
      if(!is.null(ylim)){
        graphics::plot(xgrid, Fits_mean, type="n", 
                       xlab = xnames[j], ylab = paste0("f(",xnames[j], ")"), ylim = ylim)
      } else {
        tmpylim = range(c(Fits_mean, Fits_ub, Fits_lb))
        tmpylim[1] = tmpylim[1] - 0.2; tmpylim[2] = tmpylim[2] + 0.2;
        graphics::plot(xgrid, Fits_mean, type="n", xlab = xnames[j], ylab = paste0("f(",xnames[j], ")"), ylim = tmpylim)
      }
      graphics::polygon(x = c(xgrid, rev(xgrid)), y = c(Fits_lb, rev(Fits_ub)), col = "grey", border = NA)
      graphics::lines(xgrid, Fits_mean, col = "blue", lty=2)
      if(!is.null(flist)) graphics::lines(xgrid, ytrue, col = "red")
      graphics::rug(x$X[,j])
      idx = idx + np
      if(!is.null(flist)){
        rmse = sqrt(mean((Fits_mean - ytrue)^2))
        title(main = paste0("logrmse ", signif(log(rmse), 3)))
      }
    }
    # hist(x$g_sampled, nclass=30, probability = T, 
    #      xlab= "g", col = "grey", border = "white", main ="")
    # if(x$family=="gaussian"){
    #   hist(sqrt(x$sigma2), nclass=30, probability = T, 
    #        xlab= "resid. std. dev.", col = "grey", border = "white", main ="")
    # }
    title_ = paste0(family, ", ", knotConfig, ", maxk ", max(maxk))
    if(!is.null(elapsed)) title_ = paste0(title_, ", elapsed ", paste0(signif(elapsed, 3), attr(elapsed, "units")))
    # rmse = sqrt(mean((Fits_mean - ytrue)^2))
    # title_ = paste0(title_,  ", logrmse ", signif(log(rmse), 3))
    if(!is.null(x$acceptProb)) title_ = paste0(title_,  ", Accept ", signif(x$acceptProb, 3))
    if(show_title == T) mtext(title_, outer= T, font=2)
    par(opar)
    
  }
  return(invisible(PREDtable))
}

# plothighs = function(x, ...) UseMethod("plothighs")
# plothighs.gambvs = function(x, xp, xrange = c(0, 1), smoo_index = 1, cred_int = 0.95, 
#                             np = 100, fit_col = "blue", shade_col = "gray75", fits_col = "gray75",
#                             highs_n = 1000, plot_highs_n = 100, adj.intercept = F, ...)


# plotknot = function(x){
#   par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
#   ##############################################################################
#   if(x$knotConfig == "VS"){
#     if(is.null(dim(x$knots))) x$knots = matrix(x$knots, ncol=1)
#     zmat = t(apply(x$knots, 1, function(y) idx_to_z(y, x$maxk, x$p)))
#     incprob = colMeans(zmat)
#     start = end = 1
#     if(x$p < 4){plotarray = c(1, p)} else {plotarray = c(2, ceiling(p/2))}
#     par(mfrow=plotarray)
#     for(pp in 1:x$p){
#       xgrid = seq(range(x$Xgrid[,pp])[1], range(x$Xgrid[,pp])[2], length = x$maxk + 2)
#       xgrid = xgrid[-c(1, length(xgrid))]
#       end = start + x$maxk - 1
#       plot(xgrid, incprob[start:end], type="h", ylim = c(0, 1), lwd=3, col = "steelblue",
#            main = paste0("Marginal knot inc. prob., f", pp),
#            xlab = paste0("x", pp), ylab = "")
#       start = end + 1
#     }
#     return(invisible(0))
#   }
#   ##############################################################################
#   if(x$knotConfig == "FREE"){
#     knots = unlist(x$knotLocs)
#     knots_idx = unlist(x$knotLocsIdx)
# 
#     if(x$p < 4){plotarray = c(1, p)} else {plotarray = c(2, ceiling(p/2))}
#     par(mfrow=plotarray)
#     for(pp in 1:x$p){
#       tmphist = hist(knots[knots_idx == pp], breaks = 20, plot=F)
#       counts = tmphist$counts / length(x$knotLocs)
#       plot(tmphist$breaks[-1], counts, type="h", ylim = c(0, 1), lwd=3, col = "steelblue",
#            main = paste0("Marginal (discretized) knot inc. prob., f", pp),
#            xlab = paste0("x", pp), ylab = "")
#     }
#     return(invisible(0))
#   }
#   ##############################################################################
#   if(x$knotConfig == "series"){
#     ############################################################################
#     if(x$whichmethod == "Enumerate"){
#       prob = exp(x$lpydiagnosis$lpy - max(x$lpydiagnosis$lpy))
#       prob = prob/sum(prob)
#       plot(prob)
# 
#       if(x$p < 4){plotarray = c(1, p)} else {plotarray = c(2, ceiling(p/2))}
#       par(mfrow=plotarray)
#       for(pp in 1:x$p){
#         tmpdf = data.frame(
#           knotnums = x$knotnumsAll[, pp],
#           prob = prob
#         )
#         plot(aggregate(prob ~ knotnums, data = tmpdf, sum), type="b",
#              ylim=c(0,1), lwd = 2, cex = 1.5,
#              col = "steelblue",
#              main = paste0("Marginal prob. of num. of knots, f", pp),
#              xlab = paste0("Num. of knots, f", pp),
#              ylab = "")
#       }
#       return(invisible(0))
#     }
#     ############################################################################
#     if(x$whichmethod == "MH"){
#       knotmat = matrix(0, nrow= length(x$knots), ncol = x$p)
#       for(i in 1:length(x$knots)){
#         knotmat[i,] = idx_to_knotnums(x$knots[i], x$p, x$maxk)
#       }
#       if(x$p < 4){plotarray = c(1, p)} else {plotarray = c(2, ceiling(p/2))}
#       par(mfrow=plotarray)
#       for(pp in 1:x$p){
#         tmpvec = rep(0,x$maxk+1)
#         names(tmpvec) = 0:x$maxk
#         tmpknot = table(knotmat[,pp])
#         tmpvec[names(tmpvec) %in% names(tmpknot)] = tmpknot
#         tmpvec = tmpvec / sum(tmpvec)
#         plot(tmpvec, type="b",
#              ylim=c(0,1), lwd = 2, cex = 1.5,
#              col = "steelblue",
#              main = paste0("Marginal prob. of num. of knots, f", pp),
#              xlab = paste0("Num. of knots, f", pp),
#              ylab = "")
#       }
#       return(invisible(0))
#     }
#   }
# }

# diagnosis = function(x){
#   par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
#   if(x$p > 1) stop("diagnosis only for p = 1")
#   type = "p"
#   if(x$knotConfig == "EVEN"){
#     if(x$enumerate){
#       numknot = c(x$knotnumsAll)
#       type = "b"
#     }
#   }
# 
#   if(x$knotConfig == "VS"){
#     # if(is.null(dim(x$knots))) x$knots = matrix(x$knots, ncol=1)
#     # zmat = t(apply(x$knots, 1, function(y) idx_to_z(y, x$maxk, x$p)))
#     # numknot = c(rowSums(zmat))
#     numknot= x$numknots
#   }
# 
#   loglik = x$lpydiagnosis$comp1
#   if(x$family == "gaussian") loglik = x$lpydiagnosis$r2Qm
#   pen2 = x$lpydiagnosis$comp2 
#   pen3 = x$lpydiagnosis$comp3
#   lpy = x$lpydiagnosis$lpy
# 
#   center = function(x) x - mean(x)
#   par(mfrow=c(1,4))
#   ylim = range(c(center(loglik),center(pen2), center(pen3), center(lpy)))
#   plot(numknot, center(loglik), ylim = ylim,
#        type=type, lwd = 2, cex = 1.5,
#        col = "steelblue",
#        main = "Max loglikelihood",
#        xlab = "Num. of knots",
#        ylab = "")
#   plot(numknot, center(pen2), ylim = ylim,
#        type=type, lwd = 2, cex = 1.5,
#        col = "steelblue",
#        main = "Penalty",
#        xlab = "Num. of knots",
#        ylab = "")
#   plot(numknot, center(pen3), ylim = ylim,
#        type=type, lwd = 2, cex = 1.5,
#        col = "steelblue",
#        main = "Model prior",
#        xlab = "Num. of knots",
#        ylab = "")
#   plot(numknot, center(lpy), ylim = ylim,
#        type=type, lwd = 2, cex = 1.5,
#        col = "steelblue",
#        main = "Model evidence",
#        xlab = "Num. of knots",
#        ylab = "")
#   mtext(x$prior, outer = T, font =2)
#   return(invisible(0))
# }


plotnumknot = function(x, n_row = NULL){
  old_pars <- par(no.readonly = TRUE)
  on.exit(par(old_pars), add = TRUE)
  
  if(x$p < 4){plotarray = c(1, x$p)} else {if(is.null(n_row)) n_row = 2 ;plotarray = c(n_row, ceiling(x$p/n_row))}
  par(mfrow=plotarray)
  # par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
  if(is.null(x$whichmethod)){
    x$whichmethod = "fuck"
  }

  if(is.null(colnames(x$X))) x$X = paste0("sm(x", 1:p, ")")
  substr_xnames = function(x, n){
    substr(x, n, nchar(x) - 1)
  }
  xnames = colnames(x$X)

  if(x$whichmethod == "Enumerate"){
    lpytoprob = function(x){
      x = x - max(x)
      x = exp(x)
      x = x / sum(x)
      x
    }
    TMP = as.data.frame(cbind(x$knotnumsAll, lpytoprob(x$lpydiagnosis$lpy)))
    colnames(TMP) = c(paste0("sm", 1:x$p), "py")
    for(i in 1:x$p){
      tmp = aggregate(py~., data = TMP[, c(i, x$p+1)], sum)
      tmp$py = tmp$py / sum(tmp$py)
      plot(tmp, type="h", ylab="")
      # title(main = "Marginal Posterior of Number of Knots")
    }
  } else {
    for(i in 1:x$p){
      plot(table(x$numknots[,i])/ sum(table(x$numknots[,i])),xlab = xnames[i], ylab = "", ylim = c(0, 1),
           lwd = 3, col = "skyblue")
      # title(main = "Marginal Posterior of Number of Knots")
    }
  }
  
  return(invisible(0))
}


## posterior predictive
bayesResiduals = function(fit, 
                          nReplications= 10, 
                          nBins = NULL){
  old_pars <- par(no.readonly = TRUE)
  on.exit(par(old_pars), add = TRUE)
  
  par(mfrow=c(1,1))
  # par(mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
  
  ## get posterior predictives
  y = fit$y
  XLin = fit$XLin
  linears = fit$linears
  fitted = fit$fitted
  N = length(y); M = nrow(linears); PLin = fit$pLin; PSmo = fit$p
  if(is.null(M)) M = length(linears)
  fittedSums = matrix(0, nrow = M, ncol = N)
  start = end = 1;
  for(i in 1:PSmo){
    end = start + N - 1
    fittedSums = fittedSums + fitted[, start:end] 
    start = end + 1
  }
  postPred = linears %*% t(XLin) + fittedSums
  if(is.infinite(nReplications)){
    showIdx = 1:M
  } else showIdx = sample(M, nReplications)
  
  if(fit$family == "gaussian"){
    ## identity link
    etaReps = postPred
    cutoff = max(abs(sweep(etaReps, 2, y)))
    xlim = range(y)
    plot(y-etaReps[1,], y, xlim = xlim, ylim = c(-cutoff, cutoff), type="n",
         xlab = "Observations", ylab = "Residuals",
         main = paste0("Residuals Plot, ", nReplications, " replications" ))
    for(i in showIdx){
      points(y, y-etaReps[i,], pch = ".")
    }
    abline(h=0, col="black")
    
  } else {
    ## count data residuals: do binning
    if(is.null(nBins)) nBins= ceiling(sqrt(N))
    binsIdx = rep(1:nBins, each = N%/%nBins)
    if(length(binsIdx) < N)binsIdx = c(binsIdx, rep(tail(binsIdx,1), N - length(binsIdx)))

    sortIdx = etaSort = postPredSort = ySort = rep(0, N)
    avgStd = avgExp = avgRes = matrix(0, nrow = nReplications, ncol = nBins)
    
    if(fit$family == "binomial"){
      if(fit$link == "logit"){
        linkinv = binomial(link="logit")$linkinv
        etaReps = linkinv(postPred)
      }
      for(i in seq_along(showIdx)){
        ## do binning
        sortIdx = order(postPred[showIdx[i], ])
        etaSort = etaReps[i, sortIdx]
        postPredSort = postPred[i, sortIdx]
        ySort = y[sortIdx]
        for(j in 1:nBins){
          avgExp[i, j] = mean(etaSort[binsIdx == j])
          avgRes[i, j] = mean(ySort[binsIdx ==j] - etaSort[binsIdx == j])
          # avgStd[i, j] = sd(ySort[binsIdx ==j] - etaSort[binsIdx == j])
          avgRes[i, j] = mean(ySort[binsIdx ==j] - etaSort[binsIdx == j])/
            sd(ySort[binsIdx ==j] - etaSort[binsIdx == j])
        }
      }
      
    } else if(fit$family == "poisson"){
      if(fit$link == "log"){
        linkinv = poisson(link = "log")$linkinv
        etaReps = linkinv(postPred)
      }
      
      for(i in seq_along(showIdx)){
        ## do binning
        sortIdx = order(postPred[showIdx[i], ])
        etaSort = etaReps[i, sortIdx]
        postPredSort = postPred[i, sortIdx]
        ySort = y[sortIdx]
        for(j in 1:nBins){
          avgExp[i, j] = mean(etaSort[binsIdx == j])
          avgRes[i, j] = mean(ySort[binsIdx ==j] - etaSort[binsIdx == j]) / sqrt(avgExp[i, j])
        }
      }
    }
    
    cutoff = max(3, max(abs(avgRes)))
    plot(avgExp, avgRes, ylim = c(-cutoff, cutoff),
         type="n",
         xlab = "Average Predicteds per bins", ylab = "Average Scaled Residuals per bins",
         main = paste0("Binnded Residuals Plot, ", nReplications, " replications" ))
    for(i in seq_along(showIdx)){
      points(avgExp[i, ], avgRes[i, ], pch=1)
      # lines(avgExp[i, ], avgStd[i,]*1.96, lty=2)
      # lines(avgExp[i, ], -avgStd[i,]*1.96, lty=2)
    }
    abline(h=1.96, lty=2)
    abline(h=-1.96, lty=2)
    abline(h = 0) 
  }
}

plotTrace = function(fit){
  old_pars <- par(no.readonly = TRUE)
  on.exit(par(old_pars), add = TRUE)
  
  
  
  
  par(mfrow=c(1,2))
  tmp = tryCatch(acf(fit$lpydiagnosis$lpy),
                 error = function(cnd) cnd)
  if(inherits(tmp, "error")){
    tmp = rep(0, 20)
    plot(tmp)
    ess = 1
  } else {
    ess = length(fit$lpydiagnosis$lpy)/(1+2*sum(tmp$acf))
  }
  
  plot(fit$lpydiagnosis$lpy, type="l",
       main = paste0("N ", length(fit$lpydiagnosis$lpy),
                     ", ESS ", signif(ess, 4)))
  return(invisible(ess))
}
































