library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("src/gambma.cpp")

gpriors = 1:8
PRIOR = c("Unit", "Uniform","Hyper-g", "Hyper-g/n", "Beta-prime", "ZS-adapted","Robust", "Intrinsic")
priorCode = Vectorize(function(prior){
  # c("Unit", "Hyper-g", "Uniform", "Hyper-g/n", "Beta-prime", "ZS-adapted", "Robust", "Intrinsic", "constant", "CH")
  switch(prior, 
         "Unit" = "unit", "Hyper-g" = "hyperg", "Uniform" = "uniform", "Hyper-g/n" = "hypergn", "Beta-prime" = "betaprime", 
         "ZS-adapted" = "zsadapted", "Robust" = "robust", "Intrinsic" = "intrinsic", "constant" = "constant", "CH" = "ch",
         NULL)
}, vectorize.args = "prior")
cols = c("darkgrey", "#238b45", "#66c2a4", "#74a9cf", "#d94701", "#fdbe85", "#6a51a3", "#9e9ac8")
ltys = c(1, rep(3,2),rep(2,5))

## tCCH prior
dtCCH = function(uu, a, b, r, s, nu, kap, logValue=F){
  # if(min(uu) < 0 || max(uu) > 1/nu) stop("support is [0, 1/nu]")
  if(!(a>0)) stop("must be a > 0")
  if(!(b>0)) stop("must be b > 0")
  # if(!(nu>=1)) stop("must be nu >= 1")
  if(!(nu>0)) stop("must be nu > 0")
  if(!(kap>0)) stop("must be kap > 0")
  
  ## support
  suppidx = (uu < 1/nu && uu > 0)
  logdens = rep(NA_real_, length(uu))
  uu_in = uu[suppidx]
  
  
  ## normalizing constant
  lognormconst = a/2*log(nu) + s/(2*nu) - base::lbeta(a/2, b/2) - logPhi1_cpp(b/2, r, (a+b)/2, s/(2*nu), 1-kap)
  ## density
  logdens[suppidx] = 
    (a/2-1)*log(uu[suppidx]) + (b/2-1)*log1p(1-nu*uu[suppidx]) -
    s/2*uu[suppidx] - r*log(kap + (1-kap)*nu*uu) + lognormconst
  
  if(logValue) return(logdens)
  return(exp(logdens))
}

##############################################################################################################
mdEvid = function(likelihood, gprior,
                  n, p, r2, 
                  a, b, r, s, nu, kap){
  if(hasArg(gprior)){
    if(gprior == 1) {a=b=r=s=nu=kap=0L; g=n} # fixed
    if(gprior == 2) {a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1} # uniform
    if(gprior == 3) {a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1} # hyperg
    if(gprior == 4) {a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n} # hypergn
    if(gprior == 5) {a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1} # betaprime
    if(gprior == 6) {a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1} # zsadapted
    if(gprior == 7) {a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1} # robust
    if(gprior == 8) {a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n} # intrinsic
  } 
  if(gprior != 1){
    if(!(n>0)) stop("must be n > 0")
    if(!(p>0)) stop("must be p > 0")
    if(!(a>0)) stop("must be a > 0")
    if(!(b>0)) stop("must be b > 0")
    if(!(nu>0)) stop("must be nu > 0")
    if(!(kap>0)) stop("must be kap > 0")
  }
  
  res = NA_real_
  
  if(likelihood == "gaussian"){
    if(!(r2>=0 && r2<=1)) stop("must be 0 <= R2 <= 0")
    if(!(kap==1L || s == 0L)) stop("must be either kap = 1 or s = 0")
    
    if(gprior==1L){
      res = (n-p-1)/2 *log1p(g) - (n-1)/2*log1p(g*(1-r2))
    } else {
      part_betas = base::lbeta((a+p)/2, b/2) - base::lbeta(a/2, b/2) 
      if(gprior == 5L){
        res = part_betas + (p-n+1.5)*log1p(-r2)/2
      } else if(kap == 1L){
        if(s == 0L){
          part_log2F1 = log2F1_cpp((n-1)/2, b/2, (a+b+p)/2, r2/(nu-(nu-1)*r2))
          res = part_betas + part_log2F1 - p/2*log(nu) - (n-1)/2 * log1p(- (1-1/nu)*r2)
        } else {
          part_log1F1 = log1F1_cpp(b/2, (a+b)/2, s/(2*nu))
          part_logPhi1 = logPhi1_cpp(b/2, (n-1)/2, (a+b+p)/2, s/(2*nu), r2 / (nu - (nu-1)*r2))
          res = part_betas + part_logPhi1 - part_log1F1 -p/2*log(nu) - (n-1)/2 * log1p(- (1-1/nu)*r2)
        }
      } else if (s == 0L) {
        part_log2F1 = log2F1_cpp(r, b/2, (a+b)/2, 1-kap)
        part_logF1 = logF1_cpp((a+p)/2, (a+b+p+1-n-2*r)/2, (n-1)/2, (a+b+p)/2, 1-kap, 1-kap-r2*kap/((1-r2)*nu))
        res = part_betas + part_logF1 - part_log2F1 -p/2*log(nu) - (n-1)/2 * log1p(- r2) + (a+p-2*r)/2 *log(kap)
      } else stop("must be either kap = 1 or s = 0")
    }
    
  } else if (likelihood == "exponentials"){
    Qm = -n*log(1-r2)
    if(!(Qm>0)) stop("must be Qm > 0")
    
    if(gprior==1L){
      res = -0.5*Qm/(g+1) - 0.5*p*log(1+g)
    } else {
      part_betas = base::lbeta((a+p)/2, b/2) - base::lbeta(a/2, b/2) 
      if(kap==1L){
        part_log1F1_1 = log1F1_cpp((a+p)/2, (a+b+p)/2, -(s+Qm)/(2*nu))
        if(s==0){
          res = part_betas + part_log1F1_1 - (p/2)*log(nu) 
        } else {
          part_log1F1_2 = log1F1_cpp((a)/2, (a+b)/2, -(s)/(2*nu))
          res = part_betas + part_log1F1_1 - (p/2)*log(nu) -part_log1F1_2
        }
      } else {
        part_logPhi1_1 = logPhi1_cpp(b/2, r, (a+b+p)/2, (s+Qm)/(2*nu), 1-kap)
        if(s==0L){
          part_log2F1 = log2F1_cpp(r, b/2, (a+b)/2, 1-kap)
          res = part_betas + part_logPhi1_1 - part_log2F1 - p*log(nu)/2 - Qm/(2*nu)
        } else {
          part_logPhi1_2 = logPhi1_cpp(b/2, r, (a+b)/2, (s)/(2*nu), 1-kap)
          res = part_betas + part_logPhi1_1 - part_logPhi1_2 - p*log(nu)/2 - Qm/(2*nu)
        }
      }
    }
    
  } else stop("unsupported likelihood")
  
  return(res)
}

##############################################################################################################
##############################################################################################################
## PLOT JOINT (CONTOUR PLOT)
##############################################################################################################
plotJoint = function(likelihood = "exponentials", gpriors = 1:8, plotwhat = "logBF",
                     n, p, totalp,
                     nplot_p, nplot_r2, nplot_r2_range = c(0.001, 0.999)){
  plotxlab = function(likelihood){
    switch(likelihood,
           "exponentials" = "R2pseudo",
           "gaussian" = "R2")
  }
  plotlist= list()
  pp = 1:nplot_p
  R2 = seq(nplot_r2_range[1], nplot_r2_range[2], len = nplot_r2)
  for(k in seq_along(gpriors)){
    plotmat = matrix(NA, nrow=nplot_p, ncol=nplot_r2)
    for(i in 1:nplot_p){
      for(j in 1:nplot_r2){
        if(plotwhat == "logBF"){
          tmp = tryCatch(mdEvid(likelihood, gprior=gpriors[k], n=n, p=pp[i], r2 = R2[j])-
                      mdEvid(likelihood, gprior=gpriors[k], n=n, p=pp[i]+1, r2 = R2[j]),
                      error = function(e) e)
        } else if (plotwhat == "logME"){
          tmp = tryCatch(mdEvid(likelihood, gprior=gpriors[k], n=n, p=pp[i], r2 = R2[j]),
                         error = function(e) e)
        }
        if(!is.finite(tmp)) tmp = NA
        if(!inherits(tmp, "error")) plotmat[i,j] = tmp
      }
    }
    plotlist[[k]] = plotmat
    print(gpriors[k])
  }
  
  return(plotlist)
}

##############################################################################################################
##############################################################################################################
## exponentials
n=200; p=5;
nplot_p = 50
totalp = 50;
nplot_r2 = 100; 
nplot_r2_range = c(0.001, 0.999)
pp = 1:nplot_p
R2 = seq(nplot_r2_range[1], nplot_r2_range[2], len = nplot_r2)

gridPR2 = expand.grid(pp, R2)
logBF = plotJoint("exponentials", gpriors, "logBF",
                n, p, totalp, 
                nplot_p, nplot_r2, nplot_r2_range)
range(logBF)
logBFdf = cbind(sapply(logBF, c), gridPR2)
colnames(logBFdf) = c(priorCode(PRIOR), "p", "R2")
head(logBFdf)

logME = plotJoint("exponentials", gpriors, "logME",
                  n, p, totalp, 
                  nplot_p, nplot_r2, nplot_r2_range)
range(logME)
logMEdf = cbind(sapply(logME, c), gridPR2)
colnames(logMEdf) = c(priorCode(PRIOR), "p", "R2")
head(logMEdf)

##############################################################################################################
## PLOT
nlevels = 10
filledContour = function (x = seq(0, 1, length.out = nrow(z)), 
                          y = seq(0, 1, length.out = ncol(z)), 
                          z, 
                          xlim = range(x, finite = TRUE), 
                          ylim = range(y, finite = TRUE), 
                          zlim = range(z, finite = TRUE), 
                          levels = pretty(zlim, nlevels), 
                          nlevels = 20, 
                          color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE), 
                          col = color.palette(length(levels) - 1), 
                          plot.title, plot.axes, key.title, key.axes, asp = NA, 
                          xaxs = "i", yaxs = "i", las = 1, axes = TRUE, 
                          frame.plot = axes, ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  
  
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, 
              asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}
library(latex2exp)
##############################################################################################################
plotlist = logBF
zlim = range(plotlist)
levels = pretty(zlim, nlevels)
color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
col = color.palette(length(levels) - 1)

length(plotlist)

pdf("logBFExp.pdf", width = 12, height = 6)
par(mar=c(2.5, 2.5, 2, 1/2),mgp=c(1.5,0.5,0), oma = c(0,0,1.5,0), cex.axis = 1.0, cex.main= 1)
layout(matrix(c(1,2,3,4, 9,
                5,6,7,8, 10), byrow=T, nrow=2), widths = c(rep(1,4), lcm(w)))
par(mar=c(2.5, 2.5, 2.5, 1.5), mgp=c(1.5,0.5,0), oma = c(0,0,2,0), cex.axis = 1.0, cex.main= 1)
for(i in seq_along(plotlist)){
  filledContour(R2, pp, t(plotlist[[i]]), zlim=zlim,
                main = PRIOR[i], xlab= latex2exp::TeX("$R_{\\xi_1, pseudo}^2$"), ylab = latex2exp::TeX("$J_{\\xi_1}$"))
}
par(mar=c(2, 1/2, 2, 3), mgp=c(1.5,1,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
axis(4)
mtext("log Bayes Factor", outer= T, font=2)
dev.off()

##############################################################################################################
plotlist = logME
zlim = range(plotlist)
levels = pretty(zlim, nlevels)
color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
col = color.palette(length(levels) - 1)

length(plotlist)

pdf("logMEExp.pdf", width = 12, height = 6)
par(mar=c(2.5, 2.5, 2, 1/2),mgp=c(1.5,0.5,0), oma = c(0,0,1.5,0), cex.axis = 1.0, cex.main= 1)
layout(matrix(c(1,2,3,4, 9,
                5,6,7,8, 10), byrow=T, nrow=2), widths = c(rep(1,4), lcm(w)))
for(i in seq_along(plotlist)){
  filledContour(R2, pp, t(plotlist[[i]]), zlim=zlim,
                main = PRIOR[i], xlab= latex2exp::TeX("$R_{\\xi_1, pseudo}^2$"), ylab = latex2exp::TeX("$J_{\\xi_1}$"))
}
par(mar=c(2, 1/2, 2, 3), mgp=c(1.5,1,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
axis(4)
mtext("log Model Evidence", outer= T, font=2)
dev.off()

##############################################################################################################
##############################################################################################################
## gaussian
n=200; p=5;
nplot_p = 50
totalp = 50;
nplot_r2 = 100; 
nplot_r2_range = c(0.01, 0.98)
pp = 1:nplot_p
R2 = seq(nplot_r2_range[1], nplot_r2_range[2], len = nplot_r2)

gridPR2 = expand.grid(pp, R2)
logBF = plotJoint("gaussian", gpriors, "logBF",
                  n, p, totalp, 
                  nplot_p, nplot_r2, nplot_r2_range)
range(logBF)
logBFdf = cbind(sapply(logBF, c), gridPR2)
colnames(logBFdf) = c(priorCode(PRIOR), "p", "R2")
head(logBFdf)

logME = plotJoint("gaussian", gpriors, "logME",
                  n, p, totalp, 
                  nplot_p, nplot_r2, nplot_r2_range)
range(logME)
logMEdf = cbind(sapply(logME, c), gridPR2)
colnames(logMEdf) = c(priorCode(PRIOR), "p", "R2")
head(logMEdf)

##############################################################################################################
## PLOT
nlevels = 10
##############################################################################################################
plotlist = logBF
zlim = range(plotlist)
levels = pretty(zlim, nlevels)
color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
col = color.palette(length(levels) - 1)

length(plotlist)

pdf("logBFGauss.pdf", width = 12, height = 6)
par(mar=c(2.5, 2.5, 2, 1/2),mgp=c(1.5,0.5,0), oma = c(0,0,1.5,0), cex.axis = 1.0, cex.main= 1)
layout(matrix(c(1,2,3,4, 9,
                5,6,7,8, 10), byrow=T, nrow=2), widths = c(rep(1,4), lcm(w)))
for(i in seq_along(plotlist)){
  filledContour(R2, pp, t(plotlist[[i]]), zlim=zlim,
                main = PRIOR[i], xlab= latex2exp::TeX("$R_{\\xi_1, pseudo}^2$"), ylab = latex2exp::TeX("$J_{\\xi_1}$"))
}
par(mar=c(2, 1/2, 2, 3), mgp=c(1.5,1,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
axis(4)
mtext("log Bayes Factor", outer= T, font=2)
dev.off()
##############################################################################################################
plotlist = logME
zlim = range(plotlist)
levels = pretty(zlim, nlevels)
color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE)
col = color.palette(length(levels) - 1)

length(plotlist)

pdf("logMEGauss.pdf", width = 12, height = 6)
par(mar=c(2.5, 2.5, 2, 1/2),mgp=c(1.5,0.5,0), oma = c(0,0,1.5,0), cex.axis = 1.0, cex.main= 1)
layout(matrix(c(1,2,3,4, 9,
                5,6,7,8, 10), byrow=T, nrow=2), widths = c(rep(1,4), lcm(w)))
for(i in seq_along(plotlist)){
  filledContour(R2, pp, t(plotlist[[i]]), zlim=zlim,
                main = PRIOR[i], xlab= latex2exp::TeX("$R_{\\xi_1}^2$"), ylab = latex2exp::TeX("$J_{\\xi_1}$"))
}
par(mar=c(2, 1/2, 2, 3), mgp=c(1.5,1,0))
plot.new()
plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
axis(4)
mtext("log Model Evidence", outer= T, font=2)
dev.off()

