setwd("D:/BVSGAM_internal")
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp("src/gambma.cpp")
library(latex2exp)
plotdir = paste0("")
# a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1 # uniform
# a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1 # hyperg
# a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
# a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1 # betaprime
# a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1 # zsadapted
# a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1 # robust
# a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic

gpriors = 0:7
cols = c("darkgrey", "#238b45", "#66c2a4", "#74a9cf", "#d94701", "#fdbe85", "#6a51a3", "#9e9ac8")
ltys = c(1, 
         rep(3,2),rep(2,5))


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
    (a/2-1)*log(uu[suppidx]) + (b/2-1)*log1p(1-nu*uu[suppidx]) - s/2*uu[suppidx] - r*log(kap + (1-kap)*nu*uu) + lognormconst
  
  if(logValue) return(logdens)
  return(exp(logdens))
}

##############################################################################################################
mdEvid = function(likelihood, gprior,
                      n, p, r2, 
                      a, b, r, s, nu, kap){
  if(hasArg(gprior)){
    if(gprior == 1) {a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1} # uniform
    if(gprior == 2) {a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1} # hyperg
    if(gprior == 3) {a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n} # hypergn
    if(gprior == 4) {a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1} # betaprime
    if(gprior == 5) {a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1} # zsadapted
    if(gprior == 6) {a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1} # robust
    if(gprior == 7) {a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n} # intrinsic
    if(gprior == 0) {a=b=r=s=nu=kap=0L; g=n} # fixed
  } 
  if(gprior != 0){
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
    
    if(gprior==0L){
      res = (n-p-1)/2 *log1p(g) - (n-1)/2*log1p(g*(1-r2))
    } else {
      part_betas = base::lbeta((a+p)/2, b/2) - base::lbeta(a/2, b/2) 
      if(gprior == 4L){
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
    
    if(gprior==0L){
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
## plot Bayes Factor in marginal
plotgpriors = Vectorize(function(gpriors){
  switch(as.character(gpriors), 
         "0" = "Unit (g = n)",
         "1" = "Uniform",
         "2" = "Hyper-g",
         "3" = "Hyper-g/n",
         "4" = "Beta-prime",
         "5" = "ZS-adapted",
         "6" = "Robust",
         "7" = "Intrinsic")
}, vectorize.args = "gpriors")
plotBFMarg = function(likelihood, gpriors, cols, ltys,
                    xaxis, xlab, yaxis, ylab, diffplots, plotname, 
                    ylim, legendpos, xtickloc, xticklabel){
  
  
  plotxlab = function(xlab){
    switch(xlab,
           "R2pseudo" = latex2exp::TeX("$R_{\\xi_1, pseudo}^2$"),
           "R2" = latex2exp::TeX("$R_{\\xi_1}^2$"),
           "p" = latex2exp::TeX("$J_{\\xi_1}$"))
  }
  
  plot_nrow = length(yaxis)
  plot_ncol = length(diffplots)
  par(mfcol=c(plot_nrow, plot_ncol))
  par(mar=c(5, 2.2, 1.8, 2),
      mgp=c(3.5,1,0), 
      oma = c(0,0,0.1,0), 
      cex.axis = 1.6, cex.main= 1.6, cex.lab=1.6)
  plotpos = 0
  plotlist = vector("list", plot_nrow*plot_ncol)
  for(i in seq_along(diffplots)){
    for(j in seq_along(yaxis)){
      plotmat = matrix(NA, nrow= length(gpriors), ncol = length(xaxis))
      for(k in seq_along(gpriors)){
        for(l in seq_along(xaxis)){
          if(likelihood == "exponentials"){
            if(xlab == "R2pseudo"){
              tmp = tryCatch(mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = yaxis[j], r2 = xaxis[l])-
                               mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = yaxis[j]+1, r2 = xaxis[l]),
                             error = function(e) e)
              if(!inherits(tmp, "error")){plotmat[k, l] = tmp}
            } else if (xlab == "p"){
              tmp = tryCatch(mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = xaxis[l], r2 = yaxis[j])-
                               mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = xaxis[l]+1, r2 = yaxis[j]),
                             error = function(e) e)
              if(!inherits(tmp, "error")){plotmat[k, l] = tmp}
            }
          } else if(likelihood == "gaussian"){
            if(xlab=="R2"){
              tmp = tryCatch(mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = yaxis[j], r2 = xaxis[l])-
                               mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = yaxis[j]+1, r2 = xaxis[l]),
                             error = function(e) e)
              if(!inherits(tmp, "error")){plotmat[k, l] = tmp}
            } else if (xlab == "p"){
              tmp = tryCatch(mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = xaxis[l], r2 = yaxis[j])-
                               mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = xaxis[l]+1, r2 = yaxis[j]),
                             error = function(e) e)
              if(!inherits(tmp, "error")){plotmat[k, l] = tmp}
            }
          }
        }
      }
      if(likelihood == "exponentials"){
        if(xlab == "R2pseudo"){
          title_ = latex2exp::TeX(paste0("$n=", diffplots[i], "$, ", "$J_{\\xi_1}=", yaxis[j], "$"))
        } else if (xlab == "p"){
          title_ = latex2exp::TeX(paste0("$n=", diffplots[i], "$, ", "$R_{\\xi_1, pseudo}^2=", yaxis[j], "$"))
        }
      } else if(likelihood == "gaussian"){
        if(xlab=="R2"){
          title_ = latex2exp::TeX(paste0("$n=", diffplots[i], "$, ", "$J_{\\xi_1}=", yaxis[j], "$"))
        } else if (xlab == "p"){
          title_ = latex2exp::TeX(paste0("$n=", diffplots[i], "$, ", "$R_{\\xi_1}^2=", yaxis[j], "$"))
        }
      }
      plotpos = plotpos+1
      matplot(t(plotmat), type="l", lty=(ltys), col= (cols), lwd= 2.5, ylim=ylim, xlab = plotxlab(xlab), 
              ylab = "",
              xaxt="n",
              main = title_)
      if(!missing(legendpos)){
        if(plotpos == legendpos)
          legend("bottomright",cex=0.9,
                 legend = plotgpriors(gpriors), ncol=1,
                 lty = ltys, lwd = 3, col = cols, seg.len=2.7)
      }
      if(missing(xtickloc)) xtickloc = seq(1, length(xaxis), by=10)
      if(missing(xticklabel)) xticklabel = round(xaxis[xtickloc],2)
      axis(1, at = xtickloc, label = xticklabel)
      plotlist[[(i-1)*length(yaxis) + j]] = plotmat
    }
  }
  return(invisible(plotlist))
}

gpriors = 0:7
cols = c("darkgrey", "#238b45", "#66c2a4", "#74a9cf", "#d94701", "#fdbe85", "#6a51a3", "#9e9ac8")
ltys = c(1,
         1,2,
         1,2,1,2,1)

nplot_p = 50
nplot_r2 = 200
N = c(200, 1000)
pp_yaxis = c(3, 30)
pp_xaxis = 1:nplot_p
R2_yaxis = c(0.2, 0.8)
R2_xaxis = seq(0, 1-1e-5, len = nplot_r2)


r2tickloc = c(0, 0.2, 0.4, 0.6, 0.8, 1)
ptickloc = c(1, 10, 20, 30, 40, 50)
likelihood = "exponentials"
tmp1 = plotBFMarg(likelihood, gpriors, cols, ltys,
               xaxis = R2_xaxis, xlab = "R2pseudo", yaxis = pp_yaxis, ylab = "p", diffplots = N, plotname = "N",ylim= c(0,4),3, 
               xtickloc = nplot_r2*r2tickloc, xticklabel = r2tickloc)
tmp1plot = recordPlot()
tmp2 = plotBFMarg(likelihood, gpriors, cols, ltys,
                xaxis = pp_xaxis, xlab = "p", yaxis = R2_yaxis, ylab = "Qm", diffplots = N, plotname = "N",ylim= c(0,4),
                xtickloc = ptickloc, xticklabel = ptickloc)
tmp2plot = recordPlot()

likelihood = "gaussian"
tmp1 = plotBFMarg(likelihood, gpriors, cols, ltys,
                xaxis = R2_xaxis, xlab = "R2", yaxis = pp_yaxis, ylab = "p", diffplots = N, plotname = "N",ylim= c(0,5), 
                xtickloc = nplot_r2*r2tickloc, xticklabel = r2tickloc)
tmp1plot_gau = recordPlot()
tmp2 = plotBFMarg(likelihood, gpriors, cols, ltys,
                xaxis = pp_xaxis, xlab = "p", yaxis = R2_yaxis, ylab = "R2", diffplots = N, plotname = "N",ylim= c(0,5), 1,
                xtickloc = ptickloc, xticklabel = ptickloc)
tmp2plot_gau = recordPlot()

pdf(paste0(plotdir,"log_BF_exp_r2pseudo.pdf"),  width=7, height = 7)
# tmp1plot
tmp1 = plotBFMarg(likelihood, gpriors, cols, ltys,
                  xaxis = R2_xaxis, xlab = "R2pseudo", yaxis = pp_yaxis, ylab = "p", diffplots = N, plotname = "N",ylim= c(0,4),3, 
                  xtickloc = nplot_r2*r2tickloc, xticklabel = r2tickloc)
dev.off()
pdf(paste0(plotdir,"log_BF_exp_p.pdf"),  width=7, height = 7)
tmp2plot
# legend("bottomleft",cex=1.3,
#        legend = plotgpriors(gpriors)[1:4], ncol=1,
#        lty = ltys[1:4], lwd = 3, col = cols[1:4], seg.len=2.5)
dev.off()
pdf(paste0(plotdir,"log_BF_norm_r2.pdf"), width=7, height = 7)
tmp1plot_gau
dev.off()
pdf(paste0(plotdir,"log_BF_norm_p.pdf"),  width=7, height = 7)
tmp2plot_gau
dev.off()







