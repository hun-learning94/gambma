library(latex2exp)
source("tests/specialfunctions.R")
plotdir = "C:/Users/USER/Dropbox/Kang/GAM6/plots/"
# a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1 # uniform
# a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1 # hyperg
# a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
# a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1 # betaprime
# a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1 # zsadapted
# a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1 # robust
# a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic

gpriors = 0:7
cols = c("darkgrey", "#238b45", "#66c2a4", "#74a9cf", "#d94701", "#fdbe85", "#6a51a3", "#9e9ac8")
ltys = c(4, 1,1,2,2,2,3,3)

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
  lognormconst = a/2*log(nu) + s/(2*nu) - base::lbeta(a/2, b/2) - logPhi1(b/2, r, (a+b)/2, s/(2*nu), 1-kap)
  ## density
  logdens[suppidx] = 
    (a/2-1)*log(uu[suppidx]) + (b/2-1)*log1p(1-nu*uu[suppidx]) - s/2*uu[suppidx] - r*log(kap + (1-kap)*nu*uu) + lognormconst
  
  if(logValue) return(logdens)
  return(exp(logdens))
}

#################################################################################
mdEvid = function(likelihood, gprior,
                  n, p, modelFit, 
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
    r2 = modelFit
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
          part_log2F1 = log2F1((n-1)/2, b/2, (a+b+p)/2, r2/(nu-(nu-1)*r2))
          res = part_betas + part_log2F1 - p/2*log(nu) - (n-1)/2 * log1p(- (1-1/nu)*r2)
        } else {
          part_log1F1 = log1F1(b/2, (a+b)/2, s/(2*nu))
          part_logPhi1 = logPhi1(b/2, (n-1)/2, (a+b+p)/2, s/(2*nu), r2 / (nu - (nu-1)*r2))
          res = part_betas + part_logPhi1 - part_log1F1 -p/2*log(nu) - (n-1)/2 * log1p(- (1-1/nu)*r2)
        }
      } else if (s == 0L) {
        part_log2F1 = log2F1(r, b/2, (a+b)/2, 1-kap)
        part_logF1 = logF1((a+p)/2, (a+b+p+1-n-2*r)/2, (n-1)/2, (a+b+p)/2, 1-kap, 1-kap-r2*kap/((1-r2)*nu))
        res = part_betas + part_logF1 - part_log2F1 -p/2*log(nu) - (n-1)/2 * log1p(- r2) + (a+p-2*r)/2 *log(kap)
      } else stop("must be either kap = 1 or s = 0")
    }
    
  } else if (likelihood == "exponentials"){
    Qm = modelFit
    if(!(Qm>0)) stop("must be Qm > 0")
    
    if(gprior==0L){
      res = -0.5*Qm/(g+1) - 0.5*p*log(1+g)
    } else {
      part_betas = base::lbeta((a+p)/2, b/2) - base::lbeta(a/2, b/2) 
      if(kap==1L){
        part_log1F1_1 = log1F1((a+p)/2, (a+b+p)/2, -(s+Qm)/(2*nu))
        if(s==0){
          res = part_betas + part_log1F1_1 - (p/2)*log(nu) 
        } else {
          part_log1F1_2 = log1F1((a)/2, (a+b)/2, -(s)/(2*nu))
          res = part_betas + part_log1F1_1 - (p/2)*log(nu) -part_log1F1_2
        }
      } else {
        part_logPhi1_1 = logPhi1(b/2, r, (a+b+p)/2, (s+Qm)/(2*nu), 1-kap)
        if(s==0L){
          part_log2F1 = log2F1(r, b/2, (a+b)/2, 1-kap)
          res = part_betas + part_logPhi1_1 - part_log2F1 - p*log(nu)/2 - Qm/(2*nu)
        } else {
          part_logPhi1_2 = logPhi1(b/2, r, (a+b)/2, (s)/(2*nu), 1-kap)
          res = part_betas + part_logPhi1_1 - part_logPhi1_2 - p*log(nu)/2 - Qm/(2*nu)
        }
      }
    }
    
  } else stop("unsupported likelihood")
  
  return(res)
}
# mdEvid = Vectorize(mdEvid,
#                    vectorize.args = c("p", "modelFit", "b", "nu", "kap"))

plotlist = list()
likelihood = "exponentials"; n=200; p=20;
nplot_p = 50
nplot_qm = 100
pp = 1:nplot_p
Qm = seq(0.1, 10, len = nplot_qm)*n
for(k in seq_along(gpriors)){
  plotmat = matrix(0, nrow=nplot_p, ncol=nplot_qm)
  for(i in 1:nplot_p){
    for(j in 1:nplot_qm){
      plotmat[i,j] = 
        mdEvid(likelihood, gprior=gpriors[k], n=n, p=pp[i], modelFit = Qm[j])-
        mdEvid(likelihood, gprior=gpriors[k], n=n, p=pp[i]+1, modelFit = Qm[j])
    }
  }
  plotlist[[k]] = plotmat
  print(gpriors[k])
}
# persp(pp, R2, plotmat)
library(plotly)
fig = plot_ly(showscale=F)
for(i in seq_along(gpriors)){
  fig = fig %>% add_surface(x = Qm, y= pp,z = plotlist[[i]], colorscale = list(c(0, 1), c(cols[i], cols[i])))
}
fig

plotlist= list()
likelihood = "gaussian"; n=100; p=20;
nplot_p = 50
nplot_qm = 100
pp = 1:nplot_p
R2 = seq(0.2, 0.9, len = nplot_qm)
for(k in seq_along(gpriors)){
  plotmat = matrix(0, nrow=nplot_p, ncol=nplot_qm)
  for(i in 1:nplot_p){
    for(j in 1:nplot_qm){
      plotmat[i,j] = 
        mdEvid(likelihood, gprior=gpriors[k], n=n, p=pp[i], modelFit = R2[j]) -
        mdEvid(likelihood, gprior=gpriors[k], n=n, p=pp[i]+1, modelFit = R2[j])
    }
  }
  plotlist[[k]] = plotmat
  print(gpriors[k])
}
# persp(pp, R2, plotmat)
library(plotly)
fig = plot_ly(showscale=F)
for(i in seq_along(gpriors)){
  fig = fig %>% add_surface(x = Qm, y= pp,z = plotlist[[i]], colorscale = list(c(0, 1), c(cols[i], cols[i])))
}
fig

##############################################################################################################

plotBFMarg = function(likelihood, gpriors, cols, ltys,
                    xaxis, xlab, yaxis, ylab, diffplots, plotname, 
                    ylim, legendpos){
  plotgpriors = Vectorize(function(gpriors){
    switch(as.character(gpriors), 
           "0" = "g = n",
           "1" = "Uniform",
           "2" = "Hyper-g",
           "3" = "Hyper-g/n",
           "4" = "Beta-prime",
           "5" = "ZS-adapted",
           "6" = "Robust",
           "7" = "Intrinsic")
  }, vectorize.args = "gpriors")
  
  plotxlab = function(xlab){
    switch(xlab,
           "Qm" = latex2exp::TeX("$Q_{\\xi_1}/n$"),
           "R2" = latex2exp::TeX("$R_{\\xi_1}^2$"),
           "p" = latex2exp::TeX("$J_{\\xi_1}$"))
  }
  
  plot_nrow = length(yaxis)
  plot_ncol = length(diffplots)
  par(mfcol=c(plot_nrow, plot_ncol))
  par(mar=c(3, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)
  plotlist = vector("list", plot_nrow*plot_ncol)
  for(i in seq_along(diffplots)){
    for(j in seq_along(yaxis)){
      plotmat = matrix(0, nrow= length(gpriors), ncol = length(xaxis))
      for(k in seq_along(gpriors)){
        for(l in seq_along(xaxis)){
          if(likelihood == "exponentials"){
            if(xlab == "Qm"){
              plotmat[k, l] = mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = yaxis[j], modelFit = xaxis[l]*diffplots[i])
              plotmat[k,l] = plotmat[k,l]-mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = yaxis[j]+1, modelFit = xaxis[l]*diffplots[i])
            } else if (xlab == "p"){
              plotmat[k, l] = mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = xaxis[l], modelFit = yaxis[j]*diffplots[i])
              plotmat[k,l]=plotmat[k,l]-mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = xaxis[l]+1, modelFit = yaxis[j]*diffplots[i])
            }
          } else if(likelihood == "gaussian"){
            if(xlab=="R2"){
              plotmat[k, l] = mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = yaxis[j], modelFit = xaxis[l])
              plotmat[k,l] = plotmat[k,l]-mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = yaxis[j]+1, modelFit = xaxis[l])
            } else if (xlab == "p"){
              plotmat[k, l] = mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = xaxis[l], modelFit = yaxis[j])
              plotmat[k,l]=plotmat[k,l]-mdEvid(likelihood, gprior = gpriors[k], n = diffplots[i], p = xaxis[l]+1, modelFit = yaxis[j])
            }
          }
        }
      }
      if(likelihood == "exponentials"){
        if(xlab == "Qm"){
          title_ = latex2exp::TeX(paste0("$n=", diffplots[i], "$, ", "$J_{\\xi_1}=", yaxis[j], "$"))
        } else if (xlab == "p"){
          title_ = latex2exp::TeX(paste0("$n=", diffplots[i], "$, ", "$Q_{\\xi_1}/n=", yaxis[j], "$"))
        }
      } else if(likelihood == "gaussian"){
        if(xlab=="R2"){
          title_ = latex2exp::TeX(paste0("$n=", diffplots[i], "$, ", "$J_{\\xi_1}=", yaxis[j], "$"))
        } else if (xlab == "p"){
          title_ = latex2exp::TeX(paste0("$n=", diffplots[i], "$, ", "$R_{\\xi_1}^2=", yaxis[j], "$"))
        }
      }
      matplot(t(plotmat), type="l", lty=(ltys), col= (cols), lwd= 1.5, ylim=ylim, xlab = plotxlab(xlab), ylab = "log Bayes factor", xaxt="n",
              main = title_)
      xidx = seq(1, ncol(plotmat), by=10)
      if(xlab=="Qm"){
        axis(1, at = xidx, label = round(xaxis[xidx],2))
      } else axis(1, at = xidx, label = round(xaxis[xidx],2))
      plotlist[[(i-1)*length(yaxis) + j]] = plotmat
    }
  }
  legend(legendpos, legend = plotgpriors(gpriors), lty = ltys, lwd = 1.5, col = cols, ncol=3,
         bty="n", cex=0.8, y.intersp = 3, text.width = 12)
  plotlist
}

cols = c("darkgrey", "#238b45", "#66c2a4", "#74a9cf", "#d94701", "#fdbe85", "#6a51a3", "#9e9ac8")
ltys = c(4, 1,1,2,2,2,3,3)


likelihood = "exponentials"
nplot_p = 50
nplot_qm = 100
N = c(200, 500, 1000)
pp_yaxis = c(5, 30)
pp_xaxis = 1:nplot_p
Qm_yaxis = c(1/10, 1)
Qm_xaxis = seq(0.01, 1, len = nplot_qm)
# pdf(paste0(plotdir,"log_BF_exp_Qm.pdf"), width=8, height = 4.2)
tmp1 = plotBFMarg(likelihood, gpriors, cols, ltys,
               xaxis = Qm_xaxis, xlab = "Qm", yaxis = pp_yaxis, ylab = "p", diffplots = N, plotname = "N",ylim= c(0,4), "")
tmp1plot = recordPlot()
tmp1plot
# dev.off()
# pdf(paste0(plotdir,"log_BF_exp_p.pdf"), width=8, height = 4.2)
tmp2 = plotBFMarg(likelihood, gpriors, cols, ltys,
                xaxis = pp_xaxis, xlab = "p", yaxis = Qm_yaxis, ylab = "Qm", diffplots = N, plotname = "N",ylim= c(0,4), "bottomright")
tmp2plot = recordPlot()
tmp2plot
# dev.off()

likelihood = "gaussian"
nplot_p = 50
nplot_r2 = 100
N = c(200, 500, 1000)
pp_yaxis = c(5, 30)
pp_xaxis = 1:nplot_p
R2_yaxis = c(0.2, 0.8)
R2_xaxis = seq(0.05, 0.95, len = nplot_r2)
# pdf(paste0(plotdir,"log_BF_norm_r2.pdf"), width=8, height = 4.2)
tmp1 = plotBFMarg(likelihood, gpriors, cols, ltys,
                xaxis = R2_xaxis, xlab = "R2", yaxis = pp_yaxis, ylab = "p", diffplots = N, plotname = "N",ylim= c(0,4), "")
tmp1plot_gau = recordPlot()
# dev.off()
# pdf(paste0(plotdir,"log_BF_norm_p.pdf"), width=8, height = 4.2)
tmp2 = plotBFMarg(likelihood, gpriors, cols, ltys,
                xaxis = pp_xaxis, xlab = "p", yaxis = R2_yaxis, ylab = "R2", diffplots = N, plotname = "N",ylim= c(0,4), "bottomright")
tmp2plot_gau = recordPlot()
# dev.off()

pdf(paste0(plotdir,"log_BF_exp_Qm.pdf"), width=8, height = 4.2)
tmp1plot
dev.off()
pdf(paste0(plotdir,"log_BF_exp_p.pdf"), width=8, height = 4.2)
tmp2plot
dev.off()
pdf(paste0(plotdir,"log_BF_norm_r2.pdf"), width=8, height = 4.2)
tmp1plot_gau
dev.off()
pdf(paste0(plotdir,"log_BF_norm_p.pdf"), width=8, height = 4.2)
tmp2plot_gau
dev.off()







