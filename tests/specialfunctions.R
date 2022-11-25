# library(microbenchmark)
# library(ggplot2)
##################################################################################################
## 1F1 function
## Confluent hypergeometric function (Abramowitz and Stegun 1970)
## for r > a > 0,
## 1F1(a, r, x) = 1/Beta(r-a, a) \int_0^1 u^(a-1) (1-u)^(r-a-1) exp(xu) du

log1F1 = function(aa, rr, xx, logValue=T, ...){
  if(!(rr > aa & aa >0)) stop("must be r > a > 0")
  
  tmpf = function(uu, aa, rr, xx){
    (aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) + xx*uu - base::lbeta(rr-aa, aa)
  }
  
  uu = seq(0.001, 0.999, len = 200)
  adjV = max(tmpf(uu, aa, rr, xx))- log(.Machine$double.xmax)/2
  tmpf = function(uu, aa, rr, xx, adjV){
    exp((aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) + xx*uu - base::lbeta(rr-aa, aa) - adjV)
  }
  # plot(uu, tmpf(uu, aa, rr, xx, adjV), type="l")
  tmp = log(stats::integrate(tmpf, lower=0, upper=1, 
                             aa=aa, rr=rr, xx=xx, adjV = adjV,
                             ...)$value) + adjV
  if(!logValue) tmp = exp(tmp)
  return(tmp)
}

# n = 1000
# a=2; b= 2; p = 100; Qm= n*2
# log1F1(b/2, (a+b+p)/2, Qm/2)
# log1F1((a+p)/2, (a+b+p)/2, -Qm/2)+Qm/2
# 
# n = 1000; p = 10; Qm = n*10
# aa = (1+p)/2; rr = (1+2+p)/2; xx= (n+3+Qm)/2
# BAS::hypergeometric1F1(aa, rr, xx, log=T, laplace=T)
# log1F1(aa, rr, xx)
# Rcpp::sourceCpp(paste0("src/gambma_", "grid", "_diag.cpp"))
# log_hyp1f1(aa,rr,xx)
# 
# comp = microbenchmark(
#   BAS::hypergeometric1F1(aa, rr, xx, log=T, laplace=T),
#   log_hyp1f1(aa,rr,xx),
#   log1F1(aa, rr, xx),
#   times = 1000
# )
# comp
# autoplot(comp)


##################################################################################################
## 2F1 function
## Hypergeometric function (Abramowitz and Stegun 1970)
## for r > a > 0,
## 2F1(b, a, r, x) = 1/Beta(r-a, a) \int_0^1 u^(a-1) (1-u)^(r-a-1) (1-xu)^-b du
log2F1 = function(bb, aa, rr, xx, logValue=T, ...){
  if(!(rr > aa & aa >0)) stop("must be r > a > 0")
  
  tmpf = function(uu, bb, aa, rr, xx){
    (aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-xx*uu) - base::lbeta(rr-aa, aa)
  }
  
  uu = seq(0.001, 0.999, len = 200)
  adjV = max(tmpf(uu, bb, aa, rr, xx))- log(.Machine$double.xmax)/2
  tmpf = function(uu, bb, aa, rr, xx, adjV){
    exp((aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-xx*uu) - base::lbeta(rr-aa, aa) - adjV)
  }
  tmp = log(stats::integrate(tmpf, lower=0, upper=1, 
                             bb=bb, aa=aa, rr=rr, xx=xx, adjV = adjV,
                             ...)$value) + adjV
  if(!logValue) tmp = exp(tmp)
  return(tmp)
}

# n = 300; p = 20; r2 = 0.5
# bb = (n-1)/2; aa = 1; rr = (3+p)/2; xx = r2
# BAS::hypergeometric2F1(bb,aa, rr, xx, log=T)
# log2F1(bb, aa, rr, xx)
# Rcpp::sourceCpp(paste0("src/gambma_", "grid", "_diag.cpp"))
# log_hyp2f1(bb,aa,rr,xx)
# 
# bb = (n-1)/2; aa = 1; rr = (3+p)/2; xx = ((p+1)*r2/(n+1)) / (1 - (n-p)*r2/(n+1))
# log2F1(bb, aa, rr, xx) - (n-1)/2 * log(1 - (n-p)*r2/(n+1))
# bb = (n-1)/2; aa = (1+p)/2; rr = (3+p)/2; xx = r2/(r2-1)*(p+1)/(n+1)
# log2F1(bb, aa, rr, xx) - (n-1)/2 * (log1p(-r2))
# 
# comp = microbenchmark(
#   BAS::hypergeometric2F1(bb,aa, rr, xx, log=T),
#   log_hyp2f1(bb,aa,rr,xx),
#   log2F1(bb, aa, rr, xx),
#   times = 1000
# )
# comp
# autoplot(comp)

##################################################################################################
## Phi1 function
## Confluent hypergeometric function of two variables (Gordy, 1998)
## for r > a > 0 and y < 1,
## Phi1(a, b, r, x, y) = 1/Beta(r-a, a) \int_0^1 u^(a-1) (1-u)^(r-a-1) (1-yu)^-b exp(xu) du
# Phi1_integrand = function(uu, aa, bb, rr, xx, yy){
#   return(exp((aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-yy*uu) + xx*uu- base::lbeta(rr-aa, aa)))
# }
# uu = seq(0.001, 0.999, len=200)
# plot(uu, Phi1_integrand(uu, 1, 0.75, 5, 0, 0.99), type="l")
# plot((uu), Phi1_integrand(uu, 1, 2, 1.5, 1000, 0), type="l")
# plot(exp(uu), Phi1_integrand(uu, 1, 2, 1.5, 1000, 0), type="l")
# integrate(Phi1_integrand, lower = 0, upper = 1, 
#           aa=1, bb=2, rr=5, xx=1, yy=0)

logPhi1 = function(aa, bb, rr, xx, yy, logValue = T, ...){
  # ... : arguments to R's integrate
  if(!(yy<=1)) stop("must be y < 1")
  if(!(rr > aa & aa >0)) stop("must be r > a > 0")
  
  tmpf = function(uu, aa, bb, rr, xx, yy){
    (aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-yy*uu) + xx*uu - base::lbeta(rr-aa, aa)
  }
  
  uu = seq(0.001, 0.999, len = 200)
  adjV = max(tmpf(uu, aa, bb, rr, xx, yy))- log(.Machine$double.xmax)/2
  tmpf = function(uu, aa, bb, rr, xx, yy, adjV){
    exp((aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-yy*uu) + xx*uu - base::lbeta(rr-aa, aa) - adjV)
  }
  # plot(uu, tmpf(uu, aa, bb, rr, xx, yy, adjV), type="l")
  
  tmp = log(stats::integrate(tmpf, lower=0, upper=1, 
                             aa=aa, bb=bb, rr=rr, xx=xx, yy=yy, adjV = adjV,
                             ...)$value) + adjV
  if(!logValue) tmp = exp(tmp)
  return(tmp)
}

## Gauss Example
n =40000; p = 30; r2 = 0.5
a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1 # zsadapted

logPhi1((a+p)/2, (n-1)/2, (a+b+p)/2, -s/(2*nu), r2*nu/(r2-1))
logPhi1(b/2, (n-1)/2, (a+b+p)/2, s/(2*nu), r2 / (nu - (nu-1)*r2))

BAS::phi1(1, 2, 1.5, 0, 1 / 100, log=T)
logPhi1(1, 2, 1.5, 0, 1/100)
BAS::phi1(1, 10, 1.5, 17000, 0, log=T)
logPhi1(1, 10, 1.5, 17000, 0)
BAS::phi1(1, 10, 1.5, 18000, 0, log=T)
# 
# 
# comp = microbenchmark(
#   BAS::phi1(1, 10, 1.5, 1000, 0, log=T),
#   logPhi1(1, 10, 1.5, 1000, 0),
#   times = 1000
# )
# comp
# autoplot(comp)


##################################################################################################
## F1 function
## Hypergeometric function of two variables (Appell function) (Weisstein 2019)
## F1(a,b,bp,r,x,y) = 1/Beta(r-a, a) \int_0^1 u^(a-1) (1-u)^(r-a-1) (1-xu)^-b (1-yu)^-bp du
# F1_integrand = function(uu, aa, bb, bbp, rr, xx, yy){
#   return(exp((aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-xx*uu) - bbp*log(1-yy*uu)- base::lbeta(rr-aa, aa)))
# }
# uu = seq(0.001, 0.999, len=200)
# ## hypergn setting
# n = 200; p = 20; r2 = 0.7
# aa = (1+p)/2; bb = (1+p-n)/2; bbp = (n-1)/2; rr = (3+p)/2; xx = 1-1/n; yy = 1- 1/(n*(1-r2))
# plot(uu, F1_integrand(uu, aa, bb, bbp, rr, xx, yy), type="l")
# integrate(F1_integrand, lower = 0, upper = 1, 
#           aa=aa, bb=bb, bbp=bbp, rr=rr, xx=xx, yy=yy)

logF1 = function(aa, bb, bbp, rr, xx, yy, logValue = T, ...){
  # ... : arguments to R's integrate
  if(!(rr > aa & aa >0)) stop("must be r > a > 0")
  
  tmpf = function(uu, aa, bb, bbp, rr, xx, yy){
    (aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-xx*uu) - bbp*log(1-yy*uu)- base::lbeta(rr-aa, aa)
  }
  
  uu = seq(0.001, 0.999, len = 200)
  adjV = max(tmpf(uu, aa, bb, bbp, rr, xx, yy))- log(.Machine$double.xmax)/2
  tmpf = function(uu, aa, bb, bbp, rr, xx, yy, adjV){
    exp((aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-xx*uu) - bbp*log(1-yy*uu)- base::lbeta(rr-aa, aa) - adjV)
  }
  
  # plot(uu, tmpf(uu, aa, bb, bbp, rr, xx, yy, adjV), type="l")
  tmp = log(stats::integrate(tmpf, lower=0, upper=1, 
                             aa=aa, bb=bb, bbp=bbp, rr=rr, xx=xx, yy=yy, adjV = adjV,
                             ...)$value) + adjV
  if(!logValue) tmp = exp(tmp)
  return(tmp)
}

# ## Gauss Example
# n =1000; p = 30; r2 = 0.8
# a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic
# # a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
# 
# logF1((a+p)/2, r, (n-1)/2, (a+b+p)/2, 1-1/kap, r2*nu/(r2-1))
# logF1((a+p)/2, (a+b+p+1-n-2*r)/2, (n-1)/2, (a+b+p)/2, 1-kap, 1-kap-r2*kap/((1-r2)*nu))

## F1(a,b,bp,r,0,y) = 2F1(bp, a, r, y)
## F1(a,b,bp,r,x,0) = 2F1(b, a, r, x)
# BAS::hypergeometric2F1(12, 1, 2, .65, log=T)
# logF1(1,1,12,2,0,0.65)
# logF1(1,12,2,2,0.65,0)
# 
# dGH = function(u, a, b, x, z, y, w){
#   exp((a-1)*log(u) + (b-1)*log1p(-u) - z*log1p(x*u) - w*log1p(y*u) - lbeta(a,b) - logF1(a,z,w,a+b,-x,-y))
# }
# momentf = function(r=1, a, b, x, z, y, w){
#   exp(lbeta(r+a, b) - lbeta(a,b) + logF1(a+r,z,w,a+b+r,-x,-y) - logF1(a,z,w,a+b,-x,-y))
# }
# 
# uu = seq(0.001, 0.999, len=200)
# plot(uu, dGH(uu, 2,2, 100.1, 2, 100, 2), type="l")
# n = 200; p = 30; r2 = 0.9
# 
# a_pos = (a+p)/2; b_pos = b/2; x_pos = r2/nu/(1-r2); z_pos = (n-1)/2; y_pos = (1-kap)/kap; w_pos = r
# 
# 












