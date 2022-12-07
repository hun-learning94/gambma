library(ars)
# ################################################################################
# ## tCCH distribution
# ## u = 1/(1+g) ~ tCCH(a/2, b/2, r, s/2, nu, kap)
# ## u = 1/(1+g) \in (0, 1/nu)
# 
# ## we use x = -log(u) = log(g+1)
# x_from_g = function(g) { log(g+1) }
# ## g = exp(y) - 1
# g_from_x = function(x){ exp(x) - 1 }
# 
# logdens_tCCH = function(x, a_, b_, r_, s_, nu_, kap_){
#   -a_*x/2 + (b_/2-1)*log1p(-nu_*exp(-x)) - s_*exp(-x)/2 - r_*log(kap_ + (1-kap_)*nu_*exp(-x))
# }
# dlogdens_tCCH = function(x, a_, b_, r_, s_, nu_, kap_){
#   -a_/2 + (b_/2-1) * nu_ / (exp(x)-nu_) + s_*exp(-x)/2 + r_*(1-kap_)*nu/(kap_*exp(x) + (1-kap_)*nu_)
# }
# ddlogdens_tCCH = function(x, a_, b_, r_, s_, nu_, kap_){
#   -nu_*(b_/2-1)*exp(x)/(exp(x)-nu_)^2 - s_*exp(-x)/2 - r_*(1-kap_)*nu_*kap_*exp(x)/(kap_*exp(x) + (1-kap_)*nu)^2
# }
# 
# ##
# rtCCH = function(nsamp, start, 
#                  a_, b_, r_, s_, nu_, kap_){
#   if(!(a_>0)) stop("must be a > 0")
#   if(!(b_>0)) stop("must be b > 0")
#   if(!(nu_ >= 1)) stop("must be nu >= 1")
#   if(!(kap_ > 0)) stop("must be kap > 0")
#   
#   logdens = function(x, a_, b_, r_, s_, nu_, kap_){
#     -a_*x/2 + (b_/2-1)*log1p(-nu_*exp(-x)) - s_*exp(-x)/2 - r_*log(kap_ + (1-kap_)*nu_*exp(-x))
#   }
#   dlogdens = function(x, a_, b_, r_, s_, nu_, kap_){
#     -a_/2 + (b_/2-1) * nu_ / (exp(x)-nu_) + s_*exp(-x)/2 + r_*(1-kap_)*nu/(kap_*exp(x) + (1-kap_)*nu_)
#   }
#   xlb = 0
#   if(log(nu_)>0) xlb = log(nu_)+0.1
#   return(ars::ars(nsamp,logdens,dlogdens,
#                   x=start,m=1, lb=TRUE,xlb=xlb,
#                   a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_))
# }
# 
# ## Gauss Example
# n = 1000; p = 200; r2 = 0.5
# # a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1 # zsadapted
# # a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1 # uniform
# a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1 # hyperg
# a_ = (a+p); b_ = b; r_ = (n-1)/2; s_ = s; nu_ = nu; kap_ = (1-r2)/(1-(1-1/nu)*r2)
# x = seq(log(nu_)+0.1, x_from_g(5*n), len = 200)
# par(mfrow=c(2,3))
# plot(x, 
#      exp(logdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_)),
#      type= "l")
# plot(x, 
#      logdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l")
# plot(x, 
#      dlogdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l"); abline(h=0, lty=2)
# plot(x, 
#      ddlogdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l"); abline(h=0, lty=2)
# hist(rtCCH(1000, x_from_g(n), a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_))
# 
# ## Exp example
# n = 1000; p = 30; Qm = n
# a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic, b=1 is problematic
# # a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
# # a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1 # betaprime
# a_ = (a+p); b_ = b; r_ = r; s_ = s+Qm; nu_ = nu; kap_ = kap
# x = seq(log(nu_)+0.01, x_from_g(n/2), len = 200)
# par(mfrow=c(2,3))
# plot(x, 
#      exp(logdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_)),
#      type= "l")
# plot(x, 
#      logdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l")
# plot(x, 
#      dlogdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l"); abline(h=0, lty=2)
# plot(x, 
#      ddlogdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l"); abline(h=0, lty=2)
# hist(rtCCH(1000, x_from_g(n), a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_))

################################################################################
## tCCH distribution
## u = 1/(1+g) ~ tCCH(a/2, b/2, r, s/2, nu, kap)
## u = 1/(1+g) \in (0, 1/nu)

## we use x = log(nu u / (1-nu u)) = log(nu / (g+1-nu))
x_from_g = function(g, nu) { log(nu) - log1p(g-nu)}
## g = nu exp(-x) - 1 + nu
g_from_x = function(x, nu){ nu*exp(-x)+nu-1 }

logdens_tCCH = function(x, a_, b_, r_, s_, nu_, kap_){
  sigmoidx = 1/(1+exp(-x))
  tmp = a_*x/2 - (a_+b_)/2 *log1p(exp(x)) - 0.5*s_/nu_*sigmoidx - r_*log(kap_+(1-kap_)*sigmoidx)
  tmp - median(tmp)
}
dlogdens_tCCH = function(x, a_, b_, r_, s_, nu_, kap_){
  sigmoidx = 1/(1+exp(-x))
  dsigmoidx = sigmoidx * (1-sigmoidx)
  
  a_/2 - (a_+b_)/2 *sigmoidx - 0.5*s_/nu_*dsigmoidx - 
    r_*(1-kap_)*dsigmoidx / (kap_+(1-kap_)*sigmoidx)
}
ddlogdens_tCCH = function(x, a_, b_, r_, s_, nu_, kap_){
  sigmoidx = 1/(1+exp(-x))
  dsigmoidx = sigmoidx * (1-sigmoidx)
  ddsigmoidx = sigmoidx * (1-sigmoidx) * (1-2*sigmoidx)
  
  -(a_+b_)/2 *dsigmoidx - 0.5*s_/nu_*ddsigmoidx - 
    r_*(1-kap_)* (ddsigmoidx*(kap_ + (1-kap_)*sigmoidx) - (1-kap_)*dsigmoidx^2)/(kap_+(1-kap_)*sigmoidx)^2
}

##
rtCCH = function(nsamp, 
                 a_, b_, r_, s_, nu_, kap_){
  if(!(a_>0)) stop("must be a > 0")
  if(!(b_>0)) stop("must be b > 0")
  if(!(nu_ >= 1)) stop("must be nu >= 1")
  if(!(kap_ > 0)) stop("must be kap > 0")
  
  logdens = function(x, a_, b_, r_, s_, nu_, kap_){
    sigmoidx = 1/(1+exp(-x))
    tmp = a_*x/2 - (a_+b_)/2 *log1p(exp(x)) - 0.5*s_/nu_*sigmoidx - r_*log(kap_+(1-kap_)*sigmoidx)
    tmp - median(tmp)
  }
  dlogdens = function(x, a_, b_, r_, s_, nu_, kap_){
    sigmoidx = 1/(1+exp(-x))
    dsigmoidx = sigmoidx * (1-sigmoidx)
    
    a_/2 - (a_+b_)/2 *sigmoidx - 0.5*s_/nu_*dsigmoidx - 
      r_*(1-kap_)*dsigmoidx / (kap_+(1-kap_)*sigmoidx)
  }
  # xlb = 0
  # if(log(nu_)>0) xlb = log(nu_)+0.1
  return(ars::ars(nsamp,logdens,dlogdens,m=1, 
                  a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_))
}

## Exp example
n = 1000; p = 30; Qm = n*3
# a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic, b=1 is problematic
# a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1 # betaprime
a_ = (a+p); b_ = b; r_ = r; s_ = s+Qm; nu_ = nu; kap_ = kap
x = seq(-20, 20, len = 200)
par(mfrow=c(2,2))
# plot(x, 
#      exp(logdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_)),
#      type= "l")
# plot(g_from_x(x, nu_), 
#      exp(logdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_)) / (1+g_from_x(x, nu_)-nu_), xlim = c(nu_+0.1, n),
#      type= "l")
plot(x, 
     logdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
     type= "l")
plot(x, 
     dlogdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
     type= "l"); abline(h=0, lty=2)
plot(x, 
     ddlogdens_tCCH(x, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
     type= "l"); abline(h=0, lty=2)
hist(rtCCH(1000, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_))



################################################################################
## Robust - Gauss

logdens = function(w, a_, b_, bp_, r_, x_, y_){
  (a_-1)*log(w) + (r_-a_-1)*log1p(-w) - b_*log1p(-x_*w) - bp_*log1p(-y_*w)
}
dlogdens = function(w, a_, b_, bp_, r_, x_, y_){
  (a_-1)/w + (r_-a_-1)/w + x_*b_/(1-x_*w) + y_*bp_/(1-y_*w)
}
ddlogdens = function(w, a_, b_, bp_, r_, x_, y_){
  -(a_-1)/w^2 - (r_-a_-1)/w^2 + x_^2*b_/(1-x_*w)^2 + y_^2*bp_/(1-y_*w)^2
}




################################################################################
## Intrinsic-Gauss combination (Womack, 2014)
## w = (n/(p+1))/g, g = (n/(p+1))/w
## w = sin(pi*theta/2)^2 = (n/(p+1))/g
## g = (n/(p+1)) * 1/sin(pi*theta/2)^2
par(mfrow=c(2,4))
for(r2 in c(.3,.5, .7,.8)){
  n = 1000; p = 20; r2 = r2
  theta = seq(0.01, 0.99, len=1000)
  g_from_theta = function(theta, n_, p_) (n_/(p_+1)) / sin(pi*theta)^2
  logdens = function(theta, n_, p_, r2_){
    theta = sin(pi*theta/2)^2
    (n_-p_)/2 * (log1p(theta*(p_+1)/n_) - log1p(theta*(p_+1)/(n_*(1-r2_)))) +
      (p_-1)/2 * ( log(theta*(p_+1)/n) - log1p((theta*(p_+1)/(n_*(1-r2)))))
  }
  dlogdens = function(theta, n_, p_, r2_){
    theta = sin(pi*theta/2)^2
    tmp = (n_-p_)/2 * ( (p_+1)/(n_+theta*(p_+1)) - (p_+1)/(n_*(1-r2_)+theta*(p_+1)) ) +
      (p_-1)/2 * ( 1/theta - (p_+1)/(n_*(1-r2_)+theta*(p_+1)) )
    tmp * pi * sin(pi*theta)
  }
  plot(theta, (logdens(theta,n,p,r2)), type="l", main = r2)
  plot(theta, (dlogdens(theta,n,p,r2)), type="l"); abline(h=0)
}


################################################################################
## Hypergn - Gauss combination (Liang, 2008)
## tau = log(g)
## g = exp(tau)
par(mfrow=c(4,4))
for(r2 in c(.1,.3, .7,.8)){
  n = 1000; p = 20; r2 = r2
  tau = seq(-log(n),log(n)+3, len=1000)
  logdens = function(tau, n_, p_, r2_){
    tmp = tau + (n_-p_-1)/2 * log(1+exp(tau)) - (n_-1)/2*log(1+exp(tau)*(1-r2_)) - a/2*log1p(exp(tau)/n)
    tmp - max(tmp)
  }
  dlogdens = function(tau, n_, p_, r2_){
    1 + (n_-p_-1)/2 * exp(tau)/(1+exp(tau)) - (n_-1)/2* (1-r2_)*exp(tau)/(1+(1-r2_)*exp(tau)) - a/2*exp(tau)/(n+exp(tau))
  }
  plot(exp(tau), exp(logdens(tau,n,p,r2))/exp(tau), type="l", main = r2, xlim=c(0, n))
  plot(tau, exp(logdens(tau,n,p,r2)), type="l", main = "in log(g)")
  plot(tau, (logdens(tau,n,p,r2)), type="l", main = "log dens")
  plot(tau, (dlogdens(tau,n,p,r2)), type="l", main = "first deriv"); abline(h=0)
  # hist(rfuck(1000, log(n), n_=n, p_=p, r2_=r2))
}

################################################################################
## Truncated Gamma
n = 1000; p = 20; Qm = 3000
w = seq(0.0000001, 0.999, len=10000)
par(mfrow=c(2,2))
plot(w, dgamma(w, shape = (p+2)/2, rate = Qm/2), type="l") # Uniform
plot(w, dgamma(w, shape = (p+1)/2, rate = Qm/2), type="l") # Hyper-g
plot(w, dgamma(w, shape = (p+1)/2, rate = (Qm+n+3)/2), type="l") # ZS adapted
plot(w, dgamma(w, shape = (p+1)/2, rate = Qm/2), type="l", xlim=c(0, (p+1)/(n+1))) # Robust

g = seq(0, n, len = 10000)
plot(g, dgamma(1/(g+1), shape = (p+2)/2, rate = Qm/2)/(g+1)^2, type="l") # Uniform
plot(g, dgamma(1/(g+1), shape = (p+1)/2, rate = Qm/2)/(g+1)^2, type="l") # Hyper-g
plot(g, dgamma(1/(g+1), shape = (p+1)/2, rate = (Qm+n+3)/2)/(g+1)^2, type="l") # ZS adapted
plot(g, dgamma(1/(g+1), shape = (p+1)/2, rate = Qm/2)/(g+1)^2, type="l", xlim=c((n-p)/(p+1), n)) # Robust


# ################################################################################
# ## tCCH distribution
# ## u = 1/(1+g) ~ tCCH(a/2, b/2, r, s/2, nu, kap)
# ## u = 1/(1+g) \in (0, 1/nu)
# 
# ## we use w = 1/u-1 = g
# w_from_g = function(g) { g }
# g_from_w = function(w){ w }
# 
# logdens_tCCH = function(w, a_, b_, r_, s_, nu_, kap_){
#   (r_-(a_+b_)/2)*log(1+w) + (b_/2-1)*log(w-nu_+1) - s_/(2*(w+1)) - r_*log(kap_*(w-nu_+1)+nu_)
# }
# dlogdens_tCCH = function(w, a_, b_, r_, s_, nu_, kap_){
#   (r_-(a_+b_)/2)/(1+w) + (b_/2-1)/(w-nu_+1) + 0.5*s_/(w+1)^2 - kap_*r_/(kap_*(w-nu_+1)+nu_)
# }
# ddlogdens_tCCH = function(w, a_, b_, r_, s_, nu_, kap_){
#   -(r_-(a_+b_)/2)/(1+w)^2 - (b_/2-1)/(w-nu_+1)^2 + s_/(w+1)^3 - kap_^2*r_/(kap_*(w-nu_+1)+nu_)^2
# }
# 
# ##
# rtCCH = function(nsamp, start, 
#                  a_, b_, r_, s_, nu_, kap_){
#   if(!(a_>0)) stop("must be a > 0")
#   if(!(b_>0)) stop("must be b > 0")
#   if(!(nu_ >= 1)) stop("must be nu >= 1")
#   if(!(kap_ > 0)) stop("must be kap > 0")
#   
#   logdens = function(w, a_, b_, r_, s_, nu_, kap_){
#     (r_-(a_+b_)/2)*log(1+w) + (b_/2-1)*log(w-nu_+1) - s_/(2*(w+1)) - r_*log(kap_*(w-nu_+1)+nu_)
#   }
#   dlogdens = function(w, a_, b_, r_, s_, nu_, kap_){
#     (r_-(a_+b_)/2)/(1+w) + (b_/2-1)/(w-nu_+1) + 0.5*s_/(w+1)^2 - kap_*r_/(kap_*(w-nu_+1)+nu_)
#   }
#   xlb = 0
#   if(nu_>0) xlb = nu_-1 + 0.01
#   return(ars::ars(nsamp,logdens,dlogdens,
#                   x=start,m=1, lb=TRUE,xlb=xlb,
#                   a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_))
# }
# 
# ## Gauss Example
# n = 1000; p = 20; r2 = 0.8
# # a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1 # zsadapted
# # a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1 # uniform
# a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1 # hyperg
# a_ = (a+p); b_ = b; r_ = (n-1)/2; s_ = s; nu_ = nu; kap_ = (1-r2)/(1-(1-1/nu)*r2)
# w = seq(nu_+0.01, w_from_g(n), len = 200)
# par(mfrow=c(2,3))
# plot(w, 
#      exp(logdens_tCCH(w, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_)),
#      type= "l")
# plot(w, 
#      logdens_tCCH(w, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l")
# plot(w, 
#      dlogdens_tCCH(w, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l"); abline(h=0, lty=2)
# plot(w, 
#      ddlogdens_tCCH(w, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l"); abline(h=0, lty=2)
# hist(rtCCH(1000, w_from_g(n), a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_))
# 
# ## Exp example
# n = 1000; p = 20; Qm = n*3
# # a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic, b=1 is problematic
# a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1 # robust
# # a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
# # a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1 # zsadapted
# # a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1 # betaprime
# # a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1 # uniform
# # a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1 # hyperg
# a_ = (a+p); b_ = b; r_ = r; s_ = s+Qm; nu_ = nu; kap_ = kap
# w = seq(nu_+0.01, w_from_g(n), len = 200)
# par(mfrow=c(2,3))
# plot(w, 
#      exp(logdens_tCCH(w, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_)),
#      type= "l")
# plot(w, 
#      logdens_tCCH(w, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l")
# plot(w, 
#      dlogdens_tCCH(w, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l"); abline(h=0, lty=2)
# plot(w, 
#      ddlogdens_tCCH(w, a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_),
#      type= "l"); abline(h=0, lty=2)
# hist(rtCCH(1000, w_from_g(n), a_= a_, b_=b_, r_=r_, s_=s_, nu_=nu_, kap_=kap_))


# ################################################################################
# ## Appell distribution
# 
# logF1 = function(aa, bb, bbp, rr, xx, yy, logValue = T, ...){
#   # ... : arguments to R's integrate
#   if(!(rr > aa & aa >0)) stop("must be r > a > 0")
#   
#   tmpf = function(uu, aa, bb, bbp, rr, xx, yy){
#     (aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-xx*uu) - bbp*log1p(-yy*uu)- base::lbeta(rr-aa, aa)
#   }
#   
#   uu = seq(0.001, 0.999, len = 200)
#   adjV = max(tmpf(uu, aa, bb, bbp, rr, xx, yy))- log(.Machine$double.xmax)/2
#   tmpf = function(uu, aa, bb, bbp, rr, xx, yy, adjV){
#     exp((aa-1)*log(uu) + (rr-aa-1)*log1p(-uu) - bb*log1p(-xx*uu) - bbp*log1p(-yy*uu)- base::lbeta(rr-aa, aa) - adjV)
#   }
# 
#   tmp = log(stats::integrate(tmpf, lower=0, upper=1,
#                              aa=aa, bb=bb, bbp=bbp, rr=rr, xx=xx, yy=yy, adjV = adjV,
#                              ...)$value) + adjV
# 
#   if(!logValue) tmp = exp(tmp)
#   return(tmp)
# }
# logdens_Appell = function(w, a_, b_, bp_, r_, x_, y_){
#   (a_-1)*log(w) + (r_-a_-1)*log1p(-w) - b_*log1p(-x_*w) - bp_*log1p(-y_*w)
# }
# dlogdens_Appell = function(w, a_, b_, bp_, r_, x_, y_){
#   (a_-1)/w + (r_-a_-1)/w + x_*b_/(1-x_*w) + y_*bp_/(1-y_*w)
# }
# ddlogdens_Appell = function(w, a_, b_, bp_, r_, x_, y_){
#   -(a_-1)/w^2 - (r_-a_-1)/w^2 + x_^2*b_/(1-x_*w)^2 + y_^2*bp_/(1-y_*w)^2
# }
# 
# ## Gauss Example
# par(mfrow=c(4,3))
# for(r22 in c(0.5, 0.6, 0.7, 0.8)){
#   n = 1000; p = 100; r2 = r22
#   # a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n  # intrinsic
#   a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1 # robust
#   # a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
#   a_ = (a+p)/2; b_ = (a+b+p+1-n-2*r)/2; bp_ = (n-1)/2; r_ = (a+b+p)/2; x_ = 1-kap; y_ = ((1-r2)*nu*(1-kap) - r2*kap)/((1-r2)*nu)
#   x = seq(0.01, 0.99, len = 1000)
#   plot(x,
#        exp(logdens_Appell(x, a_, b_, bp_, r_, x_, y_)),
#        type= "l", main = paste0("Gauss-Intrinsic (n ", n, ", p ",p, ", R2 ",r2, ")"),
#        ylab= "Density", xlab = "1/(g+1)")
#   sum(exp(logdens_Appell(x, a_, b_, bp_, r_, x_, y_) - base::lbeta(r_-a_, a_) - logF1(a_,b_,bp_,r_,x_,y_))) * (x[2]-x[1]) # roughly 1
#   Ex = exp(logF1(a_+1,b_, bp_,r_+1,x_,y_) - logF1(a_,b_, bp_,r_,x_,y_))
#   Ex2 = exp(logF1(a_+2,b_, bp_,r_+2,x_,y_) - logF1(a_,b_, bp_,r_,x_,y_))
#   Ex3 = exp(logF1(a_+3,b_, bp_,r_+3,x_,y_) - logF1(a_,b_, bp_,r_,x_,y_))
#   (Mean = Ex)
#   (Sd = sqrt(Ex2 - Ex^2))
#   plot(x,
#        logdens_Appell(x, a_, b_, bp_, r_, x_, y_),
#        type= "l", main = paste0("Gauss-Intrinsic (n ", n, ", p ",p, ", R2 ",r2, ")"),
#        ylab= "Log Density", xlab = "1/(g+1)")
#   # plot(x,
#   #      dlogdens_Appell(x, a_, b_, bp_, r_, x_, y_),
#   #      type= "l"); abline(h=0, lty=2)
#   # plot(x,
#   #      ddlogdens_Appell(x, a_, b_, bp_, r_, x_, y_),
#   #      type= "l"); abline(h=0, lty=2)
#   
#   rAppell = function(nsamp, start,
#                      a_, b_, bp_, r_, x_, y_){
#     if(!(r_ > a_ & a_ >0)) stop("must be r > a > 0")
#     
#     logdens = function(w, a_, b_, bp_, r_, x_, y_){
#       (a_-1)*log(w) + (r_-a_-1)*log(w) - b_*log1p(x_*w) - bp_*log(1-y_*w)
#     }
#     dlogdens = function(w, a_, b_, bp_, r_, x_, y_){
#       (a_-1)/w + (r_-a_-1)/w + x_*b_/(1-x_*w) + y_*bp_/(1-y_*w)
#     }
#     
#     return(ars::ars(nsamp,logdens,dlogdens,
#                     x=start,m=1, lb=TRUE,xlb=0, ub=TRUE, xub = 1,
#                     a_=a_, b_=b_, bp_=bp_, r_=r_, x_=x_, y_=y_))
#   }
#   hist(rAppell(1000, 0.3, a_=a_, b_=b_, bp_=bp_, r_=r_, x_=x_, y_=y_))
# }

# ## we use x = -log(u) = log(g+1)
# x_from_g = function(g) { log(g+1) }
# ## g = exp(y) - 1
# g_from_x = function(y){ exp(y) - 1 }
# 
# logdens_Appell = function(w, a_, b_, bp_, r_, x_, y_){
#   -a_*w + (r_-a_-1)*log(1-exp(-w)) - b_*log(1-x_*exp(-w)) - bp_*log(1-y_*exp(-w))
# }
# dlogdens_Appell = function(w, a_, b_, bp_, r_, x_, y_){
#   -a_ - (r_-a_-1)/(exp(w)-1) + b_*x_/(exp(w)-x_) + bp_*y_/(exp(w)-y_)
# }
# ddlogdens_Appell = function(w, a_, b_, bp_, r_, x_, y_){
#   - (r_-a_-1)*exp(w)/(exp(w)-1)^2 + b_*x_*exp(w)/(exp(w)-x_)^2 + bp_*y_*exp(w)/(exp(w)-y_)^2
# }
# 
# ## Gauss Example
# n = 1000; p = 20; r2 = 0.5
# # a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n  # intrinsic
# a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1 # robust
# # a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
# a_ = (a+p)/2; b_ = (a+b+p+1-n-2*r)/2; bp_ = (n-1)/2; r_ = (a+b+p)/2; x_ = 1-kap; y_ = ((1-r2)*nu*(1-kap) - r2*kap)/((1-r2)*nu)
# x = seq(0+0.001, x_from_g(n), len = 200)
# par(mfrow=c(2,3))
# plot(x,
#      exp(logdens_Appell(x, a_, b_, bp_, r_, x_, y_)),
#      type= "l")
# plot(x,
#      logdens_Appell(x, a_, b_, bp_, r_, x_, y_),
#      type= "l")
# plot(x,
#      dlogdens_Appell(x, a_, b_, bp_, r_, x_, y_),
#      type= "l"); abline(h=0, lty=2)
# plot(x,
#      ddlogdens_Appell(x, a_, b_, bp_, r_, x_, y_),
#      type= "l"); abline(h=0, lty=2)
# 
# 
# rAppell = function(nsamp, start,
#                    a_, b_, bp_, r_, x_, y_){
#   if(!(r_ > a_ & a_ >0)) stop("must be r > a > 0")
# 
#   logdens = function(w, a_, b_, bp_, r_, x_, y_){
#     -a_*w + (r_-a_-1)*log(1-exp(-w)) - b_*log(1-x_*exp(-w)) - bp_*log(1-y_*exp(-w))
#   }
#   dlogdens = function(w, a_, b_, bp_, r_, x_, y_){
#     -a_ - (r_-a_-1)/(exp(w)-1) + b_*x_/(exp(w)-x_) + bp_*y_/(exp(w)-y_)
#   }
# 
#   return(ars::ars(nsamp,logdens,dlogdens,
#                   x=start,m=1, lb=TRUE,xlb=0,
#                   a_=a_, b_=b_, bp_=bp_, r_=r_, x_=x_, y_=y_))
# }
# hist(rAppell(1000, x_from_g(n), a_=a_, b_=b_, bp_=bp_, r_=r_, x_=x_, y_=y_))
# 
# 
# ## we use w = log(u/(1-u)) = -log(g)
# w_from_g = function(g) { -log(g) }
# ## g = exp(-w)
# g_from_w = function(w){ exp(-w) }
# 
# logdens_Appell = function(w, a_, b_, bp_, r_, x_, y_){
#   a_*w + (b_+bp_-r_)*log1p(exp(w)) - b_*log1p(exp(w)*(1-x_)) - bp_*log1p(exp(w)*(1-y_))
# }
# dlogdens_Appell = function(w, a_, b_, bp_, r_, x_, y_){
#   a_ + (b_+bp_-r_)/(1+exp(-w)) - b_*(1-x_)/(exp(-w)+1-x_) - bp_*(1-y_)/(exp(-w)+1-y_)
# }
# ddlogdens_Appell = function(w, a_, b_, bp_, r_, x_, y_){
#   (b_+bp_-r_)*exp(w)/(1+exp(w))^2 - b_*exp(w)*(1-x_)/(1+exp(w)*(1-x_))^2 - bp_*exp(w)*(1-y_)/(1+exp(w)*(1-y_))^2
# }

## Gauss Example
# n = 1000; p = 30; r2 = 0.7
# # a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n  # intrinsic
# a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1 # robust
# # a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
# a_ = (a+p)/2; b_ = (a+b+p+1-n-2*r)/2; bp_ = (n-1)/2; r_ = (a+b+p)/2; x_ = 1-kap; y_ = ((1-r2)*nu*(1-kap) - r2*kap)/((1-r2)*nu)
# w = seq(w_from_g(n/2), -w_from_g(n), len = 200)
# par(mfrow=c(2,3))
# plot(w,
#      exp(logdens_Appell(w, a_, b_, bp_, r_, x_, y_)),
#      type= "l")
# plot(w,
#      logdens_Appell(w, a_, b_, bp_, r_, x_, y_),
#      type= "l")
# plot(w,
#      dlogdens_Appell(w, a_, b_, bp_, r_, x_, y_),
#      type= "l"); abline(h=0, lty=2)
# plot(w,
#      ddlogdens_Appell(w, a_, b_, bp_, r_, x_, y_),
#      type= "l"); abline(h=0, lty=2)
# 
# 
# rAppell = function(nsamp, start,
#                    a_, b_, bp_, r_, x_, y_){
#   if(!(r_ > a_ & a_ >0)) stop("must be r > a > 0")
# 
#   logdens = function(w, a_, b_, bp_, r_, x_, y_){
#     (a_-1)*log(w) + (r_-a_-1)*log(w) - b_*log1p(x_*w) - bp_*log(1-y_*w)
#   }
#   dlogdens = function(w, a_, b_, bp_, r_, x_, y_){
#     (a_-1)/w + (r_-a_-1)/w + x_*b_/(1-x_*w) + y_*bp_/(1-y_*w)
#   }
# 
#   return(ars::ars(nsamp,logdens,dlogdens,
#                   x=start,m=1,
#                   a_=a_, b_=b_, bp_=bp_, r_=r_, x_=x_, y_=y_))
# }
# hist(rAppell(1000, 0.7, a_=a_, b_=b_, bp_=bp_, r_=r_, x_=x_, y_=y_))














