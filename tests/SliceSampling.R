# source("tests/specialfunctions.R")
library(Rcpp)
library(RcppArmadillo)
sourceCpp("src/SliceSampling.cpp")
sourceCpp("src/specialFunctions.cpp")
# a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1 # uniform
# a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1 # hyperg
# a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
# a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1 # betaprime
# a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1 # zsadapted
# a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1 # robust
# a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic

##################################################################################################
## CONFLUENT HYPERGEOMETRIC DIST

dCH = function(uu, aa, bb, ss){
  exp((aa-1)*log(uu) + (bb-1)*log1p(-uu) - ss*uu - base::lbeta(aa, bb) - log1F1_cpp(aa, aa+bb, -ss))
}
momentf = function(r=1, a, b, x){
  exp(lbeta(r+a, b) - lbeta(a,b) + log1F1(a+r, a+b+r, -x) - log1F1(a, a+b, -x))
}
## Exp example
n = 200; p = 30; Qm = n*2
a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1 # betaprime
a_pos = (a+p)/2; b_pos = b/2; s_pos = (s+Qm)/(2*nu)
# a_pos = 2; b_pos = 2; s_pos = 10;
uu = seq(0.01, 0.99, len=200)
plot(uu, dCH(uu, a_pos, b_pos, s_pos), type="l"); abline(h=0)

# rCH = function(nsamp, burnin, a_pos, b_pos, s_pos){
#   rCH_atom = function(v_, l_, a_pos, b_pos, s_pos){
#     l_ = exp(log(runif(1))-s_pos*v_)
#     if (s_pos > 0) {
#       lb = 0
#       ub = min(-log(l_)/s_pos, 1)
#     } else {
#       lb = max(0, -log(l_)/s_pos)
#       ub = 1
#     }
#     v_ = TruncatedDistributions::rtbeta(1, a_pos, b_pos, lb,ub)
#     return(c(v_,l_))
#   }
#   
#   S = nsamp + burnin
#   V = L = rep(NA, S)
#   v_ = l_ = 0.1
#   for(i in 1:S){
#     tmp = rCH_atom(v_, l_, a_pos, b_pos, s_pos)
#     v_ = tmp[1]
#     l_ = tmp[2]
#     V[i]=v_
#   }
#   return(V[(burnin+1):S])
# }
# 
# V = rCH(1e5, 0, a_pos, b_pos, s_pos)
V = rCH_cpp(1e5, 0, a_pos, b_pos, s_pos)
# plot(V, type="l", ylim=c(0, 0.2)); abline(h = 0.1)
par(mfrow=c(1,2))
hist(V, breaks = 100, prob=T)
uu = seq(0.0001, 0.99, len=1e4)
lines(uu, dCH(uu, a_pos,b_pos, s_pos))
acf(V)

# comp = microbenchmark(
#   rCH(1e1, 0, a_pos, b_pos, s_pos),
#   rCH_cpp_vec(1e1, 0, a_pos, b_pos, s_pos),
#   times = 2000
# )
# comp
# autoplot(comp)

##################################################################################################
## GAUSSIAN HYPERGEOMETRIC DIST

dGH = function(u, a, b, x, z){
  exp((a-1)*log(u) + (b-1)*log1p(-u) - z*log1p(x*u) - lbeta(a,b) - log2F1_cpp(z,a, a+b, -x))
}
momentf = function(r=1, a, b, x, z){
  exp(lbeta(r+a, b) - lbeta(a,b) + log2F1(z, a+r, a+b+r, -x) - log2F1(z, a, a+b, -x))
}
## Gauss Example
n =200; p = 30; r2 = 0.7
# a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1 # uniform
# a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1 # hyperg
a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1 # robust
# a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
a_pos = (a+p)/2; b_pos = b/2; x_pos = (r2/nu)/(1-r2); z_pos = (n-1)/2
uu = seq(0.01, 0.99, len=200)
plot(uu, dGH(uu, a_pos,b_pos, x_pos, z_pos), type="l");abline(h=0)

S = 1e5
U = rep(NA, S)
l_ = t_ = u = 0.5
# for(i in 1:S){
#   l_ = exp(log(runif(1))-x_pos*u*t_)
#   if(x_pos>0){
#     ubt = -log(l_)/(x_pos*u)
#     lbt = 0
#   } else {
#     ubt = Inf
#     lbt = max(0, -log(l_)/(x_pos*u))
#   }
#   t_ = TruncatedDistributions::rtgamma(1, shape = z_pos, scale=1, ubt, lbt)
#   # t_ = rgamma(1, shape = z_pos, rate = 1+x_pos*u)
#   if (x_pos > 0) {
#     ub = min(-log(l_)/(x_pos*t_), 1)
#     lb = 0
#   } else {
#     lb = max(-log(l_)/(x_pos*t_), 0)
#     ub = 1
#   }
#   u = TruncatedDistributions::rtbeta(1, a_pos, b_pos, lb,ub)
#   U[i]=u
# }
U = rGH_cpp(S, 0, a_pos, b_pos, x_pos, z_pos)
hist(U, breaks = 100, prob=T)
uu = seq(0.0001, 0.99, len=1e4)
lines(uu, dGH(uu, a_pos,b_pos, x_pos, z_pos))
acf(U)


hist(U/nu, breaks = 100, prob=T,
     main = paste0("Gauss-Hypergn, n ", n, ", p ", p, ", R2 ", r2, " MC size ", S))

##################################################################################################
dCCH = function(uu, aa, bb, zz, ss, nunu, xx){
  exp((aa-1)*log(uu) + (bb-1)*log1p(-uu) - zz*log1p(xx*uu) - ss*uu/nunu + ss/nunu + zz*log1p(xx) -
        lbeta(aa,bb) - logPhi1_cpp(bb,zz,aa+bb, ss/nunu, xx/(xx+1)))
}
# momentf = function(r=1,aa, bb, zz, ss, nunu, ){
#   exp(lbeta(r+aa, bb) - lbeta(aa,bb) + logPhi1(bb,zz,aa+bb+r, ss, 1-theta) - logPhi1(bb,zz,aa+bb, ss, 1-theta))
# }

# Exp example
n = 500; p = 30; Qm = n
# a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic
a_pos = (a+p)/2; b_pos = b/2; z_pos = r; s_pos = (s+Qm)/(2); nu_pos = nu; x_pos = (1-kap)/kap

## Gauss Example
n =200; p = 30; r2 = 0.7
a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1 # zsadapted
a_pos = (a+p)/2; b_pos = b/2; z_pos = (n-1)/2; s_pos = (s)/(2); nu_pos = nu; x_pos = (r2/nu)/(1-r2)

uu = seq(0.01, 0.99, len=1000)
plot(uu, dCCH(uu, a_pos, b_pos,z_pos, s_pos, nu_pos, x_pos), type="l"); abline(h=0)
# momentf(r=1, a_pos, b_pos, z_pos, s_pos, theta_pos)

S = 1e5
L_=M_=T_ = U = rep(NA, S)
l_ = m_ = t_ = u = 0.5
for(i in 1:S){
  l_ = exp(log(runif(1))-x_pos*u*t_)
  m_ = exp(log(runif(1))-(s_pos/nu_pos)*u)
  
  if(x_pos>0){
    ubt = -log(l_)/(x_pos*u)
    lbt = 0
  } else {
    ubt = Inf
    lbt = max(0, -log(l_)/(x_pos*u))
  }
  t_ = TruncatedDistributions::rtgamma(1, shape = z_pos, scale=1, ubt, lbt)

  if (x_pos > 0 & (s_pos/nu_pos) > 0) {
    ub = min(-log(l_)/(x_pos*t_), -log(m_)/(s_pos/nu_pos))
    lb = 0
  } else if (x_pos > 0 & (s_pos/nu_pos) <= 0){
    ub = -log(l_)/(x_pos*t_)
    lb = -log(m_)/(s_pos/nu_pos)
  } else if (x_pos <= 0 & (s_pos/nu_pos) > 0){
    ub = -log(m_)/(s_pos/nu_pos)
    lb = -log(l_)/(x_pos*t_)
  } else {
    ub = 1
    lb = max(-log(l_)/(x_pos*t_), -log(m_)/(s_pos/nu_pos))
  }
  
  u = TruncatedDistributions::rtbeta(1, a_pos, b_pos, max(0, lb), min(1, ub))
  U[i]=u
  L_[i] = l_; M_[i] = m_; T_[i]=t_;
}
hist(U, breaks = 200, prob=T)

hist(U/nu, breaks=100, prob=T)
par(mfrow=c(1,2))
Unu = rtCCH_cpp(1e5, 0, a_pos, b_pos, z_pos, s_pos, nu_pos, 1/(x_pos+1))
hist(Unu/nu, breaks=100, prob=T)
uu = seq(0.0001, 0.9999, len=1e4)
lines(uu, dCCH(uu, a_pos, b_pos,z_pos, s_pos, nu_pos, x_pos))
acf(Unu)

##################################################################################################
dAPL = function(u, a, b, z, w, x, y){
  exp((a-1)*log(u) + (b-1)*log1p(-u) - z*log1p(x*u) - w*log1p(y*u) - lbeta(a,b) - logF1(a,z,w,a+b,-x,-y))
}
momentf = function(r=1, a, b, z,w,x,y){
  exp(lbeta(r+a, b) - lbeta(a,b) + logF1(a+r,z,w,a+b+r,-x,-y) - logF1(a,z,w,a+b,-x,-y))
}

## Gauss Example
n =200; p = 30; r2 = 0.7
# a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic
a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
a_pos = (a+p)/2; b_pos = b/2; z_pos = (n-1)/2; w_pos = r; x_pos = (r2/nu)/(1-r2); y_pos = (1-kap)/kap; 
uu = seq(0.01, 0.99, len=200)
plot(uu, dAPL(uu, a_pos, b_pos, z_pos, w_pos,x_pos, y_pos), type="l");abline(h=0)

S = 1e5
U = rep(NA, S)
t_=q_=l_=m_=u=0.1
for(i in 1:S){
  l_ = exp(log(runif(1))-x_pos*u*t_)
  m_ = exp(log(runif(1))-y_pos*u*q_)

  if(x_pos>0){
    ubt1 = -log(l_)/(x_pos*u)
    lbt1 = 0
  } else {
    ubt1 = Inf
    lbt1 = max(0, -log(l_)/(x_pos*u))
  }
  t_ = TruncatedDistributions::rtgamma(1, shape = z_pos, scale=1, ubt1, lbt1)

  if(y_pos>0){
    ubt2 = -log(m_)/(y_pos*u)
    lbt2 = 0
  } else {
    ubt2 = Inf
    lbt2 = max(0, -log(m_)/(y_pos*u))
  }
  q_ = TruncatedDistributions::rtgamma(1, shape = w_pos, scale=1, ubt2, lbt2)

  if (x_pos >0 & y_pos > 0) {
    lb = 0
    ub = min(-log(l_)/(x_pos*t_), -log(m_)/(y_pos*q_))
  } else if (x_pos <=0 & y_pos <=0) {
    lb = max(-log(l_)/(x_pos*t_), -log(m_)/(y_pos*q_))
    ub = 1
  } else if (x_pos > 0 & y_pos <=0) {
    lb = max(0,  -log(m_)/(y_pos*q_))
    ub = min(1, -log(l_)/(x_pos*t_))
  } else  {
    lb = max(0, -log(l_)/(x_pos*t_))
    ub = min(1,  -log(m_)/(y_pos*q_))
  }

  if(ub<lb){
    print(c(ub, lb))
    break
  }
  u = TruncatedDistributions::rtbeta(1, a_pos, b_pos, max(0, lb), min(ub,1))
  # cat(paste0("w=", w, " u=", u, " x=", x, "\n"))
  U[i]=u
}
hist(U, breaks = 100, prob=T)
lines(uu, dAPL(uu, a_pos, b_pos, z_pos, w_pos,x_pos, y_pos))
hist(U/nu, breaks=100, prob=T)

par(mfrow=c(1,2))
Unu = rAPL_cpp(1e5, 0, a_pos, b_pos, z_pos, w_pos, x_pos, y_pos)
hist(Unu, breaks=100, prob=T)
acf(Unu)














