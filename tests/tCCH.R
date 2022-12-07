## tCCH prior

tCCH = function(u, a, b, r, s, nu, k, log = F){
  if(any(u<0 || u >= 1)) stop("must be 0<u<1")
  res = rep(NA_real_, length(c(u)))
  res[u > 1/nu] = 0.0
  
  res_pos = res[u <= 1/nu]
  res_pos = 0.5*a*log(nu) + 0.5*s/nu + (0.5*a-1.0)*log(u) + 
    (0.5*b-1.0)*log1p(-nu*u) - 0.5*s*u - r*log(k + (1.0 - k)*nu*u)
  res[u <= 1/nu] = res_pos
  res = res-max(res)
  res = exp(res)
  res = res/sum(res)
  if(log == T){
    res = log(res)
  }
  return(res)
}

u = seq(0.000001, 1-0.000001, length = 1e5)

##################################################################################
maxcenter = function(x) x - max(x)
meancenter = function(x) x - mean(x)

png("tests/tCCH.png", width=1000, height = 600)
par(mfrow=c(2,3))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
N = c(100, 200, 500, 1000)
n = tail(N,1)
p = 20

##################################################################################
par(mfrow=c(2,3))
ylim = c(-50, 0)
plot(u, u,
     type="n", ylim = ylim,
     xlab = "u=1/g-1", ylab = "log density", main = "betaprime")
lty=1
for(n in N){
  a=2; b= n -p-2; r=0; s=0; nu=1; k=1;
  lines(u, tCCH(u, a, b, r, s, nu, k, log = T), lty=lty, lwd = 1.5)
  lty = lty+1
}
legend("topright", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)

####
plot(u,u,
     type="n", ylim = ylim,
     xlab = "u=1/g-1", ylab = "log density", main = "zsadpated")
lty=1
for(n in N){
  a=1; b= 2; r=0; s=n+3; nu=1; k=1;
  lines(u, tCCH(u, a, b, r, s, nu, k, log = T), lty=lty, lwd = 1.5)
  lty = lty+1
}
legend("topright", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)

####
plot(u,u,
     type="n", ylim = ylim,
     xlab = "u=1/g-1", ylab = "log density", main = "benchmark")
lty=1
for(n in N){
  a= 0.02; b= 0.02 * max(n, p^2); r=0; s = 0; nu = 1; k = 1;
  lines(u, tCCH(u, a, b, r, s, nu, k, log = T), lty=lty, lwd = 1.5)
  lty = lty+1
}
legend("topright", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)

####
plot(u,u,
     type="n", ylim = ylim,
     xlab = "u=1/g-1", ylab = "log density", main = "uniform")
lty=1
for(n in N){
  a= 2; b= 2; r=0; s = 0; nu = 1; k = 1;
  lines(u, tCCH(u, a, b, r, s, nu, k, log = T), lty=lty, lwd = 1.5)
  lty = lty+1
}
legend("bottomleft", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)


####
plot(u,u,
     type="n", ylim = ylim,
     xlab = "u=1/g-1", ylab = "log density", main = "hyperg")
lty=1
for(n in N){
  a= 1; b= 2; r=0; s = 0; nu = 1; k = 1;
  lines(u, tCCH(u, a, b, r, s, nu, k, log = T), lty=lty, lwd = 1.5)
  lty = lty+1
}
legend("bottomleft", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)


####
plot(u,u,
     type="n", ylim = ylim,
     xlab = "u=1/g-1", ylab = "log density", main = "hypergn")
lty=1
for(n in N){
  a= 1; b= 2; r=1.5; s = 0; nu = 1; k = 1/n;
  lines(u, tCCH(u, a, b, r, s, nu, k, log = T), lty=lty, lwd = 1.5)
  lty = lty+1
}
legend("bottomleft", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)

mtext(paste0("Prior on u=1/g-1, logarized (p = ", p, ")"), outer= T, font=2)

dev.off()