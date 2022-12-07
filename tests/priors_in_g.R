library(latex2exp)
plotdir = "C:/Users/USER/Dropbox/Kang/GAM8/plots/"

## hyper g prior
hyperg = function(g, log=F ){
  tmp = -3/2 * log1p(g) - log(2)
  if(log) return(tmp)
  return(exp(tmp))
}
uniform = function(g,log=F ){
  a = 2; b = 0.5*(n-p)-1
  tmp = -lbeta(a,b) - (a+b)*log1p(g) + (b-1)*log(g)
  if(log) return(tmp)
  return(exp(tmp))
}
## hyper g/n prior
hypergn = function(g, n, log=F ){
  tmp = (-3/2 * log1p(g/n) -log(2*n))
  if(log) return(tmp)
  return(exp(tmp))
}
## truncated ZS
truncZS = function(g, n, log=F ){
  a = 0.5; b = (n+3)/2
  log_incgamma = lgamma(a) + pgamma(b, shape = a, rate = 1, lower.tail = T, log= T)
  tmp = ifelse(g < 1,
               NA_real_,
               a*log(b) - log_incgamma - (a+1)*log1p(g) - b / (g+1))
  if(log) return(tmp)
  return(exp(tmp))
}
## Beta prior (Betaprime)
betaprime = function(g, n, p, log=F ){
  a = 1/4; b = 0.5*(n-p-1.5)
  tmp = -lbeta(a,b) - (a+b)*log1p(g) + (b-1)*log(g)
  if(log) return(tmp)
  return(exp(tmp))
}
 
## robust prior 
robust = function(g, n, p, log=F ){
  tmp = ifelse(g+1 < (1+n) / (1+p),
               NA_real_,
         -log(2) + 0.5* (log1p(n) - log1p(p)) - 1.5 * log1p(g))
  if(log) return(tmp)
  return(exp(tmp))
}
## intrinsic prior
intrinsic = function(g, n, p, log=F ){
  tmp = ifelse(g < n/(p+1), NA_real_,
               -lbeta(0.5, 0.5) - log(g) - 0.5*log(g*(p+1)/n - 1))
  if(log) return(tmp)
  return(exp(tmp))
}

n = 200
total_p = 30
# p = 0:total_p
p = 10
g = seq(0.0000001, 1000, len = 1e3)
N = c(200, 500, 1000)

pdf(paste0(plotdir, "gpriors.pdf"), width=8, height = 4)
par(mfrow = c(2,4))
ylim = c(0, 0.01)
# par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0),cex.axis = 1., cex.lab = 1.)
par(mar=c(1.5, 2.5, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)

plot(g, hyperg(g), type = "n", main = "g=n", 
     ylab = TeX("$p(g)$"), xlab = TeX("$g$"), ylim = ylim)
for(i in seq_along(N)){ abline(v = N[i], lty=i, col = i+1)}
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N))+1, lwd = 1.5, col = 1:length(N)+1)

plot(g, hyperg(g), type = "n", main = "Hyper-g", 
      ylab = TeX("$p(g)$"), xlab = TeX("$g$"), ylim = ylim)
for(i in seq_along(N)){ lines(g, hyperg(g), lty=i, col = i+1); }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, uniform(g), type = "n", main = "Uniform",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") , ylim = ylim)
for(i in seq_along(N)){ lines(g, uniform(g), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, hypergn(g, n), type = "n", main = "Hyper-g/n",  ylab = TeX("$p(g)$"), xlab = TeX("$g$"), ylim = ylim)
for(i in seq_along(N)){ lines(g, hypergn(g, N[i]), lty=i, col = i+1);  }
legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, truncZS(g, n), type = "n", main = "ZS-adapted",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") , ylim = ylim)
for(i in seq_along(N)){ lines(g, truncZS(g, N[i]), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, betaprime(g, n, p), type = "n", main = "Beta-prime",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") , ylim = ylim)
for(i in seq_along(N)){ lines(g, betaprime(g, N[i], p), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5,col = 1:length(N)+1)

plot(g, robust(g, n, p), type="n", main = "Robust",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") , ylim = ylim);
for(i in seq_along(N)){ lines(g, robust(g, N[i], p), lty=i, col = i+1)};
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, intrinsic(g, n, p), type="n", main = "Intrinsic",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") , ylim = ylim); 
for(i in seq_along(N)){ lines(g, intrinsic(g, N[i], p), lty=i, col = i+1)}; 
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

# mtext(paste0("Mixtures of g-priors"), outer= T, font=2)
dev.off()


pdf(paste0(plotdir, "gpriors_log.pdf"), width=6, height = 3)
par(mfrow = c(2,4))
# ylim = c(-16, -1)
# par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0),cex.axis = 1., cex.lab = 1.)
par(mar=c(3, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)

plot(g, hyperg(g, log=T), type = "n", main = "g=n", 
     ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ abline(v = N[i], lty=i, col = i+1)}
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N))+1, lwd = 1.5, col = 1:length(N)+1)

plot(g, hyperg(g, log=T), type = "n", main = "Hyper-g", 
     ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, hyperg(g, log=T), lty=i, col = i+1); }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, uniform(g, log=T), type = "n", main = "Uniform",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, uniform(g, log=T), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, hypergn(g, n, log=T), type = "n", main = "Hyper-g/n",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, hypergn(g, N[i], log=T), lty=i, col = i+1);  }
legend("bottomleft", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, truncZS(g, n, log=T), type = "n", main = "ZS-adapted",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, truncZS(g, N[i], log=T), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, betaprime(g, n, p, log=T), type = "n", main = "Beta-prime",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, betaprime(g, N[i], p, log=T), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5,col = 1:length(N)+1)

plot(g, robust(g, n, p, log=T), type="n", main = "Robust",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, robust(g, N[i], p, log=T), lty=i, col = i+1)};
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, intrinsic(g, n, p, log=T), type="n", main = "Intrinsic",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, intrinsic(g, N[i], p, log=T), lty=i, col = i+1)}; 
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

# mtext(paste0("Mixtures of g-priors"), outer= T, font=2)
dev.off()

