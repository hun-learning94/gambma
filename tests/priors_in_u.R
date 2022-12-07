library(latex2exp)
source("tests/specialfunctions.R")
plotdir = "C:/Users/USER/Dropbox/Kang/GAM5/plots/"
# a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1 # uniform
# a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1 # hyperg
# a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n # hypergn
# a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1 # betaprime
# a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1 # zsadapted
# a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1 # robust
# a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n # intrinsic

## tCCH density in g
dtCCH = function(uu, a, b, r, s, nu, kap){
  # if(min(uu) < 0 || max(uu) > 1/nu) stop("support is [0, 1/nu]")
  if(!(a>0)) stop("must be a > 0")
  if(!(b>0)) stop("must be b > 0")
  # if(!(nu>=1)) stop("must be nu >= 1")
  if(!(nu>0)) stop("must be nu > 0")
  if(!(kap>=0)) stop("must be kap > 0")
  
  # ## support
  # suppidx = (gg > nu-1)
  # uu = 1/(1+gg[suppidx])
  # logdens = rep(NA_real_, length(uu))
  # 
  # 
  # ## normalizing constant
  # lognormconst = a/2*log(nu) + s/(2*nu) - base::lbeta(a/2, b/2) - logPhi1(b/2, r, (a+b)/2, s/(2*nu), 1-kap)
  # ## density
  # logdens[suppidx] = 
  #   (a/2-1)*log(uu) + (b/2-1)*log1p(1-nu*uu) - s/2*uu - r*log(kap + (1-kap)*nu*uu) + lognormconst
  # 
  # return(logdens-0.5*log(gg+1))
  ## support
  suppidx = (uu < 1/nu && uu > 0)
  logdens = rep(NA_real_, length(uu))
  uu_in = uu[suppidx]
  
  
  ## normalizing constant
  lognormconst = a/2*log(nu) + s/(2*nu) - base::lbeta(a/2, b/2) - logPhi1(b/2, r, (a+b)/2, s/(2*nu), 1-kap)
  ## density
  logdens[suppidx] = 
    (a/2-1)*log(uu[suppidx]) + (b/2-1)*log1p(1-nu*uu[suppidx]) - s/2*uu[suppidx] - r*log(kap + (1-kap)*nu*uu) + lognormconst
  
  return(logdens)
}

momentf = function(r=1, a, b, z, s, nu, theta){
  exp(lbeta(r+a, b) - lbeta(a,b) + logPhi1(b,z,a+b+r, s, 1-theta) - logPhi1(b,z,a+b, s, 1-theta))
}

par(mfrow=c(2,4))
uu = seq(.001, 1-.001, len = 1e3)
N = c(200, 500, 1000)
n = 200; p = 10; Qm = n
for(gprior in 1:7){
  if(gprior == 1) {a = 2; b = 2; r = 0; s = 0; nu= 1; kap = 1} # uniform
  if(gprior == 2) {a = 1; b = 2; r = 0; s = 0; nu= 1; kap = 1} # hyperg
  if(gprior == 3) {a = 1; b = 2; r = 1.5; s = 0; nu= 1; kap = 1/n} # hypergn
  if(gprior == 4) {a = 1/2; b = n-p-1.5; r = 0; s = 0; nu= 1; kap = 1} # betaprime
  if(gprior == 5) {a = 1; b = 2; r = 0; s = n+3; nu= 1; kap = 1} # zsadapted
  if(gprior == 6) {a = 1; b = 2; r = 1.5; s = 0; nu= (n+1)/(p+1); kap = 1} # robust
  if(gprior == 7) {a = 1; b = 1; r = 1; s = 0; nu= (n+p+1)/(p+1); kap = (n+p+1)/n} # intrinsic
  if(gprior == 0) {a=b=r=s=nu=kap=0L; g=n} # fixed
  a_pos = (a+p)/2; b_pos = b/2; r_pos = r; s_pos = (s+Qm)/2; nu_pos = nu; kap_pos = kap
  
  plot(uu, exp(dtCCH(uu, a_pos, b_pos, r_pos, s_pos, nu_pos, kap_pos)), type = "l"); 
  abline(v=1/(1+n), lty=2); abline(v = 1/momentf(1, a_pos, b_pos,r_pos, s_pos, nu_pos,kap_pos)-1, lty=1)
}

# pdf(paste0(plotdir, "gpriors.pdf"), width=8, height = 4)
par(mfrow = c(2,4))
# par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0),cex.axis = 1., cex.lab = 1.)
par(mar=c(1.5, 2.5, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)

plot(g, hyperg(g), type = "n", main = "g=n", 
     ylab = TeX("$p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ abline(v = N[i], lty=i, col = i+1)}
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N))+1, lwd = 1.5, col = 1:length(N)+1)

plot(g, hyperg(g), type = "n", main = "Hyper-g", 
      ylab = TeX("$p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, hyperg(g), lty=i, col = i+1); }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, uniform(g), type = "n", main = "Uniform",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") )
for(i in seq_along(N)){ lines(g, uniform(g), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, hypergn(g, n), type = "n", main = "Hyper-g/n",  ylab = TeX("$p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, hypergn(g, N[i]), lty=i, col = i+1);  }
legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

# plot(g, benchmark(g, n, total_p), type = "l", main = "Benchmark (Ley, 2012)",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") )
# for(i in seq_along(N)){ lines(g, benchmark(g, N[i], total_p), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, truncZS(g, n), type = "n", main = "ZS-adapted",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") )
for(i in seq_along(N)){ lines(g, truncZS(g, N[i]), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, betaprime(g, n, p), type = "n", main = "Beta-prime",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") )
for(i in seq_along(N)){ lines(g, betaprime(g, N[i], p), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5,col = 1:length(N)+1)

plot(g, robust(g, n, p), type="n", main = "Robust",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") );
for(i in seq_along(N)){ lines(g, robust(g, N[i], p), lty=i, col = i+1)};
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, intrinsic(g, n, p), type="n", main = "Intrinsic",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") ); 
for(i in seq_along(N)){ lines(g, intrinsic(g, N[i], p), lty=i, col = i+1)}; 
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

# mtext(paste0("Mixtures of g-priors"), outer= T, font=2)
# dev.off()


pdf(paste0(plotdir, "gpriors_log.pdf"), width=8, height = 4)
par(mfrow = c(2,4))
# par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0),cex.axis = 1., cex.lab = 1.)
par(mar=c(1.5, 2.5, 1.5, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)

plot(g, hyperg(g, log=T), type = "n", main = "g=n", 
     ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ abline(v = N[i], lty=i, col = i+1)}
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N))+1, lwd = 1.5, col = 1:length(N)+1)

plot(g, hyperg(g, log=T), type = "n", main = "Hyper-g", 
     ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, hyperg(g, log=T), lty=i, col = i+1); }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, uniform(g, log=T), type = "n", main = "Uniform",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$") )
for(i in seq_along(N)){ lines(g, uniform(g, log=T), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, hypergn(g, n, log=T), type = "n", main = "Hyper-g/n",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$"))
for(i in seq_along(N)){ lines(g, hypergn(g, N[i], log=T), lty=i, col = i+1);  }
legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

# plot(g, benchmark(g, n, total_p), type = "l", main = "Benchmark (Ley, 2012)",  ylab = TeX("$p(g)$"), xlab = TeX("$g$") )
# for(i in seq_along(N)){ lines(g, benchmark(g, N[i], total_p), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, truncZS(g, n, log=T), type = "n", main = "ZS-adapted",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$") )
for(i in seq_along(N)){ lines(g, truncZS(g, N[i], log=T), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, betaprime(g, n, p, log=T), type = "n", main = "Beta-prime",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$") )
for(i in seq_along(N)){ lines(g, betaprime(g, N[i], p, log=T), lty=i, col = i+1);  }
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5,col = 1:length(N)+1)

plot(g, robust(g, n, p, log=T), type="n", main = "Robust",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$") );
for(i in seq_along(N)){ lines(g, robust(g, N[i], p, log=T), lty=i, col = i+1)};
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

plot(g, intrinsic(g, n, p, log=T), type="n", main = "Intrinsic",  ylab = TeX("$\\log\\; p(g)$"), xlab = TeX("$g$") ); 
for(i in seq_along(N)){ lines(g, intrinsic(g, N[i], p, log=T), lty=i, col = i+1)}; 
# legend("topright", legend = paste("n", c(N)), lty = 1:(length(N)), lwd = 1.5, col = 1:length(N)+1)

# mtext(paste0("Mixtures of g-priors"), outer= T, font=2)
dev.off()
