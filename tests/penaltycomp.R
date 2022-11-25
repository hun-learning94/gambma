Rcpp::sourceCpp("src/gambma_series_diag.cpp")
n = 1000
p = 1:50
a = 2
b = n - p - 2
s = 0
g = n



##################################################################################

f1 = function(p, a, b){
  lbeta((a+p)/2, b/2)
}
f2 = function(p, g){
  -p/2 * log1p(g)
}
f3 = Vectorize(function(a, b, s, Qm, p){
  log_hyp1f1((a+p)/2, (a+b+p)/2, -(Qm+s)/2.0)
}, vectorize.args = c("p", "Qm", "b"))
f4 = function(g, Qm) -Qm/(2*g + 2)

maxcenter = function(x) x - max(x)
meancenter = function(x) x - mean(x)

png("tests/penaltycomp.png", width=1000, height = 600)
par(mfrow=c(2,3))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
# N = (10^(seq(2, 4, length = 7))) %/% 100 * 100
N = c(100, 200, 500, 1000)
n = tail(N,1)
Qm_mod = 2
Qm = Qm_mod*n

##################################################################################
n = tail(N, 1)
a = 2
b = n 
s = 0
g = N[length(N)]
total = c(maxcenter(f1(p, a, b) + f3(2, b, s, Qm, p)),
          maxcenter(f2(p, g) + f4(g, Qm)))
ylim = range(total)
plot(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), type = "n", ylim = ylim,
     xlab = "model size", ylab = "log penalty", main = "fixed g=n")
lty = 1
for(n in N){
  Qm = Qm_mod*n
  a = 2
  b = n 
  s = 2
  g = n
  lines(p, maxcenter(f2(p, g)+ f4(g, Qm)), lty=lty, lwd = 1.5)
  lty = lty + 1
}
abline(h = seq(0, -300,by=-50), col = "grey", lty=2)
legend("bottomleft", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)


n = tail(N, 1)
a = 2
b = n 
s = n
g = N[length(N)]
plot(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), type = "n", ylim = ylim,
     xlab = "model size", ylab = "log penalty", main = "betaprime")
lty = 1
for(n in N){
  Qm = Qm_mod*n
  a = 2
  b = n - p -2
  s = 0
  lines(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), lty=lty, lwd = 1.5)
  lty = lty + 1
}
abline(h = seq(0, -300,by=-50), col = "grey", lty=2)
legend("bottomleft", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)

n = tail(N, 1)
a = 2
b = n 
s = n
g = N[length(N)]
plot(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), type = "n", ylim = ylim,
     xlab = "model size", ylab = "log penalty", main = "zsadapted")
lty = 1
for(n in N){
  Qm = Qm_mod*n
  a = 1
  b = 2 
  s = n+3
  lines(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), lty=lty, lwd = 1.5)
  lty = lty + 1
}
abline(h = seq(0, -300,by=-50), col = "grey", lty=2)
legend("bottomleft", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)

n = tail(N, 1)
a = 2
b = n 
s = n
g = N[length(N)]
plot(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), type = "n", ylim = ylim,
     xlab = "model size", ylab = "log penalty", main = "uniform")
lty = 1
for(n in N){
  Qm = Qm_mod*n
  a = 2
  b = 2 
  s = 0
  lines(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), lty=lty, lwd = 1.5)
  lty = lty + 1
}
abline(h = seq(0, -300,by=-50), col = "grey", lty=2)
legend("bottomleft", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)

n = tail(N, 1)
a = 2
b = n 
s = n
g = N[length(N)]
plot(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), type = "n", ylim = ylim,
     xlab = "model size", ylab = "log penalty", main = "hyperg")
lty = 1
for(n in N){
  Qm = Qm_mod*n
  a = 1
  b = 2 
  s = 0
  lines(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), lty=lty, lwd = 1.5)
  lty = lty + 1
}
abline(h = seq(0, -300,by=-50), col = "grey", lty=2)
legend("bottomleft", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)


n = tail(N, 1)
a = 2
b = n 
s = n
g = N[length(N)]
plot(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), type = "n", ylim = ylim,
     xlab = "model size", ylab = "log penalty", main = "benchmark")
lty = 1
for(n in N){
  Qm = Qm_mod*n
  a = 0.02
  b = 0.02 * max(n, p^2) 
  s = 0
  lines(p, maxcenter(f1(p, 2, b) + f3(2, b, s, Qm, p)), lty=lty, lwd = 1.5)
  lty = lty + 1
}
abline(h = seq(0, -300,by=-50), col = "grey", lty=2)
legend("bottomleft", legend = paste("n", N), lty = seq_along(N), lwd = 1.5)

# mtext(paste0("CH(2,O(n),2) ","Qm = ", Qm_mod, "*n"), outer= T, font=2)

mtext(paste0("Penalty term, logarized (Qm = n x ", Qm_mod, ")"), outer= T, font=2)


dev.off()