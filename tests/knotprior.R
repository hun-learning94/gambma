p = 0:20

pri_norm = function(p, mu){
  tmp = (-0.5/4^2 * (p - mu)^2)
  prob = exp(tmp - max(tmp))
  prob = prob/sum(prob)
  log(prob)
}

pri_poi = function(p, mu){
  tmp = -mu + p*log(mu) - lgamma(p+1)
  prob = exp(tmp - max(tmp))
  prob = prob/sum(prob)
  log(prob)
}

png("tests/knotprior.png", width=1000, height = 400)
par(mfrow=c(1,2))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
plot(p,pri_poi(p, 5),type="l", lty=1, lwd = 1.5,
     main = "knot prior",
     xlab = "number of knots",
     ylab = "normalized log probability")
lines(p, pri_poi(p, 10), type="l", lty=2, lwd = 1.5)
lines(p, pri_norm(p, 10), type="l", lty=3, lwd = 1.5)
legend("bottomleft", legend = c("poi5", "poi10", "N(10, 16)"), lty = 1:3, lwd = 1.5)

plot(p,exp(pri_poi(p, 5)),type="l", lty=1, lwd = 1.5,
     main = "knot prior",
     xlab = "number of knots",
     ylab = "normalized probability")
lines(p, exp(pri_poi(p, 10)), type="l", lty=2, lwd = 1.5)
lines(p, exp(pri_norm(p, 10)), type="l", lty=3, lwd = 1.5)
legend("topright", legend = c("poi5", "poi10", "N(10, 16)"), lty = 1:3, lwd = 1.5)
dev.off()