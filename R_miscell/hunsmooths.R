xrange = c(-1, 1)
ylim = c(-2, 2)

gressani1 = function(x){
  0.5 * (2*x^5 + 3*x^2 + cos(3*pi*x) - 1)
}
gressani2 = function(x){
  1.3 * x^5 + sin(4*x) + 0.75*x^2 - 0.25
}

gressani3 = function(x){
  # x = x * 0.85
  sin(4*pi*x)
}

gressani4 = function(x){
  exp(-x^3) * sin(2*pi*x^2) - 0.1
}

gressani5 = function(x){
  x = x * 1.05
  0.8 * x^2 * (x^3 + 2*exp(-3*x^4 + log(2*x + pi))) - 0.65
}

# hun1 = function(x){
#   x = (x-0.05) *0.72
#   -0.8*(0.5*x^4 + 1/2*x^2 + sin(3*pi*x * (log1p(abs(x*2)-0.3)))* exp(x)/1.2 - 1) - 0.7
# }

hun2 = function(x){
  x = x * 0.9
  x = (x-0.1) * 0.9
  2 * sin(5*pi*x * cos(x)) * log1p(abs(x)) + x^5/4
}

gressani6 = function(x){
  # 1.3 * (0.1 * sin(2*pi*x)^2 + 0.2 * cos(2*pi*x)^2 + 0.3 * sin(2*pi*x)^3 +
  #          0.4 * cos(2*pi*x)^3 + 0.5 * sin(2*pi*x)^3) - 0.1
  1.5 * (0.1 * sin(2*pi*x) + 0.2 * cos(2*pi*x) + 0.3 * sin(2*pi*x)^3 +
           0.4 * cos(2*pi*x)^3 + 0.5 * sin(2*pi*x)^3) - 0.22
}

wood1 = function(x){
  x = (x+1)/2
  tmp = x^11 * (10*(1-x))^6/5 + 10^4 * x^3 * (1-x)^10
  tmp / 5
}

jeong1 = function(x){
  x = (x+1)/2
  3*exp(-200 * (x - 0.2)^2) + 0.5 * exp(-50*(x-0.6)^2)
}

francom3 = function(x){
  0.75*(0.0035 * (x*3 + 1.5)^3 + (x > -0.5 & x < 0.85) * 0.07 *sin(1.7*pi*(x*3 + 1.5)^2 / 3.2)*(x*3 -2.5)^2 * exp(x*3 + 1.5))
}

dimatteo1 = function(x){
  x = 2*x
  tmp = sin(x) + 2*exp(-30*x^2)
  tmp-0.3
}

hun3 = function(x){
  x
}

hun4 = function(x){
  x = (x+1.1)
  log(x) + 0.4
}

functions = list(f1 = function(x)gressani1(x)+0,
                 f2 = function(x)dimatteo1(x)+0,
                 f3 = function(x)hun3(x)+0,
                 f4 = function(x)francom3(x) + 0,
                 f5 = function(x)gressani2(x)+0,
                 f6 = function(x)hun4(x)+0)
Linadj = as.list(c(F,F,T,F,F,F,F))
par(mfrow=c(2,3))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
xp = seq(-1, 1, length = 200)
for(i in 1:length(functions)){
  plot(xp, functions[[i]](xp), type="l", ylim = c(-2, 2)+0,
       xlab = paste0("x", i), ylab = paste0("f", i))
}
mtext("Hun 2022 smooth functions", outer=T, font=2)


pdf("results/sim1fs.pdf", width = 10, height = 2.3)
fuck = c(1,5,4)
fuck2 = c(1,2,3)
plotfs = functions[fuck]
par(mfrow= c(1,3))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)
for(i in 1:length(plotfs)){
  plot(xp, plotfs[[i]](xp), type="l", ylim = c(-2, 2)+0,lwd =1.5,
       xlab = paste0("x", fuck2[i]), ylab = paste0("f", fuck2[i]))
}
dev.off()

pdf("results/sim2fs.pdf", width = 10, height = 2.3)
fuck = c(2,3,6)
fuck2 = c(4,5,6)
plotfs = functions[fuck]
par(mfrow= c(1,3))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)
for(i in 1:length(plotfs)){
  plot(xp, plotfs[[i]](xp), type="l", ylim = c(-2, 2)+0,lwd =1.5,
       xlab = paste0("x", fuck2[i]), ylab = paste0("f", fuck2[i]))
}
dev.off()











