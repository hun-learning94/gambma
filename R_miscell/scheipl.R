xrange = c(0, 1)
ylim = c(-2, 4)

scheipl1 = function(x){
  sqrt(x*(1-x)) * sin(18*pi / (8*x+1))
}
scheipl2 = function(x){
  1.5 * exp(-600*(x-0.4)^2) + 3 * exp(-500 * (x - 0.55)^2 ) + 4.5*exp(-500 * (x - 0.7)^2)
}
scheipl3 = function(x){
  hvec = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  tvec = c(.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81)
  out = 0
  for(i in seq_along(hvec)){
    out = out + hvec[i] * (1 + sign(x - tvec[i]))
  }
  # out = out - mean(out)
  # out = out / (max(out) - min(out))
  # out = 4*out
  return(out/2)
}
scheipl4 = function(x){
  4*sin(4*pi*x) - sign(x - 0.3) - sign(0.72 - x)
}

scheipl5 = function(x){
  tmp = 3*exp(-500 * (x - 0.15)^2) + 
    1.5 * exp(-400*(x-0.3)^2) + 
    3 * exp(-500 * (x - 0.6)^2 ) + 
    2*exp(-500 * (x - 0.8)^2)
  tmp * 1.2
}

functions = list(f1 = scheipl1,
              f2 = scheipl2,
              f3 = scheipl3,
              f4 = scheipl4,
              f5 = scheipl5)
par(mfrow=c(2,3))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
xp = seq(0, 1, length = 1e3)
for(i in 1:length(functions)){
  plot(xp, functions[[i]](xp), type="l", ylim = ylim,
       xlab = paste0("x", i), ylab = paste0("f", i))
}
mtext("scheipl 2009 functions", outer=T, font=2)