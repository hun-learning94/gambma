xrange = c(0, 1)
ylim = c(-3, 3)

francom1 = function(x){
  x = 2*x - 1
  x = x*4
  x = x + 0.66
  tmp = 0.1 * x^3 + (x>0 & x<4) * 2*sin(pi*x^2 /1.9)*(x-4)^2 * exp(x)
  tmp = (tmp)/ (max(tmp) - min(tmp))
  tmp = 5*tmp
  tmp*2
}

francom2 = function(x){
  x = 2*x - 1
  x = x*4
  x = x + 1.
  tmp = 0.1 * x^3 + (x>0 & x<4) * 2*sin(pi*x^2 /1.2)*(x-4)^2 * exp(x)
  tmp = (tmp)/ (max(tmp) - min(tmp))
  tmp = 7.5*tmp
  tmp
}
francom3 = function(x){
  x = 2*x - 1
  x = x*3
  x = x + 2
  tmp = 0.1 * x^3 + (x>0 & x<4) * 2*sin(pi*x^2 / 1.)*(x-4)^2 * exp(x)
  tmp = (tmp)/ (max(tmp) - min(tmp))
  tmp = 4*tmp
  tmp + 1
}
scheipl1 = function(x){
  x = x + 0.07
  x = x / max(x)
  sqrt(x*(1-x)) * sin(15*pi / (9*x+1)) * 5
}

hundoppler = function(x){
  gressani6 = function(x){
    x = 2*x
    1.3 * (0.1 * sin(2*pi*x)^2 + 0.2 * cos(2*pi*x)^2 + 0.3 * sin(2*pi*x)^3 +
             0.4 * cos(2*pi*x)^3 + 0.5 * sin(2*pi*x)^3) - 0.1
  }
  cutoff = -0.24
  hun3 = function(x){
    x = (x - cutoff)
    x = x/2.5
    sin(2*pi*x)^2 * x^(1/5)
  }
  x = 2*x - 1
  ifelse(x<cutoff, gressani6(x), hun3(x) + gressani6(cutoff) - hun3(cutoff))*2.5
}

gressani2 = function(x){
  x = 2*x - 1
  res = 1.3 * x^5 + sin(4*x) + 0.75*x^2 - 0.25
  res * 2
}

gressani1 = function(x){
  x = 2*x - 1
  res = 0.5 * (2*x^5 + 3*x^2 + cos(3*pi*x) - 1)
  res * 1.5
}

gressani4 = function(x){
  x = 2*x - 1
  x = x * 0.95
  res = exp(-x^3) * sin(2*pi*x^2)
  res*1.5
}


xp = seq(xrange[1], xrange[2], length = 200)

functions = list(f1 = francom1,
               f2 = scheipl1,
               f3 = francom3,
               f4 = gressani2,
               f5 = gressani4,
               f6 = francom2)

par(mfrow=c(2,3))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
for(i in 1:length(functions)){
  plot(xp, functions[[i]](xp), type="l", 
       xlab = paste0("x", i), ylab = paste0("f", i), ylim = ylim)
}
if(i < 6) plot.new()
mtext("Francom 2020 functions", outer=T, font=2)