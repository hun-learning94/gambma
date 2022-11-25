xrange = c(-1, 1)
ylim = c(-2, 2)

gressani1 = function(x){
  0.5 * (2*x^5 + 3*x^2 + cos(3*pi*x) - 1)
}
gressani2 = function(x){
  1.3 * x^5 + sin(4*x) + 0.75*x^2 - 0.25
}

gressani4 = function(x){
  exp(-x^3) * sin(2*pi*x^2) - 0.1
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
  1.3 * (0.1 * sin(2*pi*x)^2 + 0.2 * cos(2*pi*x)^2 + 0.3 * sin(2*pi*x)^3 +
           0.4 * cos(2*pi*x)^3 + 0.5 * sin(2*pi*x)^3) - 0.1
}
hun3 = function(x){
  x = x + 0.1
  sin(2*pi*x)
}

functions = list(f1 = function(x)gressani1(x)+1/2,
                 f2 = function(x)gressani2(x)+1/2,
                 f3 = function(x)gressani4(x)+1/2,
                 f4 = function(x)gressani6(x)+1/2,
                 # f5 = function(x)hun1(x)+1,
                 f6 = function(x)hun2(x)+1/2,
                 f7 = function(x) hun3(x) + 1/2)

par(mfrow=c(2,3))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
xp = seq(-1, 1, length = 200)
for(i in 1:length(functions)){
  plot(xp, functions[[i]](xp), type="l", ylim = c(-2, 2)+1/2,
       xlab = paste0("x", i), ylab = paste0("f", i))
}
mtext("Hun 2022 smooth functions", outer=T, font=2)




##
# par(mfrow=c(1,3))
# par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
# y = x = seq(-1, 1, length = 100)
# z = outer(x, y, function(x, y) gressani2(x) + hun3(y))
# plot(x, gressani2(x), type = "l", main = "f1", ylab="", ylim = ylim)
# plot(y, hun3(y), type = "l", main = "f2", ylab= "", ylim = ylim)
# persp(x, y, z, theta= 45, phi =40, expand = 1, col = "lightblue", border = NA, shade = 0.5, ltheta = 30,
#       main = "f1 + f2")


















