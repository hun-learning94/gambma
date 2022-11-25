xrange = c(-1, 1)
ylim = c(-2, 2)

gressani1 = function(x){
  0.5 * (2*x^5 + 3*x^2 + cos(3*pi*x) - 1)
}

gressani3 = function(x){
  # x = x * 0.85
  sin(2*pi*x)*1.2
}

gressani6 = function(x){
  1.5 * (0.1 * sin(2*pi*x) + 0.2 * cos(2*pi*x) + 0.3 * sin(2*pi*x)^3 +
           0.4 * cos(2*pi*x)^3 + 0.5 * sin(2*pi*x)^3)
}

functions = list(f1 = function(x)gressani1(x)+0,
                 f2 = function(x)gressani3(x)+0,
                 f2 = function(x)gressani6(x)+0)
Linadj = as.list(c(F,F,F))
par(mfrow=c(1,3))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
xp = seq(-1, 1, length = 200)
for(i in 1:length(functions)){
  plot(xp, functions[[i]](xp), type="l", ylim = c(-2, 2)+0,
       xlab = paste0("x", i), ylab = paste0("f", i))
}
mtext("Hun 2022 smooth functions", outer=T, font=2)




#### 
# par(mfrow=c(1,3))
# par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0))
# y = x = seq(-1, 1, length = 100)
# z = outer(x, y, function(x, y) gressani2(x) + hun3(y))
# plot(x, gressani2(x), type = "l", main = "f1", ylab="", ylim = ylim)
# plot(y, hun3(y), type = "l", main = "f2", ylab= "", ylim = ylim)
# persp(x, y, z, theta= 45, phi =40, expand = 1, col = "lightblue", border = NA, shade = 0.5, ltheta = 30,
#       main = "f1 + f2")




