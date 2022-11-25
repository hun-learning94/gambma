
changes1 = function(x){
  hvec = c(1/2, 1/2, 1/2, 1/2,
           -1/2, -1/2, -1/2,
           0, 0,
           -1/2, -1/2, -1/2,
           0.8, 0.8)
  binnum = length(hvec)
  tvec = seq(0, 1, length = binnum + 1)
  out = 0
  for(i in seq_along(hvec)){
    out = out + hvec[i] * as.numeric((x >= tvec[i]) & (x < tvec[i+1]))
  }
  out
}

changes2 = function(x){
  hvec = c(1/2, 1/2, 1/2, 
           -1/2, -1/2, 
           0, -1/2, -1/2, 0,
           -1/2, -1/2, -1/2,
           0.8, 0.8)
  binnum = length(hvec)
  tvec = seq(0, 1, length = binnum + 1)
  out = 0
  for(i in seq_along(hvec)){
    out = out + hvec[i] * as.numeric((x >= tvec[i]) & (x < tvec[i+1]))
  }
  out
}

changes3 = function(x){
  binnum = 16
  tvec = seq(0, 1, length = binnum + 1)
  hvec = c(1, 1, .5, 1,
           -1, -1, -.5, -1,
           1, 1, .5, -1,
           -1, 0, -.5, 0)
  out = 0
  for(i in seq_along(hvec)){
    out = out + hvec[i] * as.numeric((x >= tvec[i]) & (x < tvec[i+1]))
  }
  out / 1.5
}


# changes2 = function(x){
#   delta = 2 / 3
#   hvec = rep(c(1- delta, rep(1, 3), 1 - delta, 
#                -1+delta, rep(-1, 3), -1 +delta), 3) * rep(c(1, 0.5, 1), each = 10)
#   tvec = seq(0, 1, length = 3*10 + 1)
#   out = 0
#   for(i in seq_along(hvec)){
#     out = out + hvec[i] * as.numeric((x >= tvec[i]) & (x < tvec[i+1]))
#   }
#   return(out*1.2)
# }

changes4 = function(x){
  binnum = 16
  tick = 0.2
  tvec = seq(0, 1, length = binnum + 1)
  hvec = c(1, 1 - tick, 1 - tick, 1,
           -1, -1 + tick, -1 + tick, -1,
           1, 1 - tick, 1 - tick, 1,
           -1, -1, 0, -tick)
  out = 0
  for(i in seq_along(hvec)){
    out = out + hvec[i] * as.numeric((x >= tvec[i]) & (x < tvec[i+1]))
  }
  out / 1.5
}


functions = list(f1 = changes1,
                 f2 = changes2,
                 f3 = changes3,
                 f4 = changes4)
par(mfrow=c(2, 2))
par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
xp = seq(1e-5, 1-1e-5, length = 1e4)
for(i in 1:length(functions)){
  plot(xp, functions[[i]](xp), type="l", ylim = c(-1.5, 1.5)+0,
       xlab = paste0("x", i), ylab = paste0("f", i))
}
mtext("Hun 2022 steps functions", outer=T, font=2)