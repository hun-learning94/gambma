x = fit3

plot(x, flist = f_list, show_plot = T, elapsed = elapsed3, ylim = ylim)


nknot = 5
IDX = seq_along(x$g_sampled)[sapply(x$knotlocs[-(1:500)], nrow) == nknot]
idx = sample(IDX, 1)
plot(x$Xgrid, x$predicted[idx,], type="l", lwd = 2, lty=2, col = "blue",
     xlab= "", ylab = "", axes = F)
tmp = seq(0, 1, len = nknot+2)
tmp = tmp[2:(length(tmp)-1)]
qts = quantile(x$data$x1, tmp)
x$predicted[idx, findInterval(qts, x$Xgrid)]
points(x = qts, y = x$predicted[idx, findInterval(qts, x$Xgrid)],
       cex = 2, pch = 19, col = "blue")

for(i in 1:4){
  idx = sample(IDX, 1)
  lines(x$Xgrid, x$predicted[idx,], type="l", lwd = 2, lty=2, col = "blue")
  tmp = seq(0, 1, len = nknot+2)
  tmp = tmp[2:(length(tmp)-1)]
  qts = quantile(x$data$x1, tmp)
  x$predicted[idx, findInterval(qts, x$Xgrid)]
  points(x = qts, y = x$predicted[idx, findInterval(qts, x$Xgrid)],
         cex = 2, pch = 19, col = "blue")
}

abline(v = qts, col = "lightgrey")

