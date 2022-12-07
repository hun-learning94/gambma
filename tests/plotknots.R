plotdir = "C:/Users/USER/Dropbox/Kang/GAM6/plots/"
set.seed(2)
ypts = function(pts) rep(0, length(pts))
xlim = c(-5, 5); ylim = c(-2, 2)
pts = c(-1, 1)


plotknot = function(pts, grid = F){
  plot(pts, ypts(pts), pch=16, xlim = xlim, ylim = ylim, cex = 1.5,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  if(grid)points(-5:5,  ypts(-5:5), pch="|")
  abline(h=0)
  points(pts, ypts(pts), pch = 16, col = "orange", cex = 1.1)
}

pdf(paste0(plotdir,"plotknots.pdf"), width=6, height = 3.5)
par(mfcol=c(3,3))
# par(mar=c(3, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,0,0), cex.axis = 1.0, cex.main= 1)
par(mar=c(0, 1/3, 1/2, 1/3), mgp=c(1.5,0.5,0), oma = c(1,0,2,0), cex.axis = 1.0, cex.main= 1)
## even knots
plotknot(c(-1, 1))
plotknot(c(-2, 0, 2))
plotknot(c(-3, -1, 1,3))
mtext("Even-knot", side=3, outer=T, at = 0.165, cex=0.7)
## vs knots
plotknot(c(-2, 0),T)
plotknot(c(-3,-1, 2),T)
plotknot(c(-3, 1, 4, 5),T)
mtext("VS-knot", side=3, outer=T, at = 0.5, cex=0.7)
## free knots
plotknot(c(-2-0.5, 0.3))
plotknot(c(-3,-1, 2) + c(0.2, -0.2, 0.3))
plotknot(c(-3, 1, 4, 5)+ c(-0.5, -0.1, 0.2, -0.1))
mtext("Free-knot", side=3, outer=T, at = 0.83, cex=0.7)
dev.off()