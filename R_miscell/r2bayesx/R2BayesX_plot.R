library(zoo)

plot_bayesx = function(fit, Xgrid, flist = NULL, cred.int = 0.95, maxk, ylim = NULL, elapsed = NULL) {
  if(cred.int != 0.95) stop("only 0.95 availabe")
  
  p = ncol(Xgrid)
  
  opar = par(no.readonly = T)
  if(p < 4){plotarray = c(1, p)} else {plotarray = c(2, ceiling(p/2))}
  par(mfrow=plotarray)
  par(mar=c(2.5, 3, 2, 1), mgp=c(1.5,0.5,0), oma = c(0,0,2,0))
  
  PREDtable = data.frame(
    f = integer(),
    grid = double(),
    mean = double(),
    lb = double(),
    ub = double()
  )
  if(!is.null(flist)){
    PREDtable$true = double()
    PREDtable$cover = integer()
  }
  
  for(j in 1:p){
    # from data
    tmp = fit$effects[[j]][,c(1,2,3,7)]
    xdata = tmp[,1]
    Fits_mean = tmp[,2]
    Fits_ub = tmp[,4]
    Fits_lb = tmp[,3]
    predtable = data.frame(grid = xdata, mean = Fits_mean, lb = Fits_lb, ub = Fits_ub)
    # fill in grid with na.approx
    xgrid = Xgrid[,j]
    predtable = merge(predtable, data.frame(grid = xgrid), by = "grid", all=T)
    predtable = as.data.frame(zoo::na.approx(predtable))
    predtable = as.data.frame(zoo::na.locf(predtable, fromLast = T, na.rm=F))
    predtable = as.data.frame(zoo::na.locf(predtable, fromLast = F, na.rm=F))
    predtable = predtable[match(xgrid, predtable$grid), ]
    if(!is.null(flist)){
      ytrue = flist[[j]](xgrid); ytrue = ytrue - mean(ytrue)
      cover = ifelse((ytrue > predtable$lb) & (ytrue < predtable$ub) , 1, 0)
      predtable$true = ytrue
      predtable$cover = cover
    }
    predtable$f = j
    PREDtable = rbind(PREDtable, predtable)
  }
  PREDtable = PREDtable[, c(7, 1:6)]
  
  if(is.null(ylim))ylim = c(min(PREDtable$mean), max(PREDtable$mean))
  for(j in 1:p){
    predtable = PREDtable[PREDtable$f ==j ,]
    xgrid = predtable$grid
    Fits_mean = predtable$mean
    Fits_ub = predtable$ub
    Fits_lb = predtable$lb
    ytrue = predtable$true
    graphics::plot(xgrid, Fits_mean, type="n", xlab = paste0("x", j), ylab = paste0("f",j), ylim = ylim)
    graphics::polygon(x = c(xgrid, rev(xgrid)), y = c(Fits_lb, rev(Fits_ub)), col = "grey", border = NA)
    graphics::lines(xgrid, Fits_mean, col = "blue", lty=2)
    if(!is.null(flist)) graphics::lines(xgrid, ytrue, col = "red")
    if(!is.null(flist)){
      rmse = sqrt(mean((Fits_mean - ytrue)^2))
      title(main = paste0("logrmse", signif(log(rmse), 3)))
    }
  }
  title_ = paste0(fit$call$family, ", ", "bayesx", ", K ", maxk)
  if(!is.null(elapsed)) title_ = paste0(title_, ", elapsed ", paste0(signif(elapsed, 3), attr(elapsed, "units")))
  mtext(title_, outer= T, font=2)
  par(opar)
  
  return(invisible(PREDtable))
  
  
}







