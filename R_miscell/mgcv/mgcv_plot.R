library(zoo)

plot_mgcv = function(fit, Xgrid, flist = NULL, cred.int = 0.95, maxk, ylim = NULL, elapsed = NULL,
                     Linadj = NULL, n_row=NULL) {
  if(cred.int != 0.95) stop("only 0.95 availabe")
  p = ncol(Xgrid)
  if(is.null(Linadj)) Linadj = as.list(rep(F, length=p))
  
  opar = par(no.readonly = T)
  if(p < 4){plotarray = c(1, p)} else {if(is.null(n_row)) n_row = 2 ;plotarray = c(n_row, ceiling(p/n_row))}
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
  # png(filename = tempfile())
  TMP = mgcv::plot.gam(fit, se = qnorm(1 - (1 - cred.int)/2), n = 300)
  # dev.off()
  for(j in 1:p){
    # from data
    tmp = TMP[[j]]
    xdata = tmp$x
    Fits_mean = tmp$fit
    Fits_ub = Fits_mean + tmp$se
    Fits_lb = Fits_mean - tmp$se
    
    predtable = data.frame(grid = xdata, mean = Fits_mean, lb = Fits_lb, ub = Fits_ub)
    # fill in grid with na.approx
    xgrid = Xgrid[,j]
    predtable = merge(predtable, data.frame(grid = xgrid), by = "grid", all=T)
    predtable = as.data.frame(zoo::na.approx(predtable))
    predtable = as.data.frame(zoo::na.locf(predtable, fromLast = T, na.rm=F))
    predtable = as.data.frame(zoo::na.locf(predtable, fromLast = F, na.rm=F))
    predtable = predtable[match(xgrid, predtable$grid), ]
    if(!is.null(flist)){
      ytrue = flist[[j]](xgrid); ytrue = ytrue - mean(flist[[j]](xdata))
      if(Linadj[[j]] & abs(summary(fit)$edf[j] - 1) < 0.01){ ## adjustment for linear fit
        print("mgcv adjusted for linear fit")
        adjidx = (abs(diff(predtable$mean >0)) == 1)
        ytrue = ytrue - ytrue[adjidx]
        predtable$mean = predtable$mean - predtable$mean[adjidx]
        predtable$ub = predtable$ub - predtable$mean[adjidx]
        predtable$lb = predtable$lb - predtable$mean[adjidx]
      }
      
      cover = ifelse((ytrue > predtable$lb) & (ytrue < predtable$ub) , 1, 0)
      if(Linadj[[j]] & abs(summary(fit)$edf[j] - 1) < 0.01){ ## adjustment for linear fit
        if(all(cover[c(1, length(cover))] == 1)){
          cover = rep(1, length(cover))
        }else{
          cover = rep(0, length(cover))
        }
      }
      predtable$true = ytrue
      predtable$cover = cover
    }
    predtable$f = j
    PREDtable = rbind(PREDtable, predtable)
  }
  # PREDtable = PREDtable[, c(7, 1:6)]
  
  for(j in 1:p){
    predtable = PREDtable[PREDtable$f ==j ,]
    xgrid = predtable$grid
    Fits_mean = predtable$mean
    Fits_ub = predtable$ub
    Fits_lb = predtable$lb
    ytrue = predtable$true
    if(!is.null(ylim)){
      graphics::plot(xgrid, Fits_mean, type="n", xlab = paste0("x", j), ylab = paste0("f",j), ylim = ylim)
    } else {
      tmpylim = range(c(Fits_mean, Fits_ub, Fits_lb))
      tmpylim[1] = tmpylim[1] - 0.2; tmpylim[2] = tmpylim[2] + 0.2;
      graphics::plot(xgrid, Fits_mean, type="n", xlab = paste0("x", j), ylab = paste0("f",j), ylim =tmpylim)
    }
    graphics::polygon(x = c(xgrid, rev(xgrid)), y = c(Fits_lb, rev(Fits_ub)), col = "grey", border = NA)
    graphics::lines(xgrid, Fits_mean, col = "blue", lty=2)
    # if(fit$family$family == "gaussian") points(fit$model[,1+p], fit$model[,1], pch = 20)
    # abline(h = 0, lty=2, col = "skyblue")
    if(!is.null(flist)) graphics::lines(xgrid, ytrue, col = "red")
    if(!is.null(flist)){
      rmse = sqrt(mean((Fits_mean - ytrue)^2))
      title(main = paste0("logrmse", signif(log(rmse), 3)))
    }
  }
  title_ = paste0(fit$family$family, ", ", "mgcv", ", K ", maxk)
  if(!is.null(elapsed)) title_ = paste0(title_, ", elapsed ", paste0(signif(elapsed, 3), attr(elapsed, "units")))
  mtext(title_, outer= T, font=2)
  par(opar)
  
  return(invisible(PREDtable))
  
  
}







