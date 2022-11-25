plot.gamlps = function (x, xp, smoo.index, cred.int = 0.95, plot.cred = TRUE, 
                        np = 100, fit.col = "blue", shade.col = "gray75", 
                        show.plot = TRUE, show.info = TRUE, ...) 
{
  if (smoo.index < 1 || smoo.index > x$q) 
    stop("smoo.index wrongly specified")
  smoo.index <- as.integer(smoo.index)
  if (!is.vector(cred.int, mode = "numeric") || length(cred.int) > 
      1 || is.na(cred.int) || is.infinite(cred.int) || cred.int <= 
      0 || cred.int >= 1) 
    stop("cred.int must be between 0 and 1")
  if (!is.logical(plot.cred)) 
    stop("plot.cred must be either TRUE or FALSE")
  if (np < 20 || np > 200) 
    stop("choose np between 20 and 200")
  if (is.null(x$data)) {
    mf <- stats::model.frame(x$formula)
    X <- stats::model.matrix(mf)
    smterms <- grepl("sm(", colnames(X), fixed = TRUE)
    X <- cbind(X[, as.logical(1 - smterms)], X[, smterms])
  }
  else {
    mf <- stats::model.frame(x$formula, data = x$data)
    X <- stats::model.matrix(mf, data = x$data)
    smterms <- grepl("sm(", colnames(X), fixed = TRUE)
    X <- cbind(X[, as.logical(1 - smterms)], X[, smterms])
  }
  q <- x$q
  p <- ncol(X) - q
  n <- x$n
  K <- x$K
  splines <- x$spline.estim
  j <- smoo.index
  xj <- as.numeric(X[, p + j])
  min.xgrid <- min(xj)
  max.xgrid <- max(xj)
  if (!missing(xp)) {
    if (any(xp < min.xgrid) || any(xp > max.xgrid)) 
      stop("values in xp not in the range of observed covariates")
  }
  xgrid <- seq(min.xgrid, max.xgrid, length = np)
  xj.fine <- seq(min.xgrid, max.xgrid, length = 1000)
  Bj.fine <- cubicbs(xj.fine, lower = min.xgrid, upper = max.xgrid, 
                     K = K)$Bmatrix
  Bj.fine.mean <- colMeans(Bj.fine)
  Bxgrid <- cubicbs(xgrid, lower = min.xgrid, upper = max.xgrid, 
                    K = K)$Bmatrix
  Bxgrid.centered <- Bxgrid - matrix(rep(Bj.fine.mean, np), 
                                     nrow = np, byrow = TRUE)
  Bx <- Bxgrid.centered[, -K]
  fhat <- as.numeric(Bx %*% splines[[j]])
  alpha <- 1 - cred.int
  latmaximum <- x$latmaximum
  Covmaximum <- x$Covmaximum
  thetaj.max <- latmaximum[((p + 1) + (j - 1) * (K - 1)):(p + 
                                                            j * (K - 1))]
  Sigj.max <- Covmaximum[((p + 1) + (j - 1) * (K - 1)):(p + 
                                                          j * (K - 1)), ((p + 1) + (j - 1) * (K - 1)):(p + j * 
                                                                                                         (K - 1))]
  postfj.mean <- as.numeric(Bx %*% thetaj.max)
  postfj.sd <- sqrt(diag(Bx %*% Sigj.max %*% t(Bx)))
  fj.lb <- fhat - stats::qnorm((1 - (alpha * 0.5))) * postfj.sd
  fj.ub <- fhat + stats::qnorm((1 - (alpha * 0.5))) * postfj.sd
  minf <- min(fj.lb)
  maxf <- max(fj.ub)
  shift.frame <- 0.3 * (maxf - minf)
  covariate.name <- colnames(X)[(p + 1):(p + q)][smoo.index]
  nchar.covariate <- nchar(covariate.name)
  covariate.name <- substr(covariate.name, 4, nchar.covariate - 
                             1)
  if (show.plot == TRUE) {
    graphics::plot(xgrid, fhat, type = "l", col = fit.col, xlab = covariate.name, 
                   ylab = paste0("sm(", covariate.name, ",", format(round(x$EDf[smoo.index], 2), nsmall = 2),
                                 ")"), ...)
    if (plot.cred == TRUE) {
      graphics::polygon(x = c(xgrid, rev(xgrid)), y = c(fj.lb, rev(fj.ub)), col = shade.col, border = NA)
    }
    graphics::lines(xgrid, fhat, type = "l", lwd = 2, col = fit.col)
    graphics::rug(X[, p + smoo.index])
  }
  
  if (!missing(xp)) {
    xxgrid <- xp
    xlen <- length(xp)
    Bxxgrid <- cubicbs(xp, lower = min.xgrid, upper = max.xgrid, 
                       K = K)$Bmatrix
    Bxxgrid.centered <- Bxxgrid - matrix(rep(Bj.fine.mean, 
                                             xlen), nrow = xlen, byrow = TRUE)
    Bxx <- Bxxgrid.centered[, -K]
    fhatxx <- as.numeric(Bxx %*% splines[[j]])
    postfjxx.mean <- as.numeric(Bxx %*% thetaj.max)
    postfjxx.sd <- sqrt(diag(Bxx %*% Sigj.max %*% t(Bxx)))
    fjxx.lb <- fhatxx - stats::qnorm((1 - (alpha * 0.5))) * 
      postfjxx.sd
    fjxx.ub <- fhatxx + stats::qnorm((1 - (alpha * 0.5))) * 
      postfjxx.sd
    ftable <- matrix(0, nrow = xlen, ncol = 4)
    colnames(ftable) <- c("xp", "sm", "sm.low", 
                          "sm.up")
    ftable[, 1] <- xp
    ftable[, 2] <- fhatxx
    ftable[, 3] <- fjxx.lb
    ftable[, 4] <- fjxx.ub
    ftable <- round(ftable, 4)
    if (show.info == TRUE) {
      cat("Estimated smooth function", paste0("sm(", 
                                              covariate.name, ")"), "at specified grid points (*): \n")
      cat("\n")
      print.table(format(ftable, nsmall = 4), right = TRUE)
      cat("--- \n")
      cat("* Bounds correspond to a", paste(format(round(cred.int * 
                                                           100, 2), nsmall = 2), "%", sep = ""), 
          "credible interval. \n")
    }
    
    listout <- list(xp = xp, sm.xp = fhatxx, sm.low = fjxx.lb, 
                    sm.up = fjxx.ub, cred.int = cred.int, smoo.index = smoo.index, ftable = ftable)
    return(invisible(listout))
  }
}


library(zoo)
plot_gamlps = function(fit, Xgrid, flist = NULL, cred.int = 0.95, ylim = NULL, elapsed = NULL, Linadj= NULL) {
  
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
  
  idx = 1:np
  for(j in 1:p){
    xgrid = Xgrid[,j]
    tmp = plot.gamlps(fit, xp = xgrid, smoo.index = j, cred.int=0.95, 
                      show.info = F, show.plot = F)$ftable
    xdata = tmp[,1]
    Fits_mean = tmp[,2]
    Fits_ub = tmp[,4]
    Fits_lb = tmp[,3]
    predtable = data.frame(f = j, grid = xgrid, mean = Fits_mean, lb = Fits_lb, ub = Fits_ub)
    # fill in grid with na.approx
    xgrid = Xgrid[,j]
    predtable = merge(predtable, data.frame(grid = xgrid), by = "grid", all=T)
    predtable = as.data.frame(zoo::na.approx(predtable))
    predtable = as.data.frame(zoo::na.locf(predtable, fromLast = T, na.rm=F))
    predtable = as.data.frame(zoo::na.locf(predtable, fromLast = F, na.rm=F))
    predtable = predtable[match(xgrid, predtable$grid), ]
    if(!is.null(flist)){
      ytrue = flist[[j]](xgrid); ytrue = ytrue - mean(flist[[j]](xdata))
      if(Linadj[[j]] & abs(fit$EDf[j] - 1) < 0.01){ ## adjustment for linear fit
        print("blapsr adjusted for linear fit")
        adjidx = (abs(diff(predtable$mean >0)) == 1)
        ytrue = ytrue - ytrue[adjidx]
        predtable$mean = predtable$mean - predtable$mean[adjidx]
        predtable$ub = predtable$ub - predtable$mean[adjidx]
        predtable$lb = predtable$lb - predtable$mean[adjidx]
      }
      
      cover = ifelse((ytrue > predtable$lb) & (ytrue < predtable$ub) , 1, 0)
      if(Linadj[[j]] & abs(fit$EDf[j] - 1) < 0.01){ ## adjustment for linear fit
        if(all(cover[c(1, length(cover))] == 1)){
          cover = rep(1, length(cover))
        }else{
          cover = rep(0, length(cover))
        }
      }
      predtable$true = ytrue
      predtable$cover = cover
    }
    PREDtable = rbind(PREDtable, predtable)
    idx = idx + np
  }
  
  if(is.null(ylim))ylim = c(min(PREDtable$mean), max(PREDtable$mean))
  idx = 1:np
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
    idx = idx + np
    if(!is.null(flist)){
      rmse = sqrt(mean((Fits_mean - ytrue)^2))
      title(main = paste0("logrmse", signif(log(rmse), 3)))
    }
  }
  
  title_ = paste0(fit$family, ", ", "blapsr", ", K ", fit$K)
  if(!is.null(elapsed)) title_ = paste0(title_, ", elapsed ", paste0(signif(elapsed, 3), attr(elapsed, "units")))
  mtext(title_, outer= T, font=2)
  par(opar)
  
  return(invisible(PREDtable))
  
  
}







