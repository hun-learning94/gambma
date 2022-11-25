scaleto01 = function(X, mins=NULL, maxs=NULL){
  if(is.null(mins)) stop("need mins")
  if(is.null(maxs)) stop("need maxs")
  X = sweep(X, 2, mins, "-")
  X = sweep(X, 2, maxs-mins, "/")
  return(X)
}

scalefrom01 = function(X, mins=NULL, maxs = NULL){
  if(is.null(mins)) stop("need mins")
  if(is.null(maxs)) stop("need maxs")
  X = sweep(X, 2, maxs-mins, "*")
  X = sweep(X, 2, mins, "+")
  return(X)
}

# n = 200; p =10
# X = matrix(rnorm(n*p), nrow=n)
# mins = apply(X, 2, "min")
# maxs = apply(X, 2, "max")
# np = 200
# Xgrid = matrix(rep(seq(0, 1, length=np), p), nrow=np)
# Xgrid = scalefrom01(Xgrid, mins, maxs)
