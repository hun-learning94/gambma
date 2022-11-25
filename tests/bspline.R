# https://cran.r-project.org/web/packages/crs/vignettes/spline_primer.pdf
bs = function(x, deg, i, knots){
  if(deg == 0){
    B = ifelse((x >= knots[i]) & (x < knots[i+1]), 1, 0)
  } else {
    if(knots[deg+i] - knots[i] == 0){ w1 = 0 # for clamped boundary
    } else { w1 = (x - knots[i]) / (knots[deg+i] - knots[i]) }
    
    if(knots[deg+i+1] - knots[i+1] == 0){ w2 = 0 # for clamped boundary
    } else { w2 = (knots[deg+i+1] - x) / (knots[deg+i+1] - knots[i+1]) }
    
    B = w1 * bs(x, deg-1, i, knots) + w2 * bs(x, deg-1, i+1, knots)
  }
  return(B)
}

BS = function(x, deg=3, intKnots = NULL, intercept = F, bdKnots = c(0,1)){
  # input check
  if(missing(x)) stop("nothing fed in for x")
  if(deg < 1) stop("the spline deg must be >= 1")
  
  bdKnots = sort(bdKnots)
  intKnotsSorted = NULL
  if(!is.null(intKnots)) intKnotsSorted = sort(intKnots)
  
  knots = c(rep(bdKnots[1], deg + 1),
            intKnotsSorted,
            rep(bdKnots[2], deg + 1))
  K = length(intKnots) + deg + 1 # num of columns
  Bmat = matrix(0, length(x), K)
  for(j in 1:K) Bmat[,j] = bs(x, deg, j, knots)
  
  if(any(x == bdKnots[2])) Bmat[x == bdKnots[2], K] = 1
  
  if(intercept == F){ return(Bmat[,-1])
  } else { return(Bmat)}
}

n = 1000
x = seq(0, 1, length = n)
intKnots = seq(0, 1, length = 4+2)
intKnots = intKnots[-c(1, length(intKnots))]
B = BS(x, deg=3, intKnots = intKnots, intercept = T, bdKnots = c(0,1))
matplot(x, B, type="l", lwd = 2)
