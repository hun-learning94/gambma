# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.CRAD <- function(X, X_lin, knots, knotsidx, NS, bdmargin) {
    .Call(`_gambms_CRAD_cpp`, X, X_lin, knots, knotsidx, NS, bdmargin)
}

.CRAD_test <- function(testX, X_lin, CRADlist) {
    .Call(`_gambms_CRAD_test_cpp`, testX, X_lin, CRADlist)
}

.CRAD2 <- function(X, X_lin, knots, knotsidx, NS, bdmargin) {
    .Call(`_gambms_CRAD_cpp2`, X, X_lin, knots, knotsidx, NS, bdmargin)
}

.CRAD_test2 <- function(testX, X_lin, CRADlist) {
    .Call(`_gambms_CRAD_test_cpp2`, testX, X_lin, CRADlist)
}

.gambmsEVEN <- function(y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, enumerate, numMCcandidate, MCiter, MCMCiter, Rglm, nearPDres, storeFit, forceLin, linProb, printiter) {
    .Call(`_gambms_gambmsEVEN`, y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, enumerate, numMCcandidate, MCiter, MCMCiter, Rglm, nearPDres, storeFit, forceLin, linProb, printiter)
}

.gambmsFREE <- function(y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, initS, MCMCiter, thin, bir_p, dea_p, nu, Rglm, nearPDres, storeFit, forceLin, linProb, printiter) {
    .Call(`_gambms_gambmsFREE`, y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, initS, MCMCiter, thin, bir_p, dea_p, nu, Rglm, nearPDres, storeFit, forceLin, linProb, printiter)
}

.gambmsFREE2 <- function(y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, initS, MCMCiter, thin, bir_p, dea_p, nu, Rglm, nearPDres, storeFit, forceLin, linProb, printiter) {
    .Call(`_gambms_gambmsFREE2`, y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, initS, MCMCiter, thin, bir_p, dea_p, nu, Rglm, nearPDres, storeFit, forceLin, linProb, printiter)
}

.gambmsVS <- function(y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, MCMCiter, Rglm, nearPDres, getmeMAP, storeFit, forceLin, linProb, printiter) {
    .Call(`_gambms_gambmsVS`, y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, MCMCiter, Rglm, nearPDres, getmeMAP, storeFit, forceLin, linProb, printiter)
}

.gambmsVS2 <- function(y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, MCMCiter, Rglm, nearPDres, getmeMAP, storeFit, forceLin, linProb, printiter) {
    .Call(`_gambms_gambmsVS2`, y, glmWeight, X, X_pr, XLin, offset, maxk, Lambda, familyLink, gprior, aa, bb, ss, gg, MCMCiter, Rglm, nearPDres, getmeMAP, storeFit, forceLin, linProb, printiter)
}

.log1F1 <- function(aa, rr, xx) {
    .Call(`_gambms_log1F1_cpp`, aa, rr, xx)
}

.log2F1 <- function(bb, aa, rr, xx) {
    .Call(`_gambms_log2F1_cpp`, bb, aa, rr, xx)
}

.logPhi1 <- function(aa, bb, rr, xx, yy) {
    .Call(`_gambms_logPhi1_cpp`, aa, bb, rr, xx, yy)
}

.logF1 <- function(aa, bb, bbp, rr, xx, yy) {
    .Call(`_gambms_logF1_cpp`, aa, bb, bbp, rr, xx, yy)
}

