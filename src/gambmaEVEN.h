#ifndef __gambmaEVEN__
#define __gambmaEVEN__
#include <RcppArmadillo.h>

Rcpp::List gambmaEVEN(const arma::vec &y,
                      const double &glmWeight,
                      const arma::mat &X,
                      const arma::mat &X_pr,
                      const arma::mat &XLin,
                      const arma::vec &offset,
                      const arma::uvec &maxk,
                      const arma::vec &Lambda,
                      const unsigned& familyLink,
                      const unsigned& gprior,
                      const double &aa, 
                      const double &bb, 
                      const double &ss, 
                      const double &gg,
                      const bool& enumerate,
                      unsigned numMCcandidate,
                      const unsigned &MCiter,
                      const unsigned &MCMCiter,
                      const Rcpp::Function &Rglm,
                      const Rcpp::Function &nearPDres,
                      const bool& storeFit,
                      unsigned printiter);

#endif // __gambmaEVEN__


















