#ifndef __gambmsFREE__
#define __gambmsFREE__
#include <RcppArmadillo.h>

//[[Rcpp::export(.gambmsFREE)]]
Rcpp::List gambmsFREE(const arma::vec &y,
                      const double &glmWeight,
                      const arma::mat &X,
                      const arma::mat &X_pr,
                      const arma::mat &XLin,
                      const arma::vec &offset,
                      const arma::uvec maxk,
                      const arma::vec &Lambda,
                      const unsigned& familyLink,
                      const unsigned& gprior,
                      const double &aa, 
                      const double &bb, 
                      const double &ss, 
                      const double &gg,
                      const int &initS,
                      const int &MCMCiter,
                      const int &thin,
                      double bir_p, double dea_p, double nu,
                      const Rcpp::Function &Rglm,
                      const Rcpp::Function &nearPDres,
                      const bool& storeFit,
                      const bool& forceLin,
                      double & linProb,
                      unsigned printiter);

#endif // __gambmsFREE__