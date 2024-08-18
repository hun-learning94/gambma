#ifndef __gambmsVS2__
#define __gambmsVS2__
#include <RcppArmadillo.h>

//[[Rcpp::export(.gambmsVS2)]]
Rcpp::List gambmsVS2(const arma::vec &y,
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
                    const unsigned &MCMCiter,
                    const Rcpp::Function &Rglm,
                    const Rcpp::Function &nearPDres,
                    bool getmeMAP,
                    const bool& storeFit,
                    const bool& forceLin,
                    double& linProb,
                    unsigned printiter);

#endif // __gambmsVS2__