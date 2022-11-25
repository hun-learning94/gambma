#ifndef __MLE__
#define __MLE__

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

double deviance(const double &ysize, 
                const double &scale,
                const arma::mat& y, 
                const arma::vec &eta,
                const Rcpp::String &family);

arma::vec Lpy1st(const arma::mat &y, 
                 const arma::vec &eta,
                 double ysize,
                 double yscale,
                 Rcpp::String family, 
                 Rcpp::String link);

arma::vec Wvec(const arma::mat &y, 
               const arma::vec &eta,
               const double &ysize,
               const double &yscale,
               const Rcpp::String &family, 
               const Rcpp::String &link);

arma::mat XWX(const arma::mat &y, 
              arma::mat X, 
              const arma::vec &Wvec,
              const double &ysize,
              const double &yscale,
              const Rcpp::String &family, 
              const Rcpp::String &link,
              const Rcpp::Function &crossprod);

//[[Rcpp::export]]
arma::vec etastart(const double &ysize, 
                   const double &scale,
                   const arma::mat& y,
                   const Rcpp::String &family);

//[[Rcpp::export]]
arma::vec Rglmcpp(const arma::mat &y, 
                  const arma::mat &X, 
                  const arma::vec &Etastart,
                  const bool &prevstart,
                  const arma::mat &A, 
                  const double &ysize,
                  const double &yscale,
                  const Rcpp::String &family, 
                  const Rcpp::String &link);


#endif // __MLE__



































