#ifndef __CRAD__
#define __CRAD__

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::vec armapmax(const arma::vec &x,
                   const double &y);

arma::mat CRADNS_1d_cpp(const arma::vec& x, const arma::vec &knot, 
                        bool knotalive, double bdmargin);

arma::mat CRAD_1d_cpp(const arma::vec& x, const arma::vec &knot, 
                      bool knotalive);

//[[Rcpp::export(.CRAD)]]
Rcpp::List CRAD_cpp(const arma::mat &X, 
                    const arma::mat &X_lin, 
                    const arma::vec &knots, 
                    const arma::uvec &knotsidx,
                    bool NS,
                    double bdmargin);

//[[Rcpp::export(.CRAD_test)]]
arma::mat CRAD_test_cpp(const arma::mat &testX, 
                        const arma::mat &X_lin, 
                        const Rcpp::List &CRADlist);

#endif // __CRAD__