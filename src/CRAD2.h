#ifndef __CRAD2__
#define __CRAD2__

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

arma::vec armapmax2(const arma::vec &x,
                   const double &y);

arma::mat CRADNS_1d_cpp2(const arma::vec& x, const arma::vec &knot, 
                        bool knotalive, double bdmargin);

arma::mat CRAD_1d_cpp2(const arma::vec& x, const arma::vec &knot, 
                      bool knotalive);

Rcpp::List CRAD_cpp2(const arma::mat &X, 
                    const arma::mat &X_lin, 
                    const arma::vec &knots, 
                    const arma::uvec &knotsidx,
                    bool NS,
                    double bdmargin);

arma::mat CRAD_test_cpp2(const arma::mat &testX, 
                        const arma::mat &X_lin, 
                        const Rcpp::List &CRADlist);

#endif // __CRAD2__