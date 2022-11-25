#include "MLE.h"
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double deviance(const double &ysize, 
                const double &scale,
                const arma::mat& y, 
                const arma::vec &eta,
                const Rcpp::String &family){
  arma::vec loglik(y.n_rows);
  if(family == "binomial"){
    loglik = y % eta - ysize * arma::log1p(exp(eta));
  }
  if(family == "poisson"){
    loglik = y % eta - exp(eta);
  }
  return(-2.0 * sum(loglik));
}

arma::vec Lpy1st(const arma::mat &y, 
                 const arma::vec &eta,
                 double ysize,
                 double yscale,
                 Rcpp::String family, 
                 Rcpp::String link){
  arma::vec res(eta.n_elem);
  if(family == "binomial"){
    if(link == "logit") res = y - ysize / (1.0 + arma::exp(-eta));
  }
  if(family == "poisson"){
    if(link == "log") res = y - arma::exp(eta);
  }
  return res;
}

arma::vec Wvec(const arma::mat &y, 
               const arma::vec &eta,
               const double &ysize,
               const double &yscale,
               const Rcpp::String &family, 
               const Rcpp::String &link){
  arma::vec res(y.n_rows, fill::zeros);
  if(family == "binomial"){
    if(link == "logit") res = ysize * arma::exp(eta) / arma::pow(1.0 + arma::exp(eta), 2);
  }
  if(family == "poisson"){
    if(link == "log") res = arma::exp(eta);
  }
  return res;
}

arma::mat XWX(const arma::mat &y, 
              arma::mat X, 
              const arma::vec &Wvec,
              const double &ysize,
              const double &yscale,
              const Rcpp::String &family, 
              const Rcpp::String &link,
              const Rcpp::Function &crossprod)
{
  vec res = arma::sqrt(Wvec);
  for(unsigned i=0; i<X.n_rows; ++i){
    X.row(i) = X.row(i) * res(i);
  }
  return Rcpp::as<arma::mat>(crossprod(X));
}

arma::vec etastart(const double &ysize, 
                   const double &scale,
                   const arma::mat& y,
                   const Rcpp::String &family){
  arma::vec mu(y.n_rows);
  arma::vec eta(y.n_rows);
  if(family == "binomial"){
    mu = (y / ysize + 0.5) / (1.0/ysize + 1.0);
    eta = arma::log(1.0 / (1.0 - mu) - 1.0);
  }
  if(family == "poisson"){
    mu = y + 0.1;
    eta= arma::log(mu);
  }
  return(eta);
}

arma::vec Rglmcpp(const arma::mat &y, 
                  const arma::mat &X, 
                  const arma::vec &Etastart,
                  const bool &prevstart,
                  const arma::mat &A, 
                  const double &ysize,
                  const double &yscale,
                  const Rcpp::String &family, 
                  const Rcpp::String &link)
{
  Rcpp::Environment base("package:base");
  Rcpp::Function crossprod = base["crossprod"];
  Rcpp::Function Solve = base["solve"];
  
  arma::vec Avec = A.as_col();
  int maxiter = 25;
  int p = X.n_cols;
  int n = y.n_rows;
  vec Beta(p, fill::zeros);
  double err{0.0};
  int iter{0};
  
  if(family == "normal"){
    Beta = Rcpp::as<arma::vec>(Solve(Rcpp::as<arma::mat>(crossprod(X)), X.t() * (y - Avec)));
  } 
  else {
    arma::vec Eta(n, arma::fill::zeros);
    arma::vec wvec(n, arma::fill::zeros);
    arma::mat xwx(p, p, arma::fill::zeros);
    arma::vec delta(n, arma::fill::zeros);
    arma::vec z(n, arma::fill::zeros);
    arma::vec Beta_old = Beta;
    
    if(prevstart){ Eta = Etastart;
    } else { Eta = etastart(ysize, yscale, y, family) + Avec;}
    wvec = Wvec(y, Eta, ysize, yscale, family, link);
    xwx = XWX(y, X, wvec, ysize, yscale, family, link, crossprod);
    delta = Lpy1st(y, Eta,  ysize, yscale, family, link);
    z = wvec % (Eta - Avec) + delta;
    
    double dev{0.0}, dev_new{0.0};
    double dev_const = deviance(ysize, yscale, y, Eta, family);
    
    for(int i=0; i < maxiter; ++i){
      iter++;
      Beta = Rcpp::as<arma::vec>(Solve(xwx, X.t() * z));
      if(Beta.has_nan()) Rcpp::stop("nan produced");
      Eta = X * Beta+ A.as_col();
      dev_new = (deviance(ysize, yscale, y, Eta, family) - dev_const);
      err = abs(dev_new - dev) / (abs(dev_new) + 0.1);
      
      if(err < 1e-8){
        break;
      } else {
        wvec = Wvec(y, Eta, ysize, yscale, family, link);
        xwx = XWX(y, X, wvec, ysize, yscale, family, link, crossprod);
        delta = Lpy1st(y, Eta,  ysize, yscale, family, link);
        z = wvec % (Eta - A.as_col()) + delta;
        dev = dev_new;
        Beta_old = Beta;
      }
    }
  }
  
  return Beta;
}






































