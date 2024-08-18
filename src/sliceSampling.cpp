#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <RcppArmadillo.h>
#include "stdLinspace.h"
#include "sliceSampling.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

//[[Rcpp::export]]
double rtbeta_cpp(const double &alpha, const double &beta, const double &a, const double &b){
  double x = R::runif(0.0, 1.0);
  double Fa = R::pbeta(a, alpha, beta, true, false);
  double Fb = R::pbeta(b, alpha, beta, true, false);
  x = (1.0-x)*Fa + x*Fb;
  return R::qbeta(x, alpha, beta, true, false);
}

//[[Rcpp::export]]
double rtgamma_cpp(const double &shape, const double &rate, const double &a, const double &b){
  // mean = shape * scale = shape / rate
  double x = R::runif(0.0, 1.0);
  double scale = 1.0/rate;
  double Fa = R::pgamma(a, shape, scale, true, false);
  double Fb = R::pgamma(b, shape, scale, true, false);
  x = (1.0-x)*Fa + x*Fb;
  return R::qgamma(x, shape, scale, true, false);
}

void rCH_void(double& v_, double& l_, 
              const double& aa, const double& bb, const double& ss){
  l_ = std::exp(std::log(R::runif(0.0, 1.0)) - ss*v_);
  double lb, ub;
  if(ss>0){
    lb = 0.0;
    ub = std::min(-std::log(l_)/ss, 1.0);
  } else {
    ub = std::max(-std::log(l_)/ss, 0.0);
    lb = 0.0;
  }
  v_ = rtbeta_cpp(aa, bb, lb, ub);
}

//[[Rcpp::export]]
arma::vec rCH_cpp(int nsamp, int burnin,
                  const double& aa, const double& bb, const double& ss){
  if(!(aa>0)) Rcpp::stop("Must be a > 0");
  if(!(bb>0)) Rcpp::stop("Must be b > 0");
  
  int S = nsamp + burnin;
  arma::vec RES(S);
  double v_{0.1}, l_{0.1};
  for(int i{0};i<S;i++){
    rCH_void(v_, l_, aa, bb, ss);
    RES(i) = v_;
  }
  return RES.subvec(burnin, S-1);
}

void rGH_void(double& v_, double& t_, double& l_,
              const double& aa, const double& bb, const double& xx, const double& zz){
  l_ = std::exp(std::log(R::runif(0.0, 1.0)) - xx*v_*t_);
  double lb, ub;
  if(xx>0){
    ub = -std::log(l_)/(xx*v_);
    lb = 0.0;
  } else {
    ub = std::numeric_limits<double>::infinity();
    lb = std::max(0.0, -std::log(l_)/(xx*v_));
  }
  t_ = rtgamma_cpp(zz, 1.0, lb, ub);
  if(xx>0){
    ub = std::min(-std::log(l_)/(xx*t_), 1.0);
    lb = 0.0;
  } else {
    ub = 1.0;
    lb = std::max(-std::log(l_)/(xx*t_), 0.0);
  }
  v_ = rtbeta_cpp(aa, bb, lb, ub);
} 

//[[Rcpp::export]]
arma::vec rGH_cpp(int nsamp, int burnin,
                  const double& aa, const double& bb, const double& xx, const double& zz){
  if(!(aa>0)) Rcpp::stop("Must be a > 0");
  if(!(bb>0)) Rcpp::stop("Must be b > 0");
  if(!(zz>0)) Rcpp::stop("Must be z > 0");
  
  int S = nsamp + burnin;
  arma::vec RES(S);
  double v_{0.1}, l_{0.1}, t_{0.1};
  for(int i{0};i<S;i++){
    rGH_void(v_, t_, l_, aa, bb, xx, zz);
    RES(i) = v_;
  }
  return RES.subvec(burnin, S-1);
}

void rCCH_void(double& v_, double& t_, double& l_,double& m_,
               const double& aa, const double& bb, const double& zz, const double& ss, const double& xx){
  l_ = std::exp(std::log(R::runif(0.0, 1.0)) - xx*v_*t_);
  m_ = std::exp(std::log(R::runif(0.0, 1.0)) - ss*v_);
  double lb, ub;
  if(xx>0){
    ub = -std::log(l_)/(xx*v_);
    lb = 0.0;
  } else {
    ub = std::numeric_limits<double>::infinity();
    lb = std::max(0.0, -std::log(l_)/(xx*v_));
  }
  t_ = rtgamma_cpp(zz, 1.0, lb, ub);
  if((xx>0) & (ss >0)){
    ub = std::min(1.0, std::min(-std::log(l_)/(xx*t_), -std::log(m_)/ss));
    lb = 0.0;
  } else if ((xx>0) & (ss<=0)){
    lb = std::max(0.0, -std::log(m_)/ss);
    ub = std::min(1.0, -std::log(l_)/(xx*t_));
  } else if ((xx<=0) & (ss>0)){
    ub = std::min(1.0, -std::log(m_)/ss);
    lb = std::max(0.0, -std::log(l_)/(xx*t_));
  } else {
    lb = std::max(0.0, std::max(-std::log(l_)/(xx*t_), -std::log(m_)/ss));
    ub = 1.0;
  }
  v_ = rtbeta_cpp(aa, bb, lb, ub);
} 

//[[Rcpp::export]]
arma::vec rtCCH_cpp(int nsamp, int burnin,
                    const double& aa, const double& bb, const double& zz, const double& ss, const double &nu, const double& theta){
  if(!(aa>0)) Rcpp::stop("Must be a > 0");
  if(!(bb>0)) Rcpp::stop("Must be b > 0");
  if(!(theta>0)) Rcpp::stop("Must be theta > 0");
  if(!(zz>0)) Rcpp::stop("Must be z > 0");
  
  int S = nsamp + burnin;
  arma::vec RES(S);
  double v_{0.1}, l_{0.1}, t_{0.1},m_{0.1};
  double tilde_s{ss/nu}, xx{1.0/theta - 1.0};
  for(int i{0};i<S;i++){
    rCCH_void(v_, t_, l_, m_, aa, bb, zz, tilde_s, xx);
    RES(i) = v_/nu;
  }
  return RES.subvec(burnin, S-1);
}

void rAPL_void(double& v_, double& t_, double& q_, double& l_,double& m_,
               const double& aa, const double& bb, const double& zz, const double& ww, const double& xx, const double& yy){
  l_ = std::exp(std::log(R::runif(0.0, 1.0)) - xx*v_*t_);
  m_ = std::exp(std::log(R::runif(0.0, 1.0)) - yy*v_*q_);
  double lb, ub;
  if(xx>0){
    ub = -std::log(l_)/(xx*v_);
    lb = 0.0;
  } else {
    ub = std::numeric_limits<double>::infinity();
    lb = std::max(0.0, -std::log(l_)/(xx*v_));
  }
  t_ = rtgamma_cpp(zz, 1.0, lb, ub);
  
  if(xx>0){
    ub = -std::log(m_)/(yy*v_);
    lb = 0.0;
  } else {
    ub = std::numeric_limits<double>::infinity();
    lb = std::max(0.0, -std::log(m_)/(yy*v_));
  }
  q_ = rtgamma_cpp(ww, 1.0, lb, ub);
  
  if((xx>0) & (yy >0)){
    ub = std::min(1.0, std::min(-std::log(l_)/(xx*t_), -std::log(m_)/(yy*q_)));
    lb = 0.0;
  } else if ((xx<=0) & (yy<=0)){
    lb = std::max(0.0, std::max(-std::log(l_)/(xx*t_), -std::log(m_)/(yy*q_)));
    ub = 1.0;
  } else if ((xx>0) & (yy<=0)){
    lb = std::max(0.0, -std::log(m_)/(yy*q_));
    ub = std::min(1.0, -std::log(l_)/(xx*t_));
  } else {
    ub = std::min(1.0, -std::log(m_)/(yy*q_));
    lb = std::max(0.0, -std::log(l_)/(xx*t_));
  }
  v_ = rtbeta_cpp(aa, bb, lb, ub);
} 

//[[Rcpp::export]]
arma::vec rAPL_cpp(int nsamp, int burnin,
                   const double& aa, const double& bb, const double& zz, const double& ww, const double& xx, const double& yy){
  if(!(aa>0)) Rcpp::stop("Must be a > 0");
  if(!(bb>0)) Rcpp::stop("Must be b > 0");
  if(!(xx>0)) Rcpp::stop("Must be x > -1");
  if(!(yy>0)) Rcpp::stop("Must be y > -1");
  if(!(zz>0)) Rcpp::stop("Must be z > 0");
  if(!(ww>0)) Rcpp::stop("Must be w > 0");
  
  int S = nsamp + burnin;
  arma::vec RES(S);
  double v_{0.1}, l_{0.1}, t_{0.1},m_{0.1}, q_{0.1};
  for(int i{0};i<S;i++){
    rAPL_void(v_, t_, q_, l_, m_, aa,bb,zz,ww,xx,yy);
    RES(i) = v_;
  }
  return RES.subvec(burnin, S-1);
}


