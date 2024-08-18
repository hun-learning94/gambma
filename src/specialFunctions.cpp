#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <RcppArmadillo.h>
#include <limits>
#include "stdLinspace.h"
#include "specialFunctions.h"
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]


double log1F1_cpp(const double &aa, const double &rr, const double &xx){
  double lb{0.0001}, ub{1.0-0.0001};
  std::vector<double> tt = std_linspace(lb, ub, 50);
  auto f1 = [&](double t) { 
    return (aa-1.0)*std::log(t) + (rr-aa-1.0)*std::log1p(-t) + xx*t - R::lbeta(rr-aa, aa); 
  };
  double at_max = *std::max_element(tt.begin(), tt.end(),
                                  [&](double a, double b){return f1(a)<f1(b);});
  double adjV = f1(at_max) - 0.5 * std::log(std::numeric_limits<double>::max());
  auto f2 = [&](double t) { 
    return std::exp((aa-1.0)*std::log(t) + (rr-aa-1.0)*std::log1p(-t) + xx*t - R::lbeta(rr-aa, aa) - adjV); 
  };
  
  // Rcpp::Rcout << adjV << "\n";
  double error;
  double Q = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f2, lb, ub, 5, 1e-14, &error);
  Q = std::log(Q) + adjV;
  return Q;
}

double log2F1_cpp(const double &bb, const double &aa, const double &rr, const double &xx){
  double lb{0.0001}, ub{1.0-0.0001};
  std::vector<double> tt = std_linspace(lb, ub, 50);
  auto f1 = [&](double t) { 
    return (aa-1.0)*std::log(t) + (rr-aa-1.0)*std::log1p(-t) - bb*std::log1p(-xx*t) - R::lbeta(rr-aa, aa); 
  };
  double at_max = *std::max_element(tt.begin(), tt.end(),
                                    [&](double a, double b){return f1(a)<f1(b);});
  double adjV = f1(at_max) - 0.5 * std::log(std::numeric_limits<double>::max());
  auto f2 = [&](double t) { 
    return std::exp((aa-1.0)*std::log(t) + (rr-aa-1.0)*std::log1p(-t) - bb*std::log1p(-xx*t) - R::lbeta(rr-aa, aa) - adjV); 
  };
  
  // Rcpp::Rcout << adjV << "\n";
  double error;
  double Q = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f2, lb, ub, 5, 1e-14, &error);
  Q = std::log(Q) + adjV;
  return Q;
}

double logPhi1_cpp(const double &aa, const double &bb, const double &rr, const double &xx, const double &yy){
  double lb{0.0001}, ub{1.0-0.0001};
  std::vector<double> tt = std_linspace(lb, ub, 50);
  auto f1 = [&](double t) { 
    return (aa-1.0)*std::log(t) + (rr-aa-1.0)*std::log1p(-t) - bb*std::log1p(-yy*t) + xx*t - R::lbeta(rr-aa, aa); 
  };
  double at_max = *std::max_element(tt.begin(), tt.end(),
                                    [&](double a, double b){return f1(a)<f1(b);});
  double adjV = f1(at_max) - 0.5 * std::log(std::numeric_limits<double>::max());
  auto f2 = [&](double t) { 
    return std::exp((aa-1.0)*std::log(t) + (rr-aa-1.0)*std::log1p(-t) - bb*std::log1p(-yy*t) + xx*t - R::lbeta(rr-aa, aa) - adjV); 
  };
  
  // Rcpp::Rcout << adjV << "\n";
  double error;
  return (std::log(boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f2, lb, ub, 5, 1e-14, &error)) + adjV);
}

double logF1_cpp(const double &aa, const double &bb, const double &bbp, const double &rr, const double &xx, const double &yy){
  double lb{0.0001}, ub{1.0-0.0001};
  std::vector<double> tt = std_linspace(lb, ub, 50);
  auto f1 = [&](double t) { 
    return (aa-1.0)*std::log(t) + (rr-aa-1.0)*std::log1p(-t) - bb*std::log1p(-xx*t) - bbp*std::log1p(-yy*t) - R::lbeta(rr-aa, aa); 
  };
  double at_max = *std::max_element(tt.begin(), tt.end(),
                                    [&](double a, double b){return f1(a)<f1(b);});
  double adjV = f1(at_max) - 0.5 * std::log(std::numeric_limits<double>::max());
  auto f2 = [&](double t) { 
    return std::exp((aa-1.0)*std::log(t) + (rr-aa-1.0)*std::log1p(-t) - bb*std::log1p(-xx*t) - bbp*std::log1p(-yy*t) - R::lbeta(rr-aa, aa) - R::lbeta(rr-aa, aa) - adjV); 
  };
  
  // Rcpp::Rcout << adjV << "\n";
  double error;
  return (std::log(boost::math::quadrature::gauss_kronrod<double, 15>::integrate(f2, lb, ub, 5, 1e-14, &error)) + adjV);
}

