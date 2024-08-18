#ifndef __specialFunctions__
#define __specialFunctions__
// Gauss-Kronrod quadrature powered by C++ Boost boost::math::quadrature::gauss_kronrod
// https://www.boost.org/doc/libs/1_69_0/libs/math/doc/html/math_toolkit/gauss_kronrod.html
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//[[Rcpp::export(.log1F1)]]
double log1F1_cpp(const double &aa, const double &rr, const double &xx);
  
  //[[Rcpp::export(.log2F1)]]
double log2F1_cpp(const double &bb, const double &aa, const double &rr, const double &xx);

//[[Rcpp::export(.logPhi1)]]
double logPhi1_cpp(const double &aa, const double &bb, const double &rr, const double &xx, const double &yy);

//[[Rcpp::export(.logF1)]]
double logF1_cpp(const double &aa, const double &bb, const double &bbp, const double &rr, const double &xx, const double &yy);

#endif // __specialFunctions__