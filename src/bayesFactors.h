#ifndef __bayesFactors__
#define __bayesFactors__

#include <RcppArmadillo.h>
#include "specialFunctions.h"
#include "sliceSampling.h"

////////////////////////////////////////////////////////////////////////////////
// Second log derivative calculation of binomial - probit link
Rcpp::NumericVector phi(const Rcpp::NumericVector &x);
Rcpp::NumericVector Phi(const Rcpp::NumericVector &x);
Rcpp::NumericVector PHI(const Rcpp::NumericVector &x);
Rcpp::NumericVector f1(const Rcpp::NumericVector &x);
Rcpp::NumericVector f2(const Rcpp::NumericVector &x);
template <typename T>
Rcpp::NumericVector to_NumVec(const T &x);
arma::vec lpy2nd_probit(double glmWeight, const arma::vec &x, const arma::vec &y);

////////////////////////////////////////////////////////////////////////////////
// calculation of J(eta), J(alpha), J(beta), r2, Qm
// J(Eta) = - d^2 l(eta) / d eta^2
arma::vec JEta(const double &glmWeight, 
               const arma::vec &y, 
               const arma::vec& EtaHat, 
               const unsigned& familyLink);

// J(Beta) = X_cen^T J(Eta) X_cen (X_cen = (I-P_1)X, P_1 = ortho. proj. onto 1 under info. inp. defined by J(Eta))
// rootJBeta = sqrt(J)*(I-P_1)X
arma::mat rootJBeta(const arma::vec &BetaHat, 
                    const arma::mat &Xm, 
                    const arma::vec &JEtaHat, 
                    const unsigned& familyLink);

// J(alpha)
double JAlpha(const arma::vec &JEta);

// r2 = ||betahat * X ||^2 / ||y-meany||^2
// Qm = beta^T J(beta) beta = sum((Etahat[i]-alphahat)^*J(EtaHat)[i,i])
double getR2QM(const unsigned &familyLink,
               const arma::vec &y, 
               const arma::vec &offset,
               const double& AlphaHat,
               const arma::vec &EtaHat,
               const arma::vec &JEtaHat
);

////////////////////////////////////////////////////////////////////////////////
// initial starting point of IRLS algorithm
//[[Rcpp::export]]
arma::vec etastart(const double &glmWeight, 
                   const arma::mat& y,
                   const unsigned& familyLink);

////////////////////////////////////////////////////////////////////////////////
// Bayes factors based on the null model
// we disregard any component in the model evidence that
// does not depend on the model in consideration

// GLM
void glmLoglik(double& glmLoglik_out,
               const double& glmWeight,
               const arma::vec& y, 
               const arma::mat& X,
               const arma::vec& offset,
               const arma::vec& EtaHat,
               const unsigned& familyLink);
void Lpy(double& Lpy_out,
         double& comp1,
         double& comp2,
         const double& glmWeight,
         const arma::vec& y, 
         const arma::mat& X,
         const arma::vec& offset,
         const arma::vec& EtaHat,
         const double &JAlphaHat,
         const double& r2Qm,
         const unsigned& familyLink,
         const unsigned& gprior,
         double aa, double bb, double ss, double gg);

////////////////////////////////////////////////////////////////////////////////
// gPosteriors
void gPosteriors(double& g_,
                 double& v_, 
                 double& t_, 
                 double& q_, 
                 double& l_,
                 double& m_,
                 const double &n,
                 const double &p,
                 const double &r2Qm,
                 const unsigned& familyLink,
                 const unsigned& gprior,
                 double aa, double bb, double ss, double gg);

////////////////////////////////////////////////////////////////////////////////
// posterior samples of alpha, beta
double AlphaPost(double AlphaHat, double JAlphaHat, double phi);

void basechol(arma::mat &toBeCholed,
              const Rcpp::Function &nearPDres, 
              const Rcpp::Function &Rbasechol);

arma::vec BetaPost(const arma::vec &BetaHat, 
                   const arma::mat &rootJBetaHat, 
                   const double& g, 
                   const double& phi,
                   const Rcpp::Function &nearPDres,
                   const Rcpp::Function &Rbasechol,
                   const Rcpp::Function &crossproduct);

arma::vec BetaPost_z(const arma::vec &BetaHat, 
                     const arma::mat &rootJBetaHat,
                     const double &g, 
                     const double &phi,
                     const arma::ivec &z,
                     const Rcpp::Function &nearPDres,
                     const Rcpp::Function &Rbasechol,
                     const Rcpp::Function &crossproduct);

#endif // __bayesFactors__



































