#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include "specialFunctions.h"
#include "sliceSampling.h"
#include "bayesFactors.h"

////////////////////////////////////////////////////////////////////////////////
// Second log derivative calculation of binomial - probit link
Rcpp::NumericVector phi(const Rcpp::NumericVector &x){
  return -x * Rcpp::dnorm(x);
}
Rcpp::NumericVector Phi(const Rcpp::NumericVector &x){
  return Rcpp::dnorm(x, 0.0, 1.0, true);
}
Rcpp::NumericVector PHI(const Rcpp::NumericVector &x){
  return Rcpp::pnorm(x, 0.0, 1.0, true, true);
}
Rcpp::NumericVector f1(const Rcpp::NumericVector &x){
  Rcpp::NumericVector tmp = Rcpp::exp(Phi(x) - PHI(x));
  return tmp * (tmp + x);
}
Rcpp::NumericVector f2(const Rcpp::NumericVector &x){
  Rcpp::NumericVector tmp = Rcpp::exp(Phi(x) - PHI(-x));
  return tmp * (tmp - x);
}

template <typename T>
Rcpp::NumericVector to_NumVec(const T &x){
  return Rcpp::NumericVector(x.begin(), x.end());
}

arma::vec lpy2nd_probit(double glmWeight, const arma::vec &x, const arma::vec &y)
{
  Rcpp::NumericVector tmp;
  if(glmWeight == 1.0){
    Rcpp::LogicalVector onezero{(to_NumVec(y) == 1.0)};
    tmp = Rcpp::ifelse(onezero, f1(to_NumVec(x)), f2(to_NumVec(x)));
  } else {
    tmp = to_NumVec(y) * f1(to_NumVec(x)) + (glmWeight - to_NumVec(y)) * f2(to_NumVec(x));
  }
  return Rcpp::as<arma::vec>(tmp);
}

////////////////////////////////////////////////////////////////////////////////
// calculation of J(eta), J(alpha), J(beta), r2, Qm
// J(Eta) = - d^2 l(eta) / d eta^2
arma::vec JEta(const double &glmWeight, 
               const arma::vec &y, 
               const arma::vec& EtaHat, 
               const unsigned& familyLink){
  arma::vec res;
  if(familyLink == 11){ // gaussian-identity
    res = arma::vec(y.n_rows, arma::fill::ones);
  } else if(familyLink == 21){ // binomial-logit
    res = glmWeight * arma::exp(EtaHat) / arma::pow(1.0 + arma::exp(EtaHat), 2);
  } else if(familyLink == 22){ // binomial-probit
    res = lpy2nd_probit(glmWeight, y, EtaHat);
  } else if(familyLink == 31){ // poisson-log
    res = arma::exp(EtaHat);
  }
  if(any(res < 0)){
    // Rcpp::Rcout << res.elem(find(res <0)) << '\n';
    Rcpp::stop("second derivative of the log likelihood is non-negative");
  }
  return res;
}

// J(Beta) = X_cen^T J(Eta) X_cen (X_cen = (I-P_1)X, P_1 = ortho. proj. onto 1 under info. inp. defined by J(Eta))
// rootJBeta = sqrt(J)*(I-P_1)X 
arma::mat rootJBeta(const arma::vec &BetaHat, 
                    const arma::mat &Xm, 
                    const arma::vec &JEtaHat, 
                    const unsigned& familyLink)
{
  arma::mat XmOut = Xm.cols(1, Xm.n_cols - 1);
  arma::mat JBetaHat(BetaHat.n_elem, BetaHat.n_elem, arma::fill::zeros);
  if(familyLink == 11){ // Gaussian-Identity
    XmOut.each_row() -= arma::mean(XmOut, 0);
  } else { // weighted average centering
    arma::vec weights = JEtaHat / arma::accu(JEtaHat);
    for(unsigned i{0}; i < XmOut.n_cols; i++){
      XmOut.col(i) -= arma::dot(XmOut.col(i), weights);
    }
    for(unsigned i=0; i < XmOut.n_rows; i++){
      XmOut.row(i) *= sqrt(JEtaHat(i));
    }
  }
  return XmOut;
}

// J(alpha)
double JAlpha(const arma::vec &JEta){
  return arma::accu(JEta);
}

// r2 = ||betahat * X ||^2 / ||y-meany||^2
// Qm = beta^T J(beta) beta = sum((Etahat[i]-alphahat)^*J(EtaHat)[i,i])
double getR2QM(const unsigned &familyLink,
               const arma::vec &y, 
               const arma::vec &offset,
               const double& AlphaHat,
               const arma::vec &EtaHat,
               const arma::vec &JEtaHat
               ){
  double tmp;
  if(familyLink == 11){
    tmp = std::pow(arma::norm(EtaHat - AlphaHat, 2) / arma::norm(y - arma::mean(y), 2), 2.0);
  } else {
    tmp = arma::dot(arma::pow(EtaHat - offset - AlphaHat, 2), JEtaHat);
  }
  // if(std::is_nan(tmp)){
  //   Rcpp::Rcout << AlphaHat << "\n";
  // }
  return tmp;
}

////////////////////////////////////////////////////////////////////////////////
// initial starting point of IRLS algorithm
arma::vec etastart(const double &glmWeight, 
                   const arma::mat& y,
                   const unsigned& familyLink){
  arma::vec mu(y.n_rows);
  arma::vec eta(y.n_rows);
  if((familyLink / 10)== 2){ // binomials
    mu = (y / glmWeight + 0.5) / (1.0/glmWeight + 1.0);
    eta = arma::log(1.0 / (1.0 - mu) - 1.0);
  } else if((familyLink/10)==3){ // poissons
    mu = y + 0.1;
    eta= arma::log(mu);
  }
  return(eta.zeros());
}


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
               const unsigned& familyLink){
  arma::vec EtaPlusOffset, glmLoglikVec;
  EtaPlusOffset = EtaHat + offset;
  if((familyLink / 10)== 2){ // binomials
    glmLoglikVec = y % EtaPlusOffset - glmWeight * arma::log1p(arma::exp(EtaPlusOffset));
  } else if((familyLink/10)==3){ // poissons
    glmLoglikVec = y % EtaPlusOffset - arma::exp(EtaPlusOffset);
  }
  // Rcpp::Rcout << EtaHat.t() << "\n";
  // Rcpp::Rcout <<  arma::log1p(arma::exp(EtaPlusOffset)).t() << "\n";
  // Rcpp::Rcout << glmLoglik_out << "\n";
  glmLoglik_out = arma::sum(glmLoglikVec);
  // if(std::isnan(glmLoglik_out) != 0){ 
  //   Rcpp::Rcout << "nan occured \n"; 
  //   Rcpp::Rcout << EtaHat.t() << "\n";
  //   Rcpp::Rcout <<  arma::log1p(arma::exp(EtaPlusOffset)).t() << "\n";
  //   // Rcpp::stop("\n");
  // }
}

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
         double aa, double bb, double ss, double gg){
  double n = y.n_elem; 
  double p = X.n_cols;
  double rr, nu, kappa;
  
  // comp1: related to model fit (r2 for Gaussian)
  // comp2: related to penalty, differ by gprior
  
  if((familyLink / 10) > 1){ // glm
    glmLoglik(comp1, 
              glmWeight, y, X, offset, EtaHat, familyLink);
    comp1 -= 0.5 * std::sqrt(JAlphaHat);
    if(gprior == 0){ // unit information (g=n)
      gg = n;
      comp2 = -r2Qm/(2.0+2.0*gg) - 0.5*p * std::log(gg+1.0);
    } else if(gprior== 1){ // Hyper-g
      aa = 1.0; bb = 2.0; nu = 1.0;
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + log1F1_cpp((aa+p)/2.0, (aa+bb+p)/2.0, -r2Qm/2.0);
    } else if(gprior == 2){ // Uniform
      aa = 2.0; bb = 2.0; nu = 1.0;
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + log1F1_cpp((aa+p)/2.0, (aa+bb+p)/2.0, -r2Qm/2.0);
    } else if(gprior == 3){ // Hyper-g/n
      aa = 1.0; bb = 2.0; rr = 1.5; ss =.0; nu = 1.0; kappa = 1.0/n;
      comp2 = - 0.5*p*std::log(nu) 
        - r2Qm/(2.0*nu) 
        + R::lbeta((aa+p)/2.0, bb/2.0) 
        + logPhi1_cpp(0.5*bb, rr, 0.5*(aa+bb+p), r2Qm/(2.0*nu), 1.0-kappa);
    } else if(gprior == 4){ // Beta-prime
      aa = 0.5; bb = n-p-1.5; nu = 1.0;
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0)
        + log1F1_cpp((aa+p)/2.0, (aa+bb+p)/2.0, -r2Qm/2.0);
    } else if(gprior == 5){ // ZS-adapted
      aa = 1.0; bb = 2.0; ss = n+3.0; nu = 1.0;
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + log1F1_cpp((aa+p)/2.0, (aa+bb+p)/2.0, -(ss+r2Qm)/2.0);
    } else if(gprior == 6){ // Robust
      aa = 1.0; bb = 2.0; rr = 1.5; ss =.0; nu = (n+1.0)/(p+1.0); kappa = 1.0;
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        - 0.5*p*std::log(nu) 
        + log1F1_cpp((aa+p)/2.0, (aa+bb+p)/2.0, -r2Qm/(2.0*nu)); 
    } else if(gprior == 7){ // Intrinsic
      aa = 1.0; bb = 1.0; rr = 1.0; ss =.0; nu = (n+p+1.0)/(p+1.0); kappa = (n+p+1.0)/n;
      comp2 = - 0.5*p*std::log(nu) 
        - r2Qm/(2.0*nu) 
        + R::lbeta((aa+p)/2.0, bb/2.0) 
        + logPhi1_cpp(0.5*bb, rr, 0.5*(aa+bb+p), r2Qm/(2.0*nu), 1.0-kappa);
    } else if(gprior == 8){ // g = user input
      comp2 = -r2Qm/(2.0+2.0*gg)
        - 0.5*p * std::log(gg+1.0);
    } else if(gprior == 9){ // CH(a,b,s)
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0)
        + log1F1_cpp((aa+p)/2.0, (aa+bb+p)/2.0, -(ss+r2Qm)/2.0);
    }
    // Rcpp::Rcout << comp2 << "\n";
    Lpy_out = comp1 + comp2;
    
  } else if (familyLink == 11) { // Gaussian
    comp1 = r2Qm;
    if(gprior == 0){ // unit information (g=n)
      gg = n;
      comp2 = 0.5*(n-p-1.0)*std::log1p(gg) 
        - 0.5*(n-1.0)*std::log1p(gg*(1.0-r2Qm));
    } else if(gprior== 1){ // Hyper-g
      aa = 1.0; bb = 2.0; 
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + log2F1_cpp(0.5*(n-1), bb/2.0, (aa+bb+p)/2.0, r2Qm);
    } else if(gprior == 2){ // Uniform
      aa = 2.0; bb = 2.0; 
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + log2F1_cpp(0.5*(n-1), bb/2.0, (aa+bb+p)/2.0, r2Qm);
    } else if(gprior == 3){ // Hyper-g/n
      aa = 1.0; bb = 2.0; rr = 1.5; nu = 1.0; kappa = 1.0/n;
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + logF1_cpp((aa+p)/2.0, (aa+bb+p+1.0-n-2.0*rr)/2.0, (n-1.0)/2.0, (aa+bb+p)/2.0, 1.0-kappa, 1.0-kappa-r2Qm*kappa/((1.0-r2Qm)*nu)) 
        + (aa+p-2.0*rr)*std::log(kappa)/2.0
        - 0.5*p*std::log(nu) 
        - 0.5*(n-1.0)*std::log1p(-r2Qm) 
        - log2F1_cpp(rr, 0.5*bb, 0.5*(aa+bb), 1.0 - kappa);
    } else if(gprior == 4){ // Beta-prime
      aa = 0.5; bb = n-p-1.5; 
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + 0.5*(p+1.5-n)*std::log1p(-r2Qm);
      // Rcpp::Rcout << r2Qm << "\n";
      // Rcpp::Rcout << R::lbeta((aa+p)/2.0, bb/2.0) << "\n";
      // Rcpp::Rcout << 0.5*(p+1.5-n)*std::log1p(-r2Qm) << "\n";
    } else if(gprior == 5){ // ZS-adapted
      aa = 1.0; bb = 2.0; ss = n+3.0; 
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + logPhi1_cpp(0.5*bb, 0.5*(n-1.0), 0.5*(aa+bb+p), 0.5*ss, r2Qm);
    } else if(gprior == 6){ // Robust
      aa = 1.0; bb = 2.0; rr = 1.5; nu = (n+1.0)/(p+1.0); kappa = 1.0;
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + log2F1_cpp(0.5*(n-1.0), 0.5*bb, 0.5*(aa+bb+p), r2Qm / (nu - (nu-1.0)*r2Qm)) 
        - 0.5*p*std::log(nu) 
        - 0.5*(n-1.0)*std::log(1.0-(1.0-1.0/nu)*r2Qm);
    } else if(gprior == 7){ // Intrinsic
      aa = 1.0; bb = 1.0; rr = 1.0; nu = (n+p+1.0)/(p+1.0); kappa = (n+p+1.0)/n;
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
        + logF1_cpp((aa+p)/2.0, (aa+bb+p+1.0-n-2.0*rr)/2.0, (n-1.0)/2.0, (aa+bb+p)/2.0, 1.0-kappa, 1.0-kappa-r2Qm*kappa/((1.0-r2Qm)*nu)) 
        + (aa+p-2.0*rr)*std::log(kappa)/2.0
        - 0.5*p*std::log(nu) 
        - 0.5*(n-1.0)*std::log1p(-r2Qm) 
        - log2F1_cpp(rr, 0.5*bb, 0.5*(aa+bb), 1.0 - kappa);
    } else if(gprior == 8){ // g = user input
      comp2 = 0.5*(n-p-1.0)*std::log1p(gg) 
      - 0.5*(n-1.0)*std::log1p(gg*(1.0-r2Qm));
    } else if(gprior == 9){ // CH(a,b,s)
      comp2 = R::lbeta((aa+p)/2.0, bb/2.0) 
      + logPhi1_cpp(0.5*bb, 0.5*(n-1.0), 0.5*(aa+bb+p), 0.5*ss, r2Qm);
    }
    Lpy_out = comp2;
    
  } else {
    Rcpp::stop("Lpy says, What da fuck?");
  }
  if(std::isnan(Lpy_out) != 0){
    Rcpp::Rcout << "nan occured \n";
    // Rcpp::Rcout << "comp1 : " << comp1 << '\n';
    // Rcpp::Rcout << "comp2 : " << comp2 << '\n';
    Rcpp::stop("\n");
  }
  
}

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
                 double aa, double bb, double ss, double gg){
  // g_: the actual g in the g-prior we want to sample
  // v_: nu*u_ where u_ = 1/(g_+1)
  // t_, q_, l_, m_ : auxiliary variables
  double rr{0.0}, nu{1.0}, kappa{1.0};
  double aaPos, bbPos, ssPos, xxPos, zzPos, wwPos, yyPos;
  if((familyLink / 10) > 1){ // glm
    if(gprior == 0){ // unit information (g=n)
      gg = n; g_ = gg;
      // Rcpp::Rcout << gg << "\n";
    } else if(gprior== 1){ // Hyper-g
      aa = 1.0; bb = 2.0; nu = 1.0; ss = .0;
      aaPos = 0.5*(aa+p); ssPos = (ss+r2Qm)/(2.0*nu);
      v_ = rtgamma_cpp(aaPos, ssPos, 0.0, 1.0);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 2){ // Uniform
      aa = 2.0; bb = 2.0; nu = 1.0;ss = .0;
      aaPos = 0.5*(aa+p); ssPos = (ss+r2Qm)/(2.0*nu);
      v_ = rtgamma_cpp(aaPos, ssPos, 0.0, 1.0);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 3){ // Hyper-g/n
      aa = 1.0; bb = 2.0; rr = 1.5; ss=0.0; nu = 1.0; kappa = 1.0/n;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; zzPos= rr; 
      ssPos = (ss+r2Qm)/(2.0*nu); xxPos = 1.0/kappa - 1.0; 
      rCCH_void(v_, t_, l_, m_, aaPos, bbPos, zzPos, ssPos, xxPos);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 4){ // Beta-prime
      aa = 0.5; bb = n-p-1.5; ss = 0.0; nu = 1.0;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; ssPos = (ss+r2Qm)/(2.0*nu);
      rCH_void(v_, l_, aaPos, bbPos, ssPos);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 5){ // ZS-adapted
      aa = 1.0; bb = 2.0; ss = n+3.0; nu = 1.0;
      aaPos = 0.5*(aa+p); ssPos = (ss+r2Qm)/(2.0*nu);
      v_ = rtgamma_cpp(aaPos, ssPos, 0.0, 1.0);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 6){ // Robust
      aa = 1.0; bb = 2.0; rr = 1.5; nu = (n+1.0)/(p+1.0); kappa = 1.0; ss = .0;
      aaPos = 0.5*(aa+p); ssPos = (ss+r2Qm)/(2.0*nu);
      v_ = rtgamma_cpp(aaPos, ssPos, 0.0, 1.0);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 7){ // Intrinsic
      aa = 1.0; bb = 1.0; rr = 1.0; ss=0.0; nu = (n+p+1.0)/(p+1.0); kappa = (n+p+1.0)/n;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; zzPos= rr; ssPos = (ss+r2Qm)/(2.0*nu); 
      xxPos = 1.0/kappa - 1.0; 
      rCCH_void(v_, t_, l_, m_, aaPos, bbPos, zzPos, ssPos, xxPos);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 8){ // g = user input
      g_ = gg;
    } else if(gprior == 9){ // CH(a,b,s)
      nu = 1.0; 
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; ssPos = (ss+r2Qm)/(2.0*nu);
      rCH_void(v_, l_, aaPos, bbPos, ssPos);
      g_ = 1.0/(v_/nu) - 1.0;
    }
  } else if (familyLink == 11) { // Gaussian
    if(gprior == 0){ // unit information (g=n)
      gg = n; g_ = gg;
    } else if(gprior== 1){ // Hyper-g
      aa = 1.0; bb = 2.0; nu = 1.0;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; xxPos = (r2Qm/nu)/(1.0-r2Qm); zzPos = 0.5*(n-1.0);
      rGH_void(v_, t_, l_, aaPos, bbPos, xxPos, zzPos);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 2){ // Uniform
      aa = 2.0; bb = 2.0; nu = 1.0;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; xxPos = (r2Qm/nu)/(1.0-r2Qm); zzPos = 0.5*(n-1.0);
      rGH_void(v_, t_, l_, aaPos, bbPos, xxPos, zzPos);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 3){ // Hyper-g/n
      aa = 1.0; bb = 2.0; rr = 1.5; nu = 1.0; kappa = 1.0/n;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; zzPos = 0.5*(n-1.0); wwPos = rr; xxPos = (r2Qm/nu)/(1.0-r2Qm); yyPos = 1.0/kappa-1.0;
      rAPL_void(v_, t_,q_,l_,m_, aaPos, bbPos, zzPos, wwPos, xxPos, yyPos);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 4){ // Beta-prime
      aa = 0.5; bb = n-p-1.5; 
      g_ = (1.0 / R::rbeta(0.5*(0.5+p), 0.5*(n-p-1.5)) - 1.0) / (1.0 - r2Qm);
    } else if(gprior == 5){ // ZS-adapted
      aa = 1.0; bb = 2.0; ss = n+3.0; nu = 1.0;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; zzPos = (n-1.0)/2.0; ssPos = ss/(2.0*nu); xxPos =(r2Qm/nu)/(1.0-r2Qm); 
      rCCH_void(v_, t_, l_, m_, aaPos, bbPos, zzPos, ssPos, xxPos);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 6){ // Robust
      aa = 1.0; bb = 2.0; rr = 1.5; nu = (n+1.0)/(p+1.0); kappa = 1.0;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; xxPos = (r2Qm/nu)/(1.0-r2Qm); zzPos = 0.5*(n-1.0);
      rGH_void(v_, t_, l_, aaPos, bbPos, xxPos, zzPos);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 7){ // Intrinsic
      aa = 1.0; bb = 1.0; rr = 1.0; nu = (n+p+1.0)/(p+1.0); kappa = (n+p+1.0)/n;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; zzPos = 0.5*(n-1.0); 
      wwPos = rr; xxPos = (r2Qm/nu)/(1.0-r2Qm); yyPos = 1.0/kappa-1.0;
      rAPL_void(v_, t_,q_,l_,m_, aaPos, bbPos, zzPos, wwPos, xxPos, yyPos);
      g_ = 1.0/(v_/nu) - 1.0;
    } else if(gprior == 8){ // g = user input
      g_ = gg;
    } else if(gprior == 9){ // CH(a,b,s)
      nu = 1.0;
      aaPos = 0.5*(aa+p); bbPos = bb/2.0; zzPos = (n-1.0)/2.0; ssPos = ss/(2.0*nu); xxPos =(r2Qm/nu)/(1.0-r2Qm); 
      rCCH_void(v_, t_, l_, m_, aaPos, bbPos, zzPos, ssPos, xxPos);
      g_ = 1.0/(v_/nu) - 1.0;
    }
  } else {
    Rcpp::stop("gPosterior says, What da fuck?");
  }
  // 
  // Rcpp::Rcout << "v " << v_ << " t " << t_ << " q " << q_ << " l " << l_ << " m " << m_ <<"\n";
  // Rcpp::Rcout << "aaPos " << aaPos << " ssPos " << ssPos <<"\n";
  // Rcpp::Rcout << "r2Qm " << r2Qm << "\n";
   
}

////////////////////////////////////////////////////////////////////////////////
// posterior samples of alpha, beta
double AlphaPost(double AlphaHat, double JAlphaHat, double phi){
  return AlphaHat + arma::randn<double>() / std::sqrt(JAlphaHat*phi);
}

void basechol(arma::mat &toBeCholed,
              const Rcpp::Function &nearPDres, 
              const Rcpp::Function &Rbasechol){
  try{
    toBeCholed = Rcpp::as<arma::mat>(Rbasechol(toBeCholed));
  } catch(...) {
    toBeCholed = Rcpp::as<arma::mat>(Rbasechol(Rcpp::as<arma::mat>(nearPDres(toBeCholed))));
  }
}

arma::vec BetaPost(const arma::vec &BetaHat, 
                   const arma::mat &rootJBetaHat, 
                   const double& g, 
                   const double& phi,
                   const Rcpp::Function &nearPDres,
                   const Rcpp::Function &Rbasechol,
                   const Rcpp::Function &crossproduct)
{
  // convert sqrt(J)X to XJX and do chol to it
  arma::mat cholJBetaHat = Rcpp::as<arma::mat>(crossproduct(rootJBetaHat));
  // Rcpp::Rcout << "toBeCholed dim " << arma::size(cholJBetaHat) << "\n";
  try{basechol(cholJBetaHat, nearPDres, Rbasechol);} catch(...) {Rcpp::stop("cholesky failed in sampling from beta posterior");}
  arma::vec betaM(BetaHat.n_elem, arma::fill::randn);
  betaM = (g/(g+1.0)) * BetaHat + std::sqrt((g/((g+1.0) * phi))) * arma::solve(arma::trimatu(cholJBetaHat), betaM);
  return betaM;
}

arma::vec BetaPost_z(const arma::vec &BetaHat, 
                     const arma::mat &rootJBetaHat,
                     const double &g, 
                     const double &phi,
                     const arma::ivec &z,
                     const Rcpp::Function &nearPDres,
                     const Rcpp::Function &Rbasechol,
                     const Rcpp::Function &crossproduct)
{
  // Rcpp::Rcout << "doing BetaPost \n";
  // convert sqrt(J)X to XJX and to chol on it
  arma::mat cholJBetaHat = Rcpp::as<arma::mat>(crossproduct(rootJBetaHat));
  try{basechol(cholJBetaHat, nearPDres, Rbasechol);} catch(...) {Rcpp::stop("cholesky failed in sampling from beta posterior");}
  arma::vec betaM(BetaHat.n_elem, arma::fill::randn);
  betaM = (g/(g+1.0)) * BetaHat + std::sqrt((g/((g+1.0) * phi))) * arma::solve(arma::trimatu(cholJBetaHat), betaM);
  // augment betaM to befit z
  arma::vec beta(z.n_elem - 1, arma::fill::zeros);
  arma::uvec idx = arma::find(z == 1);
  beta.elem(idx.subvec(1, idx.n_elem - 1) - 1) = betaM;
  return beta;
}

				

























