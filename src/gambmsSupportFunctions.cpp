#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include "CRAD.h"
#include "bayesFactors.h"
#include "gambmsSupportFunctions.h"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// repPermu_arma (EVENKNOT)
// a recursion formula to list out all the possible knot combinations
void repPermu_void(const unsigned &maxk, 
                   const unsigned &P, 
                   arma::umat &OUT, 
                   unsigned depth, 
                   unsigned &rowidx){
  if(depth == P){
    rowidx +=1;
    if(rowidx < OUT.n_rows) OUT.row(rowidx) = OUT.row(rowidx-1);
    return;
  }
  for(unsigned i{0}; i<(maxk+1); i++){ // 0, 1, 2, ..., maxk (including zero)
    OUT(rowidx, depth) = i;
    repPermu_void(maxk, P, OUT, depth+1, rowidx);
  }
}

arma::umat repPermu_arma(const arma::uvec &maxk, 
                         const unsigned &P){
  
  unsigned msdim = std::pow(arma::max(maxk)+1, P);
  arma::umat OUT(msdim, P, arma::fill::ones);
  unsigned depth{0}, rowidx{0};
  repPermu_void(arma::max(maxk), P, OUT, depth, rowidx);
  
  // applying component-wise maxk restriction
  arma::uvec idx(msdim, arma::fill::ones);
  for(unsigned i{0}; i<msdim; i++){
    if(arma::any(OUT.row(i).t() > maxk)) { idx(i) -= 1; }
  }
  return OUT.rows(arma::find(idx > 0));
}


////////////////////////////////////////////////////////////////////////////////
// knotnums_to_idx (EVENKNOT)
// (3, 4, 1, 2) -> number in base maxk
int64_t knotnums_to_idx(const arma::ivec &knotnums, 
                        const unsigned &maxk){
  int64_t idx{0};
  for(unsigned j{0}; j < knotnums.n_elem; j++){
    idx += std::pow(maxk + 1, j) * knotnums(j); // include 0
  }
  return idx;
}

////////////////////////////////////////////////////////////////////////////////
// idx_to_knotnums (EVENKNOT)
// inverse of knotnums_to_idx
arma::ivec idx_to_knotnums(int64_t idx, 
                           const unsigned &knotnums_nelem, 
                           const unsigned &maxk)
{
  arma::ivec knotnums(knotnums_nelem, arma::fill::zeros);
  for(unsigned j{0}; j < knotnums_nelem; j++){
    knotnums(knotnums_nelem - 1 - j) = idx / std::pow(maxk+1, knotnums_nelem - 1 - j);
    idx -= knotnums(knotnums_nelem - 1 - j) * std::pow(maxk+1, knotnums_nelem - 1 - j);
  }
  return knotnums;
}

////////////////////////////////////////////////////////////////////////////////
// deprecated: will not store for VSKNOT
// int64_t z_to_longidx(const arma::ivec &z){
//   int64_t longidx{0};
//   for(unsigned i{0}; i<z.n_elem; i++){
//     longidx += std::pow(2, i) * z(i);
//   }
//   if(longidx < 0){
//     Rcpp::Rcout<< "Numerical overflow for __int64 for z.n_elem = " << z.n_elem << "\n";
//     Rcpp::stop("");
//   }
//   // Rcpp::Rcout<< "longidx " << longidx << "\n";
//   return longidx;
// }
// 
// arma::ivec longidx_to_z(int64_t longidx,
//                         const unsigned &z_nelem){
//   arma::ivec z(z_nelem, arma::fill::zeros);
//   int64_t tmp{0};
//   for(unsigned i{0}; i < z.n_elem; i++){
//     tmp = std::pow(2, z.n_elem - 1 - i);
//     z(z.n_elem - 1 - i) = longidx / tmp;
//     longidx -= z(z.n_elem - 1 - i) * tmp;
//   }
//   return z;
// }

////////////////////////////////////////////////////////////////////////////////
// knotsandidx (EVENKNOT)
// knotnums to location and index of knots
void knotsandidx(arma::vec &knots,
                 arma::uvec &knotsidx,
                 const arma::uvec &knotnums,
                 const arma::mat &X){
  // Rcpp::Rcout << "knotnums " << knotnums.t() << "\n";
  
  knots.set_size(arma::accu(knotnums));
  knots.zeros();
  knotsidx.set_size(arma::accu(knotnums));
  knotsidx.zeros();
  double bdmargin{.0};
  
  arma::vec qts;
  unsigned knotnum;
  unsigned start{0}, end{0};
  for(unsigned k{0}; k < knotnums.n_elem; k++){
    // Rcpp::Rcout << k << "\n";
    knotnum = knotnums(k);
    if(knotnum == 0){ // no knot selected
      continue;
    } else {
      qts = arma::linspace(bdmargin, 1.0 - bdmargin, knotnum + 2); // all predictors rescaled to [0,1] in the first place
      qts = qts.subvec(1, knotnum); // n_elem = k
      end = start + knotnum - 1;
      knots.subvec(start, end) = arma::quantile(arma::unique(X.col(k)), qts);
      knotsidx.subvec(start, end) += (k + 1);
      start = end + 1;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// knotprior 
// componentwise independent Poisson prior for the number of knots
void knotPrior(double &comp3,
               const arma::uvec &betaidx,
               const arma::uvec &maxk,
               const arma::vec &Lambda,
               const bool& isgrid){
  comp3 = 0.0; 
  double knotnums{.0};
  bool LambdaPos = (arma::accu(Lambda) > 0);
  for(unsigned i{0}; i<Lambda.n_elem; i++){
    knotnums = (arma::accu(betaidx == (i+1)) - 1);
    if(isgrid) {comp3 -= (std::lgamma(maxk(i) + 1.0) - std::lgamma(maxk(i) - knotnums + 1.0) - std::lgamma(knotnums + 1.0));}
    if(LambdaPos) {comp3 += knotnums * std::log(Lambda(i)) - Lambda(i) - std::lgamma(knotnums + 1);}
  }
}



////////////////////////////////////////////////////////////////////////////////
// MATX_TO_LPY
// from B_tr the design matrix, yield model evidence, mle and rootJBetaHat 
// input: y, B_tr, offset, glmWeight, 
//        knotnums,
//        Rglm(y, B_tr, offset, glmWeight, EtaHat)
//        Lambda, familyLink, gprior, aa, bb, ss, gg
// output: lpy, comp1, comp2, comp3, r2Qm, mle, rootJBetaHat
// comp1: related to model fit (r2 for Gaussian)
// comp2: related to penalty, differ by gprior
// comp3: related to prior on the number of knots
void MATX_TO_LPY(double &lpy,
                 double &comp1,
                 double &comp2,
                 double &comp3,
                 double &r2Qm,
                 arma::vec &mle,
                 arma::vec &EtaHat,
                 arma::mat &rootJBetaHat,
                 const arma::vec &y,
                 const double &glmWeight,
                 const arma::mat &B_tr,
                 const arma::uvec &betaidx,
                 const arma::vec &offset,
                 const arma::vec &Lambda,
                 const unsigned& familyLink,
                 const unsigned& gprior,
                 const double &aa, 
                 const double &bb, 
                 const double &ss, 
                 const double &gg,
                 const Rcpp::Function &Rglm,
                 const bool& isgrid,
                 const arma::uvec &maxk){
  
  // arma::vec knots;
  // arma::uvec knotsidx;
  // knotsandidx(knots, knotsidx, knotnums, X);
  // Rcpp::List CRAD_tr = CRAD_cpp(X, XLin, knots, knotsidx, true, 0.0);
  // arma::mat B_tr = Rcpp::as<arma::mat>(CRAD_tr["X"]); // expanded for train, intercept term included
  // Rcout << "trying mle \n";
  // Rcout << "y.n_elem " << y.n_elem << "\n";
  // Rcout << "size(B_tr) " << arma::size(B_tr) << "\n";
  try{
    // Rcpp::Rcout << "colSums(B_tr) " << arma::sum(B_tr) << "\n";
    mle = Rcpp::as<arma::vec>(Rglm(y, B_tr, offset, glmWeight, EtaHat));
    // Rcpp::Rcout << mle.t() << "\n";
  } catch(...) {
    // Rcpp::Rcout << "mle (base::glm.fit) error \n";
    // Rcpp::Rcout << EtaHat.t() << "\n";
    // Rcpp::Rcout << glmWeight << "\n";
    // Rcpp::Rcout << familyLink << "\n";
    // Rcpp::stop("");
  }

  
  double AlphaHat = mle(0);
  EtaHat = AlphaHat + B_tr.cols(1, B_tr.n_cols-1) * mle.subvec(1, mle.n_elem-1) + offset;
  arma::vec JEtaHat = JEta(glmWeight, y, EtaHat, familyLink);
  double JAlphaHat = JAlpha(JEtaHat);
  rootJBetaHat = rootJBeta(mle.subvec(1, mle.n_elem-1), B_tr, JEtaHat, familyLink);
  
  // Rcpp::Rcout << "trying r2Qm \n";
  r2Qm = getR2QM(familyLink, y, offset, AlphaHat,EtaHat, JEtaHat);
  // Rcpp::Rcout << "trying Lpy \n";
  Lpy(lpy, comp1, comp2, glmWeight, 
      y, B_tr.cols(1, B_tr.n_cols-1), offset, EtaHat, JAlphaHat, r2Qm,
      familyLink, gprior, aa, bb, ss, gg);
  // Rcpp::Rcout << "trying knotPrior \n";
  knotPrior(comp3, betaidx, maxk, Lambda, isgrid);
  lpy += comp3;
  // Rcpp::Rcout << "lpy done \n";
}

void KNOT_TO_LPY(const arma::uvec &knotnums,
                 arma::vec &knots,
                 arma::uvec &knotsidx,
                 double &lpy,
                 double &comp1,
                 double &comp2,
                 double &comp3,
                 double &r2Qm,
                 arma::vec &mle,
                 arma::vec &EtaHat,
                 arma::mat &rootJBetaHat,
                 const arma::vec &y,
                 const double &glmWeight,
                 const arma::mat &X,
                 const arma::mat &XLin,
                 const arma::vec &offset,
                 const arma::vec &Lambda,
                 const unsigned& familyLink,
                 const unsigned& gprior,
                 const double &aa, 
                 const double &bb, 
                 const double &ss, 
                 const double &gg,
                 const Rcpp::Function &Rglm,
                 const bool& isgrid,
                 const arma::uvec &maxk){
  
  bool NS{true};
  knotsandidx(knots, knotsidx, knotnums, X);
  Rcpp::List CRAD_tr = CRAD_cpp(X, XLin, knots, knotsidx, NS, 0.0);
  arma::mat B_tr = Rcpp::as<arma::mat>(CRAD_tr["X"]); // expanded for train, intercept term included
  arma::uvec betaidx = Rcpp::as<arma::uvec>(CRAD_tr["betaidx"]);
  // betaidx.shed_row(0);
  
  // Rcpp::Rcout << knots.t() << "\n";
  // Rcpp::Rcout << B_tr << "\n";
  // Rcpp::Rcout << "KNOT_TO_LPY basis expansion done \n";
  MATX_TO_LPY(lpy, comp1, comp2, comp3, r2Qm, mle, EtaHat, rootJBetaHat,
              y, glmWeight, B_tr, betaidx, offset, Lambda, 
              familyLink, gprior, aa, bb, ss, gg,
              Rglm, isgrid, maxk);
}

////////////////////////////////////////////////////////////////////////////////
// MATX_to_SAMPLE
// from (mle, EtaHat, rootJBetaHat, r2Qm), produce new posterior sample
// input: y, glmWeight, mle, EtaHat, rootJBetaHat, r2Qm, 
//        B_tr, B_pr, betaidx, 
//        familyLink, gprior, aa, bb, ss, gg
// output: PredSmooths, PredLinears, phi, g
// we do reuse stored mle rootJBetaHat but recreate CRAD basis as
// it is too costly memory-wise to store B_tr for every knot configuration
void MATX_TO_SAMPLE(arma::vec &FittedSmooths,
                    arma::vec &PredSmooths,
                    arma::vec &PredLinears,
                    double &phi,
                    double &g_,
                    double &v_,
                    double &t_,
                    double &q_,
                    double &l_,
                    double &m_,
                    const arma::vec &y,
                    const double &glmWeight,
                    const arma::vec &mle,
                    const arma::vec &EtaHat,
                    const arma::mat &rootJBetaHat,
                    const double &r2Qm,
                    const arma::mat &B_tr,
                    const arma::mat &B_pr,
                    const arma::uvec &betaidx,
                    const double &P,
                    const double &PLin,
                    const unsigned& familyLink,
                    const unsigned& gprior,
                    const double &aa, 
                    const double &bb, 
                    const double &ss, 
                    const double &gg,
                    const Rcpp::Function &nearPDres,
                    const Rcpp::Function &crossproduct,
                    const Rcpp::Function &Rbasechol,
                    const bool &storeFit){
  
  unsigned N{B_tr.n_rows};
  unsigned np{B_pr.n_rows};
  
  // PredSmooths.set_size(P*np);
  FittedSmooths.zeros(); PredSmooths.zeros();
  // PredLinears.set_size(PLin);
  PredLinears.zeros();
  
  // arma::vec knots;
  // arma::uvec knotsidx;
  // knotsandidx(knots, knotsidx, knotnums, X);
  // Rcpp::List CRAD_tr = CRAD_cpp(X, XLin, knots, knotsidx, true, 0.0);
  // arma::mat B_tr = Rcpp::as<arma::mat>(CRAD_tr["X"]); // expanded for train, intercept term included
  // arma::mat B_pr = CRAD_test_cpp(X_pr, XLin, CRAD_tr);
  
  // arma::vec EtaHat = mle(0) + B_tr.cols(1, B_tr.n_cols-1) * mle.subvec(1, mle.n_elem-1) + offset;
  arma::vec JEtaHat = JEta(glmWeight, y, EtaHat, familyLink);
  double JAlphaHat = JAlpha(JEtaHat);
  
  // 1. g posterior
  gPosteriors(g_, v_, t_, q_, l_, m_,
              static_cast<double>(N), static_cast<double>(P), 
              r2Qm, familyLink, gprior, 
              aa, bb, ss, gg);
  
  // 2. phi, alpha, beta posteriors
  if(familyLink == 11){
    phi = R::rgamma(static_cast<double>(N) / 2.0, 
                    2.0 / (std::pow(arma::norm(y-arma::as_scalar(arma::mean(y)),2),2.0) * (1.0-r2Qm+r2Qm/(g_+1.0))));
  } else {phi = 1.0;}
  double alpha = AlphaPost(mle(0), JAlphaHat, phi);
  // Rcpp::Rcout << "getting betaposterior \n";
  // Rcpp::Rcout << "mle dim " << arma::size(mle) << "\n";
  // Rcpp::Rcout << "rootJBetaHat dim " << arma::size(rootJBetaHat) << "\n";
  arma::vec beta = BetaPost(mle.subvec(1, mle.n_elem-1),
                            rootJBetaHat,
                            g_, phi, nearPDres, Rbasechol, crossproduct);
  // arma::uvec betaidx = Rcpp::as<arma::uvec>(CRAD_tr["betaidx"]); // betaidx includes intercept term
  // betaidx.shed_row(0); // subtract the first to match length
  // Rcpp::Rcout << "got betaposterior \n";
  // 3. make prediction, store results
  // PredLinears.fill(alpha);
  PredLinears(0) = alpha;
  if(PLin > 1){PredLinears.subvec(1, PLin-1) = beta.subvec(0, PLin-1-1);}
  unsigned start_fit = 0, end_fit = 0, start_prd = 0, end_prd = 0;
  arma::vec beta_sub;
  // Rcpp::Rcout << "got PredLinears \n";
  for(unsigned i{0}; i < P; i++){
    // Rcpp::Rcout << "i " << i << "\n";
    beta_sub = beta;
    beta_sub.elem(arma::find(betaidx.subvec(1, betaidx.n_elem - 1) != (i+1))).zeros(); 
    
    end_fit = start_fit + N - 1; end_prd = start_prd + np - 1;
    if(storeFit) FittedSmooths(arma::span(start_fit, end_fit)) = (B_tr.cols(1, B_tr.n_cols - 1) * beta_sub);
    PredSmooths(arma::span(start_prd, end_prd)) = (B_pr.cols(1, B_pr.n_cols - 1) * beta_sub);
    start_prd = end_prd + 1; start_fit = end_fit + 1;
  }
  // Rcpp::Rcout << "got predicted \n";
}

void KNOT_TO_SAMPLE(const arma::uvec &knotnums,
                    arma::vec &knots,
                    arma::uvec &knotsidx,
                    arma::vec &FittedSmooths,
                    arma::vec &PredSmooths,
                    arma::vec &PredLinears,
                    double &phi,
                    double &g_,
                    double &v_,
                    double &t_,
                    double &q_,
                    double &l_,
                    double &m_,
                    const arma::vec &y,
                    const double &glmWeight,
                    const arma::mat &X,
                    const arma::mat &X_pr,
                    const arma::mat &XLin,
                    const arma::vec &mle,
                    const arma::vec &EtaHat,
                    const arma::mat &rootJBetaHat,
                    const double &r2Qm,
                    const unsigned& familyLink,
                    const unsigned& gprior,
                    const double &aa, 
                    const double &bb, 
                    const double &ss, 
                    const double &gg,
                    const Rcpp::Function &nearPDres,
                    const Rcpp::Function &crossproduct,
                    const Rcpp::Function &Rbasechol,
                    const bool& storeFit){
  
  
  knotsandidx(knots, knotsidx, knotnums, X);
  // Rcpp::Rcout << "knotsandidx done\n";
  Rcpp::List CRAD_tr = CRAD_cpp(X, XLin, knots, knotsidx, true, 0.0);
  // Rcpp::Rcout << "CRAD_tr done\n";
  arma::uvec betaidx = Rcpp::as<arma::uvec>(CRAD_tr["betaidx"]);
  arma::mat B_tr = Rcpp::as<arma::mat>(CRAD_tr["X"]); 
  arma::mat B_pr = CRAD_test_cpp(X_pr, XLin, CRAD_tr);
  // Rcpp::Rcout << "B_pr done\n";
  
  // Rcpp::Rcout << "getting into MATX_TO_SAMPLE \n";
  MATX_TO_SAMPLE(FittedSmooths, PredSmooths, PredLinears,
                 phi, g_, v_, t_, q_, l_, m_,
                 y, glmWeight, mle, EtaHat, rootJBetaHat, r2Qm,
                 B_tr, B_pr, betaidx, X.n_cols, XLin.n_cols,
                 familyLink, gprior, 
                 aa, bb, ss, gg,
                 nearPDres, crossproduct, Rbasechol, storeFit);
}


////////////////////////////////////////////////////////////////////////////////
// MC_POSTERIOR (EVENKNOT)
// from list of (knotnums, mle, EtaHat, rootJBetaHat, r2Qm) and WEIGHT, produce MCiter posterior samples 
// input: y, glmWeight, WEIGHT,
//        list of (knotnums, mle, EtaHat, rootJBetaHat, r2Qm), 
//        familyLink, gprior, aa, bb, ss, gg
// output: list of (PredSmooths, PredLinears, phi, g), KNOTLOCS, KNOTLOCSIDX
// we do reuse stored mle rootJBetaHat but recreate CRAD basis as
// it is too costly memory-wise to store B_tr for every knot configuration

void MC_POSTERIOR(arma::mat &FITTEDSMOOTHS,
                  arma::mat &PREDSMOOTHS,
                  arma::mat &PREDLINEARS,
                  Rcpp::List &MCKNOTLOCS,
                  Rcpp::List &MCKNOTLOCSIDX,
                  arma::vec &PHI,
                  arma::vec &G,
                  const unsigned &numMCcandidate,
                  const unsigned &MCiter,
                  const arma::vec &WEIGHT,
                  const arma::vec &y,
                  const double &glmWeight,
                  const arma::mat &X,
                  const arma::mat &X_pr,
                  const arma::mat &XLin,
                  const arma::umat &KNOTNUMS_MCcandidate,
                  const Rcpp::List &MLE,
                  const arma::mat &ETAHAT,
                  const Rcpp::List &ROOTJBETAHAT,
                  const arma::vec &R2QM,
                  const unsigned& familyLink,
                  const unsigned& gprior,
                  const double &aa, 
                  const double &bb, 
                  const double &ss, 
                  const double &gg,
                  const Rcpp::Function &nearPDres,
                  const Rcpp::Function &crossproduct,
                  const Rcpp::Function &Rbasechol,
                  const bool& storeFit){
  
  unsigned P{X.n_cols}; // num of predictors
  unsigned N{X.n_rows};
  unsigned PLin{XLin.n_cols};
  unsigned np{X_pr.n_rows};
  
  arma::vec knots;
  arma::uvec knotsidx, knotnums, betaidx;
  Rcpp::List CRAD_tr;
  arma::mat B_tr, B_pr;
  arma::vec mle, EtaHat;
  arma::mat rootJBetaHat;
  // double r2Qm;
  
  double phi{1.0}, g_{static_cast<double>(N)};
  double v_{0.1}, t_{0.1}, q_{0.1}, l_{0.1}, m_{0.1};
  arma::vec FittedSmooths, PredSmooths, PredLinears;
  PredSmooths.set_size(P*np); PredLinears.set_size(PLin);FittedSmooths.set_size(P*N);
  
  arma::uvec orders = Rcpp::RcppArmadillo::sample(
    arma::regspace<arma::uvec>(0, numMCcandidate - 1), 
    MCiter, true, WEIGHT);
  unsigned idx;
  
  // Rcpp::Rcout << "got orders \n";
  for(unsigned s{0}; s<MCiter; s++){
    idx = orders(s);
    
    // knotnums = KNOTNUMS_MCcandidate.row(idx).t();
    // knotsandidx(knots, knotsidx, knotnums, X);
    // CRAD_tr = CRAD_cpp(X, XLin, knots, knotsidx, true, 0.0);
    // betaidx = Rcpp::as<arma::uvec>(CRAD_tr["betaidx"]);
    // B_tr = Rcpp::as<arma::mat>(CRAD_tr["X"]); 
    // B_pr = CRAD_test_cpp(X_pr, XLin, CRAD_tr);
    // mle = Rcpp::as<arma::vec>(MLE(idx));
    // EtaHat = ETAHAT.col(idx);
    // rootJBetahat = Rcpp::as<arma::mat>(ROOTJBETAHAT(idx));
    // r2Qm = R2QM(idx);
    
    // MATX_TO_SAMPLE(PredSmooths, PredLinears,
    //                phi, g_, v_, t_, q_, l_, m_,
    //                y, glmWeight, mle, EtaHat, rootJBetaHat, r2Qm,
    //                B_tr, B_pr, betaidx, X.n_cols, XLin.n_cols,
    //                familyLink, gprior, 
    //                aa, bb, ss, gg,
    //                nearPDres, crossproduct, Rbasechol);
    
    // Rcpp::Rcout << "getting into KNOT_TO_SAMPLE \n";
    // Rcpp::Rcout << KNOTNUMS_MCcandidate.col(idx) << "\n";
    KNOT_TO_SAMPLE(KNOTNUMS_MCcandidate.col(idx), 
                   knots, knotsidx,
                   FittedSmooths, PredSmooths, PredLinears,
                   phi, g_, v_, t_, q_, l_, m_, 
                   y, glmWeight, X, X_pr, XLin, 
                   Rcpp::as<arma::vec>(MLE(idx)), 
                   ETAHAT.col(idx),
                   Rcpp::as<arma::mat>(ROOTJBETAHAT(idx)),
                   R2QM(idx),
                   familyLink, gprior, aa, bb, ss, gg,
                   nearPDres, crossproduct, Rbasechol,
                   storeFit);
    
    FITTEDSMOOTHS.row(s) = FittedSmooths.t();
    PREDSMOOTHS.row(s) = PredSmooths.t();
    PREDLINEARS.row(s) = PredLinears.t();
    MCKNOTLOCS(s) = knots;
    MCKNOTLOCSIDX(s) = knotsidx;
    G(s) = g_; PHI(s) = phi;
  }
}

////////////////////////////////////////////////////////////////////////////////
// mark non-knot basis terms (linear etc) (VSKNOT)
arma::ivec markIsKont(const arma::uvec &betaidx, 
                      const unsigned &howmanybasisterms){
  arma::ivec isknot = arma::ivec(betaidx.n_elem, arma::fill::ones);
  isknot.elem(arma::find(betaidx < 1)) -= 1;
  for(unsigned i{1}; i< isknot.n_elem; i++){
    if(betaidx(i) > 0){
      if((betaidx(i) - betaidx(i-1)) > 0){
        for(unsigned j{0}; j<howmanybasisterms; j++) isknot(i+j) = 0;
      } 
    }
  }
  return isknot;
}

////////////////////////////////////////////////////////////////////////////////
// z to knots (VSKNOT, FREEKNOT)
void zToKnots(arma::vec &knotsZ, 
              arma::uvec &knotsidxZ,
              arma::uvec &betaidxZ,
              arma::mat &B_trZ,
              arma::mat &B_prZ,
              const arma::ivec &z,
              const arma::uvec &betaidx,
              const arma::mat &B_tr,
              const arma::mat &B_pr,
              const arma::ivec &isknot,
              const arma::vec &knots,
              const arma::uvec &knotsidx){
  
  knotsZ.set_size(isknot.n_elem);
  knotsZ.zeros();
  knotsidxZ.set_size(isknot.n_elem);
  knotsidxZ.zeros();
  knotsZ.elem(arma::find(isknot > 0)) = knots;
  knotsZ = knotsZ.elem(arma::find(z > 0));
  knotsidxZ.elem(arma::find(isknot > 0)) = knotsidx;
  knotsidxZ = knotsidxZ.elem(arma::find(z > 0));
  
  betaidxZ = betaidx.elem(arma::find(z > 0));
  B_trZ = B_tr.cols(arma::find(z > 0));
  B_prZ = B_pr.cols(arma::find(z > 0));
}















































